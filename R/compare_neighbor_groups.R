#' Compare neighbour composition at sample or FOV level
#'
#' Provides three explicitly labelled analysis strategies for grouped FOV
#' neighbour summaries. `"sample_wilcox"` first averages FOVs within each
#' biological sample. `"fov_wilcox"` treats FOVs as analysis units and is
#' intended as an exploratory reproduction of the original workflow.
#' `"cluster_permutation"` retains the FOV input but permutes group labels for
#' whole biological samples; all FOVs from a sample always move together.
#'
#' @param data Long FOV-level neighbour summary containing sample, FOV, group,
#'   neighbour cell type, and numeric value columns.
#' @param analysis_level Either `"sample"` or `"FOV"`.
#' @param method One of `"sample_wilcox"`, `"fov_wilcox"`, or
#'   `"cluster_permutation"`. The method must match `analysis_level`.
#' @param value_col Numeric neighbour statistic, usually `"mean_prop"`.
#' @param group_col,sample_col,fov_col,celltype_col Column names.
#' @param group_order Optional two-element vector fixing comparison direction;
#'   the reported effect is first group minus second group.
#' @param p_adjust Multiple-testing correction passed to [stats::p.adjust()].
#'   Use `"none"` to retain raw P values only.
#' @param wilcox_exact Passed to [stats::wilcox.test()] for Wilcoxon methods.
#' @param permutations Number of Monte Carlo permutations when the number of
#'   possible sample-label allocations exceeds `max_exact`.
#' @param max_exact Maximum number of allocations to enumerate exactly.
#' @param seed Random seed used only for Monte Carlo permutations.
#'
#' @return A data frame with one row per neighbour cell type, group estimates,
#'   effect, P values, sample/FOV counts, method labels, and permutation count.
#'
#' @details
#' Ordinary FOV-level Wilcoxon ignores dependence between FOVs from the same
#' sample and can underestimate uncertainty. The cluster-permutation method
#' uses an equal-sample-weighted difference in mean FOV values and permutes
#' labels only between biological samples. It requires samples to be
#' independent and exchangeable under the null hypothesis.
#'
#' Methodological references:
#'
#' - Winkler AM et al. (2014). Permutation inference for the general linear
#'   model. NeuroImage 92:381-397. doi:10.1016/j.neuroimage.2014.01.060.
#' - Rosner B, Glynn RJ, Lee MLT (2003). Incorporation of clustering effects
#'   for the Wilcoxon rank sum test. Biometrics 59:1089-1098.
#' - Datta S, Satten GA (2005). Rank-sum tests for clustered data. JASA
#'   100:908-915. doi:10.1198/016214504000001583.
#'
#' @export
#'
#' @examples
#' d <- expand.grid(
#'   Sample = paste0("S", 1:6), FOV = paste0("F", 1:3),
#'   CellType = c("cDC", "CD8Tn"), stringsAsFactors = FALSE
#' )
#' d$Group <- ifelse(d$Sample %in% paste0("S", 1:3), "R", "NR")
#' d$mean_prop <- seq_len(nrow(d)) / nrow(d)
#' compare_neighbor_groups(d, "sample", "sample_wilcox")
#' compare_neighbor_groups(d, "FOV", "cluster_permutation")
compare_neighbors <- function(
    data,
    analysis_level = c("sample", "FOV"),
    method = NULL,
    value_col = "mean_prop",
    group_col = "Group",
    sample_col = "Sample",
    fov_col = "FOV",
    celltype_col = "CellType",
    group_order = NULL,
    p_adjust = "BH",
    wilcox_exact = FALSE,
    permutations = 9999L,
    max_exact = 100000L,
    seed = 2608L) {
  analysis_level <- match.arg(analysis_level)
  if (is.null(method)) {
    method <- if (analysis_level == "sample") {
      "sample_wilcox"
    } else {
      "cluster_permutation"
    }
  }
  method <- match.arg(
    method, c("sample_wilcox", "fov_wilcox", "cluster_permutation")
  )
  if (analysis_level == "sample" && method != "sample_wilcox") {
    stop("Sample-level analysis requires `method = 'sample_wilcox'`.",
         call. = FALSE)
  }
  if (analysis_level == "FOV" && method == "sample_wilcox") {
    stop("FOV-level analysis requires `fov_wilcox` or `cluster_permutation`.",
         call. = FALSE)
  }

  vc <- .find_col(data, value_col)
  gc <- .find_col(data, group_col)
  sc <- .find_col(data, sample_col)
  cc <- .find_col(data, celltype_col)
  fc <- .find_col(data, fov_col)
  if (any(vapply(list(vc, gc, sc, cc), is.null, logical(1)))) {
    stop("Could not find value, group, sample, and/or cell-type columns.",
         call. = FALSE)
  }
  if (!is.numeric(data[[vc]]) || !nrow(data)) {
    stop("`value_col` must be numeric and `data` must not be empty.",
         call. = FALSE)
  }
  if (analysis_level == "FOV" && is.null(fc)) {
    stop("FOV-level analysis requires an FOV column.", call. = FALSE)
  }

  d <- data[
    is.finite(data[[vc]]) & !is.na(data[[gc]]) & !is.na(data[[sc]]) &
      !is.na(data[[cc]]), , drop = FALSE
  ]
  sample_groups <- unique(d[, c(sc, gc), drop = FALSE])
  if (anyDuplicated(sample_groups[[sc]])) {
    stop("Every sample must belong to exactly one group.", call. = FALSE)
  }
  groups <- if (is.null(group_order)) unique(as.character(d[[gc]])) else
    as.character(group_order)
  if (length(groups) != 2L || !all(groups %in% d[[gc]])) {
    stop("Exactly two valid groups are required.", call. = FALSE)
  }

  d$.__group <- factor(as.character(d[[gc]]), levels = groups)
  d$.__sample <- as.character(d[[sc]])
  d$.__celltype <- as.character(d[[cc]])
  d$.__value <- d[[vc]]
  if (!is.null(fc)) d$.__fov <- as.character(d[[fc]])

  sample_values <- d |>
    dplyr::group_by(.data$.__sample, .data$.__group, .data$.__celltype) |>
    dplyr::summarise(.__value = mean(.data$.__value), .groups = "drop")
  celltypes <- unique(d$.__celltype)

  compare_one <- function(celltype) {
    raw <- d[d$.__celltype == celltype, , drop = FALSE]
    sv <- sample_values[sample_values$.__celltype == celltype, , drop = FALSE]
    n_samples <- table(factor(sv$.__group, levels = groups))
    if (any(n_samples < 1L)) {
      stop("Every cell type must contain observations from both groups.",
           call. = FALSE)
    }
    n_fovs <- if (is.null(fc)) c(NA_integer_, NA_integer_) else {
      vapply(groups, function(g) {
        rows <- raw$.__group == g
        nrow(unique(raw[rows, c(".__sample", ".__fov"), drop = FALSE]))
      }, integer(1))
    }

    if (method == "sample_wilcox") {
      x <- sv$.__value[sv$.__group == groups[1]]
      y <- sv$.__value[sv$.__group == groups[2]]
      test <- suppressWarnings(stats::wilcox.test(x, y, exact = wilcox_exact))
      estimate <- c(mean(x), mean(y))
      statistic <- unname(test$statistic)
      p_value <- test$p.value
      n_perm <- NA_integer_
      label <- paste0(
        "Wilcoxon rank-sum; analysis unit: biological sample; exact=",
        wilcox_exact
      )
      permutation_unit <- NA_character_
    } else if (method == "fov_wilcox") {
      x <- raw$.__value[raw$.__group == groups[1]]
      y <- raw$.__value[raw$.__group == groups[2]]
      test <- suppressWarnings(stats::wilcox.test(x, y, exact = wilcox_exact))
      estimate <- c(mean(x), mean(y))
      statistic <- unname(test$statistic)
      p_value <- test$p.value
      n_perm <- NA_integer_
      label <- paste0(
        "Exploratory Wilcoxon rank-sum; analysis unit: FOV; exact=",
        wilcox_exact, "; FOV nesting not adjusted"
      )
      permutation_unit <- NA_character_
    } else {
      sample_ids <- sv$.__sample
      sample_group <- as.character(sv$.__group)
      values <- sv$.__value
      n_group1 <- sum(sample_group == groups[1])
      observed <- mean(values[sample_group == groups[1]]) -
        mean(values[sample_group == groups[2]])
      total_allocations <- choose(length(sample_ids), n_group1)

      if (total_allocations <= max_exact) {
        allocations <- utils::combn(seq_along(sample_ids), n_group1)
        perm_stats <- apply(allocations, 2, function(idx) {
          mean(values[idx]) - mean(values[-idx])
        })
        p_value <- mean(abs(perm_stats) >= abs(observed) - sqrt(.Machine$double.eps))
        n_perm <- length(perm_stats)
        label <- paste0(
          "Exact sample-block permutation; analysis input: FOV; ",
          "permutation unit: biological sample; equal sample weights"
        )
      } else {
        if (!is.numeric(permutations) || length(permutations) != 1L ||
            permutations < 1) {
          stop("`permutations` must be a positive integer.", call. = FALSE)
        }
        set.seed(seed)
        perm_stats <- replicate(as.integer(permutations), {
          idx <- sample(seq_along(sample_ids), n_group1)
          mean(values[idx]) - mean(values[-idx])
        })
        p_value <- (1 + sum(abs(perm_stats) >= abs(observed))) /
          (length(perm_stats) + 1)
        n_perm <- length(perm_stats)
        label <- paste0(
          "Monte Carlo sample-block permutation; analysis input: FOV; ",
          "permutation unit: biological sample; equal sample weights"
        )
      }
      estimate <- c(
        mean(values[sample_group == groups[1]]),
        mean(values[sample_group == groups[2]])
      )
      statistic <- observed
      permutation_unit <- "biological sample"
    }

    data.frame(
      CellType = celltype,
      group1 = groups[1], group2 = groups[2],
      estimate1 = estimate[1], estimate2 = estimate[2],
      effect = estimate[1] - estimate[2], statistic = statistic,
      p.value = p_value,
      n_sample1 = unname(n_samples[1]), n_sample2 = unname(n_samples[2]),
      n_fov1 = n_fovs[1], n_fov2 = n_fovs[2],
      analysis_level = analysis_level, method = method,
      method_label = label, permutation_unit = permutation_unit,
      n_permutations = n_perm, stringsAsFactors = FALSE
    )
  }

  out <- dplyr::bind_rows(lapply(celltypes, compare_one))
  out$p.adjust <- if (identical(p_adjust, "none")) out$p.value else
    stats::p.adjust(out$p.value, method = p_adjust)
  out$p_adjust_method <- p_adjust
  out
}

#' @rdname compare_neighbors
#' @export
compare_neighbor_groups <- compare_neighbors
