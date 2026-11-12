## Make sure we run from the package root regardless of the caller's cwd
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grepl("^--file=", args)])
if (length(script_path) && nzchar(script_path)) {
  setwd(dirname(dirname(normalizePath(script_path))))
}
## Generate the small example dataset shipped with cellfiers
## Run once from the package root:  Rscript data-raw/make_example_data.R
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

set.seed(2026)

## ---- 1. Rule table (same structure as Docu/SLN_mIF_anno_rule.txt) ----
mif_rules <- tibble::tribble(
  ~CellType,            ~CD11c, ~CD3, ~CD8, ~PDL1, ~CD45RA, ~CD69,
  "cDC",                "+",    "-",  "-",  "-",   NA,      NA,
  "cDC(Mature/Active)", "+",    "-",  "-",  "+",   NA,      NA,
  "CD8Tn(Quiescent)",   "-",    "+",  "+",  NA,    "+",     "-",
  "CD8Tn(Primed)",      "-",    "+",  "+",  NA,    "+",     "+",
  "CD8Teff(Resting)",   "-",    "+",  "+",  NA,    "-",     "-",
  "CD8Teff(Active)",    "-",    "+",  "+",  NA,    "-",     "+",
  "CD4Tn(Quiescent)",   "-",    "+",  "-",  NA,    "+",     "-",
  "CD4Tn(Primed)",      "-",    "+",  "-",  NA,    "+",     "+",
  "CD4Teff(Resting)",   "-",    "+",  "-",  NA,    "-",     "-",
  "CD4Teff(Active)",    "-",    "+",  "-",  NA,    "-",     "+"
)

## ---- 2. Per-sample cutoffs (structure of Docu/SLN_mIF_6channel_cutoff.txt) ----
mif_cutoffs <- tibble::tribble(
  ~SampleID, ~CD11c, ~CD69, ~CD3, ~PDL1, ~CD8, ~CD45RA,
  "C1001",   450,    320,   650,  480,   420,  600,
  "C1002",   430,    360,   600,  520,   400,  640,
  "C1003",   470,    300,   680,  500,   440,  580
)

## ---- 3. Group table (structure of Docu/SLN_mIF_Group_pCRpNR.txt) ----
mif_groups <- tibble::tribble(
  ~Sample, ~Group,
  "C1001", "R",
  "C1002", "NR",
  "C1003", "R"
)

## ---- 4. Simulated single-cell data ----
marker_map <- c(
  CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
  CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
  PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
)
markers <- names(marker_map)
intensity_cols <- unname(marker_map)

simulate_sample <- function(sample_id, n = 150, prop_other = 0.08) {
  types <- sample(mif_rules$CellType, n, replace = TRUE,
                  prob = c(0.06, 0.06, 0.12, 0.08, 0.10, 0.08,
                           0.12, 0.08, 0.12, 0.10))
  # a few "Other" (unassigned) cells
  other <- runif(n) < prop_other
  types[other] <- "Other"

  rule_idx <- match(types, mif_rules$CellType)
  dat <- matrix(NA_real_, nrow = n, ncol = length(markers),
                dimnames = list(NULL, markers))
  for (m in seq_along(markers)) {
    status <- rep(NA_character_, n)
    valid <- !is.na(rule_idx)
    status[valid] <- vapply(mif_rules[[markers[m]]][rule_idx[valid]],
                            as.character, character(1))
    for (i in seq_len(n)) {
      st <- status[i]
      if (is.na(st)) {
        dat[i, m] <- runif(1, 200, 500)
      } else if (st == "+") {
        dat[i, m] <- rnorm(1, mean = 750, sd = 120)
      } else {
        dat[i, m] <- rnorm(1, mean = 180, sd = 60)
      }
    }
  }
  dat[dat < 0] <- 0
  colnames(dat) <- intensity_cols

  # two well-separated spatial clusters per sample (gap > min_distance)
  cluster <- sample(c("A", "B"), n, replace = TRUE)
  x <- ifelse(cluster == "A", runif(n, 0, 300), runif(n, 900, 1200))
  y <- ifelse(cluster == "A", runif(n, 0, 300), runif(n, 900, 1200))

  # most cells inside the T-region, a few outside (for region filtering demo)
  parent <- sample(c("SLN_Tregion", "outline"),
                   n, replace = TRUE, prob = c(0.92, 0.08))

  d <- as.data.frame(dat, stringsAsFactors = FALSE)
  d[[".sample"]] <- sample_id
  d[[".type"]] <- types
  d[[".x"]] <- x
  d[[".y"]] <- y
  d[[".parent"]] <- parent
  d
}

cells <- lapply(c("C1001", "C1002", "C1003"), simulate_sample) %>%
  bind_rows()

mif_cells <- cells %>%
  mutate(
    Image = paste0("S1_", .sample, "_Scan1.er.qptiff - resolution #1"),
    Classification = NA_character_,
    `Nucleus: Area` = round(runif(n(), 20, 90), 2),
    `Nucleus/Cell area ratio` = round(runif(n(), 0.2, 0.9), 4),
    Parent = .parent,
    `Centroid X um` = .x,
    `Centroid Y um` = .y,
    sample_id = .sample,
    cell_id = paste0(.sample, "_", seq_len(n()))
  ) %>%
  select(sample_id, cell_id, Image, Classification, Parent,
         `Centroid X um`, `Centroid Y um`, `Nucleus: Area`,
         all_of(intensity_cols), `Nucleus/Cell area ratio`)

## ---- 5. Save package datasets ----
dir.create("data", showWarnings = FALSE)
save(mif_cells, file = "data/mif_cells.rda", compress = "xz")
save(mif_cutoffs, file = "data/mif_cutoffs.rda", compress = "xz")
save(mif_rules, file = "data/mif_rules.rda", compress = "xz")
save(mif_groups, file = "data/mif_groups.rda", compress = "xz")

## ---- 6. Example config files for inst/extdata ----
dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
write_tsv(mif_cutoffs, "inst/extdata/cutoffs_example.txt")
write_tsv(mif_rules, "inst/extdata/rules_example.txt")
write_tsv(mif_groups, "inst/extdata/groups_example.txt")

message("Example data written:")
message("  data/mif_cells.rda   : ", nrow(mif_cells), " cells x ", ncol(mif_cells), " cols")
message("  data/mif_cutoffs.rda : ", nrow(mif_cutoffs), " samples")
message("  data/mif_rules.rda   : ", nrow(mif_rules), " rules")
message("  data/mif_groups.rda  : ", nrow(mif_groups), " samples")
message("  inst/extdata/*_example.txt")
