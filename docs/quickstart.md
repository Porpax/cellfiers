# Quick start

[Home](../README.md) · [Installation](installation.md) · [Function reference](reference.md)

Use the bundled simulated data to select a tissue region, phenotype cells, apply FOV-level quality control, analyse cellular neighbourhoods, and draw a heat-pie map.

## 1. Select a region and assign cell phenotypes

```r
library(cellfiers)
data("mif_cells", package = "cellfiers")
data("mif_cutoffs", package = "cellfiers")
data("mif_rules", package = "cellfiers")
data("mif_groups", package = "cellfiers")

marker_map <- c(
  CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
  CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
  PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
)
cells <- run_annotation(
  mif_cells, mif_cutoffs, mif_rules, marker_map,
  regions = "SLN_Tregion", region_match = "exact"
)
summarize_celltypes(cells)
plot_proportions(cells)
```

## 2. Assign cells to FOVs and apply quality control

```r
centres <- c("cDC", "cDC(Mature/Active)")
qc <- run_fov_qc(
  cells, grid_size = 200, min_distance = 500,
  dc_types = centres, max_dc_ratio = 0.70, probs = c(0, 1)
)
plot_qc_summary(qc$qc)
clean <- qc$data_final
```

With `probs = c(0, 1)`, the quantile filter retains all simulated FOVs that pass the initial QC filters. Choose quantiles appropriate to your study. Inspect `qc$qc` for the FOV-level QC results and use `qc$data_final` to access cells in retained FOVs.

## 3. Compare cellular neighbourhood composition between groups

```r
neighbors <- analyze_neighbors(clean, radius = 20, center_types = centres)
grouped <- add_groups(neighbors$fov_summary, mif_groups, sample_col = "Sample")
tests <- compare_neighbors(grouped, analysis_level = "sample")
plot_neighbors(grouped, tests, geom = "boxplot")
```

`analyze_neighbors()` searches for neighbours within each FOV. `compare_neighbors()` supports sample-level inference, exploratory FOV-level Wilcoxon testing, and permutation of group labels at the biological-sample level. The group comparison illustrates the interface using simulated data.

## 4. Quantify target-cell neighbourhoods and draw a heat-pie map

```r
one_sample <- clean[clean$sample_id == "C1001", ]
targets <- c("CD8Tn(Primed)", "CD8Teff(Active)",
             "CD4Tn(Primed)", "CD4Teff(Active)")
result <- target_neighbors(
  one_sample, center_types = centres, target_types = targets,
  radius = 20, grid_size = 200, shared_weight = "inverse_degree"
)
head(result$fov_level)
p <- plot_heatpie(result)
p
```

`target_neighbors()` processes one biological sample at a time. Its search can cross FOV boundaries within the sample, and results are aggregated by the FOV containing each centre cell. The search radius and grid side length are in micrometres. Tile shading represents mean target cell-type richness per centre cell. Pie sectors show weighted target cell-type composition, and pie area represents mean weighted target-cell abundance per centre cell. When a target cell is shared by `m` centre cells, inverse-degree weighting assigns `1/m` to each centre–target pair.

## Use your own data

Retain the region label (by default, `Parent`), cell centroid coordinates, and fluorescence intensities in your QuPath export. Verify that the coordinates are in micrometres and map each marker in the threshold table to the corresponding measurement column. Specify one cell type per row of the rule table, with marker requirements recorded as `+`, `-`, or `NA`. The first matching rule determines the phenotype.

```r
# Set your file paths, sample ID, region name, and marker-to-channel mapping.
cells <- read_qupath("path/to/qupath_export.txt", sample_id = "S01")
cutoffs <- read_cutoffs("path/to/cutoffs.tsv")
rules <- read_rules("path/to/rules.tsv")
annotated <- run_annotation(
  cells, cutoffs, rules, marker_map, regions = "SLN_Tregion"
)
```

Use `read_qupath_batch()` for multiple TSV exports. For a CSV export, import it with `readr::read_csv()` and supply `sample_id` and a `cell_id` unique within each sample before calling `run_annotation()`. Review phenotyping and QC by sample before choosing centre and target cell types and the search radius. Use `ggplot2::ggsave()` to save figures as PNG or PDF.
