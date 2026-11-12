# Function reference

[Home](../README.md) | [Installation](installation.md) | [Quick start](quickstart.md)

Load the package with `library(cellfiers)`. The same reference is available in R through `help("function_name", package = "cellfiers")`.

- **Input and annotation** - [`read_qupath()`](#read-qupath), [`read_qupath_batch()`](#read-qupath-batch), [`read_cutoffs()`](#read-cutoffs), [`read_rules()`](#read-rules), [`read_groups()`](#read-groups), [`filter_region()`](#filter-region), [`apply_cutoffs()`](#apply-cutoffs), [`annotate_cells()`](#annotate-cells), [`run_annotation()`](#run-annotation), [`annotate_files()`](#annotate-files), [`summarize_celltypes()`](#summarize-celltypes), [`add_groups()`](#add-groups)
- **FOV and quality control** - [`assign_fov()`](#assign-fov), [`qc_fov_summary()`](#qc-fov-summary), [`filter_qc_fovs()`](#filter-qc-fovs), [`remove_extreme_fovs()`](#remove-extreme-fovs), [`remove_outliers_iqr()`](#remove-outliers-iqr), [`run_fov_qc()`](#run-fov-qc)
- **Neighbourhoods and statistics** - [`find_neighbors()`](#find-neighbors), [`neighbor_proportions()`](#neighbor-proportions), [`neighbor_fov_summary()`](#neighbor-fov-summary), [`neighbor_sample_summary()`](#neighbor-sample-summary), [`nearest_neighbor_rank()`](#nearest-neighbor-rank), [`analyze_neighbors()`](#analyze-neighbors), [`target_neighbors()`](#target-neighbors), [`compare_groups()`](#compare-groups), [`compare_neighbors()`](#compare-neighbors)
- **Plotting** - [`plot_cell_map()`](#plot-cell-map), [`plot_cells_by_sample()`](#plot-cells-by-sample), [`plot_proportions()`](#plot-proportions), [`plot_markers()`](#plot-markers), [`plot_fov()`](#plot-fov), [`plot_fovs_by_sample()`](#plot-fovs-by-sample), [`plot_qc_summary()`](#plot-qc-summary), [`plot_metric_density()`](#plot-metric-density), [`plot_dc_counts()`](#plot-dc-counts), [`plot_groups()`](#plot-groups), [`plot_neighbors()`](#plot-neighbors), [`plot_celltype_radar()`](#plot-celltype-radar), [`plot_ring_radar()`](#plot-ring-radar), [`plot_nearest_rank()`](#plot-nearest-rank), [`plot_heatpie()`](#plot-heatpie), [`theme_mif()`](#theme-mif)
- **Compatibility** - [`analyze_dc_tcell_tribe()`](#analyze-dc-tcell-tribe)

## Previous function names

The names below remain available with identical parameters and outputs.

| Previous name | Preferred name |
| --- | --- |
| `analyze_target_neighborhoods()` | [`target_neighbors()`](#target-neighbors) |
| `compare_neighbor_groups()` | [`compare_neighbors()`](#compare-neighbors) |
| `plot_cell_maps_by_sample()` | [`plot_cells_by_sample()`](#plot-cells-by-sample) |
| `plot_celltype_proportions()` | [`plot_proportions()`](#plot-proportions) |
| `plot_marker_distributions()` | [`plot_markers()`](#plot-markers) |
| `plot_fov_distribution()` | [`plot_fov()`](#plot-fov) |
| `plot_fov_distributions_by_sample()` | [`plot_fovs_by_sample()`](#plot-fovs-by-sample) |
| `plot_dc_count_distribution()` | [`plot_dc_counts()`](#plot-dc-counts) |
| `plot_group_comparison()` | [`plot_groups()`](#plot-groups) |
| `plot_neighbor_comparison()` | [`plot_neighbors()`](#plot-neighbors) |
| `plot_neighbor_ring_radar()` | [`plot_ring_radar()`](#plot-ring-radar) |
| `plot_neighborhood_heatpie()` | [`plot_heatpie()`](#plot-heatpie) |

## Input and annotation

<a id="read-qupath"></a>
### `read_qupath()`

Read a QuPath cell classification export

Reads a single tab-separated classification export produced by QuPath
(e.g. Shanxi_SLN_C2309160.txt), adds a sample_id column and a unique
cell_id column (sample_id + "_" + row number). Column names are preserved
as they appear in the QuPath file (for example "Centroid X µm" or
"Cell: Opal 780 mean").

**Usage**

```r
read_qupath(
  file,
  sample_id = NULL,
  sample_pattern = ".*_(C[0-9]+)[.]txt$",
  encoding = "UTF-8",
  ...
)
```

**Parameters**

- **`file`** - Path to a QuPath export file (tab-separated text).
- **`sample_id`** - Optional sample id. If NULL, it is extracted from the
file name using sample_pattern.
- **`sample_pattern`** - Regular expression with exactly one capture group
used to extract the sample id from the file name. Defaults to
".*_(C[0-9]+)[.]txt$" which matches names such as
Shanxi_SLN_C2309160.txt and yields C2309160.
- **`encoding`** - File encoding passed to readr::read_tsv(). Defaults to
"UTF-8".
- **`...`** - Further arguments passed to readr::read_tsv().

**Returns**

A tibble with the original QuPath columns plus sample_id and
  cell_id.

**Example**

```r
# Requires local input files
cells <- read_qupath("classfy_data/Shanxi_SLN_C2309160.txt")
```

<a id="read-qupath-batch"></a>
### `read_qupath_batch()`

Read all QuPath classification exports in a directory

Batch version of read_qupath(). Reads every file matching pattern in dir,
extracts a sample id from each file name and returns either a list of data
frames (one per sample) or a single combined data frame.

**Usage**

```r
read_qupath_batch(
  dir,
  pattern = "[.]txt$",
  sample_pattern = ".*_(C[0-9]+)[.]txt$",
  combine = TRUE,
  encoding = "UTF-8",
  skip_invalid = TRUE,
  ...
)
```

**Parameters**

- **`dir`** - Path to a directory containing QuPath export files.
- **`pattern`** - File-name pattern passed to base::list.files().
- **`sample_pattern`** - Regular expression with one capture group used to
extract sample ids (see read_qupath()).
- **`combine`** - Logical; if TRUE (default) all samples are combined into a
single data frame with a sample_id column. If FALSE, a named list is
returned instead.
- **`encoding`** - File encoding passed to readr::read_tsv().
- **`skip_invalid`** - Logical; if `TRUE` (default), files that cannot be
read or whose columns have invalid (oversized/duplicated) names are
skipped with a warning instead of aborting the batch.
- **`...`** - Further arguments passed to read_qupath().

**Returns**

A data frame (if combine = TRUE) or a named list of data frames.

**Example**

```r
# Requires local input files
cells <- read_qupath_batch("classfy_data/")
```

<a id="read-cutoffs"></a>
### `read_cutoffs()`

Read configuration tables used by the annotation pipeline

Convenience wrappers around [readr::read_tsv()] for the tab-delimited
configuration files used by the mIF workflow:

**Usage**

```r
read_cutoffs(path, ...)
```

**Parameters**

- **`path`** - Path to a tab-delimited file.
- **`...`** - Further arguments passed to [readr::read_tsv()].

**Returns**

A [tibble::tibble()].

**Details**

- `read_cutoffs()`: per-sample fluorescence channel cutoffs (first column
  sample id, remaining columns marker names).
- `read_rules()`: cell-type classification rules (first column cell type,
  remaining columns marker statuses `"+"`/`"-"`/`NA`).
- `read_groups()`: sample-to-group mapping (columns `Sample`, `Group`).

**Example**

```r
# Requires local input files
cutoffs <- read_cutoffs("Docu/SLN_mIF_6channel_cutoff.txt")
rules <- read_rules("Docu/SLN_mIF_anno_rule.txt")
groups <- read_groups("Docu/SLN_mIF_Group_pCRpNR.txt")
```

<a id="read-rules"></a>
### `read_rules()`

Read configuration tables used by the annotation pipeline

Convenience wrappers around [readr::read_tsv()] for the tab-delimited
configuration files used by the mIF workflow:

**Usage**

```r
read_rules(path, ...)
```

**Parameters**

- **`path`** - Path to a tab-delimited file.
- **`...`** - Further arguments passed to [readr::read_tsv()].

**Returns**

A [tibble::tibble()].

**Details**

- `read_cutoffs()`: per-sample fluorescence channel cutoffs (first column
  sample id, remaining columns marker names).
- `read_rules()`: cell-type classification rules (first column cell type,
  remaining columns marker statuses `"+"`/`"-"`/`NA`).
- `read_groups()`: sample-to-group mapping (columns `Sample`, `Group`).

**Example**

```r
# Requires local input files
cutoffs <- read_cutoffs("Docu/SLN_mIF_6channel_cutoff.txt")
rules <- read_rules("Docu/SLN_mIF_anno_rule.txt")
groups <- read_groups("Docu/SLN_mIF_Group_pCRpNR.txt")
```

<a id="read-groups"></a>
### `read_groups()`

Read configuration tables used by the annotation pipeline

Convenience wrappers around [readr::read_tsv()] for the tab-delimited
configuration files used by the mIF workflow:

**Usage**

```r
read_groups(path, ...)
```

**Parameters**

- **`path`** - Path to a tab-delimited file.
- **`...`** - Further arguments passed to [readr::read_tsv()].

**Returns**

A [tibble::tibble()].

**Details**

- `read_cutoffs()`: per-sample fluorescence channel cutoffs (first column
  sample id, remaining columns marker names).
- `read_rules()`: cell-type classification rules (first column cell type,
  remaining columns marker statuses `"+"`/`"-"`/`NA`).
- `read_groups()`: sample-to-group mapping (columns `Sample`, `Group`).

**Example**

```r
# Requires local input files
cutoffs <- read_cutoffs("Docu/SLN_mIF_6channel_cutoff.txt")
rules <- read_rules("Docu/SLN_mIF_anno_rule.txt")
groups <- read_groups("Docu/SLN_mIF_Group_pCRpNR.txt")
```

<a id="filter-region"></a>
### `filter_region()`

Filter cells by parent region

Keeps (or drops) cells whose parent annotation matches user-selected
regions. Region names are not restricted to `Tregion` or `SLN_Tregion`;
any QuPath annotation name can be selected by exact, fixed-substring or
regular-expression matching.

**Usage**

```r
filter_region(
  data,
  parent_col = NULL,
  regions,
  keep = TRUE,
  match = c("exact", "fixed", "regex"),
  ignore_case = FALSE
)
```

**Parameters**

- **`data`** - A data frame with QuPath cell data.
- **`parent_col`** - Name of the parent/annotation column. If NULL, the
function looks for a column named "Parent" (or a normalised match).
- **`regions`** - Character vector of parent values to keep (or drop).
- **`keep`** - Logical; if TRUE keep the matching rows, if FALSE drop them.
- **`match`** - Matching mode: `"exact"`, `"fixed"` substring, or `"regex"`.
- **`ignore_case`** - Logical; ignore case for `fixed` or `regex` matching.

**Returns**

A filtered data frame.

**Example**

```r
data(mif_cells)
tzone <- filter_region(mif_cells, regions = "SLN_Tregion")
tumour <- filter_region(mif_cells, regions = "Tumou?r|Tumor",
                        match = "regex", ignore_case = TRUE)
```

<a id="apply-cutoffs"></a>
### `apply_cutoffs()`

Apply fluorescence channel cutoffs to call markers positive/negative

For each marker in `marker_map`, cells with a channel intensity greater
than or equal to the sample-specific cutoff are labelled positive, others
negative. Cutoffs are stored in a data frame such as
`Docu/SLN_mIF_6channel_cutoff.txt` with one row per sample and one column
per marker.

**Usage**

```r
apply_cutoffs(
  data,
  cutoffs,
  marker_map,
  sample_col = "sample_id",
  cutoff_id_col = NULL,
  pos = "+",
  neg = "-"
)
```

**Parameters**

- **`data`** - A data frame with QuPath cell data (must contain a sample
column and the intensity columns named by `marker_map`).
- **`cutoffs`** - A data frame with sample ids and one numeric cutoff column
per marker. The first column (or `cutoff_id_col`) holds the sample ids.
- **`marker_map`** - A named character vector mapping marker names to
intensity column names, e.g.
`c(CD8 = "Cell: Opal 780 mean", CD3 = "Cell: Opal 480 mean")`.
- **`sample_col`** - Name of the sample column in `data`.
- **`cutoff_id_col`** - Name of the sample-id column in `cutoffs`. Defaults
to the first column of `cutoffs`.
- **`pos`** - Character used for positive calls (default `"+"`).
- **`neg`** - Character used for negative calls (default `"-"`).

**Returns**

The input data frame with one extra column per marker holding
  `pos`/`neg` (or `NA` when no cutoff is available for the sample, or when
  the intensity is missing).

**Example**

```r
data(mif_cells)
data(mif_cutoffs)
marker_map <- c(
  CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
  CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
  PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
)
cells_pos <- apply_cutoffs(mif_cells, mif_cutoffs, marker_map)
```

<a id="annotate-cells"></a>
### `annotate_cells()`

Assign cell types from a classification rule table

Each cell is matched against the rows of a rule table (such as
`Docu/SLN_mIF_anno_rule.txt`). A rule row lists, for each marker, a
required status (`"+"`, `"-"` or `NA` = "don't care"). The first rule
whose non-NA markers all agree with the cell's marker status wins. Cells
matching no rule are labelled `unmatched` (default `"Other"`).

**Usage**

```r
annotate_cells(
  data,
  rules,
  markers = NULL,
  celltype_col = "CellType",
  rule_col = NULL,
  unmatched = "Other"
)
```

**Parameters**

- **`data`** - A data frame with marker status columns (e.g. produced by
[apply_cutoffs()]).
- **`rules`** - A data frame whose first column (or `rule_col`) holds the
cell type labels and whose remaining columns hold marker statuses.
- **`markers`** - Optional character vector of marker columns to use. Defaults
to the columns of `rules` other than `rule_col` that also exist in
`data`.
- **`celltype_col`** - Name of the output column holding the assigned cell
type. Defaults to `"CellType"`.
- **`rule_col`** - Name of the cell-type column in `rules`. Defaults to the
first column of `rules`.
- **`unmatched`** - Label for cells that match no rule (default `"Other"`).

**Returns**

The input data frame with an extra `celltype_col` column.

**Details**

Marker status columns in `data` may hold `"+"`/`"-"`, `1`/`0` or
`TRUE`/`FALSE`; they are normalised internally.

**Example**

```r
data(mif_cells)
data(mif_cutoffs)
data(mif_rules)
marker_map <- c(
  CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
  CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
  PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
)
pos <- apply_cutoffs(mif_cells, mif_cutoffs, marker_map)
annotated <- annotate_cells(pos, mif_rules)
```

<a id="run-annotation"></a>
### `run_annotation()`

Run the full Phase-1 annotation pipeline

Convenience wrapper that chains [filter_region()] (optional),
[apply_cutoffs()] and [annotate_cells()] on a single data frame.

**Usage**

```r
run_annotation(
  data,
  cutoffs,
  rules,
  marker_map,
  regions = NULL,
  parent_col = NULL,
  region_match = c("exact", "fixed", "regex"),
  region_ignore_case = FALSE,
  sample_col = "sample_id",
  celltype_col = "CellType",
  unmatched = "Other",
  ...
)
```

**Parameters**

- **`data`** - A data frame with QuPath cell data (must contain a sample
column and the intensity columns named by `marker_map`).
- **`cutoffs`** - A cutoff data frame (see [apply_cutoffs()]).
- **`rules`** - A rule data frame (see [annotate_cells()]).
- **`marker_map`** - Named character vector mapping markers to intensity
columns (see [apply_cutoffs()]).
- **`regions`** - Optional character vector of parent regions to keep. If
`NULL`, no region filtering is performed.
- **`parent_col`** - Name of the parent column (see [filter_region()]).
- **`region_match`** - Region matching mode passed to [filter_region()].
- **`region_ignore_case`** - Logical; ignore case when matching region names.
- **`sample_col`** - Name of the sample column in `data`.
- **`celltype_col`** - Name of the output cell-type column.
- **`unmatched`** - Label for cells matching no rule (see [annotate_cells()]).
- **`...`** - Further arguments passed to [annotate_cells()].

**Returns**

The annotated data frame.

**Example**

```r
data(mif_cells)
data(mif_cutoffs)
data(mif_rules)
marker_map <- c(
  CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
  CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
  PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
)
annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
```

<a id="annotate-files"></a>
### `annotate_files()`

Annotate a directory of QuPath exports and save results

Batch pipeline over a directory of QuPath classification exports: reads
each file, optionally filters to a region, applies cutoffs and rules,
writes one annotated CSV per sample plus a combined cell-type summary to
`output_dir`.

**Usage**

```r
annotate_files(
  input_dir,
  output_dir,
  cutoffs,
  rules,
  marker_map,
  regions = NULL,
  pattern = "[.]txt$",
  sample_pattern = ".*_(C[0-9]+)[.]txt$",
  encoding = "UTF-8",
  celltype_col = "CellType",
  write_summary = TRUE,
  ...
)
```

**Parameters**

- **`input_dir`** - Directory containing QuPath export files.
- **`output_dir`** - Directory where annotated CSVs and the summary are
written (created if needed).
- **`cutoffs`** - A cutoff data frame (see [apply_cutoffs()]).
- **`rules`** - A rule data frame (see [annotate_cells()]).
- **`marker_map`** - Named character vector mapping markers to intensity
columns (see [apply_cutoffs()]).
- **`regions`** - Optional character vector of parent regions to keep.
- **`pattern`** - File-name pattern passed to [read_qupath_batch()].
- **`sample_pattern`** - Regex with one capture group to extract sample ids
(see [read_qupath()]).
- **`encoding`** - File encoding passed to [readr::read_tsv()].
- **`celltype_col`** - Name of the output cell-type column.
- **`write_summary`** - Logical; if `TRUE` (default) also write
`celltype_summary.csv` with per-sample cell-type counts.
- **`...`** - Further arguments passed to [annotate_cells()].

**Returns**

Invisibly, a named list of annotated data frames (one per sample).

**Example**

```r
# Requires local input files
annotate_files(
  input_dir = "classfy_data",
  output_dir = "annotated_results",
  cutoffs = cutoffs_df,
  rules = rules_df,
  marker_map = marker_map,
  regions = c("Annotation (SLN_Tregion)", "SLN_Tregion")
)
```

<a id="summarize-celltypes"></a>
### `summarize_celltypes()`

Summarise cell-type counts per sample

Builds a per-sample count table with one column per cell type and a
`Total` column. Cell-type columns are ordered by `levels` (or, by default,
by the order of appearance in `data`).

**Usage**

```r
summarize_celltypes(
  data,
  celltype_col = "CellType",
  sample_col = "sample_id",
  levels = NULL,
  include_total = TRUE
)
```

**Parameters**

- **`data`** - A data frame with a cell-type column and a sample column
(e.g. the output of [run_annotation()]).
- **`celltype_col`** - Name of the cell-type column.
- **`sample_col`** - Name of the sample column.
- **`levels`** - Optional character vector giving the order of cell types in
the output. If `NULL`, the order of appearance is used.
- **`include_total`** - Logical; if `TRUE` (default) add a `Total` column.

**Returns**

A data frame with one row per sample.

**Example**

```r
data(mif_cells)
data(mif_cutoffs)
data(mif_rules)
marker_map <- c(
  CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
  CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
  PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
)
annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
summarize_celltypes(annotated)
```

<a id="add-groups"></a>
### `add_groups()`

Add experimental group information to a data frame

Left-joins a sample-to-group mapping (such as
`Docu/SLN_mIF_Group_pCRpNR.txt`, columns `Sample` and `Group`) onto cell or
FOV level data.

**Usage**

```r
add_groups(
  data,
  groups,
  sample_col = "sample_id",
  sample_id_col = "Sample",
  group_col = "Group"
)
```

**Parameters**

- **`data`** - A data frame with a sample column.
- **`groups`** - A data frame with sample ids and group labels.
- **`sample_col`** - Name of the sample column in `data`.
- **`sample_id_col`** - Name of the sample-id column in `groups` (defaults to
`"Sample"`).
- **`group_col`** - Name of the output group column (defaults to `"Group"`).

**Returns**

The input data frame with an extra `group_col` column.

**Example**

```r
data(mif_groups)
data(mif_cells)
with_group <- add_groups(mif_cells, mif_groups)
head(with_group[, c("sample_id", "Group")])
```

## FOV and quality control

<a id="assign-fov"></a>
### `assign_fov()`

Assign fields of view (FOVs) by rasterising tissue coordinates

Reproduces the FOV rasterisation step of the original workflow
(`mIF_Tzone_02_GridArea_Process.R`). Cells are placed on a square grid of
`grid_size` micrometres. Grid cells are grouped into spatially separated
tissue regions by hierarchical clustering of the grid centres (single
linkage, cut at `min_distance`). Within each tissue region, grids are
re-numbered from the minimum x/y grid origin (row-major) and each cell is labelled
`"<sample>_TR<tregion>_FOV<fov>"`.

**Usage**

```r
assign_fov(
  data,
  sample_col = "sample_id",
  x = NULL,
  y = NULL,
  grid_size = 200,
  min_distance = 500,
  fov_col = "FOV"
)
```

**Parameters**

- **`data`** - A data frame with sample, x and y columns (e.g. annotated cell
data from [run_annotation()]).
- **`sample_col`** - Name of the sample column.
- **`x`, `y`** - Names of the centroid columns. Defaults to QuPath-style
`"Centroid X µm"` / `"Centroid Y µm"` (or normalised matches).
- **`grid_size`** - Grid side length in micrometres.
- **`min_distance`** - Minimum distance (µm) between two grids for them to be
considered separate tissue regions.
- **`fov_col`** - Name of the output FOV column.

**Returns**

The input data frame with an extra `fov_col` column. Row order is
  preserved.

**Example**

```r
data(mif_cells)
data(mif_cutoffs)
data(mif_rules)
marker_map <- c(
  CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
  CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
  PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
)
annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
annotated_fov <- assign_fov(annotated)
head(annotated_fov$FOV)
```

<a id="qc-fov-summary"></a>
### `qc_fov_summary()`

Per-FOV quality control summary

Computes per-FOV QC metrics from annotated cell data (mirroring
`mIF_Tzone_05_QC_Visual.R`): total cell number, DC and T cell counts and
their ratios, the mean cell count per sample and a low-cell-number flag
(a FOV is flagged `"LowCellNum"` when it contains fewer than
`mean_cells / low_cell_factor` cells).

**Usage**

```r
qc_fov_summary(
  data,
  fov_col = "FOV",
  celltype_col = "CellType",
  sample_col = "sample_id",
  dc_types = NULL,
  tcell_types = NULL,
  low_cell_factor = 4,
  max_dc_ratio = 0.7
)
```

**Parameters**

- **`data`** - A data frame with FOV, cell-type and sample columns (output of
[assign_fov()] applied to annotated data).
- **`fov_col`** - Name of the FOV column.
- **`celltype_col`** - Name of the cell-type column.
- **`sample_col`** - Name of the sample column.
- **`dc_types`** - Character vector of cell types counted as DC. If `NULL`,
types whose name contains `"DC"` are used (e.g. `"cDC"`,
`"cDC(Mature/Active)"`).
- **`tcell_types`** - Character vector of cell types counted as T cells. If
`NULL`, types starting with `"CD4"` or `"CD8"` are used.
- **`low_cell_factor`** - Denominator of the low-cell-number threshold
(default 4, i.e. a FOV with fewer than one quarter of the sample mean is
flagged).
- **`max_dc_ratio`** - Maximum accepted fraction of DC cells in a FOV. FOVs
above this value are flagged as `"HighDCRatio"`. The default is `0.70`;
use `NULL` to disable this rule.

**Returns**

A data frame with one row per FOV and columns `Sample`, FOV,
  `TotalCells`, `DC_count`, `Tcell_count`, `DC_ratio`, `Tcell_ratio`,
  `MeanCells`, `LowCellFlag`, `HighDCRatioFlag`, and `InitialQC`.

**Example**

```r
data(mif_cells)
data(mif_cutoffs)
data(mif_rules)
marker_map <- c(
  CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
  CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
  PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
)
annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
fov <- assign_fov(annotated)
qc <- qc_fov_summary(fov)
head(qc)
```

<a id="filter-qc-fovs"></a>
### `filter_qc_fovs()`

Keep only FOVs that pass quality control

Filters cell-level (or summary-level) data down to the FOVs whose
`LowCellFlag` is one of `keep` (default `"OK"`).

**Usage**

```r
filter_qc_fovs(data, qc, keep = "OK", fov_col = "FOV")
```

**Parameters**

- **`data`** - Cell-level data (or a QC summary) with an FOV column.
- **`qc`** - The QC summary returned by [qc_fov_summary()].
- **`keep`** - Character vector of `LowCellFlag` values to keep.
- **`fov_col`** - Name of the FOV column in `data` and `qc`.

**Returns**

The filtered data frame.

**Example**

```r
data(mif_cells)
data(mif_cutoffs)
data(mif_rules)
marker_map <- c(
  CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
  CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
  PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
)
annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
fov <- assign_fov(annotated)
qc <- qc_fov_summary(fov)
clean <- filter_qc_fovs(fov, qc)
```

<a id="remove-extreme-fovs"></a>
### `remove_extreme_fovs()`

Remove FOVs with extreme values of a QC metric

Drops FOVs whose value of `metric` (e.g. `DC_count`) lies outside the
central `probs` quantile interval, mirroring the "remove extreme DC
counts" step of `mIF_Tzone_05_QC_Visual.R`. When cell-level `data` is
supplied, both the filtered QC summary and the filtered cell data are
returned.

**Usage**

```r
remove_extreme_fovs(
  qc,
  metric = "DC_count",
  probs = c(0.05, 0.95),
  data = NULL,
  fov_col = "FOV"
)
```

**Parameters**

- **`qc`** - The QC summary returned by [qc_fov_summary()].
- **`metric`** - Name of the numeric QC metric column to filter on.
- **`probs`** - Quantile probabilities defining the kept interval
(default `c(0.05, 0.95)`).
- **`data`** - Optional cell-level data with an FOV column; if given, it is
filtered to the retained FOVs.
- **`fov_col`** - Name of the FOV column.

**Returns**

If `data` is `NULL`, the filtered QC summary. Otherwise a list with
  elements `qc` (filtered summary) and `data` (filtered cell data).

**Example**

```r
data(mif_cells)
data(mif_cutoffs)
data(mif_rules)
marker_map <- c(
  CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
  CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
  PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
)
annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
fov <- assign_fov(annotated)
qc <- qc_fov_summary(fov)
qc_f <- remove_extreme_fovs(qc, metric = "DC_count")
```

<a id="remove-outliers-iqr"></a>
### `remove_outliers_iqr()`

Remove outliers with the IQR rule

Drops rows whose `value_col` lies outside
`[Q1 - k * IQR, Q3 + k * IQR]`, computed either globally or per group
(e.g. per cell type). This is the outlier cleaning step used before box
plots in `mIF_Tzone_05_QC_Visual.R`.

**Usage**

```r
remove_outliers_iqr(data, value_col, group_col = NULL, k = 1.5)
```

**Parameters**

- **`data`** - A data frame.
- **`value_col`** - Name of the numeric column to clean.
- **`group_col`** - Optional name of a grouping column; bounds are computed
within each group.
- **`k`** - IQR multiplier (default 1.5).

**Returns**

The data frame with outlier rows removed.

**Example**

```r
df <- data.frame(
  group = rep(c("A", "B"), each = 10),
  value = c(rnorm(10, 5, 1), rnorm(10, 50, 10))
)
clean <- remove_outliers_iqr(df, value_col = "value", group_col = "group")
nrow(clean)
```

<a id="run-fov-qc"></a>
### `run_fov_qc()`

Rasterise cells and run field-of-view quality control

Runs the Phase 2 workflow in its intended order: assign FOVs, flag FOVs
below the per-sample mean cell count divided by `low_cell_factor`, flag
FOVs above `max_dc_ratio`, and then remove FOVs outside the requested
quantiles of a QC metric. Quantiles are estimated only from FOVs that
passed the low-cell and DC-ratio rules.

**Usage**

```r
run_fov_qc(
  data,
  sample_col = "sample_id",
  x = NULL,
  y = NULL,
  grid_size = 200,
  min_distance = 500,
  fov_col = "FOV",
  celltype_col = "CellType",
  dc_types = NULL,
  tcell_types = NULL,
  low_cell_factor = 4,
  max_dc_ratio = 0.7,
  metric = "DC_count",
  probs = c(0.05, 0.95)
)
```

**Parameters**

- **`data`** - A data frame with sample, x and y columns (e.g. annotated cell
data from [run_annotation()]).
- **`sample_col`** - Name of the sample column.
- **`x`, `y`** - Names of the centroid columns. Defaults to QuPath-style
`"Centroid X µm"` / `"Centroid Y µm"` (or normalised matches).
- **`grid_size`** - Grid side length in micrometres.
- **`min_distance`** - Minimum distance (µm) between two grids for them to be
considered separate tissue regions.
- **`fov_col`** - Name of the output FOV column.
- **`celltype_col`** - Name of the cell-type column.
- **`dc_types`** - Character vector of cell types counted as DC. If `NULL`,
types whose name contains `"DC"` are used (e.g. `"cDC"`,
`"cDC(Mature/Active)"`).
- **`tcell_types`** - Character vector of cell types counted as T cells. If
`NULL`, types starting with `"CD4"` or `"CD8"` are used.
- **`low_cell_factor`** - Denominator of the low-cell-number threshold
(default 4, i.e. a FOV with fewer than one quarter of the sample mean is
flagged).
- **`max_dc_ratio`** - Maximum accepted fraction of DC cells in a FOV. FOVs
above this value are flagged as `"HighDCRatio"`. The default is `0.70`;
use `NULL` to disable this rule.
- **`metric`** - Numeric QC metric used for the extreme-FOV filter.
- **`probs`** - Two probabilities defining the retained quantile interval.

**Returns**

A named list containing `data_fov`, `qc`, `qc_pass`, `qc_final`,
  `data_final`, and `metric_quantiles`.

**Example**

```r
data(mif_cells)
data(mif_cutoffs)
data(mif_rules)
marker_map <- c(
  CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
  CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
  PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
)
annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
result <- run_fov_qc(annotated, probs = c(0, 1))
names(result)
```

## Neighbourhoods and statistics

<a id="find-neighbors"></a>
### `find_neighbors()`

Find neighbouring cells within a radius (per FOV)

For every cell, searches for cells of the same FOV within `radius`
micrometres using a KD-tree ([RANN::nn2()]). This mirrors the neighbour
analysis of `mIF_Tzone_03_FOV_cell_Distance.R`.

**Usage**

```r
find_neighbors(
  data,
  radius = 20,
  celltype_col = "CellType",
  fov_col = "FOV",
  x = NULL,
  y = NULL,
  id_col = NULL,
  max_k = 1000,
  include_self = FALSE
)
```

**Parameters**

- **`data`** - A data frame with FOV, cell-type and centroid columns (output
of [assign_fov()] applied to annotated data).
- **`radius`** - Search radius in micrometres.
Radii above 30 micrometres are allowed but trigger a warning because
runtime and neighbour density increase and biological interpretation may
become less local.
- **`celltype_col`** - Name of the cell-type column.
- **`fov_col`** - Name of the FOV column.
- **`x`, `y`** - Names of the centroid columns. Defaults to QuPath-style
`"Centroid X µm"` / `"Centroid Y µm"`.
- **`id_col`** - Optional column with unique cell ids (e.g. `cell_id`). If
`NULL`, a synthetic id is used internally.
- **`max_k`** - Maximum number of neighbours requested from the KD-tree
(performance cap).
- **`include_self`** - Logical; if `TRUE`, the centre cell itself counts as a
neighbour (default `FALSE`).

**Returns**

The input data frame with neighbour-count columns, `n_neighbors`
  and `nearestNGB` added. Row order is preserved.

**Details**

For each centre cell the result contains:

  
- 
 one count column per cell type (number of neighbours of that type
    within the radius);
  
- 
 `n_neighbors`: total number of neighbours within the radius;
  
- 
 `nearestNGB`: cell type(s) of the closest neighbour(s) (ties are
    comma-joined), or `"None"` when the cell has no neighbours.

**Example**

```r
data(mif_cells)
data(mif_cutoffs)
data(mif_rules)
marker_map <- c(
  CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
  CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
  PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
)
annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
fov <- assign_fov(annotated)
nb <- find_neighbors(fov, radius = 60)
head(nb[, c("cell_id", "n_neighbors", "nearestNGB")])
```

<a id="neighbor-proportions"></a>
### `neighbor_proportions()`

Per-FOV neighbour statistics for centre-cell types

Two summary functions built on top of [find_neighbors()]:

**Usage**

```r
neighbor_proportions(data, center_types = NULL, celltype_col = "CellType")
```

**Parameters**

- **`data`** - Output of [find_neighbors()].
- **`center_types`** - Cell types to use as centre cells (e.g.
`c("cDC", "cDC(Mature/Active)")`). If `NULL`, all cell types are used.
- **`celltype_col`** - Name of the cell-type column.

**Returns**

- `neighbor_proportions()`: the centre-cell data frame with
  `total_neighbors` and `Pro_<type>` columns.
- `neighbor_fov_summary()`: a data frame with columns `Sample`, FOV,
  `CellType`, `mean_count`, `mean_prop`, `n_center` (long), or one row per
  FOV with `mean_count_<type>` / `mean_prop_<type>` columns (wide).

**Details**

- `neighbor_proportions()`: keeps centre cells of the requested types and,
  for each of them, computes the total number of neighbours and the
  proportion of neighbours of each cell type (mirroring the per-DC
  computation of `mIF_Tzone_04_statistic.R`).
- `neighbor_fov_summary()`: aggregates neighbour counts and proportions per
  FOV (mean count and mean proportion per neighbour cell type, plus the
  number of centre cells per FOV), in long (default) or wide format.

**Example**

```r
data(mif_cells)
data(mif_cutoffs)
data(mif_rules)
marker_map <- c(
  CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
  CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
  PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
)
annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
fov <- assign_fov(annotated)
nb <- find_neighbors(fov, radius = 60)
dc <- neighbor_proportions(nb, center_types = c("cDC", "cDC(Mature/Active)"))
summ <- neighbor_fov_summary(nb, center_types = c("cDC", "cDC(Mature/Active)"))
head(summ)
```

<a id="neighbor-fov-summary"></a>
### `neighbor_fov_summary()`

Per-FOV neighbour statistics for centre-cell types

Two summary functions built on top of [find_neighbors()]:

**Usage**

```r
neighbor_fov_summary(
  data,
  center_types = NULL,
  celltype_col = "CellType",
  fov_col = "FOV",
  sample_col = "sample_id",
  wide = FALSE
)
```

**Parameters**

- **`data`** - Output of [find_neighbors()].
- **`center_types`** - Cell types to use as centre cells (e.g.
`c("cDC", "cDC(Mature/Active)")`). If `NULL`, all cell types are used.
- **`celltype_col`** - Name of the cell-type column.
- **`fov_col`** - Name of the FOV column.
- **`sample_col`** - Name of the sample column.
- **`wide`** - For `neighbor_fov_summary()`, return wide format (one row per
FOV) instead of long format.

**Returns**

- `neighbor_proportions()`: the centre-cell data frame with
  `total_neighbors` and `Pro_<type>` columns.
- `neighbor_fov_summary()`: a data frame with columns `Sample`, FOV,
  `CellType`, `mean_count`, `mean_prop`, `n_center` (long), or one row per
  FOV with `mean_count_<type>` / `mean_prop_<type>` columns (wide).

**Details**

- `neighbor_proportions()`: keeps centre cells of the requested types and,
  for each of them, computes the total number of neighbours and the
  proportion of neighbours of each cell type (mirroring the per-DC
  computation of `mIF_Tzone_04_statistic.R`).
- `neighbor_fov_summary()`: aggregates neighbour counts and proportions per
  FOV (mean count and mean proportion per neighbour cell type, plus the
  number of centre cells per FOV), in long (default) or wide format.

**Example**

```r
data(mif_cells)
data(mif_cutoffs)
data(mif_rules)
marker_map <- c(
  CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
  CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
  PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
)
annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
fov <- assign_fov(annotated)
nb <- find_neighbors(fov, radius = 60)
dc <- neighbor_proportions(nb, center_types = c("cDC", "cDC(Mature/Active)"))
summ <- neighbor_fov_summary(nb, center_types = c("cDC", "cDC(Mature/Active)"))
head(summ)
```

<a id="neighbor-sample-summary"></a>
### `neighbor_sample_summary()`

Summarise FOV neighbour statistics at the biological-sample level

Averages FOV-level neighbour statistics within each sample and neighbour
cell type. This creates one observation per biological sample and cell type,
which is generally more appropriate than treating FOVs as independent
replicates in between-group tests.

**Usage**

```r
neighbor_sample_summary(data, sample_col = "Sample", celltype_col = "CellType")
```

**Parameters**

- **`data`** - Long-format output of [neighbor_fov_summary()].
- **`sample_col`** - Name of the sample column.
- **`celltype_col`** - Name of the neighbour cell-type column.

**Returns**

A data frame with `Sample`, `CellType`, `mean_count`, `mean_prop`,
  `n_fov`, and `n_center` columns.

**Example**

```r
data(mif_cells)
data(mif_cutoffs)
data(mif_rules)
marker_map <- c(
  CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
  CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
  PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
)
annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
nb <- analyze_neighbors(assign_fov(annotated), radius = 60)
neighbor_sample_summary(nb$fov_summary)
```

<a id="nearest-neighbor-rank"></a>
### `nearest_neighbor_rank()`

Rank nearest-neighbour cell types per FOV

For each FOV, tabulates the `nearestNGB` values (the cell type of the
closest neighbour) and ranks them by count (mirroring the nearest-neighbour
ranking of `mIF_Tzone_04_statistic.R`).

**Usage**

```r
nearest_neighbor_rank(
  data,
  fov_col = "FOV",
  sample_col = "sample_id",
  nearest_col = "nearestNGB",
  center_types = NULL,
  celltype_col = "CellType",
  exclude = "None"
)
```

**Parameters**

- **`data`** - Output of [find_neighbors()].
- **`fov_col`** - Name of the FOV column.
- **`sample_col`** - Name of the sample column.
- **`nearest_col`** - Name of the nearest-neighbour column.
- **`center_types`** - Optional cell types to retain as centre cells. If
`NULL`, all cells are used.
- **`celltype_col`** - Name of the centre-cell type column.
- **`exclude`** - Character vector of `nearestNGB` values to drop (default
`"None"` for cells without neighbours).

**Returns**

A data frame with columns FOV, Sample, `nearest_col`, `Count`,
  `Proportion` and `Rank` (1 = most frequent), sorted by FOV and rank.

**Example**

```r
data(mif_cells)
data(mif_cutoffs)
data(mif_rules)
marker_map <- c(
  CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
  CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
  PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
)
annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
fov <- assign_fov(annotated)
nb <- find_neighbors(fov, radius = 60)
rank <- nearest_neighbor_rank(nb)
head(rank)
```

<a id="analyze-neighbors"></a>
### `analyze_neighbors()`

Run the neighbour-analysis pipeline

Convenience wrapper that runs [find_neighbors()], [neighbor_fov_summary()]
and [nearest_neighbor_rank()] in one call.

**Usage**

```r
analyze_neighbors(
  data,
  radius = 20,
  center_types = NULL,
  celltype_col = "CellType",
  fov_col = "FOV",
  sample_col = "sample_id",
  x = NULL,
  y = NULL,
  id_col = NULL,
  max_k = 1000,
  include_self = FALSE
)
```

**Parameters**

- **`data`** - A data frame with FOV, cell-type and centroid columns.
- **`radius`** - Search radius in micrometres.
- **`center_types`** - Cell types used as centre cells for the FOV summary.
- **`celltype_col`** - Name of the cell-type column.
- **`fov_col`** - Name of the FOV column.
- **`sample_col`** - Name of the sample column.
- **`x`, `y`** - Names of the centroid columns. Defaults to QuPath-style
`"Centroid X µm"` / `"Centroid Y µm"`.
- **`id_col`** - Optional column with unique cell ids (e.g. `cell_id`). If
`NULL`, a synthetic id is used internally.
- **`max_k`** - Maximum number of neighbours requested from the KD-tree
(performance cap).
- **`include_self`** - Logical; if `TRUE`, the centre cell itself counts as a
neighbour (default `FALSE`).

**Returns**

A list with elements:

  
- 
 `neighbors`: per-cell neighbour counts (output of
    [find_neighbors()]);
  
- 
 `fov_summary`: per-FOV neighbour statistics (output of
    [neighbor_fov_summary()]);
  
- 
 `nearest_rank`: per-FOV nearest-neighbour ranking (output of
    [nearest_neighbor_rank()]).

**Example**

```r
data(mif_cells)
data(mif_cutoffs)
data(mif_rules)
marker_map <- c(
  CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
  CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
  PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
)
annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
fov <- assign_fov(annotated)
res <- analyze_neighbors(fov, radius = 20,
                         center_types = c("cDC", "cDC(Mature/Active)"))
names(res)
```

<a id="target-neighbors"></a>
### `target_neighbors()`

Analyse user-defined centre--target neighbourhoods

Searches around any user-selected centre cell type and quantifies any set
of target cell types. A target shared by `m` centre cells is weighted by
`1/m` by default, preventing the same target cell from being counted in
full for every nearby centre.

**Usage**

```r
target_neighbors(
  data,
  center_types,
  target_types,
  radius = 20,
  other_types = NULL,
  shared_weight = c("inverse_degree", "none"),
  sample_col = "sample_id",
  celltype_col = "CellType",
  fov_col = "FOV",
  x = NULL,
  y = NULL,
  grid_size = 200,
  max_k = 1000L
)
```

**Parameters**

- **`data`** - One biological sample of annotated, FOV-assigned cell data.
- **`center_types`** - Character vector of cell types used as centres.
- **`target_types`** - Character vector of target cell types used to calculate
richness, raw abundance and weighted abundance.
- **`radius`** - Search radius in micrometres. Values above 30 are allowed but
trigger a warning about computation and biological interpretation.
- **`other_types`** - Optional non-target types used to label an
`Other-selected` neighbourhood when no requested target is present.
- **`shared_weight`** - Either `"inverse_degree"` (default, `1/m`) or
`"none"` (every centre--target edge has weight 1).
- **`sample_col`, `celltype_col`, `fov_col`** - Column names.
- **`x`, `y`** - Optional centroid column names.
- **`grid_size`** - FOV side length in micrometres, used for map coordinates.
- **`max_k`** - Maximum KD-tree neighbours requested per centre.

**Returns**

A list with `center_level`, `fov_level`, `edges`, `settings`, and
  the target-composition column names in `weighted_columns`.

**Example**

```r
cells <- data.frame(
  sample_id = "S1", FOV = "F1",
  CellType = c("DC", "DC", "CD4 T", "CD8 T"),
  x = c(0, 10, 5, 15), y = c(0, 0, 0, 0)
)
result <- target_neighbors(cells, "DC", c("CD4 T", "CD8 T"), radius = 20)
result$fov_level
```

<a id="compare-groups"></a>
### `compare_groups()`

Compare groups with a statistical test

Runs a Wilcoxon rank-sum (default) or Student t test between groups for a
numeric variable, either globally or separately for each level of a
`by` column (e.g. cell type). All pairwise group combinations are tested.
Observations are assumed to be independent. For imaging studies with
multiple FOVs per biological sample, first aggregate with
[neighbor_sample_summary()] to avoid treating FOVs as biological replicates.

**Usage**

```r
compare_groups(
  data,
  value_col,
  group_col,
  by = NULL,
  method = c("wilcox", "t"),
  paired = FALSE
)
```

**Parameters**

- **`data`** - A data frame with a numeric value column, a group column and
optionally a `by` column.
- **`value_col`** - Name of the numeric column to compare.
- **`group_col`** - Name of the group column.
- **`by`** - Optional column whose levels are compared separately (e.g.
`"CellType"`).
- **`method`** - Test to use: `"wilcox"` (Wilcoxon rank-sum / Mann-Whitney)
or `"t"` (Student t test).
- **`paired`** - Logical; if `TRUE`, a paired test is used.

**Returns**

A data frame with one row per comparison: `group1`, `group2`,
  `statistic`, `p.value`, `n1`, `n2`, `method`, plus the `by` column when
  provided.

**Example**

```r
df <- data.frame(
  Group = rep(c("R", "NR"), each = 10),
  value = c(rnorm(10, 10, 2), rnorm(10, 14, 2))
)
compare_groups(df, value_col = "value", group_col = "Group")
```

<a id="compare-neighbors"></a>
### `compare_neighbors()`

Compare neighbour composition at sample or FOV level

Provides three explicitly labelled analysis strategies for grouped FOV
neighbour summaries. `"sample_wilcox"` first averages FOVs within each
biological sample. `"fov_wilcox"` treats FOVs as analysis units and is
intended as an exploratory reproduction of the original workflow.
`"cluster_permutation"` retains the FOV input but permutes group labels for
whole biological samples; all FOVs from a sample always move together.

**Usage**

```r
compare_neighbors(
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
  seed = 2608L
)
```

**Parameters**

- **`data`** - Long FOV-level neighbour summary containing sample, FOV, group,
neighbour cell type, and numeric value columns.
- **`analysis_level`** - Either `"sample"` or `"FOV"`.
- **`method`** - One of `"sample_wilcox"`, `"fov_wilcox"`, or
`"cluster_permutation"`. The method must match `analysis_level`.
- **`value_col`** - Numeric neighbour statistic, usually `"mean_prop"`.
- **`group_col`, `sample_col`, `fov_col`, `celltype_col`** - Column names.
- **`group_order`** - Optional two-element vector fixing comparison direction;
the reported effect is first group minus second group.
- **`p_adjust`** - Multiple-testing correction passed to [stats::p.adjust()].
Use `"none"` to retain raw P values only.
- **`wilcox_exact`** - Passed to [stats::wilcox.test()] for Wilcoxon methods.
- **`permutations`** - Number of Monte Carlo permutations when the number of
possible sample-label allocations exceeds `max_exact`.
- **`max_exact`** - Maximum number of allocations to enumerate exactly.
- **`seed`** - Random seed used only for Monte Carlo permutations.

**Returns**

A data frame with one row per neighbour cell type, group estimates,
  effect, P values, sample/FOV counts, method labels, and permutation count.

**Details**

Ordinary FOV-level Wilcoxon ignores dependence between FOVs from the same
sample and can underestimate uncertainty. The cluster-permutation method
uses an equal-sample-weighted difference in mean FOV values and permutes
labels only between biological samples. It requires samples to be
independent and exchangeable under the null hypothesis.

Methodological references:

- Winkler AM et al. (2014). Permutation inference for the general linear
  model. NeuroImage 92:381-397. doi:10.1016/j.neuroimage.2014.01.060.
- Rosner B, Glynn RJ, Lee MLT (2003). Incorporation of clustering effects
  for the Wilcoxon rank sum test. Biometrics 59:1089-1098.
- Datta S, Satten GA (2005). Rank-sum tests for clustered data. JASA
  100:908-915. doi:10.1198/016214504000001583.

**Example**

```r
d <- expand.grid(
  Sample = paste0("S", 1:6), FOV = paste0("F", 1:3),
  CellType = c("cDC", "CD8Tn"), stringsAsFactors = FALSE
)
d$Group <- ifelse(d$Sample %in% paste0("S", 1:3), "R", "NR")
d$mean_prop <- seq_len(nrow(d)) / nrow(d)
compare_neighbor_groups(d, "sample", "sample_wilcox")
compare_neighbor_groups(d, "FOV", "cluster_permutation")
```

## Plotting

<a id="plot-cell-map"></a>
### `plot_cell_map()`

Plot a spatial map of cells coloured by cell type

Scatter plot of cell centroids, coloured by a categorical column (cell
type by default). With multiple samples the panels are faceted by sample,
or pass `sample_id` to plot a single sample. This is the "display"
companion of the annotation step and helps to visually inspect where each
cell type sits in the tissue.

**Usage**

```r
plot_cell_map(
  data,
  x = NULL,
  y = NULL,
  color = NULL,
  sample_col = "sample_id",
  sample_id = NULL,
  facet = TRUE,
  point_size = 0.6,
  alpha = 0.6,
  palette = NULL,
  title = NULL
)
```

**Parameters**

- **`data`** - A data frame with x/y coordinates and a colour column.
- **`x`, `y`** - Names of the centroid columns. Defaults to QuPath-style
`"Centroid X µm"` / `"Centroid Y µm"` (or normalised matches).
- **`color`** - Name of the column used to colour the points. Defaults to a
column named `"CellType"`, then `"Name"`, then `"Classification"`.
- **`sample_col`** - Name of the sample column used for faceting/filtering.
- **`sample_id`** - Optional sample id; if given, only that sample is plotted
(and used as the default title).
- **`facet`** - Logical; if `TRUE` (default) facet by `sample_col` when more
than one sample is present and `sample_id` is `NULL`.
- **`point_size`** - Point size passed to [ggplot2::geom_point()].
- **`alpha`** - Point transparency passed to [ggplot2::geom_point()].
- **`palette`** - Optional character vector of colours; by default a
colour-blind friendly palette is used.
- **`title`** - Optional plot title.

**Returns**

A [ggplot2::ggplot()] object.

**Example**

```r
data(mif_cells)
data(mif_cutoffs)
data(mif_rules)
marker_map <- c(
  CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
  CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
  PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
)
annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
plot_cell_map(annotated)
plot_cell_map(annotated, sample_id = "C1001")
```

<a id="plot-cells-by-sample"></a>
### `plot_cells_by_sample()`

Save one cell map per sample

Loops over samples and saves an individual [plot_cell_map()] figure per
sample into `output_dir`, named `<sample>_cell_map.<filetype>`.

**Usage**

```r
plot_cells_by_sample(
  data,
  output_dir,
  filetype = "png",
  width = 8,
  height = 6,
  dpi = 150,
  ...
)
```

**Parameters**

- **`data`** - A data frame with sample, centroid and colour columns.
- **`output_dir`** - Directory where figures are written (created if needed).
- **`filetype`** - Output file type passed to [ggplot2::ggsave()] (e.g.
`"png"` or `"pdf"`).
- **`width`, `height`, `dpi`** - Figure size in inches and resolution.
- **`...`** - Further arguments passed to [plot_cell_map()] (e.g. `color`,
`palette`, `point_size`).

**Returns**

Invisibly, the file paths written.

**Example**

```r
data(mif_cells)
data(mif_cutoffs)
data(mif_rules)
marker_map <- c(
  CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
  CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
  PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
)
annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)

# Requires local input files
plot_cell_maps_by_sample(annotated, "figures")
```

<a id="plot-proportions"></a>
### `plot_proportions()`

Plot cell-type proportions per sample (stacked bars)

Stacked bar chart showing, for each sample, either the proportion or the
absolute count of each cell type.

**Usage**

```r
plot_proportions(
  data,
  celltype_col = "CellType",
  sample_col = "sample_id",
  levels = NULL,
  proportion = TRUE,
  palette = NULL,
  title = NULL
)
```

**Parameters**

- **`data`** - A data frame with a cell-type column and a sample column.
- **`celltype_col`** - Name of the cell-type column.
- **`sample_col`** - Name of the sample column.
- **`levels`** - Optional ordering of cell types (bottom to top of the bars).
- **`proportion`** - Logical; if `TRUE` (default) the y axis shows the
proportion within each sample, otherwise absolute counts.
- **`palette`** - Optional character vector of colours.
- **`title`** - Optional plot title.

**Returns**

A [ggplot2::ggplot()] object.

**Example**

```r
data(mif_cells)
data(mif_cutoffs)
data(mif_rules)
marker_map <- c(
  CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
  CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
  PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
)
annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
plot_celltype_proportions(annotated)
```

<a id="plot-markers"></a>
### `plot_markers()`

Plot marker intensity distributions (with optional cutoffs)

Density plots of the fluorescence intensity of each marker, one panel per
marker, coloured by sample. When a cutoff table is supplied, a vertical
dashed line is drawn at each sample's cutoff. This helps to choose and
validate the channel thresholds used by [apply_cutoffs()].

**Usage**

```r
plot_markers(
  data,
  marker_map,
  sample_col = "sample_id",
  cutoffs = NULL,
  cutoff_id_col = NULL,
  facet_sample = FALSE,
  alpha = 0.4,
  palette = NULL
)
```

**Parameters**

- **`data`** - A data frame with the intensity columns named by `marker_map`.
- **`marker_map`** - Named character vector mapping marker names to intensity
columns (see [apply_cutoffs()]).
- **`sample_col`** - Name of the sample column.
- **`cutoffs`** - Optional cutoff data frame; if given, per-sample cutoff
lines are drawn.
- **`cutoff_id_col`** - Name of the sample-id column in `cutoffs` (defaults to
the first column).
- **`facet_sample`** - Logical; if `TRUE` (default) also facet by sample.
- **`alpha`** - Density fill transparency.
- **`palette`** - Optional character vector of colours.

**Returns**

A [ggplot2::ggplot()] object.

**Example**

```r
data(mif_cells)
data(mif_cutoffs)
marker_map <- c(
  CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
  CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
  PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
)
plot_marker_distributions(mif_cells, marker_map, cutoffs = mif_cutoffs)
```

<a id="plot-fov"></a>
### `plot_fov()`

Plot the FOV spatial distribution

Mirrors the FOV distribution figure of the original workflow
(`mIF_Tzone_02_GridArea_Process.R`): cell centroids coloured by FOV with
`geom_point(size = 0.4, alpha = 0.6)`, a turbo colour ramp from the viridis
package, `theme_classic()` without grid lines, fixed aspect ratio, and axis
labels "Centroid X (µm)" / "Centroid Y (µm)". The title is
"FOV distribution - <sample>".

**Usage**

```r
plot_fov(
  data,
  fov_col = "FOV",
  sample_col = "sample_id",
  sample_id = NULL,
  x = NULL,
  y = NULL,
  point_size = 0.4,
  alpha = 0.6,
  palette = NULL,
  show_legend = FALSE,
  title = NULL,
  x_lab = "Centroid X (<U+00B5>m)",
  y_lab = "Centroid Y (<U+00B5>m)"
)
```

**Parameters**

- **`data`** - A data frame with FOV labels and centroid coordinates (output
of [assign_fov()]).
- **`fov_col`** - Name of the FOV column.
- **`sample_col`** - Name of the sample column (used for faceting/filtering).
- **`sample_id`** - Optional sample id; if given, only that sample is plotted.
- **`x`, `y`** - Names of the centroid columns. Defaults to QuPath-style
`"Centroid X µm"` / `"Centroid Y µm"`.
- **`point_size`** - Point size (default 0.4, as in the original).
- **`alpha`** - Point transparency (default 0.6, as in the original).
- **`palette`** - Optional character vector of colours (one per FOV). If
`NULL`, a turbo palette is generated with `viridis::viridis(30,
option = "turbo")` (recycled and shuffled when > 30 FOVs).
- **`show_legend`** - Logical; if `TRUE` show the FOV legend (default `FALSE`
as in the original).
- **`title`** - Optional plot title (defaults to "FOV distribution -
<sample>").
- **`x_lab`, `y_lab`** - Axis labels (defaults match the original figure).

**Returns**

A [ggplot2::ggplot()] object.

**Example**

```r
data(mif_cells)
data(mif_cutoffs)
data(mif_rules)
marker_map <- c(
  CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
  CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
  PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
)
annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
fov <- assign_fov(annotated)
plot_fov_distribution(fov)
plot_fov_distribution(fov, sample_id = "C1001")
```

<a id="plot-fovs-by-sample"></a>
### `plot_fovs_by_sample()`

Save one FOV distribution plot per sample

Loops over samples and saves an individual [plot_fov_distribution()]
figure per sample into `output_dir`, named
`<sample>_FOV_distribution.<filetype>`.

**Usage**

```r
plot_fovs_by_sample(
  data,
  output_dir,
  filetype = "png",
  width = 8,
  height = 6,
  dpi = 150,
  ...
)
```

**Parameters**

- **`data`** - A data frame with sample, FOV and centroid columns (output of
[assign_fov()]).
- **`output_dir`** - Directory where figures are written (created if needed).
- **`filetype`** - Output file type passed to [ggplot2::ggsave()].
- **`width`, `height`, `dpi`** - Figure size in inches and resolution.
- **`...`** - Further arguments passed to [plot_fov_distribution()].

**Returns**

Invisibly, the file paths written.

**Example**

```r
data(mif_cells)
data(mif_cutoffs)
data(mif_rules)
marker_map <- c(
  CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
  CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
  PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
)
annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
fov <- assign_fov(annotated)

# Requires local input files
plot_fov_distributions_by_sample(fov, "figures")
```

<a id="plot-qc-summary"></a>
### `plot_qc_summary()`

Plot QC metrics by sample

Box plot (with optional jitter points) of one or more per-FOV QC metrics,
faceted by metric, in the original workflow style: `theme_classic()`,
black axis text (size 10), bold axis titles (size 11).

**Usage**

```r
plot_qc_summary(
  qc,
  metrics = c("TotalCells", "DC_count", "Tcell_count", "DC_ratio", "Tcell_ratio"),
  sample_col = "Sample",
  color_by = NULL,
  jitter = TRUE,
  palette = NULL,
  title = NULL
)
```

**Parameters**

- **`qc`** - The QC summary returned by [qc_fov_summary()].
- **`metrics`** - Character vector of numeric QC metrics to plot.
- **`sample_col`** - Name of the sample column in `qc`.
- **`color_by`** - Optional column (e.g. `LowCellFlag`) used to colour the
points.
- **`jitter`** - Logical; if `TRUE`, individual FOV values are shown as
jittered points.
- **`palette`** - Optional character vector of colours.
- **`title`** - Optional plot title.

**Returns**

A [ggplot2::ggplot()] object.

**Example**

```r
data(mif_cells)
data(mif_cutoffs)
data(mif_rules)
marker_map <- c(
  CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
  CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
  PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
)
annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
fov <- assign_fov(annotated)
qc <- qc_fov_summary(fov)
plot_qc_summary(qc)
```

<a id="plot-metric-density"></a>
### `plot_metric_density()`

Plot the density of a QC metric (original workflow style)

Density plot of a per-FOV QC metric (e.g. `DC_count`), mirroring the QC
density figures of `mIF_Tzone_05_QC_Visual.R`: steelblue fill,
`alpha = 0.5`, black outline, `theme_classic()` with black axis text and
bold axis titles. By default all samples are pooled into one density;
pass `by_sample = TRUE` to colour by sample instead.

**Usage**

```r
plot_metric_density(
  qc,
  metric = "DC_count",
  sample_col = "Sample",
  by_sample = FALSE,
  cutoffs = NULL,
  palette = NULL,
  title = NULL
)
```

**Parameters**

- **`qc`** - The QC summary returned by [qc_fov_summary()].
- **`metric`** - Numeric QC metric to plot.
- **`sample_col`** - Name of the sample column in `qc`.
- **`by_sample`** - Logical; if `TRUE`, densities are coloured by sample.
- **`cutoffs`** - Optional numeric vector; vertical dashed lines drawn at
these values (e.g. the quantile boundaries).
- **`palette`** - Optional character vector of colours (when
`by_sample = TRUE`).
- **`title`** - Optional plot title.

**Returns**

A [ggplot2::ggplot()] object.

**Example**

```r
data(mif_cells)
data(mif_cutoffs)
data(mif_rules)
marker_map <- c(
  CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
  CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
  PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
)
annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
fov <- assign_fov(annotated)
qc <- qc_fov_summary(fov)
plot_metric_density(qc, metric = "DC_count")
```

<a id="plot-dc-counts"></a>
### `plot_dc_counts()`

Plot the original-workflow DC-count density

Reproduces the density figures in `mIF_Tzone_05_QC_Visual.R` without
changing their geometry, colours, typography, or labels.

**Usage**

```r
plot_dc_counts(data, metric = "DC_count", filtered = FALSE)
```

**Parameters**

- **`data`** - A data frame containing a numeric DC-count column.
- **`metric`** - Name of the DC-count column.
- **`filtered`** - Logical; use the original title for the central 90 percent
filtered data when `TRUE`.

**Returns**

A [ggplot2::ggplot()] object. Save at 8 by 4 inches to reproduce
  the original output dimensions.

**Example**

```r
qc <- data.frame(DC_count = c(0, 1, 2, 2, 3, 5, 8))
plot_dc_count_distribution(qc)
```

<a id="plot-groups"></a>
### `plot_groups()`

Group comparison violin/box plot (original workflow style)

Reproduces the group-comparison figure of the original workflow
(`mIF_Tzone_05_QC_Visual.R`): opaque violin plots filled with the group
colours, white boxplots (`width = 0.2, alpha = 0.5`), `ggthemes::theme_few`
with a black panel border, one facet row, a percentage axis
(`values x 100`, breaks every 2.5) and pairwise Wilcoxon p-values drawn by
`ggpubr::stat_compare_means` (`label = "p.format"`, `size = 3.5`).

**Usage**

```r
plot_groups(
  data,
  value_col,
  group_col,
  by = NULL,
  method = "wilcox.test",
  comparisons = NULL,
  geom = c("both", "box", "violin"),
  percent = TRUE,
  y_lab = "% of DC 20<U+00B5>m neighbors",
  palette = NULL,
  jitter = FALSE,
  point_size = 0.8,
  alpha = 0.5,
  title = NULL,
  iqr_filter = TRUE
)
```

**Parameters**

- **`data`** - A data frame with a numeric value column, a group column and
optionally a `by` column.
- **`value_col`** - Name of the numeric column to plot.
- **`group_col`** - Name of the group column (x axis).
- **`by`** - Optional column whose levels are faceted in one row (e.g.
`"CellType"`).
- **`method`** - Test used for p-values (passed to
`ggpubr::stat_compare_means`), e.g. `"wilcox.test"` or `"t.test"`.
- **`comparisons`** - Optional list of group pairs to test, e.g.
`list(c("R", "NR"))`. If `NULL`, all pairwise combinations are used.
- **`geom`** - Plot geometry: `"both"` (violin + box, default, as in the
original), `"box"` or `"violin"`.
- **`percent`** - Logical; if `TRUE` (default, original behaviour), values are
multiplied by 100 and the y axis uses breaks every 2.5.
- **`y_lab`** - Y axis label (default `"
original).
- **`palette`** - Optional character vector of colours, one per group. Use a
named vector to fix colours to groups, e.g.
`c(R = "#E64B35", NR = "#4DBBD5")`. Defaults to the NPG red/blue pair
for two groups.
- **`jitter`** - Logical; if `TRUE`, individual points are shown (the original
figure does not show them).
- **`point_size`, `alpha`** - Jitter point size and transparency.
- **`title`** - Optional plot title. The original figure has no title.
- **`iqr_filter`** - Logical; apply the original per-cell-type 1.5-IQR filter
before plotting.

**Returns**

A [ggplot2::ggplot()] object.

**Example**

```r
df <- data.frame(
  Group = rep(c("R", "NR"), each = 30),
  CellType = rep(rep(c("DC", "CD8"), each = 15), 2),
  value = c(rnorm(15, 0.2, 0.05), rnorm(15, 0.35, 0.05),
            rnorm(15, 0.15, 0.05), rnorm(15, 0.3, 0.05))
)
plot_group_comparison(df, value_col = "value", group_col = "Group",
                      by = "CellType")
```

<a id="plot-neighbors"></a>
### `plot_neighbors()`

Plot neighbour-group comparisons with explicit test labels

Plots the data used by [compare_neighbors()] and annotates each cell
type with its P value. The subtitle records the analysis and permutation
units so exploratory FOV-level results cannot be mistaken for clustered or
sample-level inference.

**Usage**

```r
plot_neighbors(
  data,
  statistics,
  value_col = "mean_prop",
  group_col = "Group",
  sample_col = "Sample",
  celltype_col = "CellType",
  p_column = "p.value",
  percent = TRUE,
  palette = c(R = "#E64B35", NR = "#4DBBD5"),
  geom = c("both", "violin", "boxplot"),
  show_points = TRUE,
  point_size = 1.4,
  alpha = 0.65,
  violin_alpha = 0.4,
  box_color = "black",
  visual_iqr_filter = FALSE,
  visual_iqr_multiplier = 1.5,
  title = NULL,
  y_lab = NULL
)
```

**Parameters**

- **`data`** - Grouped long FOV-level neighbour summary.
- **`statistics`** - Output of [compare_neighbors()].
- **`value_col`, `group_col`, `sample_col`, `celltype_col`** - Column names in `data`.
- **`p_column`** - P-value column in `statistics`. Defaults to the raw
`"p.value"`; use `"p.adjust"` to display multiplicity-adjusted values.
- **`percent`** - Multiply proportions by 100 for display.
- **`palette`** - Named group colours.
- **`geom`** - One of `"both"`, `"violin"`, or `"boxplot"`.
- **`show_points`** - Show individual sample/FOV observations.
- **`point_size`, `alpha`** - Point size and transparency. Points use a black
outline and the group colour as fill so they remain visible over violins.
- **`violin_alpha`** - Transparency of the violin fill.
- **`box_color`** - Colour of boxplot outlines, whiskers and median lines.
- **`visual_iqr_filter`** - If `TRUE`, observations outside the pooled
cell-type-specific `Q1 - visual_iqr_multiplier * IQR` to
`Q3 + visual_iqr_multiplier * IQR` interval are omitted from the figure
only. Statistics supplied in `statistics` are not recalculated.
- **`visual_iqr_multiplier`** - IQR multiplier used only when
`visual_iqr_filter = TRUE`.
- **`title`** - Optional title.
- **`y_lab`** - Optional y-axis label.

**Returns**

A [ggplot2::ggplot()] object.

**Example**

```r
d <- expand.grid(Sample = paste0("S", 1:6), FOV = c("F1", "F2"),
                 CellType = "CD8 T", stringsAsFactors = FALSE)
d$Group <- ifelse(d$Sample %in% paste0("S", 1:3), "R", "NR")
d$mean_prop <- seq_len(nrow(d)) / nrow(d)
tests <- compare_neighbors(d, analysis_level = "sample")
plot_neighbors(d, tests)
```

<a id="plot-celltype-radar"></a>
### `plot_celltype_radar()`

Radar chart of cell-type composition by group (fmsb style)

Draws a group-comparison radar chart with [fmsb::radarchart()], following
the reference fmsb code of the original workflow: axes start at the top and
are connected clockwise, polygons are filled with a semi-transparent colour
(`alpha = 0.1`), borders are solid (`plwd = 1`), grid lines are grey
(`cglwd = 0.8`), axis labels are black, and the legend is drawn on the
right with `pch = 20`. Cell types listed in `exclude` (by default
`"Other"`) are removed before computing proportions, and (for
`value = "proportion"`) values are shown as percentages with `seg = 4`
grid segments.

**Usage**

```r
plot_celltype_radar(
  data,
  celltype_col = "CellType",
  sample_col = "sample_id",
  group_col = "Group",
  value = c("proportion", "count"),
  exclude = "Other",
  percent = TRUE,
  levels = NULL,
  palette = NULL,
  file = NULL,
  width = 15,
  height = 15,
  res = 600,
  family = "serif",
  seg = 4,
  plwd = 1,
  cglwd = 0.8,
  vlcex = 0.8,
  legend_cex = 1.2,
  title = NULL,
  ...
)
```

**Parameters**

- **`data`** - A data frame with cell-type, sample and group columns (e.g.
annotated cells after [add_groups()]).
- **`celltype_col`** - Name of the cell-type column.
- **`sample_col`** - Name of the sample column.
- **`group_col`** - Name of the group column.
- **`value`** - What to average per group: `"proportion"` (default, cell-type
proportion within each sample, in percent) or `"count"` (mean cell
count).
- **`exclude`** - Cell types to remove before analysis (default `"Other"`).
Use `NULL` to keep every cell type.
- **`percent`** - Logical; if `TRUE` (default) and `value = "proportion"`,
values are multiplied by 100 and the axis labels show percentages.
- **`levels`** - Optional ordering of cell types around the radar (defaults to
the order of appearance).
- **`palette`** - Optional character vector of colours, one per group. Use a
named vector to fix colours to groups, e.g.
`c(R = "#E64B35", NR = "#4DBBD5")`. Defaults to the reference palette
`c("#1B9E77", "#7570B3", "#D95F02")`.
- **`file`** - Optional output file (`.png`, `.jpg`, `.pdf`, `.bmp`). If
`NULL`, the chart is drawn on the current graphics device.
- **`width`, `height`** - Figure size in centimetres (for raster output).
- **`res`** - Resolution in dpi for raster output (default 600, as in the
reference).
- **`family`** - Font family (default `"serif"`, as in the reference).
- **`seg`** - Number of grid segments (default 4).
- **`plwd`** - Polygon border line width (default 1).
- **`cglwd`** - Grid line width (default 0.8).
- **`vlcex`** - Vertex (cell-type) label size (default 0.8).
- **`legend_cex`** - Legend text size (default 1.2).
- **`title`** - Optional plot title.
- **`...`** - Further arguments passed to [fmsb::radarchart()].

**Returns**

Invisibly, the data frame used by `radarchart` (rows: `max`, `min`,
  then one row per group; columns: cell types).

**Details**

Because the chart uses base graphics, it is drawn either to `file` (when
given) or to the current graphics device.

**Example**

```r
data(mif_cells)
data(mif_cutoffs)
data(mif_rules)
data(mif_groups)
marker_map <- c(
  CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
  CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
  PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
)
annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
annotated <- add_groups(annotated, mif_groups)

# Requires local input files
plot_celltype_radar(annotated)
plot_celltype_radar(annotated, file = "radar.png")
```

<a id="plot-ring-radar"></a>
### `plot_ring_radar()`

Plot grouped neighbour proportions as a ring radar chart

Draws one closed radar line per group and an outer segmented annotation
ring. The function is designed for sample-level output from
[neighbor_sample_summary()], but accepts any long data frame containing a
group, neighbour cell type, and numeric proportion.

**Usage**

```r
plot_ring_radar(
  data,
  group_col = "Group",
  category_col = "CellType",
  value_col = "mean_prop",
  group_colors = c(R = "#C44732", NR = "#087DB5"),
  ring_colors = "#D8C0EA",
  category_order = NULL,
  category_labels = NULL,
  p_values = NULL,
  p_category_col = category_col,
  p_value_col = "p.value",
  show_p = TRUE,
  breaks = NULL,
  limits = NULL,
  ring_width = 0.28,
  ring_alpha = 1,
  line_width = 1,
  point_size = 2.2,
  label_size = 3.6,
  p_size = 3,
  grid_label_size = 2.8,
  title = NULL
)
```

**Parameters**

- **`data`** - Long-format data frame.
- **`group_col`** - Name of the grouping column.
- **`category_col`** - Name of the neighbour cell-type/category column.
- **`value_col`** - Name of the numeric proportion column.
- **`group_colors`** - Colours for group lines and points. A named vector is
recommended, for example `c(R = "#C44732", NR = "#087DB5")`.
- **`ring_colors`** - One outer-ring colour, an unnamed vector in category
order, or a named vector keyed by category.
- **`category_order`** - Optional category order around the circle.
- **`category_labels`** - Optional named character vector replacing displayed
category labels without changing the data.
- **`p_values`** - Optional data frame containing category and p-value columns.
When omitted and `show_p = TRUE`, two-group Wilcoxon tests are calculated
from the rows in `data`.
- **`p_category_col`, `p_value_col`** - Column names in `p_values`.
- **`show_p`** - Logical; display one p-value per category.
- **`breaks`** - Numeric radial grid breaks. If `NULL`, pretty breaks are used.
- **`limits`** - Numeric radial limits for the data region. If `NULL`, limits
are inferred from zero and the observed values.
- **`ring_width`** - Width of the outer ring relative to the data radius.
- **`ring_alpha`** - Outer-ring opacity.
- **`line_width`, `point_size`** - Radar line width and point size.
- **`label_size`, `p_size`, `grid_label_size`** - Text sizes for ring labels,
p-values, and radial tick labels.
- **`title`** - Optional plot title.

**Returns**

A [ggplot2::ggplot()] object.

**Example**

```r
d <- expand.grid(
  Sample = paste0("S", 1:6),
  CellType = c("T cell", "DC", "Other")
)
d$Group <- ifelse(d$Sample %in% paste0("S", 1:3), "R", "NR")
d$mean_prop <- seq(0.1, 0.6, length.out = nrow(d))
plot_neighbor_ring_radar(
  d, group_col = "Group", category_col = "CellType",
  value_col = "mean_prop",
  group_colors = c(R = "#C44732", NR = "#087DB5"),
  ring_colors = "#D8C0EA"
)
```

<a id="plot-nearest-rank"></a>
### `plot_nearest_rank()`

Bar chart of nearest-neighbour cell-type ranking

Visualises the output of [nearest_neighbor_rank()]: for each nearest cell
type, the number of FOVs in which it appears among the `top_n` most
frequent nearest neighbours. Optionally faceted by sample. Uses the
[ggplot2::theme_classic()] publication style.

**Usage**

```r
plot_nearest_rank(
  data,
  top_n = 3,
  nearest_col = "nearestNGB",
  rank_col = "Rank",
  count_col = "Count",
  sample_col = "sample_id",
  by_sample = FALSE,
  fill_colour = "#4DBBD5",
  palette = NULL,
  title = NULL
)
```

**Parameters**

- **`data`** - Output of [nearest_neighbor_rank()] (or any data frame with
nearest-type, rank and count columns).
- **`top_n`** - Keep only types whose `Rank` is `<= top_n` (default 3).
- **`nearest_col`** - Name of the nearest-neighbour type column.
- **`rank_col`** - Name of the rank column.
- **`count_col`** - Name of the count column (used only for validation).
- **`sample_col`** - Name of the sample column (used when `by_sample = TRUE`).
- **`by_sample`** - Logical; if `TRUE`, bars are split/faceted by sample.
- **`fill_colour`** - Fill colour for the aggregated bars.
- **`palette`** - Optional character vector of colours (when `by_sample =
TRUE`).
- **`title`** - Optional plot title.

**Returns**

A [ggplot2::ggplot()] object.

**Example**

```r
data(mif_cells)
data(mif_cutoffs)
data(mif_rules)
marker_map <- c(
  CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
  CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
  PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
)
annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
fov <- assign_fov(annotated)
nb <- find_neighbors(fov, radius = 60)
rank <- nearest_neighbor_rank(nb)
plot_nearest_rank(rank, top_n = 3)
```

<a id="plot-heatpie"></a>
### `plot_heatpie()`

Plot an FOV heat-pie map for centre--target neighbourhoods

Tile colour represents mean target-type richness, pie sectors represent
weighted target composition, and pie area represents mean weighted target
abundance per centre. By default, target-rich FOVs use the project's
pale-pink to deep-rose palette.

**Usage**

```r
plot_heatpie(
  result,
  target_colors = NULL,
  richness_colors = c("#FECEF5", "#FECEF5", "#EC407A", "#B6066F"),
  other_selected_color = "#7ABD31",
  centre_only_color = "#FFFFCC",
  other_isolated_color = "#CBCBCB",
  no_centre_color = "white",
  border_color = "#666666",
  max_richness = NULL,
  pie_max = NULL,
  pie_radius_fraction = 0.38,
  title = NULL,
  subtitle = NULL
)
```

**Parameters**

- **`result`** - Output from [target_neighbors()].
- **`target_colors`** - Named colours for target types. Names may be target
type labels or the corresponding `Weighted_...` columns.
- **`richness_colors`** - Colours from low to high richness.
- **`other_selected_color`** - Colour for FOVs dominated by user-specified
non-target types.
- **`centre_only_color`** - Colour for centre--centre FOVs without targets.
- **`other_isolated_color`** - Colour for other/isolated centre FOVs.
- **`no_centre_color`** - Colour for FOVs without a centre cell.
- **`border_color`** - Grid-border colour.
- **`max_richness`** - Upper limit of the richness legend. Defaults to the
number of requested target types.
- **`pie_max`** - Optional common upper reference for pie size.
- **`pie_radius_fraction`** - Maximum pie radius as a fraction of grid size.
- **`title`, `subtitle`** - Optional plot text.

**Returns**

A ggplot object.

**Example**

```r
cells <- data.frame(
  sample_id = "S1", FOV = "F1", CellType = c("DC", "CD8 T"),
  x = c(0, 10), y = c(0, 0)
)
result <- target_neighbors(cells, "DC", "CD8 T", radius = 20)
plot_heatpie(result)
```

<a id="theme-mif"></a>
### `theme_mif()`

cellfiers publication theme

A clean theme in the style of the original workflow figures: classic axes,
a thin black panel border, black axis text, bold axis/strip titles and a
centred bold plot title. Use it on any ggplot, e.g.
`my_plot + theme_mif()`.

**Usage**

```r
theme_mif(base_size = 11)
```

**Parameters**

- **`base_size`** - Base font size in points.

**Returns**

A [ggplot2::theme()] object.

**Example**

```r
library(ggplot2)
ggplot(mtcars, aes(mpg, wt)) + geom_point() + theme_mif()
```

## Compatibility

<a id="analyze-dc-tcell-tribe"></a>
### `analyze_dc_tcell_tribe()`

Analyse DC-centred target-T neighbourhoods (compatibility wrapper)

This legacy name is retained for existing scripts. New analyses should use
[target_neighbors()], which accepts arbitrary centre and target
cell types.

**Usage**

```r
analyze_dc_tcell_tribe(
  data,
  radius = 20,
  center_types = c("cDC", "cDC(Mature/Active)"),
  target_t_types = c("CD8Tn(Primed)", "CD8Teff(Active)", "CD4Tn(Primed)",
    "CD4Teff(Active)"),
  other_t_types = c("CD8Tn(Quiescent)", "CD8Teff(Resting)", "CD4Tn(Quiescent)",
    "CD4Teff(Resting)"),
  shared_weight = c("inverse_degree", "none"),
  sample_col = "sample_id",
  celltype_col = "CellType",
  fov_col = "FOV",
  x = NULL,
  y = NULL,
  grid_size = 200,
  max_k = 1000L
)
```

**Parameters**

- **`data`** - One biological sample of annotated, FOV-assigned cell data.
- **`radius`** - Search radius in micrometres. Values above 30 are allowed but
trigger a warning about computation and biological interpretation.
- **`center_types`** - Centre cell types.
- **`target_t_types`** - Target cell types.
- **`other_t_types`** - Optional other cell types used for labelling.
- **`shared_weight`** - Either `"inverse_degree"` (default, `1/m`) or
`"none"` (every centre--target edge has weight 1).
- **`sample_col`, `celltype_col`, `fov_col`** - Column names.
- **`x`, `y`** - Optional centroid column names.
- **`grid_size`** - FOV side length in micrometres, used for map coordinates.
- **`max_k`** - Maximum KD-tree neighbours requested per centre.

**Returns**

A list containing the legacy `dc_level` element plus `fov_level`,
  `edges`, `settings`, and `weighted_columns`.

**Example**

```r
cells <- data.frame(
  sample_id = "S1", FOV = "F1",
  CellType = c("cDC", "CD8Tn(Primed)"), x = c(0, 10), y = c(0, 0)
)
analyze_dc_tcell_tribe(cells)$dc_level
```
