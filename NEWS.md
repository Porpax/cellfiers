# cellfiers 0.2.0

- Added shorter names for the centre–target workflow, group comparisons and
  plotting. The previous names remain exported with the same arguments and
  return values. See the [function reference](docs/reference.md).

# cellfiers 0.1.0

- Added `compare_neighbor_groups()` with selectable sample-level Wilcoxon,
  exploratory FOV-level Wilcoxon, and whole-sample cluster permutation tests.
- Added adjustable `max_dc_ratio` FOV QC (default 0.70), explicit exclusion
  flags, and final QC audit labels.

- Added QuPath import and region filtering.
- Added sample-specific marker calling and rule-based cell phenotyping.
- Added field-of-view rasterisation and quality-control summaries.
- Added `run_fov_qc()` to enforce low-cell filtering before extreme-FOV
  quantile filtering.
- Added radius-based neighbourhood analysis, group comparisons and plots.
- Restricted nearest-neighbour ranks to requested centre-cell types in
  `analyze_neighbors()` and retained distinct cells with identical centroids.
- Added sample-level neighbour summaries to avoid FOV pseudoreplication in
  group comparisons.
- Reimplemented the canonical FOV, DC-density, and R/NR violin figures
  directly from the original `mIF_Tzone` scripts with fixed visual defaults.
- Added a configurable ring radar chart for grouped centre-cell neighbour
  proportions.
