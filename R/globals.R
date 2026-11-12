# Variables created inside dplyr/tidyr verbs (NSE). Declared here to silence
# "no visible binding for global variable" notes from R CMD check.
utils::globalVariables(c(
  # plot_marker_distributions / plot_qc_summary pivots
  "cutoff", "intensity", "marker_intensity", "metric", "n", "sample_id",
  "value",
  # summarize_celltypes
  ".ct", ".ct_fct", ".prop", ".sm_fct", ".total",
  # qc_fov_summary
  ".__sample", ".__dc", ".__tcell",
  # assign_fov
  ".__fov", ".__fov_num", ".__gx", ".__gx_n", ".__gy", ".__gy_n",
  ".__max_x", ".__rowid", ".__tr",
  # remove_outliers_iqr
  ".__hi", ".__lo", ".__q1", ".__q3",
  # qc_fov_summary columns created by summarise/mutate
  "DC_count", "MeanCells", "Sample", "Tcell_count", "TotalCells",
  # plot_metric_density
  "xintercept",
  # neighbor analysis
  "CellType", "Count", "Proportion", "Rank", "count", "mean_count",
  "mean_prop", "n_center", "prop", "total_neighbors", ".__sample", ".__id", ".__rowid",
  # add_groups
  ".__group", ".__sample_key",
  # compare_groups / plot_group_comparison annotation
  "group1", "group2", "label", "p.label", "statistic", "p.value", "x",
  "x1", "x2", "xend", "y", "yend", ".__idx", ".__ymax", ".__cc", ".__total", "FOV_count"
))