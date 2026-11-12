# 快速入门

[首页](../README.md) · [安装](installation.md) · [函数参考](reference.md)

下面用包内模拟数据完成选区、注释、FOV 质控、邻域分析和热饼图。可在安装后直接运行。

## 1. 按选区和规则注释细胞

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

## 2. 栅格化和质控

```r
centres <- c("cDC", "cDC(Mature/Active)")
qc <- run_fov_qc(
  cells, grid_size = 200, min_distance = 500,
  dc_types = centres, max_dc_ratio = 0.70, probs = c(0, 1)
)
plot_qc_summary(qc$qc)
clean <- qc$data_final
```

这里用 `probs = c(0, 1)` 保留通过初筛的模拟 FOV。实际数据可按项目的质控方案设定分位数。质控结果保存在 `qc$qc`，通过质控的细胞在 `qc$data_final`。

## 3. 比较中心细胞的邻居组成

```r
neighbors <- analyze_neighbors(clean, radius = 20, center_types = centres)
grouped <- add_groups(neighbors$fov_summary, mif_groups, sample_col = "Sample")
tests <- compare_neighbors(grouped, analysis_level = "sample")
plot_neighbors(grouped, tests, geom = "boxplot")
```

`analyze_neighbors()` 在 FOV 内搜索并汇总邻居。`compare_neighbors()` 可按生物样本汇总后检验，也可做探索性的 FOV 级 Wilcoxon 检验或整样本标签置换。模拟数据的组间结果仅用于展示接口。

## 4. 计算指定目标细胞的邻域并出图

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

`target_neighbors()` 一次处理一个生物样本，样本内搜索可跨 FOV 边界，再把结果按中心细胞所属 FOV 汇总。半径和网格边长都以 µm 为单位。热饼图的底色为每个 FOV 内中心细胞的平均目标类型数，扇区为加权目标组成，饼图面积为每个中心细胞的平均加权目标丰度。一个目标细胞同时出现在 `m` 个中心细胞的邻域时，每条边贡献 `1/m`。

## 换成自己的数据

QuPath 导出时保留区域列（默认 `Parent`）、细胞质心坐标及荧光强度。确认坐标以 µm 记录，并将阈值表中的标记物逐个映射到正确的荧光测量列。规则表每行一个细胞类型，标记物列取值为 `+`、`-` 或 `NA`。规则按行匹配，首次命中的类型生效。

```r
# 改成自己的路径、样本编号、选区名称和通道映射。
cells <- read_qupath("path/to/qupath_export.txt", sample_id = "S01")
cutoffs <- read_cutoffs("path/to/cutoffs.tsv")
rules <- read_rules("path/to/rules.tsv")
annotated <- run_annotation(
  cells, cutoffs, rules, marker_map, regions = "SLN_Tregion"
)
```

批量导入 TSV 用 `read_qupath_batch()`。若 QuPath 导出为 CSV，可用 `readr::read_csv()`，补齐 `sample_id` 和样本内唯一的 `cell_id` 后传入 `run_annotation()`。按样本检查注释和 QC，再选择中心细胞、目标细胞和半径。输出的 ggplot 图可用 `ggplot2::ggsave()` 保存为 PNG 或 PDF。
