# 安装 cellfiers

[首页](../README.md) · [快速入门](quickstart.md) · [函数参考](reference.md)

## R 与依赖

包声明需要 **R ≥ 4.1**。第三方运行时依赖均来自 CRAN。

| 用途 | R 包 |
| --- | --- |
| 数据读取、整理与列操作 | `readr`, `dplyr`, `tidyr`, `stringr`, `rlang` |
| 半径邻居搜索 | `RANN` |
| 可视化 | `ggplot2`, `ggpubr`, `ggthemes`, `ggnewscale`, `scales`, `scatterpie`, `viridis`, `fmsb` |
| R 自带 | `graphics`, `grDevices`, `stats`, `utils` |
| 从 GitHub 安装 | `remotes` |
| 构建教程及运行测试 | `knitr`, `rmarkdown`, `testthat`（≥ 3.0.0） |
| 修改源码后重新生成帮助 | `roxygen2`（维护者工具） |

安装包时会解析运行依赖。需要手动预装时，可以运行：

```r
install.packages(c(
  "readr", "dplyr", "tidyr", "stringr", "rlang", "RANN",
  "ggplot2", "ggpubr", "ggthemes", "ggnewscale", "scales",
  "scatterpie", "viridis", "fmsb"
))
```

部分依赖的新版可能要求更高的 R 版本。Windows 从源码编译依赖时，可能需要与 R 版本匹配的 Rtools。

## 从 GitHub 安装

在 R 中运行以下命令，从 GitHub 安装 cellfiers。

```r
install.packages("remotes")
remotes::install_github("Porpax/cellfiers", dependencies = NA, upgrade = "never")
library(cellfiers)
packageVersion("cellfiers")
```

也可以在 RStudio 中安装本地源码：

```r
install.packages("remotes")
remotes::install_local(
  "path/to/cellfiers", dependencies = NA, upgrade = "never"
)
library(cellfiers)
help("target_neighbors", package = "cellfiers")
```

接下来用[包内模拟数据跑完整流程](quickstart.md)。
