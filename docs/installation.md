# Install cellfiers

[Home](../README.md) · [Quick start](quickstart.md) · [Function reference](reference.md)

## R version and dependencies

cellfiers requires **R ≥ 4.1**. All third-party runtime dependencies are available from CRAN.

| Purpose | R packages |
| --- | --- |
| Data import and manipulation | `readr`, `dplyr`, `tidyr`, `stringr`, `rlang` |
| Radius-based neighbour search | `RANN` |
| Visualisation | `ggplot2`, `ggpubr`, `ggthemes`, `ggnewscale`, `scales`, `scatterpie`, `viridis`, `fmsb` |
| Included with R | `graphics`, `grDevices`, `stats`, `utils` |
| Installation from GitHub | `remotes` |
| Building vignettes and running tests | `knitr`, `rmarkdown`, `testthat` (≥ 3.0.0) |
| Regenerating R help files | `roxygen2` (for maintainers) |

The installer resolves runtime dependencies. To install them separately, run the following commands.

```r
install.packages(c(
  "readr", "dplyr", "tidyr", "stringr", "rlang", "RANN",
  "ggplot2", "ggpubr", "ggthemes", "ggnewscale", "scales",
  "scatterpie", "viridis", "fmsb"
))
```

Recent versions of some dependencies may require a newer version of R. On Windows, compiling dependencies from source may require the version of Rtools that corresponds to your R installation.

## Install from GitHub

Run the following commands in R.

```r
install.packages("remotes")
remotes::install_github("Porpax/cellfiers", dependencies = NA, upgrade = "never")
library(cellfiers)
packageVersion("cellfiers")
```

To install from a local copy of the source package, run the following commands in R or RStudio.

```r
install.packages("remotes")
remotes::install_local(
  "path/to/cellfiers", dependencies = NA, upgrade = "never"
)
library(cellfiers)
help("target_neighbors", package = "cellfiers")
```

Next, follow the [quick start with bundled simulated data](quickstart.md).
