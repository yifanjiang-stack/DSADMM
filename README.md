# DSADMM

Helpers and an ADMM solver for quantile regression (double-splitting style).

## Install

```r
# install.packages("remotes")
remotes::install_github("YOUR_GITHUB/DSADMM")
```

## Usage

The intended entry point for reproducing the workflow is the single-run example script:

- `inst/examples/one_run_simulation.R`

### Run the example

Option 1 (recommended): clone the repository and run the script locally

```r
setwd("path/to/DSADMM")
source("inst/examples/one_run_simulation.R")
