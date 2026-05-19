# DSADMM

`DSADMM` packages the double-splitting ADMM implementation for penalized
rank regression used in the updated rank-ADMM simulation code.

The package exports reusable fitting, tuning, pairwise rank-transformation,
simulation, and evaluation helpers for:

- L1-penalized rank regression.
- SCAD-penalized rank regression using LLA/MM with ADMM inner solves.

## Example Settings

The scripts in `inst/examples/` preserve the simulation settings from the
updated Example 1 and Example 2 scripts:

- `n = 60`, `p = 200`
- AR(1) Gaussian design with `rho = 0.5`
- normal error with standard deviation `1`
- true beta has the first four entries equal to `1` and all others `0`
- pairwise rank transformation after centering/scaling
- `M = 2`, `G = 4`, `phi = 1 / 400`, `rho = 1.618`
- HBIC lambda selection with `max_iter = 500`
- final reported ADMM estimator with `max_iter = 400`
- raw beta L1 estimation error, with no thresholding or top-K evaluation

The SCAD example uses the updated damped LLA path with monotone SCAD objective
line search, initialized from the L1 ADMM estimator.

## Running Examples

After installation:

```r
library(DSADMM)
source(system.file("examples", "example1_L1.R", package = "DSADMM"))
source(system.file("examples", "example2_SCAD.R", package = "DSADMM"))
```

By default each example runs one seed and includes the competing methods from
the original scripts when their optional packages are installed. Missing
optional competitors are reported as `*_missing` rows. To run only DSADMM:

```powershell
$env:DSADMM_RUN_COMPETITORS = "FALSE"
```

To run more seeds:

```powershell
$env:DSADMM_EXAMPLE_REPS = "10"
Rscript -e "library(DSADMM); source(system.file('examples','example1_L1.R',package='DSADMM'))"
Rscript -e "library(DSADMM); source(system.file('examples','example2_SCAD.R',package='DSADMM'))"
```

`glmnet` is used in the examples to preserve the original cold-start rule.
If it is unavailable, set `DSADMM_ALLOW_ZERO_START=TRUE` to run a reduced
smoke check with zero initialization.

The comparison helpers are kept directly inside the two example scripts, as in
the previous simulation-script organization. They cover the previous-script
baselines: QPADM, padmmR when available, conquer, hqreg, QICD, rqPen, and the
subsampling/rqPen comparison.
Install the optional comparison packages before running the full comparison
examples if you want non-missing rows for every baseline.
