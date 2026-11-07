# TwoSampleIV

An R package for two-sample instrumental variables estimation from individual-level data (`tsiv()`) or regression summary statistics (`tsiv_sum()`).

## Install from GitHub

After uploading this package to a GitHub repository:

```r
install.packages("devtools")
devtools::install_github("stat0357/TwoSampleIV")
library(TwoSampleIV)
```

The package depends on `geigen`, also install `sandwich`.

## Documentation

- [tsiv()](docs/tsiv.md)
- [tsiv_sum()](docs/tsiv_sum.md)

## Example summary

```r
summary(results1)                       # All available estimators
summary(results1, method = "two_step")
summary(results1, method = "cue")       # If CUE=TRUE
coef(results1, method = "ts2sls")
```

The summary-statistics function labels the first estimator *One-step*, whereas the individual-data function labels it *TS2SLS*. First-stage F diagnostics and Anderson–Rubin confidence sets are implemented for one endogenous regressor.
