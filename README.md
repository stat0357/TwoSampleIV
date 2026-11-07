# TwoSampleIV

An R package for two-sample instrumental variables estimation from individual-level data (`tsiv()`) and regression summary statistics (`tsiv_sum()`).

## Install from GitHub

After uploading this package to a GitHub repository:

```r
install.packages("devtools")
devtools::install_github("YOUR_USERNAME/TwoSampleIV")
library(TwoSampleIV)
```

Replace `YOUR_USERNAME` with your GitHub username. The package depends on `geigen`. For clustered covariance examples, also install `sandwich`.

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

## Before publication

Replace the placeholder author/email in DESCRIPTION and copyright holder in LICENSE files. Confirm the intended license. Run `devtools::check()` after installing dependencies.
