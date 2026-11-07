# `tsiv()` — individual-level two-sample IV

## Usage

```r
tsiv(dfY, dfX, Y, X, exog, inst, Vhety_str, Vhetx_str,
     full=TRUE, CUE=FALSE, homosk=FALSE, alpha=0.05)
```

`dfY` and `dfX` are outcome and exposure data frames. `Y` is the outcome column name; `X` is a character vector of endogenous exposure names. `inst` is a character vector of excluded instrument names and `exog` is a character vector of exogenous control names (`""` for none). All instruments and controls must exist in both samples.

`Vhety_str` and `Vhetx_str` are character strings containing covariance-matrix expressions evaluated after the reduced-form `lmY` and first-stage `lmX` regressions are fitted. The exposure covariance argument is spelled `Vhetx_str` (lowercase `x`). `homosk=TRUE` instead uses conventional `vcov()` matrices. `full=FALSE` estimates endogenous-variable coefficients from excluded-instrument coefficients; `full=TRUE` additionally includes the intercept and controls. `CUE=TRUE` requests the continuously updated estimator and, with one endogenous variable, the Anderson–Rubin confidence set. `alpha` sets its significance level.

## Example

```r
library(sandwich)
library(geigen)
source("tsivfuncAR_with_sum_onestep_stars_fixed.R")

res <- tsiv(
  dfY, dfX, Y="outcome", X="schooling",
  exog=c("age", "year"), inst=c("z1", "z2"),
  Vhety_str='vcovCL(lmY, cluster=dfY$state, type="HC1", cadjust=TRUE)',
  Vhetx_str='vcovCL(lmX, cluster=dfX$state, type="HC1", cadjust=TRUE)',
  full=FALSE, CUE=TRUE
)
summary(res)
summary(res, method="two_step")
coef(res, method="ts2sls")
```

The data-frame and variable names are placeholders. Cluster vectors must align with the observations used by the regressions.

## Returned object

The result is a `tsiv` list with `sample`, `ts2sls`, `two_step`, `hansen`, and, when available, `first_stage`, `cue`, `cue_hansen`, and `ar_confint`. `summary(res)` prints all available methods. For `tsiv()`, the first estimator is labelled **TS2SLS**. Four-decimal printing does not round stored estimates.

The first-stage diagnostics include conventional F and, for non-homoskedastic covariance specifications, robust F and generalized effective F diagnostics. These are implemented only for one endogenous regressor. The AR confidence set is also limited to one endogenous regressor and requires `CUE=TRUE`.

Coefficient p-values use large-sample normal z-tests. Hansen J inference is meaningful as an overidentification test only with more excluded instruments than endogenous regressors. The validity of heteroskedasticity- or cluster-robust inference depends on the covariance estimator and sampling assumptions. Covariance expressions are evaluated as R code; use trusted inputs only.
