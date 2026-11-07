# `tsiv_sum()` — two-sample IV from summary statistics

## Usage

```r
tsiv_sum(piy, Vhomy, Vhety, pix, Vhomx, Vhetx,
         exog, inst, full=TRUE, CUE=FALSE,
         homosk=FALSE, alpha=0.05,
         nY=NA_integer_, nX=NA_integer_,
         dependent="Not supplied", vcov_info=NULL)
```

The six primary inputs are:

| Argument | Required content |
|---|---|
| `piy` | Named reduced-form coefficient vector, including intercept, instruments and controls. |
| `Vhomy` | Conventional covariance matrix of `piy`. |
| `Vhety` | User-supplied alternative covariance matrix of `piy`. |
| `pix` | First-stage coefficient matrix; rows in the same order as `piy`, columns named for endogenous variables. |
| `Vhomx` | Conventional covariance matrix of the vectorized first-stage coefficients. |
| `Vhetx` | User-supplied alternative covariance matrix of the vectorized first-stage coefficients. |

For `p` reduced-form coefficients and `k` endogenous regressors, the two outcome covariance matrices must be `p x p` and the two exposure covariance matrices `pk x pk`. With multiple endogenous regressors, covariance ordering must follow R's column-major `as.vector(pix)` ordering, including cross-equation covariances. `names(piy)` and `rownames(pix)` must agree exactly. Instruments must occur among those names; `colnames(pix)` identifies endogenous regressors.

`exog` lists control coefficient names and `inst` lists excluded instrument coefficient names. `full=FALSE` estimates endogenous-variable effects; `full=TRUE` additionally estimates intercept and control coefficients. `CUE=TRUE` enables CUE and (for one endogenous regressor) AR confidence sets. `homosk=TRUE` substitutes `Vhomy` and `Vhomx` for the alternative covariance matrices. `alpha` controls AR confidence level.

`nY`, `nX`, `dependent`, and `vcov_info=list(outcome="...", exposure="...")` are optional reporting metadata. Sample sizes, variable names and the provenance of covariance matrices cannot be inferred from their numeric values alone. If unspecified, sample sizes display as `NA` and robustness is labelled unspecified when matrices differ from conventional covariance estimates.

## Example from fitted regressions

```r
library(sandwich)
library(geigen)
source("tsivfuncAR_with_sum_onestep_stars_fixed.R")

inst <- c("z1", "z2")
exog <- c("age", "year")
pred <- c(inst, exog)

lmY <- lm(reformulate(pred, response="outcome"), data=dfY)
lmX <- lm(reformulate(pred, response="schooling"), data=dfX)

piy <- coef(lmY)
pix <- matrix(coef(lmX), ncol=1,
              dimnames=list(names(coef(lmX)), "schooling"))
Vhomy <- vcov(lmY)
Vhety <- vcovCL(lmY, cluster=dfY$state, type="HC1", cadjust=TRUE)
Vhomx <- vcov(lmX)
Vhetx <- vcovCL(lmX, cluster=dfX$state, type="HC1", cadjust=TRUE)

res <- tsiv_sum(
  piy, Vhomy, Vhety, pix, Vhomx, Vhetx, exog, inst,
  full=FALSE, CUE=TRUE, nY=nobs(lmY), nX=nobs(lmX),
  dependent="outcome",
  vcov_info=list(
    outcome="Cluster-robust variance estimate (clustered by state)",
    exposure="Cluster-robust variance estimate (clustered by state)")
)
summary(res)
```

`piy` is a vector: use `names(piy)`, not `colnames(piy)`. Cluster vectors must match the observations retained by each regression.

## Output and limitations

The returned object has the same `tsiv` class and result elements as `tsiv()`. Its first estimator is labelled **One-step**, but its internal method name remains `"ts2sls"`:

```r
summary(res)                        # All available estimates
summary(res, method="ts2sls")       # One-step only
summary(res, method="two_step")
summary(res, method="cue")          # Requires CUE=TRUE
coef(res, method="ts2sls")
```

Four-decimal formatting affects printing only. First-stage F diagnostics and AR confidence sets are available only for one endogenous regressor; AR requires `CUE=TRUE`. An informative Hansen J overidentification test requires more excluded instruments than endogenous regressors. Coefficient z-tests rely on large-sample approximations. Input covariance matrices must be valid for the sampling design; `tsiv_sum()` cannot verify their statistical provenance. The function does not accept a cross-sample covariance block.
