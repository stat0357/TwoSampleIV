\name{tsiv}
\alias{tsiv}
\title{Two-sample IV estimation from individual-level data}
\usage{tsiv(dfY, dfX, Y, X, exog, inst, Vhety_str, Vhetx_str,
  full = TRUE, CUE = FALSE, homosk = FALSE, alpha = 0.05)}
\arguments{
\item{dfY}{Outcome-sample data frame.}
\item{dfX}{Exposure-sample data frame.}
\item{Y}{Outcome variable name.}
\item{X}{Character vector of endogenous exposure names.}
\item{exog}{Exogenous control names; use \code{""} for none.}
\item{inst}{Excluded instrument names.}
\item{Vhety_str}{String expression for outcome covariance matrix, evaluated with \code{lmY} available.}
\item{Vhetx_str}{String expression for exposure covariance matrix, evaluated with \code{lmX} available.}
\item{full}{Include intercept and controls in the IV coefficient vector.}
\item{CUE}{Compute continuously updated estimates and, for one endogenous regressor, an AR confidence set.}
\item{homosk}{Use conventional homoskedastic covariance estimates.}
\item{alpha}{Significance level for AR confidence set.}
}
\description{Fits reduced-form and first-stage regressions on separate samples and constructs two-sample IV estimates.}
\value{An object of class \code{tsiv} with sample information, one-step/TS2SLS and two-step estimates, Hansen J results, and optional first-stage diagnostics, CUE estimates and AR confidence sets.}
\seealso{\code{\link{tsiv_sum}}}
