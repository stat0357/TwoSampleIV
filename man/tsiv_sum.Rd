\name{tsiv_sum}
\alias{tsiv_sum}
\title{Two-sample IV estimation from regression summary statistics}
\usage{tsiv_sum(piy, Vhomy, Vhety, pix, Vhomx, Vhetx,
  exog, inst, full = TRUE, CUE = FALSE, homosk = FALSE,
  alpha = 0.05, nY = NA_integer_, nX = NA_integer_,
  dependent = "Not supplied", vcov_info = NULL)}
\arguments{
\item{piy}{Named reduced-form coefficient vector.}
\item{Vhomy}{Conventional covariance matrix for \code{piy}.}
\item{Vhety}{Alternative covariance matrix for \code{piy}.}
\item{pix}{First-stage coefficient matrix with row names aligned to \code{piy} and columns named for endogenous regressors.}
\item{Vhomx}{Conventional covariance matrix for column-major vectorization of \code{pix}.}
\item{Vhetx}{Alternative covariance matrix for column-major vectorization of \code{pix}.}
\item{exog}{Exogenous control coefficient names.}
\item{inst}{Excluded instrument coefficient names.}
\item{full}{Include intercept and control coefficients.}
\item{CUE}{Compute CUE and available AR confidence sets.}
\item{homosk}{Substitute conventional covariance matrices for alternatives.}
\item{alpha}{Significance level for AR confidence set.}
\item{nY}{Optional outcome sample size.}
\item{nX}{Optional exposure sample size.}
\item{dependent}{Optional dependent-variable label.}
\item{vcov_info}{Optional list with \code{outcome} and \code{exposure} covariance-description strings.}
}
\description{Constructs two-sample IV estimates using coefficients and covariance matrices, without individual-level observations.}
\details{For p reduced-form coefficients and k endogenous regressors, the outcome covariance matrices are p by p and the exposure covariance matrices are pk by pk. For multiple endogenous regressors, exposure covariance matrices use column-major coefficient ordering.}
\value{An object of class \code{tsiv}; the first estimator is printed as One-step.}
\seealso{\code{\link{tsiv}}}
