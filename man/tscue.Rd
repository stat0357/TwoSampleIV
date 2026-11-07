\name{tscue}\alias{tscue}\title{Continuously updated two-sample IV estimation}
\usage{tscue(piy, Vhety, pix, Vhetx, bstart, alpha)}
\arguments{\item{piy}{Selected reduced-form coefficients.}\item{Vhety}{Reduced-form covariance matrix.}\item{pix}{Selected first-stage coefficients.}\item{Vhetx}{First-stage covariance matrix.}\item{bstart}{Starting parameter values.}\item{alpha}{Significance level.}}
\description{Low-level continuously updated estimation with optional AR confidence sets.}\value{A list of CUE estimates, Hansen test results and, for one endogenous regressor, AR confidence intervals.}
