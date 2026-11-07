\name{tsest}\alias{tsest}\title{Two-sample IV estimation from selected coefficients}
\usage{tsest(piy, Vhomy, Vhety, pix, Vhomx, Vhetx, indX)}
\arguments{\item{piy}{Reduced-form coefficient vector.}\item{Vhomy}{Conventional reduced-form covariance.}\item{Vhety}{Alternative reduced-form covariance.}\item{pix}{First-stage coefficient matrix.}\item{Vhomx}{Conventional first-stage covariance.}\item{Vhetx}{Alternative first-stage covariance.}\item{indX}{Indices of endogenous coefficients.}}
\description{Low-level one-step and two-step two-sample IV estimation.}\value{A list of estimates and Hansen J results.}
