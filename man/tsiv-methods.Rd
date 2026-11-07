\name{tsiv-methods}
\alias{summary.tsiv}
\alias{print.tsiv}
\alias{print.summary.tsiv}
\alias{print.summary.tsiv.all}
\alias{coef.tsiv}
\title{Methods for two-sample IV results}
\usage{\method{summary}{tsiv}(object, method = c("all", "ts2sls", "two_step", "cue"), ...)
\method{print}{tsiv}(x, ...)
\method{coef}{tsiv}(object, method = c("ts2sls", "two_step", "cue"), ...)}
\arguments{
\item{object}{A \code{tsiv} object.}
\item{x}{A \code{tsiv} or summary object.}
\item{method}{Which estimator to display or extract.}
\item{...}{Additional arguments.}
}
\description{Summarize, print and extract coefficients from two-sample IV estimation objects.}
\details{The default summary displays all available estimators. Printed numbers use four decimal places and significance stars.}
