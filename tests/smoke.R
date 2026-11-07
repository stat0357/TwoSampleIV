library(TwoSampleIV)
stopifnot(is.function(tsiv), is.function(tsiv_sum))
stopifnot("all" %in% eval(formals(getS3method("summary", "tsiv"))$method))
