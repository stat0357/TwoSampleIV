#' @export
Fcalc <- function(pix,Vhomy,Vhomx,Vhetx,kz) {

  tau <- 0.1
  F <- t(pix)%*%solve(Vhomx)%*%pix/kz
  F_rob <- t(pix)%*%solve(Vhetx)%*%pix/kz

  iVhomy <- solve(Vhomy)

  if (dim(iVhomy)[1]>1) {
    eigVhomy <- eigen(Vhomy)
    iVhomh <- eigVhomy$vectors%*%diag(1/sqrt(eigVhomy$values))%*%t(eigVhomy$vectors)
  } else {
    iVhomh <- sqrt(iVhomy)
  }

  Vhethom <- iVhomh%*%Vhetx%*%iVhomh
  tr_Vhethom <- sum(diag(Vhethom))

  F_geffhet <- (t(pix)%*%iVhomy%*%pix)/tr_Vhethom

  evalhet <- eigen(Vhethom)$values
  maxevalhet <- max(evalhet)
  minevalhet <- min(evalhet)

  B <- abs(1-2*minevalhet/tr_Vhethom)
  d <- B/tau
  k_eff <- tr_Vhethom^2*(1+2*d)/(sum(diag(t(Vhethom)%*%Vhethom))+2*d*tr_Vhethom*maxevalhet)
  cv_het <- qchisq(p=0.95,df=k_eff,ncp=k_eff*d)/k_eff

  resF <- data.frame(
    F_stat = as.numeric(F),
    rob_F_stat = as.numeric(F_rob),
    F_geff = as.numeric(F_geffhet),
    cv_F_geff = as.numeric(cv_het)
  )

  return(resF)
}

#' @export
tscue <- function(piy,Vhety,pix,Vhetx,bstart,alpha) {

  object <- function(b) {
    krbi <- kronecker(b,diag(kz))
    t((piy-pix%*%b))%*%solve(Vhety+t(krbi)%*%Vhetx%*%krbi)%*%(piy-pix%*%b)
  }

  kx <- dim(pix)[2]
  kz <- dim(pix)[1]
  dfs <- kz-kx

  if (kx<2) {
    kz_F_rob <- t(pix)%*%solve(Vhetx)%*%pix
    cvc <- qchisq(1-alpha,kz)
    m0 <- matrix(0,nrow=kz,ncol=kz)
    v0 <- matrix(0,nrow=kz,ncol=1)
    D0 <- rbind(cbind(Vhety,piy),cbind(t(piy),cvc))
    D1 <- -(rbind(cbind(m0,pix),cbind(t(pix),0)))
    D2 <- rbind(cbind(Vhetx,v0),cbind(t(v0),0))
    
    m0 <- matrix(0,nrow=(kz+1),ncol=(kz+1))
    A <- rbind(cbind(m0,diag(kz+1)),cbind(D0,D1))
    B <- rbind(cbind(diag(kz+1),m0),cbind(m0,-D2))
    
    roots <- geigen(A,B,symmetric=FALSE)
    roots <- roots$values
    roots <- roots[is.finite(roots)]
    roots <- sort(Re(roots[abs(Im(roots))<1e-9]))
    nroots = length(roots)
   
    lower <- numeric(0)
    upper <- numeric(0)
    
    if (nroots>0) {
      if (kz_F_rob>cvc) {
        for (j in 1:(nroots/2)) {
          lower <- c(lower,roots[2*j-1])
          upper <- c(upper,roots[2*j])
          } 
        } else {
          lower <- -Inf  
          upper <- Inf
          for (j in seq_len(nroots/2 - 1L)) {
            lower <- c(lower,roots[2*j])
            upper <- c(upper,roots[(2*j)+1])
          }
        }
      } else {
      if (kz_F_rob<cvc) {
        lower <- -Inf
        upper <- Inf
        }
      }

   
    AR_set <- data.frame(
      lower,upper
    )
   
    nset <- length(lower)
    resval <- 99999
    if (nset>0) {
      if (is.finite(lower[1])) {
        for (j in 1:nset) {
          bstart <- (upper[j]+lower[j])/2
          rescuej <- optim(par=bstart, fn = object, method = "Brent", lower=lower[j], upper=upper[j])
          if (rescuej$value<resval) {
            rescue <- rescuej
            resval <- rescuej$value
          }
        }
      } else {
        rescuej <- optim(par=bstart, fn = object, method = "Brent", lower=-999, upper=upper[1])
        if (rescuej$value<resval) {
          rescue <- rescuej
          resval <- rescuej$value
        }
        rescuej <- optim(par=bstart, fn = object, method = "Brent", lower=lower[nset], upper=999)
        if (rescuej$value<resval) {
          rescue <- rescuej
          resval <- rescuej$value
        }
        if (nset>2) {
          for (j in 2:(nset-1)) {
            bstart <- (upper[j]+lower[j])/2
            rescuej <- optim(par=bstart, fn = object, method = "Brent", lower=lower[j], upper=upper[j])
            if (rescuej$value<resval) {
              rescue <- rescuej
              resval <- rescuej$value
            }
          }
        }
      }
    }
        
  } else {
    rescue <- optim(par=bstart, fn = object)
  }

  if (!exists("rescue", inherits = FALSE)) rescue <- optim(par=bstart, fn=object, method="BFGS")
  bcue <- rescue$par
  Han_cue <- rescue$value
  pval_Han <- pchisq(Han_cue,df=dfs,lower.tail=FALSE)

  krbi <- kronecker(bcue,diag(kz))
  se_cue <- sqrt(diag(solve((t(pix)%*%solve(Vhety+t(krbi)%*%Vhetx%*%krbi)%*%pix))))

  rescue <- data.frame(
    Variable = colnames(pix),
    bcue = bcue,
    se_rob_cue = se_cue
  )

  res_Han_cue <- data.frame(
    Han_J_cue = as.numeric(Han_cue),
    p_val = as.numeric(pval_Han)
  )

  if (kx<2) {
    return(list(rescue,res_Han_cue,AR_set))
  } else {
    return(list(rescue,res_Han_cue))
  }
}



#' @export
tsest <- function(piy,Vhomy,Vhety,pix,Vhomx,Vhetx,indX) {

  kx <- dim(pix)[2]
  kz <- dim(pix)[1]

  dfs <- kz-kx

  W1 <- solve(Vhomy)

  c <- W1%*%pix
  den <- solve(t(c)%*%pix)
  num <- t(c)%*%piy
  bt1 <- den%*%num
  b1 <- bt1[indX]
  krbc <- kronecker(b1,c)
  krbi <- kronecker(b1,diag(kz))
  W2hom <- solve(Vhomy+t(krbi)%*%Vhomx%*%krbi)

  Varb1hom <- den%*%(t(c)%*%Vhomy%*%c+t(krbc)%*%Vhomx%*%krbc)%*%den
  Varb1het <- den%*%(t(c)%*%Vhety%*%c+t(krbc)%*%Vhetx%*%krbc)%*%den

  W2het <- solve(Vhety+t(krbi)%*%Vhetx%*%krbi)
  c <- W2het%*%pix
  den <- solve(t(c)%*%pix)
  num <- t(c)%*%piy
  bt2 <- den%*%num
  Varb2 <- den

  Han <- t(piy-pix%*%bt2)%*%W2het%*%(piy-pix%*%bt2)
  pval_Han <- pchisq(Han,df=dfs,lower.tail=FALSE)


  res_ts2sls <- data.frame(
    Variable = colnames(pix),
    ts2sls = bt1,
    se_hom = sqrt(diag(Varb1hom)),
    se_rob = sqrt(diag(Varb1het))
  )

  res_twostep <- data.frame(
    Variable = colnames(pix),
    twostep = bt2,
    se_rob_2s = sqrt(diag(Varb2))
    )

  resultsHan <- data.frame(
    Hansen_J = Han,
    p_val = pval_Han
  )


  output <- list(res_ts2sls,res_twostep,resultsHan)

  return(output)

}

#' @export
tsiv <- function(dfY,dfX,Y,X,exog,inst,Vhety_str,Vhetx_str,full=TRUE,CUE=FALSE,homosk=FALSE,alpha=0.05) {

  if (length(exog)<2) {
    if ((exog != "")) {
      pred <- c(inst,exog)
    } else {
      pred <- inst
    }
  } else {
    pred <- c(inst,exog)
  }



  formulaY <- as.formula(paste(Y, "~", paste(pred, collapse = " + ")))

  formulaX <- as.formula(paste("cbind(",paste(X, collapse = ", "),
                               ") ~ ", paste(pred, collapse = " + "), sep = ""))


  lmY <- lm(formulaY, data = dfY)
  piy <- lmY$coefficients
  Vhomy <- vcov(lmY)
  if (homosk==TRUE) {
    Vhety = Vhomy
  } else {
    Vhety <- eval(parse(text=Vhety_str))
  }

  nY <- nobs(lmY)

  lmX <- lm(formulaX, data = dfX)
  pix <- as.matrix(lmX$coefficients)
  colnames(pix) = X
  Vhomx <- vcov(lmX)
  if (homosk==TRUE) {
    Vhetx = Vhomx
  } else {
    Vhetx <- eval(parse(text=Vhetx_str))
  }

  nX <- nobs(lmX)

  resdat <- data.frame(
    Dep_Var = Y,
    nY = nY,
    nX = nX
  )

  if (full==FALSE || CUE==TRUE || length(X)<2) {

    indZ <- which(names(piy)==inst[1])
    if (length(inst)>1){
      for (j in 2:length(inst)){
        indZ <- c(indZ,which(names(piy)==inst[j]))
      }
    }

    piysel <- piy[indZ]
    Vhomysel <- Vhomy[indZ,indZ,drop=FALSE]
    Vhetysel <- Vhety[indZ,indZ,drop=FALSE]

    k_var <- dim(pix)[1]
    k_x <- dim(pix)[2]

    selm <- indZ
    if (k_x>1) {
      for (j in 2:k_x){
        selm <- c(selm,indZ+(j-1)*k_var)
      }
    }
    pixsel <- as.matrix(pix[indZ,])
    colnames(pixsel) <- X
    Vhomxsel <- Vhomx[selm,selm,drop=FALSE]
    Vhetxsel <- Vhetx[selm,selm,drop=FALSE]

    if (length(X)<2) {
      resF <- Fcalc(pixsel,Vhomysel,Vhomxsel,Vhetxsel,dim(pixsel)[1])
    }

    if (CUE==TRUE) {
      rescue <- tscue(piysel,Vhetysel,pixsel,Vhetxsel,rep(0,k_x),alpha)
    }

    if (full==FALSE) {
      res <- tsest(piysel,Vhomysel,Vhetysel,pixsel,Vhomxsel,Vhetxsel,1:k_x)
    }
  }
  if (full==TRUE){

    npred <- dim(pix)[1]
    Inp <- diag(npred)

    ss = 0

    if (length(exog)<2) {
      if (exog != ""){
        ss <- which(pred==exog[1])
      }
    } else {
      ss <- which(pred==exog[1])
    }

    if (ss>0){
      if (ss==1) {
        pix = cbind(Inp[,1:(length(exog)+1)],pix)
        colnames(pix) = c("Const", exog, X)
      } else {
        pix =cbind(Inp[,1],pix,Inp[,(npred-length(exog)+1):npred])
        colnames(pix) = c("Const", X, exog)
      }
    } else {
      pix =cbind(Inp[,1],pix)
      colnames(pix) = c("Const", X)
    }

    indX <- which(colnames(pix)==X[1])
    if (length(X)>1) {
      for (j in 2:length(X)){
        indX <- c(indX,which(colnames(pix)==X[j]))
      }
    }

    res <- tsest(piy,Vhomy,Vhety,pix,Vhomx,Vhetx,indX)
  }

  output <- res


  if (length(X) < 2) {
    output <- append(output,list(resF),after=0)
  }

  output <- append(output,list(resdat),after=0)

  if (CUE==TRUE) {
    output <- append(output,rescue)
  }

  names(output) <- c(
    "sample", if (length(X) < 2) "first_stage",
    "ts2sls", "two_step", "hansen",
    if (CUE) c("cue", "cue_hansen", if (length(X) < 2) "ar_confint")
  )
  class(output) <- c("tsiv", "list")
  attr(output, "call") <- match.call()
  attr(output, "dependent") <- Y
  attr(output, "endogenous") <- X
  attr(output, "instruments") <- inst
  attr(output, "alpha") <- alpha
  attr(output, "homosk_covariance") <- isTRUE(homosk) || (
    .tsiv_is_vcov(Vhety_str) && .tsiv_is_vcov(Vhetx_str))
  attr(output, "vcov_info") <- list(
    outcome = .tsiv_vcov_label(Vhety_str, homosk),
    exposure = .tsiv_vcov_label(Vhetx_str, homosk))
  return(output)
}

# Identify the conventional homoskedastic covariance specification.
.tsiv_is_vcov <- function(expr) {
  compact <- gsub("[[:space:]]+", "", paste(expr, collapse = ""))
  compact %in% c("vcov", "Vhomy", "Vhomx") ||
    grepl("^(stats::)?vcov\\(", compact)
}

# Hide robust first-stage diagnostics only when both samples use homoskedastic covariance.
.tsiv_first_stage <- function(dat, homosk_covariance = FALSE) {
  if (is.null(dat) || !isTRUE(homosk_covariance)) return(dat)
  dat[, setdiff(names(dat), c("rob_F_stat", "F_geff", "cv_F_geff")), drop = FALSE]
}

# Describe the covariance estimators actually used in the two samples.
.tsiv_vcov_label <- function(expr, homosk = FALSE) {
  if (isTRUE(homosk)) return("Variance estimate valid under homoskedasticity")
  expr <- paste(expr, collapse = "")
  compact <- gsub("[[:space:]]+", "", expr)
  if (grepl("^(stats::)?vcov\\(", compact) || compact %in% c("vcov", "Vhomy", "Vhomx"))
    return("Variance estimate valid under homoskedasticity")
  if (grepl("(^|::)vcovCL\\(", compact)) {
    cluster <- sub(".*cluster[[:space:]]*=[[:space:]]*([^,)]+).*", "\\1", expr)
    if (identical(cluster, expr)) return("Cluster-robust variance estimate")
    cluster <- trimws(cluster)
    cluster <- sub("^(dfY|dfX)\\$", "", cluster)
    return(paste0("Cluster-robust variance estimate (clustered by ", cluster, ")"))
  }
  if (grepl("(^|::)vcovHC\\(", compact)) return("Heteroskedasticity-robust variance estimate")
  paste0("Robust/custom covariance estimator: ", expr)
}

.tsiv_print_vcov <- function(info) {
  if (is.null(info)) return(invisible(NULL))
  cat("Outcome sample: ", info$outcome, "\n", sep = "")
  cat("Exposure sample: ", info$exposure, "\n", sep = "")
  invisible(NULL)
}

# Summary and printing methods for two-sample IV estimates.
# The z tests use a large-sample normal approximation, not lm's residual df.
.tsiv_coef_table <- function(dat, estimate, se) {
  if (is.null(dat) || !nrow(dat)) return(NULL)
  b <- as.numeric(dat[[estimate]])
  s <- as.numeric(dat[[se]])
  z <- b / s
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  tab <- cbind("Estimate" = b, "Std. Error" = s,
               "z value" = z, "Pr(>|z|)" = p)
  rownames(tab) <- as.character(dat$Variable)
  tab
}

#' @export
summary.tsiv <- function(object, method = c("all", "ts2sls", "two_step", "cue"), ...) {
  method <- match.arg(method)
  if (method == "all") {
    methods <- c("ts2sls", "two_step", if (!is.null(object$cue)) "cue")
    ans <- list(call=attr(object,"call"), dependent=attr(object,"dependent"),
                sample=object$sample, first_stage=.tsiv_first_stage(object$first_stage, attr(object,"homosk_covariance")),
                vcov_info=attr(object,"vcov_info"),
                from_summary_stats=isTRUE(attr(object,"from_summary_stats")),
                tables=setNames(lapply(methods, function(m) summary(object, method=m)), methods),
                ar_confint=object$ar_confint, alpha=attr(object,"alpha"))
    class(ans) <- "summary.tsiv.all"
    return(ans)
  }
  if (method == "cue" && is.null(object$cue))
    stop("CUE estimates were not requested; rerun tsiv(..., CUE = TRUE).")
  spec <- switch(method,
    ts2sls = list(data = object$ts2sls, b = "ts2sls", se = "se_rob"),
    two_step = list(data = object$two_step, b = "twostep", se = "se_rob_2s"),
    cue = list(data = object$cue, b = "bcue", se = "se_rob_cue"))
  ans <- list(call = attr(object, "call"), method = method,
              coefficients = .tsiv_coef_table(spec$data, spec$b, spec$se),
              sample = object$sample, first_stage = .tsiv_first_stage(object$first_stage, attr(object,"homosk_covariance")),
              vcov_info = attr(object,"vcov_info"),
              from_summary_stats = isTRUE(attr(object,"from_summary_stats")),
              hansen = if (method == "cue") object$cue_hansen else object$hansen,
              ar_confint = object$ar_confint,
              dependent = attr(object, "dependent"),
              endogenous = attr(object, "endogenous"))
  class(ans) <- "summary.tsiv"
  ans
}

# Print numerical estimates to exactly four decimal places without rounding stored results.
.tsiv_signif <- function(p) {
  out <- rep("", length(p))
  ok <- !is.na(p)
  out[ok & p <= 0.1] <- "."
  out[ok & p <= 0.05] <- "*"
  out[ok & p <= 0.01] <- "**"
  out[ok & p <= 0.001] <- "***"
  out
}

.tsiv_print_numeric <- function(dat, stars = FALSE) {
  if (is.null(dat)) return(invisible(NULL))
  if (!nrow(dat)) return(invisible(NULL))
  vals <- as.data.frame(dat, check.names = FALSE)
  for (nm in names(vals)) {
    if (is.numeric(vals[[nm]])) {
      vals[[nm]] <- ifelse(is.na(vals[[nm]]), "NA",
                           sprintf("%.4f", vals[[nm]]))
    }
  }
  if (stars) {
    pcol <- if ("Pr(>|z|)" %in% colnames(dat)) "Pr(>|z|)" else
      if ("p_val" %in% names(dat)) "p_val" else NULL
    if (!is.null(pcol)) vals[[" "]] <- .tsiv_signif(if (is.matrix(dat)) dat[, pcol] else dat[[pcol]])
  }
  print(vals, row.names = !is.null(rownames(dat)) &&
        !identical(rownames(dat), as.character(seq_len(nrow(dat)))),
        right = TRUE, quote = FALSE)
  invisible(NULL)
}

.tsiv_print_hansen <- function(dat) {
  .tsiv_print_numeric(dat, stars = TRUE)
}

#' @export
print.summary.tsiv <- function(x, digits = 4L,
                               signif.stars = getOption("show.signif.stars"), ...) {
  cat("\nTwo-sample instrumental variables estimation\n")
  cat("Dependent variable:", x$dependent, "\n")
  cat("Estimator:", switch(x$method, ts2sls = if (isTRUE(x$from_summary_stats)) "One-step" else "TS2SLS",
                             two_step = "Two-step efficient IV", cue = "CUE"), "\n\n")
  .tsiv_print_vcov(x$vcov_info)
  if (!is.null(x$vcov_info) &&
      all(grepl("valid under homoskedasticity", unlist(x$vcov_info), fixed = TRUE)))
    cat("Test statistics are valid under homoskedasticity only.\n")
  if (!is.null(x$first_stage)) {
    cat("\nFirst-stage diagnostics:\n")
    .tsiv_print_numeric(x$first_stage)
  }
  cat("\nCoefficients:\n")
  .tsiv_print_numeric(x$coefficients, stars = signif.stars)
  cat("\nObservations: outcome sample =", x$sample$nY,
      "; exposure sample =", x$sample$nX, "\n\n")
  if (!is.null(x$hansen)) {
    cat("\nOveridentification test (Hansen J):\n")
    .tsiv_print_hansen(x$hansen)
  }
  if (x$method == "cue" && !is.null(x$ar_confint)) {
    cat("\nAnderson-Rubin confidence set:\n")
    .tsiv_print_numeric(x$ar_confint)
  }
  if (isTRUE(signif.stars)) {
    cat("---\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  }
  invisible(x)
}

#' @export
print.tsiv <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' @export
coef.tsiv <- function(object, method = c("ts2sls", "two_step", "cue"), ...) {
  summary(object, method = match.arg(method))$coefficients[, "Estimate"]
}


#' @export
print.summary.tsiv.all <- function(x, digits=4L,
                                   signif.stars=getOption("show.signif.stars"), ...) {
  cat("\nTwo-sample instrumental variables estimation\n")
  cat("Dependent variable:", x$dependent, "\n")
  cat("Observations: outcome sample =", x$sample$nY,
      "; exposure sample =", x$sample$nX, "\n\n")
  .tsiv_print_vcov(x$vcov_info)
  if (!is.null(x$vcov_info) &&
      all(grepl("valid under homoskedasticity", unlist(x$vcov_info), fixed = TRUE)))
    cat("Test statistics are valid under homoskedasticity only.\n")
  if (!is.null(x$first_stage)) {
    cat("\nFirst-stage diagnostics:\n")
    .tsiv_print_numeric(x$first_stage)
  }
  for (m in names(x$tables)) {
    tab <- x$tables[[m]]
    cat("\n", switch(m, ts2sls=if (isTRUE(x$from_summary_stats)) "One-step" else "TS2SLS", two_step="Two-step efficient IV",
                       cue="CUE"), " coefficients:\n", sep="")
    .tsiv_print_numeric(tab$coefficients, stars=signif.stars)
    if (!is.null(tab$hansen)) {
      cat("\nOveridentification test (Hansen J):\n")
      .tsiv_print_hansen(tab$hansen)
    }
  }
  if (!is.null(x$ar_confint)) {
    cat("\nAnderson-Rubin ", 100*(1-x$alpha), "% confidence set:\n", sep="")
    if (nrow(x$ar_confint)) .tsiv_print_numeric(x$ar_confint)
    else cat("Empty confidence set\n")
  }
  if (isTRUE(signif.stars)) {
    cat("---\n")
    cat("Signif. codes:  0 \'***\' 0.001 \'**\' 0.01 \'*\' 0.05 \'.\' 0.1 \' \' 1\n")
  }
  invisible(x)
}


#' Two-sample IV estimation from first-stage and reduced-form summary statistics.
#'
#' piy: reduced-form coefficient vector (including intercept and controls).
#' pix: first-stage coefficient matrix, rows aligned with piy, columns with X.
#' Covariance matrices must follow the corresponding coefficient ordering;
#' for multiple endogenous regressors, Vhomx/Vhetx use column-major ordering.
#' nY/nX and dependent are optional metadata, unavailable from six moments alone.
#' @export
tsiv_sum <- function(piy, Vhomy, Vhety, pix, Vhomx, Vhetx,
                     exog, inst, full = TRUE, CUE = FALSE,
                     homosk = FALSE, alpha = 0.05,
                     nY = NA_integer_, nX = NA_integer_,
                     dependent = 'Not supplied',
                     vcov_info = NULL) {
  piy_names <- names(piy)
  if (is.matrix(piy)) {
    if (ncol(piy) != 1L) stop('piy must be a coefficient vector or one-column matrix.')
    if (is.null(piy_names)) piy_names <- rownames(piy)
  }
  piy <- as.numeric(piy)
  if (is.null(piy_names) || anyNA(piy_names) || anyDuplicated(piy_names))
    stop('piy must have unique coefficient names, including the instruments.')
  names(piy) <- piy_names
  pix <- as.matrix(pix)
  if (!is.numeric(pix) || !is.numeric(piy)) stop('piy and pix must be numeric.')
  if (nrow(pix) != length(piy)) stop('pix must have one row for each element of piy.')
  if (is.null(rownames(pix))) rownames(pix) <- names(piy)
  if (!identical(rownames(pix), names(piy)))
    stop('rownames(pix) must match names(piy) in the same order.')
  X <- colnames(pix)
  if (is.null(X) || anyDuplicated(X)) stop('pix needs unique column names for the endogenous variables.')
  k_x <- ncol(pix)
  ncoef <- length(piy)
  kz <- length(inst)
  if (!kz || anyDuplicated(inst) || !all(inst %in% names(piy)))
    stop('inst must contain distinct coefficient names present in piy.')
  if (any(inst %in% X)) stop('Instrument names and endogenous regressor names must differ.')
  if (kz < k_x) stop('Number of excluded instruments must be at least the number of endogenous regressors.')
  if (length(exog) == 1L && identical(exog, '')) exog <- character(0)
  if (length(exog) && !all(exog %in% names(piy)))
    stop('All exog names must occur in piy.')
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1)
    stop('alpha must be between 0 and 1.')
  validate_cov <- function(V, size, label) {
    V <- as.matrix(V)
    if (!is.numeric(V) || !identical(dim(V), c(size, size)))
      stop(label, ' must be a ', size, ' by ', size, ' numeric covariance matrix.')
    V
  }
  Vhomy <- validate_cov(Vhomy, ncoef, 'Vhomy')
  Vhety <- validate_cov(Vhety, ncoef, 'Vhety')
  Vhomx <- validate_cov(Vhomx, ncoef * k_x, 'Vhomx')
  Vhetx <- validate_cov(Vhetx, ncoef * k_x, 'Vhetx')
  if (isTRUE(homosk)) {
    Vhety <- Vhomy
    Vhetx <- Vhomx
  }
  is_hom_y <- isTRUE(all.equal(unname(Vhety), unname(Vhomy), tolerance = 1e-12))
  is_hom_x <- isTRUE(all.equal(unname(Vhetx), unname(Vhomx), tolerance = 1e-12))
  homosk_covariance <- is_hom_y && is_hom_x
  if (is.null(vcov_info)) {
    vcov_info <- list(
      outcome = if (is_hom_y) 'Variance estimate valid under homoskedasticity'
                else 'User-supplied variance estimate (robustness unspecified)',
      exposure = if (is_hom_x) 'Variance estimate valid under homoskedasticity'
                 else 'User-supplied variance estimate (robustness unspecified)')
  }

  indZ <- match(inst, names(piy))
  piysel <- piy[indZ]
  Vhomysel <- Vhomy[indZ, indZ, drop = FALSE]
  Vhetysel <- Vhety[indZ, indZ, drop = FALSE]
  selm <- unlist(lapply(seq_len(k_x), function(j) indZ + (j - 1L) * ncoef), use.names = FALSE)
  pixsel <- pix[indZ, , drop = FALSE]
  Vhomxsel <- Vhomx[selm, selm, drop = FALSE]
  Vhetxsel <- Vhetx[selm, selm, drop = FALSE]

  if (k_x == 1L) {
    resF <- Fcalc(pixsel, Vhomysel, Vhomxsel, Vhetxsel, kz)
  }
  if (isTRUE(CUE)) {
    rescue <- tscue(piysel, Vhetysel, pixsel, Vhetxsel, rep(0, k_x), alpha)
  }
  if (!isTRUE(full)) {
    res <- tsest(piysel, Vhomysel, Vhetysel, pixsel, Vhomxsel, Vhetxsel, seq_len(k_x))
  } else {
    # Match the full-model construction used by tsiv().
    Inp <- diag(ncoef)
    if (length(exog)) {
      ss <- match(exog[1], names(piy))
      if (ss == 1L) {
        pix_full <- cbind(Inp[, seq_len(length(exog) + 1L), drop = FALSE], pix)
        colnames(pix_full) <- c('Const', exog, X)
      } else {
        pix_full <- cbind(Inp[, 1, drop = FALSE], pix,
                          Inp[, (ncoef - length(exog) + 1L):ncoef, drop = FALSE])
        colnames(pix_full) <- c('Const', X, exog)
      }
    } else {
      pix_full <- cbind(Inp[, 1, drop = FALSE], pix)
      colnames(pix_full) <- c('Const', X)
    }
    indX <- match(X, colnames(pix_full))
    res <- tsest(piy, Vhomy, Vhety, pix_full, Vhomx, Vhetx, indX)
  }
  output <- res
  if (k_x == 1L) output <- append(output, list(resF), after = 0)
  output <- append(output, list(data.frame(Dep_Var = dependent, nY = nY, nX = nX)), after = 0)
  if (isTRUE(CUE)) output <- append(output, rescue)
  names(output) <- c('sample', if (k_x == 1L) 'first_stage',
                     'ts2sls', 'two_step', 'hansen',
                     if (CUE) c('cue', 'cue_hansen', if (k_x == 1L) 'ar_confint'))
  class(output) <- c('tsiv', 'list')
  attr(output, 'call') <- match.call()
  attr(output, 'dependent') <- dependent
  attr(output, 'endogenous') <- X
  attr(output, 'instruments') <- inst
  attr(output, 'alpha') <- alpha
  attr(output, 'homosk_covariance') <- homosk_covariance
  attr(output, 'vcov_info') <- vcov_info
  attr(output, 'from_summary_stats') <- TRUE
  output
}
