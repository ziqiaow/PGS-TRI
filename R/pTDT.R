#' @title pTDT
#' @description Polygenic transmission disequilibrium test in case-parent trios for direct PGS effect
#'
#' @param pgs_c The PGS values of the affected probands (children). A vector of length N, no missing values are allowed
#' @param pgs_m The PGS values of mothers that corresponds to the children. A vector of same length N, no missing values are allowed
#' @param pgs_f The PGS values of fathers that corresponds to the children. A vector of same length N, no missing values are allowed
#' @param side0=2 Sided of the t test, default is 2-sided.
#'
#'
#' @return Results of pTDT
#'  \item{res_beta}{Results of direct PGS effect, if GxE_int is TRUE, then the result will also include PGSxE interaction effects}
#'
#' @export
#'
#'
ptdt=function(pgs_c, #The PGS values of the affected probands (children). A vector of length N, no missing values are allowed
              pgs_m, #The PGS values of mothers that corresponds to the children. A vector of same length N, no missing values are allowed
              pgs_f, #The PGS values of fathers that corresponds to the children. A vector of same length N, no missing values are allowed
              side0=2 #Sided of the t test, default is 2-sided.
              ){
  pgs_mp=(pgs_m+pgs_f)/2
  ptdt.deviation=(pgs_c-pgs_mp)/sd(pgs_mp)
  # t.test(ptdt.deviation, mu = 0, alternative = "two.sided")
  # t.ptdt=mean(ptdt.deviation)/sd(ptdt.deviation)*sqrt(length(pgs_c))

  #Summarize the results
  res.sum.t=function (parms=beta_hat, cov=var_beta_hat, df0,sided)
  {
    if (sided != 1)
      sided <- 2
    cols <- c("Estimate", "Std.Error", "t.value", "Pvalue")
    n <- length(parms)
    ret <- matrix(data = NA, nrow = n, ncol = 4)
    pnames <- names(parms)
    rownames(ret) <- pnames
    colnames(ret) <- cols
    ret[, 1] <- parms

    if (is.null(pnames))
      pnames <- 1:n
    t.ptdt=parms/cov
    ret[, 2] <- cov
    ret[, 3] <- t.ptdt
    ret[, 4] <- sided * pt(abs(t.ptdt), df0,lower.tail = F)
    ret
  }

  res_beta=res.sum.t(parms=mean(ptdt.deviation),cov=sd(ptdt.deviation)/sqrt(length(pgs_c)),df0=length(pgs_c)-1,sided = side0)

  res=list(res_beta=res_beta)
  return(res)
}

