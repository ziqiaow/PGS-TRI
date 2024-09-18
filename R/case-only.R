#' @title Case only analysis for PGSxE interactions
#' @description Polygenic PGSxE interactions using the case-only method in unrelated individuals with the disease/condition
#'
#' @param pgs_c The PGS values of the affected individuals. A vector of length N, no missing values are allowed
#' @param envir The environmental variables of interest for interaction effect. A vector of length N for one environmental variable or a data frame/data matrix of NxP for P environmental variables are allowed.
#' @param mean_prs The mean PGS value from an external dataset (or could use the same dataset, i.e., mean(pgs_c) as input). This step only influence the PGS main effect result.
#'
#'
#' @return Results of case-only analysis for PGSxE
#'  \item{res}{Results of PGSxE interaction effects}
#'
#' @export
#'
#'
function_caseonly <- function(prs_c,
                              envir,
                              mean_prs
                              ){


  summary.caseonly=function (parms, sd, sided = 2)
  {
    if (sided != 1)
      sided <- 2
    cols <- c("Estimate", "Std.Error", "Z.value", "Pvalue")
    n <- length(parms)
    ret <- matrix(data = NA, nrow = n, ncol = 4)
    pnames <- c("prs",paste0("prs:",names(parms)[-1]))
    rownames(ret) <- pnames
    colnames(ret) <- cols
    ret[, 1] <- parms
    if (is.null(pnames))
      pnames <- 1:n
    cov <- sd
    ret[, 2] <- cov
    ret[, 3] <- parms/cov
    ret[, 4] <- sided * pnorm(abs(ret[, 3]), lower.tail = FALSE)
    ret
  }


  fit_caseonly <- lm(prs_c ~ .,data=data.frame(envir))
  beta_ge_ini=fit_caseonly$coefficients[-1]/(sd(fit_caseonly$residuals))^2


  beta=fit_caseonly$coefficients
  beta_int=fit_caseonly$coefficients[-1]/(sd(fit_caseonly$residuals))^2
  sd_int=summary(fit_caseonly)$coef[-1,2]/(sd(fit_caseonly$residuals))^2

  beta_prs=(fit_caseonly$coefficients[1]-mean_prs)/(sd(fit_caseonly$residuals))^2
  sd_prs1= sqrt(( (summary(fit_caseonly)$coef[1,2])^2 +(sd(fit_caseonly$residuals))^2/1000000) / (sd(fit_caseonly$residuals))^4)

  res=summary.caseonly(parms = c(beta_prs,beta_int),sd=c(sd_prs1,sd_int))


  return(res)
}




