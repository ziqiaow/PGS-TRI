#' @title Simulation
#' @description Simulate case-parent trio data and environmental variables in the population
#'
#' @param n_fam The total number of trios in the population, recommend values larger than 100000
#' @param cor_e_prs Whether there are correlations between PGS and E
#' @param cor_strat  Level of stratification biases for PGS main effect
#' @param rho2 Level of stratification biases for PGSxE interactions
#' @param alpha_fam Intercept term to control for population disease prevalence
#' @param betaG_normPRS True log RR for direct PGS main effect
#' @param betaE_bin True log RR for E1, a binary variable
#' @param betaE_norm True log RR for E2, a continuous variable
#' @param betaGE_normPRS_bin True log RR for PGSxE1
#' @param betaGE_normPRS_norm True log RR for PGSxE2
#' @param envir Whether there are environmental variables and PGSxE interactions in the true underlying disease model
#'
#'
#' @return A list of data from simulation
#'  \item{D_sim}{Disease indicator of the offspring}
#'  \item{E_sim}{A data frame of simulated environmental variables}
#'  \item{pgs_fam}{PGS values of trios}
#'
#'
#' @importFrom MASS mvrnorm
#' @export
#'
sim_prospective_population=function(n_fam=100000,cor_e_prs=TRUE,cor_strat=0.25,rho2=0.25,alpha_fam=-5.5,betaG_normPRS=0.4,betaE_bin=0.2, betaE_norm=-0.6,betaGE_normPRS_bin=0.2, betaGE_normPRS_norm=-0.4,envir=TRUE){
  # random generation for sigma_i^2 of the PGS values for the family, assume a mixture of gamma distributions
  rmixgamma <- function(n, pi, alpha, beta) {
    k <- sample.int(length(pi), n, replace = TRUE, prob = pi)
    rgamma(n, alpha[k], beta[k])
  }

  pi <- c(0.6,0.3,0.1)
  alpha <- c(6, 15,60)
  beta <- c(15, 30,150)

  var_sim=rmixgamma(n_fam, pi, alpha, beta)

  if(envir == TRUE){
    if(cor_e_prs == TRUE){
      rho=runif(n_fam,0,0.5)} else {
        rho=rep(0,n_fam)
      }
    mu1=rnorm(n_fam)
    e1=rnorm(n_fam)
    e2=rnorm(n_fam)
    mu_alpha=array(0,c(n_fam,3))
    for(i in 1:n_fam){
      sigma_sim=matrix(c(1,rho[i],rho[i],rho[i],1,0,rho[i],0,1),3,3)
      sigma_sim=chol(sigma_sim)
      mu_alpha[i,]=c(mu1[i],e1[i],e2[i]) %*% t(sigma_sim)
    }


    mu_var=cbind(mu_alpha[,1],var_sim)
    myfunctionpgs=function(x){
      sigma_tmp=matrix(c(1,0.5,0.5,0.5,1,0,0.5,0,1)*x[2],3,3)
      MASS::mvrnorm(n=1,mu=rep(x[1],3),Sigma=sigma_tmp)

    }
    pgs_fam=t(apply(mu_var,1,myfunctionpgs))

    colnames(pgs_fam)=c("pgs_c","pgs_m","pgs_f")
    prs <- as.numeric(as.vector(pgs_fam[,1]))

    E_sim_bin <- ifelse(mu_alpha[,2]>0,0,1)
    E_sim_norm <- mu_alpha[,3]
    strat=rnorm(n=n_fam,mean=alpha_fam+cor_strat*mu_alpha[,1]+rho2*mu_alpha[,1]*E_sim_bin+rho2*mu_alpha[,1]*E_sim_norm, sd=1)

    xb <- strat + betaG_normPRS*prs + betaE_bin*E_sim_bin + betaE_norm*E_sim_norm + betaGE_normPRS_bin*E_sim_bin*prs+betaGE_normPRS_norm*E_sim_norm*prs
  } else {
    mu_sim=c(0,alpha_fam) #this is the mean vector for family PRS mu_i and alpha_i in the disease model
    sigma_sim=matrix(c(1,cor_strat,cor_strat,1),2,2) #assume there is population stratification
    mu_alpha=MASS::mvrnorm(n=n_fam,mu=mu_sim,Sigma=sigma_sim)
    mu_var=cbind(mu_alpha[,1],var_sim)
    myfunctionpgs=function(x){
      sigma_tmp=matrix(c(1,0.5,0.5,0.5,1,0,0.5,0,1)*x[2],3,3)
      MASS::mvrnorm(n=1,mu=rep(x[1],3),Sigma=sigma_tmp)

    }
    pgs_fam=t(apply(mu_var,1,myfunctionpgs))

    colnames(pgs_fam)=c("pgs_c","pgs_m","pgs_f")
    prs <- as.numeric(as.vector(pgs_fam[,1]))

    xb <- mu_alpha[,2] + betaG_normPRS*prs}


  prob <- 1/(1 + exp(-xb)) #Even though our assumption is that the disease model follows log-linear, we simulate data using logistic regression
  D_sim <- rbinom(n = n_fam, size = 1, prob = prob)  # 1=case, 0=control
  cat("Disease prevalence: ",mean(D_sim))

  if(envir==TRUE){
    E_sim = cbind(E_sim_bin,E_sim_norm)} else {E_sim=NULL}

  #### Bind D, G_logistic,PRS_fam, and E into a list
  dat = list(D_sim=D_sim, E_sim=E_sim, pgs_fam=pgs_fam)
  return(dat)


}
