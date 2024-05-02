#------------------------------------------------------------------
#PGS in case-parent trio analysis
#1. The proposed retrospective likelihood-based method, pgs.cpt()
#2. pTDT test proposed by Weiner et al., Nat Genet 2017
#April 26, 2023
#------------------------------------------------------------------

pgs.cpt = function(pgs_offspring, #The PGS values of the affected probands (children). A vector of length N, no missing values are allowed
                   pgs_mother, #The PGS values of mothers that corresponds to the children. A vector of same length N, no missing values are allowed
                   pgs_father, #The PGS values of fathers that corresponds to the children. A vector of same length N, no missing values are allowed
                   GxE_int = FALSE, #Whether there are interaction effect between pgs and environmental variables that are of interest in the model. If FALSE, then "formula" and "E" are ignored. 
                   parental_nurture = FALSE, #Whether to estimate potential parental nurturing effect, returns an estimated difference of mother and father parental effect (delta_MF = beta_M - beta_F). Note that when this is TRUE, GxE_int will be ignored.
                   formula= ~ envir1 +envir2+factor(s1), #The environmental variables of interest for the PGSxE interaction effect
                   E, #The environmental variables of interest for interaction effect. A vector of length N for one environmental variable or a data frame/data matrix of NxP for P environmental variables are allowed.
                   side = 2, #Sided of the Wald test, default is 2-sided.
                   smalltriosize = FALSE, #Whether number of trios is small (<100), if TRUE, a t test will be used rather than a wald test.
                   numDeriv = FALSE #Whether to use the score function (FALSE) or numerically derive (TRUE) for parameter estimation, by default, using the score function (it is computationally faster and more accurate). 
){
  
  if(any(is.na(c(pgs_offspring,pgs_mother,pgs_mother))) == TRUE) {
    stop("There are missing values in the family PGS, remove them to run the analysis") }
  
  if( (length(pgs_offspring) != length(pgs_mother) | length(pgs_offspring) != length(pgs_father) ) == TRUE ) {
    stop(paste0("The number of family PGS values do not match with each other! N(Offspring) = ", length(pgs_offspring),", N(Mother) =",length(pgs_mother), ", N(Father) = ",length(pgs_father))) }
  
  
  pgs.tdt=function(pgs_c,pgs_m,pgs_f,side0=2){
    cat(paste("The complete number of trios is",length(pgs_c),"\n"))
    var_fam=1/2*(pgs_m-pgs_f)^2
    beta_hat=2*sum(pgs_c-(pgs_m+pgs_f)/2)/sum(var_fam)
    var_beta_hat=2/sum(var_fam)
    
    #Summarize the results
    res.sum=function (parms=beta_hat, cov=var_beta_hat, sided) 
    {
      if (sided != 1) 
        sided <- 2
      cols <- c("Estimate", "Std.Error", "Z.value", "Pvalue")
      n <- length(parms)
      ret <- matrix(data = NA, nrow = n, ncol = 4)
      pnames <- names(parms)
      rownames(ret) <- pnames
      colnames(ret) <- cols
      ret[, 1] <- parms
      cov <- sqrt(cov)
      if (is.null(pnames)) 
        pnames <- 1:n
      ret[, 2] <- cov
      ret[, 3] <- parms/cov
      ret[, 4] <- sided * pnorm(abs(ret[, 3]), lower.tail = FALSE)
      ret
    }
    
    res.sum.t=function (parms=beta_hat, cov=var_beta_hat, df0, sided) 
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
      cov <- sqrt(cov)
      if (is.null(pnames)) 
        pnames <- 1:n
      ret[, 2] <- cov
      ret[, 3] <- parms/cov
      ret[, 4] <- sided * pt(abs(ret[, 3]), df0,lower.tail = F)
      ret
    }
    
    if(smalltriosize == FALSE){
    res_beta=res.sum(sided = side0)} else {
      res_beta=res.sum.t(df0 = (length(pgs_c)-1) ,sided = side0)
    }
    rownames(res_beta)="PGS"
    res=list(res_beta=res_beta, var_fam=var_fam)
    return(res)
  }
  
  
  
  pgs.tdt.nurture=function(pgs_c,pgs_m,pgs_f,side0=2){
    cat(paste("The complete number of trios is",length(pgs_c),"\n"))
    n_family=length(pgs_c)
    x_bar=sum(pgs_m-pgs_f)/n_family
    
    var_fam_sum=1/2*sum((pgs_m-pgs_f)^2)-1/2/n_family*(sum(pgs_m-pgs_f))^2
    delta_mf=sum(pgs_m-pgs_f)/var_fam_sum
    var_delta=2/var_fam_sum

    beta_hat=2*sum(pgs_c-(pgs_m+pgs_f)/2)/var_fam_sum
    var_beta_hat=2/var_fam_sum
    
    #Summarize the results
    res.sum=function (parms=beta_hat, cov=var_beta_hat, sided) 
    {
      if (sided != 1) 
        sided <- 2
      cols <- c("Estimate", "Std.Error", "Z.value", "Pvalue")
      n <- length(parms)
      ret <- matrix(data = NA, nrow = n, ncol = 4)
      pnames <- names(parms)
      rownames(ret) <- pnames
      colnames(ret) <- cols
      ret[, 1] <- parms
      # cols <- colnames(cov)
      #cov <- sqrt(diag(cov))
      cov <- sqrt(cov)
      # names(cov) <- cols
      if (is.null(pnames)) 
        pnames <- 1:n
      # cov <- cov[pnames]
      ret[, 2] <- cov
      ret[, 3] <- parms/cov
      ret[, 4] <- sided * pnorm(abs(ret[, 3]), lower.tail = FALSE)
      ret
    }
    
    res.sum.t=function (parms=beta_hat, cov=var_beta_hat, df0, sided) 
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
      cov <- sqrt(cov)
      if (is.null(pnames)) 
        pnames <- 1:n
      ret[, 2] <- cov
      ret[, 3] <- parms/cov
      ret[, 4] <- sided * pt(abs(ret[, 3]), df0,lower.tail = F)
      ret
    }
    
    if(smalltriosize == FALSE){
      
      res_beta=res.sum(sided = side0)
      res_delta=res.sum(parms=delta_mf, cov=var_delta, sided= side0) } else {
        
      res_beta=res.sum.t(df0 = (length(pgs_c)-1) ,sided = side0)
      res_delta=res.sum.t(parms=delta_mf, cov=var_delta, df0 = (length(pgs_c)-1), sided= side0)
      
      }
    
   # res_beta=res.sum(sided = side0)
    rownames(res_beta)="PGS"
   # res_delta=res.sum(parms=delta_mf, cov=var_delta, sided= side0)
    rownames(res_delta)="Nurture_Diff_MF"
    
    res=list(res_beta=res_beta,res_delta=res_delta,var_fam_sum=var_fam_sum)
    return(res)
    
  }
  
  pgs.tdt.gxe=function(pgs_c,pgs_m,pgs_f,formula0,envir0,side0=2,numDeriv0=F){
    
    envir0=data.frame(envir0)
    options(na.action='na.pass')
    envir=model.matrix(formula0,data=envir0)[,-1] #remove intercept term
    
    #Remove NA values in the environmental variables
    id=complete.cases(envir)
    #if( (length(formula0[[2]]) == 1 & is.null(dim(envir)))){
      if(  is.null(dim(envir))){
      envir=envir[id]
      
    } else {
      
      envir=envir[id,]
      
    }
    
    pgs_c=pgs_c[id]
    pgs_m=pgs_m[id]
    pgs_f=pgs_f[id]
    
    cat(paste("The complete number of trios with non-missing environmental variables is",length(pgs_c),"\n"))
    
    var_fam=1/2*(pgs_m-pgs_f)^2
    mu_c=(pgs_m+pgs_f)/2
    var_c=var_fam/2
    
    #Generate initial values
    #Initial value of the main PGS effect using PGS only tdt 
    beta_ini=2*sum(pgs_c-(pgs_m+pgs_f)/2)/sum(var_fam)
    
    #Initial value of the PGSxE interaction effect using case-only analysis
    fit_caseonly <- lm(pgs_c ~ .,data=data.frame(envir))
    beta_ge_ini=fit_caseonly$coefficients[-1]/(sd(fit_caseonly$residuals))^2
    
    param0=c(beta_ini,beta_ge_ini)
    if(length(beta_ge_ini)==1){
      names(param0)=c("PGS",paste0("PGS x ",labels(terms(formula0)))) #"beta_pgsxE")
    } else {
      names(param0)=c("PGS",paste0("PGS x ",colnames(envir)))
    }
    
    #Likelihood, find MLE
    loglilke_strat <- function(param){
      
      beta=param[1]
      beta_ge=param[2:length(param)]
      
      if(length(beta_ge)==1){
        
        log_l1=-var_c*(beta+envir*beta_ge)^2/2+(pgs_c-mu_c)*beta+(pgs_c-mu_c)*envir*beta_ge
        
      } else{
        
        log_l1=-var_c*(beta+as.vector(envir%*%beta_ge))^2/2+(pgs_c-mu_c)*beta+as.vector((pgs_c-mu_c)*(envir%*%beta_ge))
        
      }
      
      ret <- sum(log_l1)
      return(ret)
      
    }
    
    fn_gr <- function(param){
      
      beta=param[1]
      beta_ge=param[2:length(param)]
      
      if(length(beta_ge)==1){
        
        beta.tmp=sum(-var_c*(beta+envir*beta_ge)+(pgs_c-mu_c))
        beta_ge.tmp=sum(-var_c*(beta+as.vector(envir*beta_ge))*envir+(pgs_c-mu_c)*envir)
        
      } else{
        
        beta.tmp=sum(-var_c*(beta+as.vector(envir%*%beta_ge))+(pgs_c-mu_c))
        beta_ge.tmp=as.vector(-var_c*(beta+as.vector(envir%*%beta_ge)))%*%envir + as.vector(pgs_c-mu_c)%*%envir
        
      }
      
      grr <- c(beta.tmp,beta_ge.tmp)
      names(grr)=names(param)
      return(grr)
      
    }
    
    control <- list(fnscale = -1,trace = TRUE,
                    REPORT = 50,maxit=20000)
    
    if(numDeriv0==T){
      
      ret <- optim(param0, loglilke_strat, method="BFGS",
                   control = control, hessian = TRUE)
      
    } else {
      
      ret <- optim(param0, loglilke_strat,gr=fn_gr, method="BFGS",
                   control = control, hessian = TRUE)
      
    }
    
    cov <- chol(-ret$hessian)
    cov <- chol2inv(cov)
    cov <- sqrt(diag(cov))
    names(cov) <- names(param0)
    
    cov1=cov 
    
    #Summarize the results
    res.sum=function (parms=ret$par, cov=cov1, sided) 
    {
      if (sided != 1) 
        sided <- 2
      cols <- c("Estimate", "Std.Error", "Z.value", "Pvalue")
      n <- length(parms)
      ret <- matrix(data = NA, nrow = n, ncol = 4)
      pnames <- names(parms)
      rownames(ret) <- pnames
      colnames(ret) <- cols
      ret[, 1] <- parms
      if (is.null(pnames)) 
        pnames <- 1:n
      cov <- cov[pnames]
      ret[, 2] <- cov
      ret[, 3] <- parms/cov
      ret[, 4] <- sided * pnorm(abs(ret[, 3]), lower.tail = FALSE)
      ret
    }
    
    res.sum.t=function (parms=ret$par, cov=cov1, df0, sided) 
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
      cov <- cov[pnames]
      ret[, 2] <- cov
      ret[, 3] <- parms/cov
      ret[, 4] <- sided * pt(abs(ret[, 3]), df0,lower.tail = F)
      ret
    }
    
    if(smalltriosize == FALSE){
      
      res_beta=res.sum(sided = side0) } else {
        
      res_beta=res.sum.t(df0 = (length(pgs_c)-1) ,sided = side0)
        
      }
    
    
    res=list(res_beta=res_beta, var_fam=var_fam)
    return(res)
    
  }
  
  
  if (parental_nurture == TRUE){
    
    pgs.tdt.nurture(pgs_c = pgs_offspring,
                    pgs_m = pgs_mother,
                    pgs_f = pgs_father,
                    side0 = side)
    
  } else if (GxE_int == FALSE){
    
    pgs.tdt(pgs_c = pgs_offspring,
            pgs_m = pgs_mother,
            pgs_f = pgs_father,
            side0 = side)
    
  } else {
    
    pgs.tdt.gxe(pgs_c = pgs_offspring,
                pgs_m = pgs_mother,
                pgs_f = pgs_father,
                formula0 = formula,
                envir0 = E,
                side0 = side,
                numDeriv0 = numDeriv)
    
  }
  
  
  
}





#pTDT test proposed by Weiner et al., Nat Genet 2017

ptdt=function(pgs_c,pgs_m,pgs_f,side0=2){
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
