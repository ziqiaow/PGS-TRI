#' @title PGS.TRI
#' @description A unified framework for the estimation of polygenic direct effect, indirect effect, and GxE interactions in case-parent trio studies
#'
#' @param pgs_offspring The PGS values of the affected probands (children). A vector of length N, no missing values are allowed
#' @param pgs_mother The PGS values of mothers that corresponds to the children. A vector of same length N, no missing values are allowed
#' @param pgs_father The PGS values of fathers that corresponds to the children. A vector of same length N, no missing values are allowed
#' @param GxE_int Whether there are interaction effect between pgs and environmental variables that are of interest in the model. If FALSE, then "formula" and "E" are ignored.
#' @param parental_indirect Whether to estimate potential parental indirent effect, returns an estimated difference of mother and father parental effect (delta_MF = beta_M - beta_F).
#' @param formula The environmental variables of interest for the PGSxE interaction effect
#' @param E The environmental variables of interest for interaction effect. A vector of length N for one environmental variable or a data frame/data matrix of NxP for P environmental variables are allowed.
#' @param side Sided of the Wald test, default is 2-sided.
#' @param smalltriosize Whether number of trios is small (<100), if TRUE, a t test will be used rather than a wald test.
#'
#' @return A list of results of PGS.TRI
#'  \item{res_beta}{Results of direct PGS effect, if GxE_int is TRUE, then the result will also include PGSxE interaction effects}
#'  \item{res_delta}{Results of indirect parental PGS effect difference: PGS_mother - PGS_father}
#'  \item{var_fam}{Within-family variances for each family}
#'  \item{var_fam_sum}{Sum of within-family variances}
#'
#' @export
#'
#'
PGS.TRI = function(pgs_offspring, #The PGS values of the affected probands (children). A vector of length N, no missing values are allowed
                   pgs_mother, #The PGS values of mothers that corresponds to the children. A vector of same length N, no missing values are allowed
                   pgs_father, #The PGS values of fathers that corresponds to the children. A vector of same length N, no missing values are allowed
                   GxE_int = FALSE, #Whether there are interaction effect between pgs and environmental variables that are of interest in the model. If FALSE, then "formula" and "E" are ignored.
                   parental_indirect = FALSE, #Whether to estimate potential parental nurturing effect, returns an estimated difference of mother and father parental effect (delta_MF = beta_M - beta_F). Note that when this is TRUE, GxE_int will be ignored.
                   formula= ~ envir1 +envir2+factor(s1), #The environmental variables of interest for the PGSxE interaction effect
                   E, #The environmental variables of interest for interaction effect. A vector of length N for one environmental variable or a data frame/data matrix of NxP for P environmental variables are allowed.
                   side = 2, #Sided of the Wald test, default is 2-sided.
                   smalltriosize = FALSE #Whether number of trios is small (<100), if TRUE, a t test will be used rather than a wald test.
){

  if(any(is.na(c(pgs_offspring,pgs_mother,pgs_mother))) == TRUE) {
    stop("There are missing values in the family PGS, remove them to run the analysis") }

  if( (length(pgs_offspring) != length(pgs_mother) | length(pgs_offspring) != length(pgs_father) ) == TRUE ) {
    stop(paste0("The number of family PGS values do not match with each other! N(Offspring) = ", length(pgs_offspring),", N(Mother) =",length(pgs_mother), ", N(Father) = ",length(pgs_father))) }

  #Summarize the results
  res.sum=function(parms=beta_hat, sd=sd_beta_hat, sided)
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
    # cov <- sqrt(cov)
    if (is.null(pnames))
      pnames <- 1:n
    ret[, 2] <- sd
    ret[, 3] <- parms/sd
    ret[, 4] <- sided * pnorm(abs(ret[, 3]), lower.tail = FALSE)
    ret
  }

  res.sum.t=function(parms=beta_hat, sd=sd_beta_hat, df0, sided)
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
    # cov <- sqrt(cov)
    if (is.null(pnames))
      pnames <- 1:n
    ret[, 2] <- sd
    ret[, 3] <- parms/sd
    ret[, 4] <- sided * pt(abs(ret[, 3]), df0,lower.tail = F)
    ret
  }




  pgs.tdt=function(pgs_c,pgs_m,pgs_f,side0=2){
    cat(paste("The complete number of trios is",length(pgs_c),"\n"))
    var_fam=1/2*(pgs_m-pgs_f)^2
    beta_hat=2*sum(pgs_c-(pgs_m+pgs_f)/2)/sum(var_fam)
    sigma_4 = 1/12*(pgs_m-pgs_f)^4 #2/3*((pgs_m-(pgs_m+pgs_f)/2)^4 + (pgs_f - (pgs_m+pgs_f)/2)^4)
    var_beta_hat=2/sum(var_fam)+2*beta_hat^2*sum(sigma_4)/(sum(var_fam))^2
    sd_beta_hat=sqrt(var_beta_hat)


    if(smalltriosize == FALSE){
      res_beta=res.sum(parms=beta_hat, sd=sd_beta_hat, sided = side0)} else {
        res_beta=res.sum.t(parms=beta_hat, sd=sd_beta_hat, df0 = (length(pgs_c)-1) ,sided = side0)
      }
    rownames(res_beta)="PGS"
    res=list(res_beta=res_beta, var_fam=var_fam)
    return(res)
  }



  pgs.tdt.nurture=function(pgs_c,pgs_m,pgs_f,side0=2){
    cat(paste("The complete number of trios is",length(pgs_c),"\n"))
    n_family=length(pgs_c)
    x_bar=sum(pgs_m-pgs_f)/n_family

    # var_fam_sum=1/(2*(n_family-1))*sum((pgs_m-pgs_f)^2)-1/(2*(n_family-1)*n_family)*(sum(pgs_m-pgs_f))^2
    var_fam_sum=n_family/(2*(n_family-1))*sum((pgs_m-pgs_f)^2)-1/(2*(n_family-1))*(sum(pgs_m-pgs_f))^2
    delta_mf=sum(pgs_m-pgs_f)/var_fam_sum
    var_delta=2/var_fam_sum #+ 8/(n_family-1)*delta_mf^2
    sd_delta = sqrt(var_delta)

    beta_hat=2*sum(pgs_c-(pgs_m+pgs_f)/2)/var_fam_sum
    #sigma_4 = (n_family-1)/((n_family+1)*n_family)*var_fam_sum^2
    sigma_4 = 1/12*(pgs_m-pgs_f - mean(pgs_m - pgs_f))^4*n_family^2/(n_family -1)^2
    var_beta_hat=2/var_fam_sum+2*beta_hat^2*sum(sigma_4)/(var_fam_sum)^2
    sd_beta_hat = sqrt(var_beta_hat)

    if(smalltriosize == FALSE){

      res_beta=res.sum(parms=beta_hat, sd=sd_beta_hat, sided = side0)
      res_delta=res.sum(parms=delta_mf, sd=sd_delta, sided= side0) } else {

        res_beta=res.sum.t(parms=beta_hat, sd=sd_beta_hat, df0 = (length(pgs_c)-1) ,sided = side0)
        res_delta=res.sum.t(parms=delta_mf, sd=sd_delta, df0 = (length(pgs_c)-1), sided= side0)

      }

    rownames(res_beta)="PGS"
    rownames(res_delta)="Indirect_Diff_MF"

    res=list(res_beta=res_beta,res_delta=res_delta,var_fam_sum=var_fam_sum)
    return(res)

  }

  pgs.tdt.gxe=function(pgs_c,pgs_m,pgs_f,formula0,envir0,side0=2,numDeriv0=F){

    envir0=data.frame(envir0)
    options(na.action='na.pass')
    envir=model.matrix(formula0,data=envir0)

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
    sigma_4 = 1/12*(pgs_m-pgs_f)^4


    m_list = list()
    for(i in 1:dim(envir)[1]){

      m_list[[i]] = var_fam[i]*envir[i,] %*% t(envir[i,])
    }

    m_w = Reduce(`+`, lapply(m_list,function(x) {x[is.na(x)] <-0;x}))
    m_inv = chol2inv(chol(m_w))
    beta_hat = 2* m_inv %*% as.vector(apply(((pgs_c - mu_c)* envir), 2, sum))

    # var_beta = 2* m_inv %*% m_w %*% t(m_inv)
    # sd_beta = sqrt(diag(var_beta))



    m_list2 = list()
    for(i in 1:dim(envir)[1]){
      m_list2[[i]] = sigma_4[i]*(m_inv %*% envir[i,] %*% t(envir[i,])  %*%  m_inv %*% (m_w %*% beta_hat)) %*% t(m_inv %*% envir[i,] %*% t(envir[i,])  %*%  m_inv %*% (m_w %*% beta_hat))

    }

    m_w2 = Reduce(`+`, lapply(m_list2,function(x) {x[is.na(x)] <-0;x}))

    cov_y = 2*m_w2

    mu_x = m_w %*% beta_hat

    var_beta_taylor = 2* m_inv %*% t(m_inv) %*% m_w  + cov_y
    sd_beta_taylor = sqrt(diag(var_beta_taylor))



    if(dim(envir)[2]==2){
      names(beta_hat)=c("PGS",paste0("PGS x ",labels(terms(formula0)))) #"beta_pgsxE")
    } else {
      names(beta_hat)=c("PGS",paste0("PGS x ",colnames(envir)[-1]))
    }
    names(sd_beta_taylor)=names(beta_hat)




    if(smalltriosize == FALSE){
      res_beta=res.sum(parms=beta_hat,sd=sd_beta_taylor,sided = side0)
    } else {
      res_beta=res.sum.t(parms=beta_hat,sd=sd_beta_taylor,df0 = (length(pgs_c)-1) ,sided = side0)
    }

    res=list(res_beta=res_beta,  var_fam=var_fam)
    return(res)

  }


  pgs.tdt.nurture.gxe=function(pgs_c,pgs_m,pgs_f,formula0,envir0,side0=2){

    envir0=data.frame(envir0)
    options(na.action='na.pass')
    envir=model.matrix(formula0,data=envir0)

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



    n_family=length(pgs_c)
    x_bar=sum(pgs_m-pgs_f)/n_family

    var_fam_sum=n_family/(2*(n_family-1))*sum((pgs_m-pgs_f)^2)-1/(2*(n_family-1))*(sum(pgs_m-pgs_f))^2
    delta_mf=sum(pgs_m-pgs_f)/var_fam_sum
    var_delta=2/var_fam_sum
    sd_delta = sqrt(var_delta)

    sigma_4 = 1/12*(pgs_m-pgs_f - mean(pgs_m - pgs_f))^4*n_family^2/(n_family -1)^2
    var_fam = 1/2*(pgs_m-pgs_f - mean(pgs_m - pgs_f))^2*n_family/(n_family-1)

    mu_c=(pgs_m+pgs_f)/2


    m_list = list()
    for(i in 1:dim(envir)[1]){

      m_list[[i]] = var_fam[i]*envir[i,] %*% t(envir[i,])
    }

    m_w = Reduce(`+`, lapply(m_list,function(x) {x[is.na(x)] <-0;x}))
    m_inv = chol2inv(chol(m_w))
    beta_hat = 2* m_inv %*% as.vector(apply(((pgs_c - mu_c)* envir), 2, sum))


    m_list2 = list()
    for(i in 1:dim(envir)[1]){
      m_list2[[i]] = sigma_4[i]*(m_inv %*% envir[i,] %*% t(envir[i,])  %*%  m_inv %*% (m_w %*% beta_hat)) %*% t(m_inv %*% envir[i,] %*% t(envir[i,])  %*%  m_inv %*% (m_w %*% beta_hat))

    }

    m_w2 = Reduce(`+`, lapply(m_list2,function(x) {x[is.na(x)] <-0;x}))

    cov_y = 2*m_w2

    mu_x = m_w %*% beta_hat

    var_beta_taylor = 2* m_inv %*% t(m_inv) %*% m_w  + cov_y
    sd_beta_taylor = sqrt(diag(var_beta_taylor))





    if(dim(envir)[2]==2){
      names(beta_hat)=c("PGS",paste0("PGS x ",labels(terms(formula0)))) #"beta_pgsxE")
    } else {
      names(beta_hat)=c("PGS",paste0("PGS x ",colnames(envir)[-1]))
    }
    names(sd_beta_taylor)=names(beta_hat)



    if(smalltriosize == FALSE){

      res_beta=res.sum(parms=beta_hat,sd=sd_beta_taylor,sided = side0)
      res_delta=res.sum(parms=delta_mf, sd=sd_delta, sided= side0)
    } else {
      res_beta=res.sum.t(parms=beta_hat,sd=sd_beta_taylor,df0 = (length(pgs_c)-1) ,sided = side0)
      res_delta=res.sum.t(parms=delta_mf, sd=sd_delta, df0 = (length(pgs_c)-1), sided= side0)

    }

    rownames(res_delta)="Indirect_Diff_MF"

    res=list(res_beta=res_beta,res_delta=res_delta,var_fam_sum=var_fam_sum)
    return(res)

  }


  if (parental_indirect == TRUE){

    if (GxE_int == FALSE){
      pgs.tdt.nurture(pgs_c = pgs_offspring,
                      pgs_m = pgs_mother,
                      pgs_f = pgs_father,
                      side0 = side)

    } else {

      pgs.tdt.nurture.gxe(pgs_c = pgs_offspring,
                          pgs_m = pgs_mother,
                          pgs_f = pgs_father,
                          formula0 = formula,
                          envir0 = E,
                          side0 = side)

    }

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
                side0 = side)

  }



}



