SCMM.gene <- function(Y, Xmu, Xp, sample_x, offsetMu=NULL, offsetP=NULL, 
  sample_id=NULL, ZINB_F=FALSE, ZINB_R=FALSE, debug=FALSE){
  
  options(stringsAsFactors=FALSE)
  #if(debug){
  #  cat(paste0("class(x)=", class(x), "\n"))
  #}
  x <- Xmu
  Y <- matrix(Y, ncol=1)
  if(ZINB_F){
    if(!is.null(sample_x)){ 
        Xmu0 <- as.matrix(model.matrix(lm(Y ~ sample_id))) # nested design 
    }else{
        Xmu0 <- as.matrix(model.matrix(lm(Y ~ Xmu + sample_id))) # one subject can have both cells in case and control 
    }
  }else{
    Xmu0 <- as.matrix(model.matrix(lm(Y ~ Xmu)))
  }

  #ifelse(is.null(Xp)==TRUE, Xp0 <- matrix(1,nrow(Y),1), Xp0 <- as.matrix(model.matrix(lm(Y~Xp))) )
  Xp0 <- Xmu0 # make logit part the same as NB part


  if(is.null(offsetMu)==TRUE){offsetMu <- matrix(rep(0, length(Y)))}
  if(is.null(offsetP)==TRUE){offsetP <- matrix(rep(0, length(Y)))}
  offsetMu <- matrix(offsetMu, ncol=1); offsetP <- matrix(offsetP, ncol=1)
  
  # better initial (remove intercept) 
  Xmu0_noInt <- Xmu0[,-1, drop=FALSE]  
  #Xp0_noInt <- Xp0[,-1, drop=FALSE]
  Xp0_noInt <- Xmu0_noInt # make logit part the same as NB part

  f.nb <- tryCatch(glm.nb(Y ~ offset(offsetMu) + Xmu0_noInt),
    error=function(e) NA)

  f.p <- tryCatch(glm((Y!=0) ~ offset(offsetP) + Xp0_noInt, family = "binomial"),
    error=function(e) NA)

  if(all(is.na(f.nb))) {
    f.nb <- list(coefficients=rep(0, ncol(Xmu0_noInt), theta=1))
  }

  if(all(is.na(f.p))){
    f.p <- list(coefficients=rep(0, ncol(Xp0_noInt)))
  }


  if(ZINB_R){
    nK <- length(unique(sample_id))
    Data <- list(y=Y, Xu=Xmu0, Xp=Xp0, osu=offsetMu, osp=offsetP, group=sample_id)
    Parameters <- list(logtheta=log(f.nb$theta+10^(-6)), 
      logsd_mu=0,logsd_p=0, b_mu=rnorm(nK),
      b_p=rnorm(nK), 
      BETAu=f.nb$coefficients, BETAp=f.p$coefficients)
    if(debug){
      print("parameters are")
      print(Parameters)
    }
    Obj <- MakeADFun(data=Data, parameters=Parameters, random=c("b_mu","b_p"), 
      DLL="SCMMRan",silent=TRUE)
    Opt <- nlminb(Obj$par, Obj$fn, Obj$gr)
    rep <- sdreport(Obj, getReportCovariance = TRUE)
    rep.fixed <- summary(rep,"fixed")
    if(all(is.nan(rep.fixed[,"Std. Error"]))){
      Parameters <- list(logtheta=0, logsd_mu=0,logsd_p=0, b_mu=rnorm(nK), 
        b_p=rnorm(nK), 
        BETAu=rep(0, ncol(Xmu0)), BETAp=rep(0, ncol(Xp0)))
      Obj <- MakeADFun(data=Data, parameters=Parameters, random=c("b_mu", "b_p"), 
        DLL="SCMMRan",silent=TRUE)
      Opt <- nlminb(Obj$par, Obj$fn, Obj$gr)
      rep <- sdreport(Obj, getReportCovariance = TRUE)
      rep.fixed <- summary(rep,"fixed")
    }
  }else{
    Data <- list(y=Y, Xu=Xmu0, Xp=Xp0, osu=offsetMu, osp=offsetP)
    Parameters <- list(logtheta=log(f.nb$theta+10^(-6)), 
      BETAu=f.nb$coefficients, BETAp=f.p$coefficients)
    Obj <- MakeADFun(data=Data, parameters=Parameters, DLL="SCMM", silent=TRUE)
    Opt <- nlminb(Obj$par, Obj$fn, Obj$gr)
    rep <- sdreport(Obj, getReportCovariance = TRUE)
    rep.fixed <- summary(rep,"fixed")
    if(all(is.nan(rep.fixed[,"Std. Error"]))){
      Parameters <- list(logtheta=0, 
            BETAu=rep(0, ncol(Xmu0)), BETAp=rep(0, ncol(Xp0)))
      Obj <- MakeADFun(data=Data, parameters=Parameters, DLL="SCMM", silent=TRUE)
      Opt <- nlminb(Obj$par, Obj$fn, Obj$gr)
      rep <- sdreport(Obj, getReportCovariance = TRUE)
      rep.fixed <- summary(rep,"fixed")
    }
  }


  nb0 <- rep.fixed[rownames(rep.fixed)=="BETAu",]; nb.z <- nb0[,"Estimate"]/nb0[,"Std. Error"]
  nb <- cbind(nb0,nb.z,2*(pnorm(-abs(nb.z))))
  colnames(nb) <- c("Estimate","Std. Error","z value","Pr(>|z|)")
  rownames(nb) <- colnames(Xmu0)
  logit0 <- rep.fixed[rownames(rep.fixed)=="BETAp",]; logit.z <- logit0[,"Estimate"]/logit0[,"Std. Error"]
  logit <- cbind(logit0,logit.z,2*(pnorm(-abs(logit.z))))
  colnames(logit) <- c("Estimate","Std. Error","z value","Pr(>|z|)")
  rownames(logit) <- colnames(Xp0)

  nb <- rbind(nb, rep(NA, ncol(nb)))
  rownames(nb)[nrow(nb)] <- "nb_beta"

  logit <- rbind(logit, rep(NA, ncol(logit)))
  rownames(logit)[nrow(logit)] <- "logit_beta"

  if(ZINB_F){
    par_cov_nb <- rep$cov.fixed[(2:(ncol(Xmu0)+1)), (2:(ncol(Xmu0)+1))] # first column/row is the logtheta (overdispersion)
    par_cov_logit <- rep$cov.fixed[((ncol(Xmu0)+2):(2*ncol(Xmu0)+1)), ((ncol(Xmu0)+2):(2*ncol(Xmu0)+1))]

    if((!is.na(mean(par_cov_nb))) & (!is.null(sample_x))) {

      n_x <- table(sample_x)
      H <- matrix(c(0, ifelse(sample_x==levels(x)[2],1/n_x[2],-1/n_x[1])[-1]), 1, length(sample_x))
      if(debug){
        print("sample_x")
        print(sample_x)
        print("H")
        print(H)
        print("levels(x)")
        print(levels(x))
        print("x[1:5]")
        print(x[1:5])
      }
      
      est_hat_nb <- H %*% nb[-nrow(nb),1]
      cov_hat_nb <- H %*% par_cov_nb %*% t(H)
      t_stat_nb <- est_hat_nb/sqrt(cov_hat_nb)
      df_nb <- nrow(Xmu0)-ncol(Xmu0)
      pval_nb <- 2*pt(abs(t_stat_nb), df_nb, lower.tail=FALSE)
      nb["nb_beta",] <- c(est_hat_nb, sqrt(cov_hat_nb), t_stat_nb, pval_nb)

      est_hat_logit <- H %*% logit[-nrow(logit),1]
      cov_hat_logit <- H %*% par_cov_logit %*% t(H)
      t_stat_logit <- est_hat_logit/sqrt(cov_hat_logit)
      df_logit <- nrow(Xp0)-ncol(Xp0)
      pval_logit <- 2*pt(abs(t_stat_logit), df_logit, lower.tail=FALSE)
      logit["logit_beta",] <- c(est_hat_logit, sqrt(cov_hat_logit), t_stat_logit, pval_logit)

    }else{
      nb["nb_beta",] <- nb[2, ]
      logit["logit_beta",] <- logit[2, ]
    }
  }

  if(ZINB_R){
    sd_mu_random <- summary(rep,"report")["sd_mu",1]
    sd_p_random <- summary(rep,"report")["sd_p",1]
  }else{
    sd_mu_random <- sd_p_random <- NA
  }

  if(debug){
    print("nb is")
    print(nb)
  }
  
  return(list(nb=nb,logit=logit,
                 theta=summary(rep,"report")["theta",1], 
                 sd_mu_random=sd_mu_random, sd_p_random=sd_p_random))  
}

