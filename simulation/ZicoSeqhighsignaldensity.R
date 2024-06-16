library(GUniFrac)
library(tools)
library(MASS)
library(mgcv)
library(fql)
library(mgcv); library(MASS)
library(VGAM)
library(MASS)
library(rmutil)
library(fql)
library(actuar)
library(boot)
library(sandwich)
library(lmtest)
library(robust)
library(robustbase)
library('RSimNorm')
library('MicrobiomeStat')
library(phyloseq)
library(ANCOMBC)

load("adenomas.RData")
com=data.obj$otu.tab
com=com[rowSums(com)!=0,]
com <- com[rowMeans(com != 0) > 0.25, ] 


#orginal code
rdirichlet.m <- function (alpha) {
  Gam <- matrix(rgamma(length(alpha), shape = alpha), nrow(alpha), ncol(alpha))
  t(t(Gam) / colSums(Gam))
}

EstPara <- function (ref.otu.tab) {
  
  if (is.null(rownames(ref.otu.tab))) {
    rownames(ref.otu.tab) <- paste0('OTU', 1 : nrow(ref.otu.tab))
  } # otu * sample
  samplenames = colnames(ref.otu.tab)
  taxnames = rownames(ref.otu.tab)
  
  dirmult.paras <- dirmult::dirmult(t(ref.otu.tab))
  
  gamma = dirmult.paras$gamma
  names(gamma) = names(dirmult.paras$pi)
  
  # Add pseduo count(each OTU add gamma estimated from dirmult)
  ref.otu.tab = sapply(1:ncol(ref.otu.tab), function (i) gamma + ref.otu.tab[,i]) # C_ij otu * sample
  
  # back to dirchlet, calculate the true proportion
  ref.otu.tab.p <- rdirichlet.m(ref.otu.tab) # P_ij nOTU*nSam
  colnames(ref.otu.tab.p) = samplenames
  rownames(ref.otu.tab.p) = taxnames
  
  # order OTUs by mean OTU proportion, for later selection
  ord = order(rowMeans(ref.otu.tab.p), decreasing = TRUE)
  ref.otu.tab.p =  ref.otu.tab.p[ord,]
  
  # apply size factor
  Si = exp(rnorm(ncol(ref.otu.tab.p)))
  ref.otu.tab0 = t(t(ref.otu.tab.p)*Si)
  colnames(ref.otu.tab0) = colnames(ref.otu.tab.p)
  return(list(mu = ref.otu.tab.p, ref.otu.tab = ref.otu.tab0))
}

SimulateMSeqUni<-function (ref.otu.tab, model.paras,nSam = 100, nOTU = 500, diff.otu.pct = 0.1, 
                           diff.otu.direct = c("balanced", "unbalanced"), diff.otu.mode = c("abundant", 
                                                                                            "rare", "mix"), covariate.type = c("binary", "continuous"), 
                           grp.ratio = 1, covariate.eff.mean = 1, covariate.eff.sd = 0, 
                           confounder.type = c("none", "binary", "continuous", "both"), 
                           conf.cov.cor = 0.6, conf.diff.otu.pct = 0, conf.nondiff.otu.pct = 0.1, 
                           confounder.eff.mean = 0, confounder.eff.sd = 0, error.sd = 0, 
                           depth.mu = 10000, depth.theta = 5, depth.conf.factor = 0) 
{
  diff.otu.direct <- match.arg(diff.otu.direct)
  diff.otu.mode <- match.arg(diff.otu.mode)
  covariate.type <- match.arg(covariate.type)
  confounder.type <- match.arg(confounder.type)
  ref.otu.tab0 <- ref.otu.tab
#  model.paras <- EstPara(ref.otu.tab = ref.otu.tab)
  sample.names <- colnames(model.paras$ref.otu.tab)
  ref.otu.tab <- model.paras$ref.otu.tab[(1:(nOTU)), ]
  dim(ref.otu.tab)
  idx.otu <- rownames(ref.otu.tab)
  idx.sample <- sample(sample.names, nSam)
  idx.nonsample <- colnames(ref.otu.tab)[!(colnames(ref.otu.tab) %in% 
                                             idx.sample)]
  ref.otu.tab = ref.otu.tab[, idx.sample]
  ref.otu.tab.unselect = ref.otu.tab0[c(1:(nOTU)), ][, idx.nonsample]
  if (confounder.type == "none") {
    confounder.type <- "continuous"
    confounder.eff.mean <- 0
    confounder.eff.sd <- 0
    Z <- NULL
  }
  if (confounder.type == "continuous") 
    Z <- cbind(rnorm(nSam))
  if (confounder.type == "binary") 
    Z <- cbind(c(rep(0, nSam%/%2), rep(1, nSam - nSam%/%2)))
  if (confounder.type == "both") 
    Z <- cbind(rnorm(nSam), c(rep(0, nSam%/%2), rep(1, nSam - 
                                                      nSam%/%2)))
  rho <- sqrt(conf.cov.cor^2/(1 - conf.cov.cor^2))
  if (covariate.type == "continuous") {
    X <- rho * scale(scale(Z) %*% rep(1, ncol(Z))) + runif(nSam)
  }
  if (covariate.type == "binary") {
    X <- rho * scale(scale(Z) %*% rep(1, ncol(Z))) + runif(nSam)
    X <- cbind(ifelse(X <= quantile(X, grp.ratio/(1 + grp.ratio)), 
                      0, 1))
  }
  names(X) <- colnames(ref.otu.tab)
  covariate.eff.mean1 = covariate.eff.mean
  covariate.eff.mean2 = covariate.eff.mean
  if (diff.otu.direct == "balanced") {
    if (diff.otu.mode == "abundant") {
      eta.diff <- sample(c(rnorm(floor(nOTU/2), mean = -covariate.eff.mean2, 
                                 sd = covariate.eff.sd), rnorm(nOTU - floor(nOTU/2), 
                                                               mean = covariate.eff.mean2, sd = covariate.eff.sd))) %*% 
        t(scale(X))
    }
    else if (diff.otu.mode == "rare") {
      eta.diff <- sample(c(rnorm(floor(nOTU/2), mean = -covariate.eff.mean2, 
                                 sd = covariate.eff.sd), rnorm(nOTU - floor(nOTU/2), 
                                                               mean = covariate.eff.mean2, sd = covariate.eff.sd))) %*% 
        t(scale(X))
    }
    else {
      eta.diff <- c(sample(c(rnorm(floor(nOTU/4), mean = -covariate.eff.mean1, 
                                   sd = covariate.eff.sd), rnorm(floor(nOTU/2) - 
                                                                   floor(nOTU/4), mean = covariate.eff.mean1, sd = covariate.eff.sd))), 
                    sample(c(rnorm(floor((nOTU - floor(nOTU/2))/2), 
                                   mean = -covariate.eff.mean2, sd = covariate.eff.sd), 
                             rnorm(nOTU - floor(nOTU/2) - floor((nOTU - 
                                                                   floor(nOTU/2))/2), mean = covariate.eff.mean2, 
                                   sd = covariate.eff.sd)))) %*% t(scale(X))
    }
  }
  if (diff.otu.direct == "unbalanced") {
    if (diff.otu.mode == "abundant") {
      eta.diff <- rnorm(nOTU, mean = covariate.eff.mean2, 
                        sd = covariate.eff.sd) %*% t(scale(X))
    }
    else if (diff.otu.mode == "rare") {
      eta.diff <- rnorm(nOTU, mean = covariate.eff.mean2, 
                        sd = covariate.eff.sd) %*% t(scale(X))
    }
    else {
      eta.diff <- c(sample(c(rnorm(floor(nOTU/2), mean = covariate.eff.mean1, 
                                   sd = covariate.eff.sd))), sample(c(rnorm(nOTU - 
                                                                              floor(nOTU/2), mean = covariate.eff.mean2, sd = covariate.eff.sd)))) %*% 
        t(scale(X))
    }
  }
  eta.conf <- sample(c(rnorm(floor(nOTU/2), mean = -confounder.eff.mean, 
                             sd = confounder.eff.sd), rnorm(nOTU - floor(nOTU/2), 
                                                            mean = confounder.eff.mean, sd = confounder.eff.sd))) %*% 
    t(scale(scale(Z) %*% rep(1, ncol(Z))))
  otu.ord <- 1:(nOTU)
  diff.otu.ind <- NULL
  diff.otu.num <- round(diff.otu.pct * nOTU)
  if (diff.otu.mode == "mix") 
    diff.otu.ind <- c(diff.otu.ind, sample(otu.ord, diff.otu.num))
  if (diff.otu.mode == "abundant") 
    diff.otu.ind <- c(diff.otu.ind, sample(otu.ord[1:round(length(otu.ord)/4)], 
                                           diff.otu.num))
  if (diff.otu.mode == "rare") 
    diff.otu.ind <- c(diff.otu.ind, sample(otu.ord[round(3 * 
                                                           length(otu.ord)/4):length(otu.ord)], diff.otu.num))
  if (length(diff.otu.ind) >= round(nOTU * conf.diff.otu.pct)) {
    conf.otu.ind1 <- sample(diff.otu.ind, round(nOTU * conf.diff.otu.pct))
  }else {
    conf.otu.ind1 <- diff.otu.ind
  }
  conf.otu.ind <- c(conf.otu.ind1, sample(setdiff(1:(nOTU), 
                                                  diff.otu.ind), round(conf.nondiff.otu.pct * nOTU)))
  eta.diff[setdiff(1:(nOTU), diff.otu.ind), ] <- 0
  eta.conf[setdiff(1:(nOTU), conf.otu.ind), ] <- 0
  eta.error <- matrix(rnorm(nOTU * nSam, 0, error.sd), nOTU, 
                      nSam)
  eta.exp <- exp(t(eta.diff + eta.conf + eta.error))
  eta.exp <- eta.exp * t(ref.otu.tab)
  ref.otu.tab.prop <- eta.exp/rowSums(eta.exp)
  ref.otu.tab.prop <- t(ref.otu.tab.prop)
  nSeq <- rnegbin(nSam, mu = depth.mu * exp(scale(X) * depth.conf.factor), 
                  theta = depth.theta)
  otu.tab.sim <- sapply(1:ncol(ref.otu.tab.prop), function(i) rmultinom(1, 
                                                                        nSeq[i], ref.otu.tab.prop[, i]))
  colnames(otu.tab.sim) <- rownames(eta.exp)
  rownames(otu.tab.sim) <- rownames(ref.otu.tab)
  diff.otu.ind = (1:nOTU) %in% diff.otu.ind
  conf.otu.ind = (1:nOTU) %in% conf.otu.ind
  return(list(otu.tab.sim = otu.tab.sim, covariate = X, confounder = Z, 
              diff.otu.ind = diff.otu.ind, otu.names = idx.otu, conf.otu.ind = conf.otu.ind))
}

MYFIT=function(sim.obj){
  
  meta.dat <- data.frame(X = sim.obj$covariate, Z1 = sim.obj$confounder[, 1]) 
  otu.tab.sim <- sim.obj$otu.tab.sim
  
  # Create a phyloseq object
  OTU <- otu_table(otu.tab.sim, taxa_are_rows = TRUE)
  rownames(meta.dat)<-colnames(otu.tab.sim)
  SAM <- sample_data(meta.dat)
  
  # Creating a phyloseq object
  physeq <- phyloseq(OTU, SAM)
  
  # Running ANCOMBC2
  ancombc2.obj <- ancombc2(
    data = physeq, tax_level = 'Species',
    #assay.type = "counts",
    fix_formula = "X",
    p_adj_method = "fdr",
    alpha = 0.05,
    neg_lb = FALSE,
    prv_cut = 0
  )
  
  #zicoseq
  # Run ZicoSeq for differential abundance analysis 
  zico.obj <- ZicoSeq(meta.dat = meta.dat, feature.dat = otu.tab.sim, grp.name = 'X',  
                      feature.dat.type = "count", 
                      # Filter to remove rare taxa 
                      prev.filter = 0, mean.abund.filter = 0, max.abund.filter = 0, min.prop = 0, 
                      # Winsorization to replace outliers 
                      is.winsor = FALSE,# outlier.pct = 0.03, winsor.end = 'top',
                      # Posterior sampling to impute zeros 
                      is.post.sample = TRUE, post.sample.no = 25, 
                      # Multiple link functions to capture diverse taxon-covariate relation 
                      link.func = list( function(x) log(x+10e-10)),
                      stats.combine.func = max, 
                      # Permutation-based multiple testing correction 
                      perm.no = 100, strata = NULL, 
                      # Reference-based multiple stage normalization 
                      ref.pct = 0.5, stage.no = 6, excl.pct = 0.2, 
                      # Family-wise error rate control 
                      is.fwer = FALSE, verbose = TRUE, return.feature.dat = FALSE) 
  iters=dim(otu.tab.sim)[1]
  estall=matrix(999, iters, 6);
  estall2=matrix(999, iters, 6);
  estall3=matrix(999, iters, 6);
  estall5=matrix(999, iters, 6);
  for (III in 1:iters) {
    set.seed(III)
    mydata=data.frame(Y=sim.obj$otu.tab.sim[III,],meta.dat)
    mydata
    
    #nb
    datas=data.frame(mydata); 
    fitnb<-glm.nb(Y~X,link="log",data=datas)
    summary(fitnb)$coefficients
    estall[III,1:6]=c(summary(fitnb)$coefficients[1,-3], summary(fitnb)$coefficients[2,-3])
    colnames(estall)<-c('intercept','se_intercept','p_intercept','beta1','se_beta1','p_beta1')
    #Estimate   Std. Error      z value     Pr(>|z|)
    
    #poi
    fitpoi<-glm(Y~X,poisson(link="log"),data=datas)
    summary(fitpoi)$coefficients
    estall2[III,1:6]=c(summary(fitpoi)$coefficients[1,-3], summary(fitpoi)$coefficients[2,-3])
    colnames(estall2)<-c('intercept','se_intercept','p_intercept','beta1','se_beta1','p_beta1')
    #Estimate   Std. Error      z value     Pr(>|z|)
    
    #poission bootstrap
    fitg3coef <- glm_coef(data= datas,formula = Y~X)
    fitg3results <- boot(data = datas, statistic = glm_coef, R = 1000)
    # Calculate bootstrap standard errors
    fitg3boot_se <- apply(fitg3results$t, 2, sd)
    Zvalue = fitg3coef/fitg3boot_se
    Pvalue = 2 * (1 - pnorm(abs(fitg3coef/fitg3boot_se)))
    estall3[III,1:6] = c(fitg3coef,fitg3boot_se,Pvalue)
    colnames(estall3)<-c('intercept','beta1','se_intercept','se_beta1','p_intercept','p_beta1')
    
    #poisson sandwich
    # Fit a Poisson GLM
    model <- glm(Y~X, family = poisson(), data = datas)
    estall5[III,1:6]  = as.numeric(coeftest(model, vcov. = vcovHC(model, type = "HC3"))[,c(1,2,4)])
    colnames(estall5) = c('intercept','beta1','se_intercept','se_beta1','p_intercept','p_beta1')
    
  }
  
  linda.obj  <- linda(feature.dat = otu.tab.sim,meta.dat = meta.dat,formula = '~X',feature.dat.type = 'count',zero.handling = c('pseudo-count'),alpha=0.05)
  
  return(list(zico.obj,estall,estall2,estall3,estall5,linda.obj$output$X,ancombc2.obj$res))
}

MytestResult<-function(myfit,sim.obj,cutoff=0.625,pvalue=0.05){
  
  clean_error <- function(df){
    df = data.frame(df)
    # Identify columns containing the value 999
    rows_with_999 <- apply(df,1, function(x) any(x == 999))
    # Drop these columns
    df_cleaned <- df[ !rows_with_999 , ]
    return(df_cleaned)
  }
  
  zico.obj=myfit[[1]]
  estall=clean_error(myfit[[2]])
  estall2=clean_error(myfit[[3]])
  estall3=clean_error(myfit[[4]])
  estall5=clean_error(myfit[[5]])
  linda.obj=myfit[[6]]
  ancombc2.obj=myfit[[7]]
  
  common.taxa=rowMeans(sim.obj$otu.tab.sim!= 0)>cutoff
  rare.taxa=!common.taxa
  dim(sim.obj$otu.tab.sim[common.taxa,])
  sim.obj.common=sim.obj
  sim.obj.common$otu.tab.sim=sim.obj$otu.tab.sim[common.taxa,]
  sim.obj.common$diff.otu.ind=sim.obj$diff.otu.ind[common.taxa]
  sim.obj.common$otu.names=sim.obj$otu.names[common.taxa]
  
  sim.obj.rare=sim.obj
  sim.obj.rare$otu.tab.sim=sim.obj$otu.tab.sim[rare.taxa,]
  sim.obj.rare$diff.otu.ind=sim.obj$diff.otu.ind[rare.taxa]
  sim.obj.rare$otu.names=sim.obj$otu.names[rare.taxa]
  
  
  ziconames=colnames(zico.obj$coef.list[[1]])
  ziconames.common=ziconames[common.taxa]
  ziconames.rare=ziconames[rare.taxa]
  diffnames=colnames(zico.obj$coef.list[[1]])[sim.obj$diff.otu.ind]
  diffnames.common=colnames(zico.obj$coef.list[[1]])[sim.obj$diff.otu.ind&common.taxa]
  diffnames.rare=colnames(zico.obj$coef.list[[1]])[sim.obj$diff.otu.ind&rare.taxa]
  
  if(sum(zico.obj$p.adj.fdr<pvalue)==0){zicoseqfdr=0
  }else{
    zicoseqfdr=1-sum(ziconames[zico.obj$p.adj.fdr<pvalue]%in%diffnames)/sum(zico.obj$p.adj.fdr<pvalue)
  }
  zicoseqfdr
  
  if(sum(zico.obj$p.adj.fdr[common.taxa]<pvalue)==0){zicoseqfdr.common=0
  }else{
    zicoseqfdr.common=1-sum(ziconames.common[zico.obj$p.adj.fdr[common.taxa]<pvalue]%in%diffnames.common)/sum(zico.obj$p.adj.fdr[common.taxa]<pvalue)
  }
  zicoseqfdr.common
  
  if(sum(zico.obj$p.adj.fdr[rare.taxa]<pvalue)==0){zicoseqfdr.rare=0
  }else{
    zicoseqfdr.rare=1-sum(ziconames.rare[zico.obj$p.adj.fdr[rare.taxa]<pvalue]%in%diffnames.rare)/sum(zico.obj$p.adj.fdr[rare.taxa]<pvalue)
  }
  zicoseqfdr.rare

  zicoseqtpr=sum(ziconames[zico.obj$p.adj.fdr<pvalue]%in%diffnames)/length(diffnames)
  zicoseqtpr
  zicoseqtpr.common=sum(ziconames.common[zico.obj$p.adj.fdr[common.taxa]<pvalue]%in%diffnames.common)/length(diffnames.common)
  zicoseqtpr.common
  zicoseqtpr.rare=sum(ziconames.rare[zico.obj$p.adj.fdr[rare.taxa]<pvalue]%in%diffnames.rare)/length(diffnames.rare)
  zicoseqtpr.rare
  
  if(sum(p.adjust(estall[,'p_beta1'],method = 'fdr')<pvalue)==0){nbfdr=0
  }else{
  nbfdr=1-sum(ziconames[p.adjust(estall[,'p_beta1'],method = 'fdr')<pvalue]%in%diffnames)/sum(p.adjust(estall[,'p_beta1'],method = 'fdr')<pvalue)
  }
  nbfdr
  if(sum(p.adjust(estall[,'p_beta1'],method = 'fdr')[common.taxa]<pvalue)==0){nbfdr.common=0
  }else{
  nbfdr.common=1-sum(ziconames.common[p.adjust(estall[,'p_beta1'],method = 'fdr')[common.taxa]<pvalue]%in%diffnames.common)/sum(p.adjust(estall[,'p_beta1'],method = 'fdr')[common.taxa]<pvalue)
  }
  nbfdr.common
  if(sum(p.adjust(estall[,'p_beta1'],method = 'fdr')[rare.taxa]<pvalue)==0){nbfdr.rare=0
  }else{
  nbfdr.rare=1-sum(ziconames.rare[p.adjust(estall[,'p_beta1'],method = 'fdr')[rare.taxa]<pvalue]%in%diffnames.rare)/sum(p.adjust(estall[,'p_beta1'],method = 'fdr')[rare.taxa]<pvalue)
  }
  nbfdr.rare
  
  nbtpr=sum(ziconames[p.adjust(estall[,'p_beta1'],method = 'fdr')<pvalue]%in%diffnames)/length(diffnames)
  nbtpr
  nbtpr.common=sum(ziconames.common[p.adjust(estall[,'p_beta1'],method = 'fdr')[common.taxa]<pvalue]%in%diffnames.common)/length(diffnames.common)
  nbtpr.common
  nbtpr.rare=sum(ziconames.rare[p.adjust(estall[,'p_beta1'],method = 'fdr')[rare.taxa]<pvalue]%in%diffnames.rare)/length(diffnames.rare)
  nbtpr.rare
  
  if(sum(p.adjust(estall2[,'p_beta1'],method = 'fdr')<pvalue)==0){poifdr=0
  }else{
  poifdr=1-sum(ziconames[p.adjust(estall2[,'p_beta1'],method = 'fdr')<pvalue]%in%diffnames)/sum(p.adjust(estall2[,'p_beta1'],method = 'fdr')<pvalue)
  }
  poifdr
  if(sum(p.adjust(estall2[,'p_beta1'],method = 'fdr')[common.taxa]<pvalue)==0){poifdr.common=0
  }else{
  poifdr.common=1-sum(ziconames.common[p.adjust(estall2[,'p_beta1'],method = 'fdr')[common.taxa]<pvalue]%in%diffnames.common)/sum(p.adjust(estall2[,'p_beta1'],method = 'fdr')[common.taxa]<pvalue)
  }
  poifdr.common
  if(sum(p.adjust(estall2[,'p_beta1'],method = 'fdr')[rare.taxa]<pvalue)==0){poifdr.rare=0
  }else{
  poifdr.rare=1-sum(ziconames.rare[p.adjust(estall2[,'p_beta1'],method = 'fdr')[rare.taxa]<pvalue]%in%diffnames.rare)/sum(p.adjust(estall2[,'p_beta1'],method = 'fdr')[rare.taxa]<pvalue)
  }
  poifdr.rare
  
  poitpr=sum(ziconames[p.adjust(estall2[,'p_beta1'],method = 'fdr')<pvalue]%in%diffnames)/length(diffnames)
  poitpr
  poitpr.common=sum(ziconames.common[p.adjust(estall2[,'p_beta1'],method = 'fdr')[common.taxa]<pvalue]%in%diffnames.common)/length(diffnames.common)
  poitpr.common
  poitpr.rare=sum(ziconames.rare[p.adjust(estall2[,'p_beta1'],method = 'fdr')[rare.taxa]<pvalue]%in%diffnames.rare)/length(diffnames.rare)
  poitpr.rare

####  
  if(sum(p.adjust(estall3[,'p_beta1'],method = 'fdr')<pvalue)==0){bootfdr=0
  }else{
    bootfdr=1-sum(ziconames[p.adjust(estall3[,'p_beta1'],method = 'fdr')<pvalue]%in%diffnames)/sum(p.adjust(estall3[,'p_beta1'],method = 'fdr')<pvalue)
  }
  bootfdr
  if(sum(p.adjust(estall3[,'p_beta1'],method = 'fdr')[common.taxa]<pvalue)==0){bootfdr.common=0
  }else{
    bootfdr.common=1-sum(ziconames.common[p.adjust(estall3[,'p_beta1'],method = 'fdr')[common.taxa]<pvalue]%in%diffnames.common)/sum(p.adjust(estall3[,'p_beta1'],method = 'fdr')[common.taxa]<pvalue)
  }
  bootfdr.common
  if(sum(p.adjust(estall3[,'p_beta1'],method = 'fdr')[rare.taxa]<pvalue)==0){bootfdr.rare=0
  }else{
    bootfdr.rare=1-sum(ziconames.rare[p.adjust(estall3[,'p_beta1'],method = 'fdr')[rare.taxa]<pvalue]%in%diffnames.rare)/sum(p.adjust(estall3[,'p_beta1'],method = 'fdr')[rare.taxa]<pvalue)
  }
  bootfdr.rare
  
  boottpr=sum(ziconames[p.adjust(estall3[,'p_beta1'],method = 'fdr')<pvalue]%in%diffnames)/length(diffnames)
  boottpr
  boottpr.common=sum(ziconames.common[p.adjust(estall3[,'p_beta1'],method = 'fdr')[common.taxa]<pvalue]%in%diffnames.common)/length(diffnames.common)
  boottpr.common
  boottpr.rare=sum(ziconames.rare[p.adjust(estall3[,'p_beta1'],method = 'fdr')[rare.taxa]<pvalue]%in%diffnames.rare)/length(diffnames.rare)
  boottpr.rare
###
  if(sum(p.adjust(estall5[,'p_beta1'],method = 'fdr')<pvalue)==0){sandwichfdr=0
  }else{
    sandwichfdr=1-sum(ziconames[p.adjust(estall5[,'p_beta1'],method = 'fdr')<pvalue]%in%diffnames)/sum(p.adjust(estall5[,'p_beta1'],method = 'fdr')<pvalue)
  }
  sandwichfdr
  if(sum(p.adjust(estall5[,'p_beta1'],method = 'fdr')[common.taxa]<pvalue)==0){sandwichfdr.common=0
  }else{
    sandwichfdr.common=1-sum(ziconames.common[p.adjust(estall5[,'p_beta1'],method = 'fdr')[common.taxa]<pvalue]%in%diffnames.common)/sum(p.adjust(estall5[,'p_beta1'],method = 'fdr')[common.taxa]<pvalue)
  }
  sandwichfdr.common
  if(sum(p.adjust(estall5[,'p_beta1'],method = 'fdr')[rare.taxa]<pvalue)==0){sandwichfdr.rare=0
  }else{
    sandwichfdr.rare=1-sum(ziconames.rare[p.adjust(estall5[,'p_beta1'],method = 'fdr')[rare.taxa]<pvalue]%in%diffnames.rare)/sum(p.adjust(estall5[,'p_beta1'],method = 'fdr')[rare.taxa]<pvalue)
  }
  sandwichfdr.rare
  
  sandwichtpr=sum(ziconames[p.adjust(estall5[,'p_beta1'],method = 'fdr')<pvalue]%in%diffnames)/length(diffnames)
  sandwichtpr
  sandwichtpr.common=sum(ziconames.common[p.adjust(estall5[,'p_beta1'],method = 'fdr')[common.taxa]<pvalue]%in%diffnames.common)/length(diffnames.common)
  sandwichtpr.common
  sandwichtpr.rare=sum(ziconames.rare[p.adjust(estall5[,'p_beta1'],method = 'fdr')[rare.taxa]<pvalue]%in%diffnames.rare)/length(diffnames.rare)
  sandwichtpr.rare  
###    
  ###    
  if(sum(linda.obj$reject)==0){lindafdr=0
  }else{
    lindafdr=1-sum(ziconames[linda.obj$reject]%in%diffnames)/sum(linda.obj$reject)
  }
  lindafdr
  
  if(sum(linda.obj$reject[common.taxa])==0){lindafdr.common=0
  }else{
    lindafdr.common=1-sum(ziconames.common[linda.obj$reject[common.taxa]]%in%diffnames.common)/sum(linda.obj$reject[common.taxa])
  }
  lindafdr.common
  
  if(sum(linda.obj$reject[rare.taxa])==0){lindafdr.rare=0
  }else{
    lindafdr.rare=1-sum(ziconames.rare[linda.obj$reject[rare.taxa]]%in%diffnames.rare)/sum(linda.obj$reject[rare.taxa])
  }
  lindafdr.rare
  
  
  lindatpr=sum(ziconames[linda.obj$reject]%in%diffnames)/length(diffnames)
  lindatpr
  
  lindatpr.common=sum(ziconames.common[linda.obj$reject[common.taxa]]%in%diffnames.common)/length(diffnames.common)
  lindatpr.common
  
  lindatpr.rare=sum(ziconames.rare[linda.obj$reject[rare.taxa]]%in%diffnames.rare)/length(diffnames.rare)
  lindatpr.rare
  
  
  ###  
  if(sum(ancombc2.obj$diff_X)==0){ancombc2fdr=0
  }else{
    ancombc2fdr=1-sum(ziconames[ancombc2.obj$diff_X]%in%diffnames)/sum(ancombc2.obj$diff_X)
  }
  ancombc2fdr
  
  if(sum(ancombc2.obj$diff_X[common.taxa])==0){ancombc2fdr.common=0
  }else{
    ancombc2fdr.common=1-sum(ziconames.common[ancombc2.obj$diff_X[common.taxa]]%in%diffnames.common)/sum(ancombc2.obj$diff_X[common.taxa])
  }
  ancombc2fdr.common
  
  if(sum(ancombc2.obj$diff_X[rare.taxa])==0){ancombc2fdr.rare=0
  }else{
    ancombc2fdr.rare=1-sum(ziconames.rare[ancombc2.obj$diff_X[rare.taxa]]%in%diffnames.rare)/sum(ancombc2.obj$diff_X[rare.taxa])
  }
  ancombc2fdr.rare
  
  
  ancombc2tpr=sum(ziconames[ancombc2.obj$diff_X]%in%diffnames)/length(diffnames)
  ancombc2tpr
  
  ancombc2tpr.common=sum(ziconames.common[ancombc2.obj$diff_X[common.taxa]]%in%diffnames.common)/length(diffnames.common)
  ancombc2tpr.common
  
  ancombc2tpr.rare=sum(ziconames.rare[ancombc2.obj$diff_X[rare.taxa]]%in%diffnames.rare)/length(diffnames.rare)
  ancombc2tpr.rare
  ###  
  fdrvalue=c(nbfdr,poifdr,zicoseqfdr,bootfdr,sandwichfdr,lindafdr,ancombc2fdr)
  fdrvalue.common=c(nbfdr.common,poifdr.common,zicoseqfdr.common,bootfdr.common,sandwichfdr.common,lindafdr.common,ancombc2fdr.common)
  fdrvalue.rare=c(nbfdr.rare,poifdr.rare,zicoseqfdr.rare,bootfdr.rare,sandwichfdr.rare,lindafdr.rare,ancombc2fdr.rare)
  tprvalue=c(nbtpr,poitpr,zicoseqtpr,boottpr,sandwichtpr,lindatpr,ancombc2tpr,length(common.taxa))
  tprvalue.common=c(nbtpr.common,poitpr.common,zicoseqtpr.common,boottpr.common,sandwichtpr.common,lindatpr.common,ancombc2tpr.common,length(diffnames.common))
  tprvalue.rare=c(nbtpr.rare,poitpr.rare,zicoseqtpr.rare,boottpr.rare,sandwichtpr.rare,lindatpr.rare,ancombc2tpr.rare,length(diffnames.rare))
  
  return(list(fdrvalue,tprvalue,fdrvalue.common,tprvalue.common,fdrvalue.rare,tprvalue.rare))
}

glm_coef <- function(data, indices,formula=Y~X) {
  fit <- glm(formula, family = poisson(), data = data[indices,])
  return(coef(fit))
}
# Simulate binary covariate, 10% signal density, abundant differential OTUs, unbalanced change 
# This setting simulates strong compositional effects 
model.paras <- EstPara(ref.otu.tab = com)
fdr=c()
tpr=c()
fdr.common=c()
tpr.common=c()
fdr.rare=c()
tpr.rare=c()
#set the differential percentage as 30% for high signal density
diffpct=0.3
suppressWarnings({
for (i in 1:100) {
  print(i)
  set.seed(i)
  sim.obj <- SimulateMSeqUni( ref.otu.tab = com,model.paras = model.paras,
                              nSam =400, nOTU = 100, 
                              # True signal setting 
                              diff.otu.pct = diffpct, diff.otu.direct = c("balanced"), 
                              diff.otu.mode = c("mix"), covariate.type = c("continuous"),
                              covariate.eff.mean = 0.2,
                              covariate.eff.sd = 0, 
                              confounder.type = 'continuous',conf.cov.cor=0,confounder.eff.mean = 0,
                              confounder.eff.sd = 0,conf.diff.otu.pct = 0,conf.nondiff.otu.pct = 0,
                              # Depth setting 
                              depth.mu = 10000, depth.theta = 5,depth.conf.factor = 0) 
  
  
  
  res1=try({myfit=MYFIT(sim.obj)},silent = F)
  errorj=i
  if(inherits(res1,"try-error")==TRUE){
    repeat {
      errorj=errorj+100
      set.seed(errorj)
      print(errorj)
      sim.obj <- SimulateMSeqUni( ref.otu.tab = com,model.paras = model.paras, nSam =400, nOTU = 100, 
                                  # True signal setting 
                                  diff.otu.pct = diffpct, diff.otu.direct = c("balanced"), 
                                  diff.otu.mode = c("mix"), covariate.type = c("continuous"),
                                  covariate.eff.mean = 0.2,
                                  covariate.eff.sd = 0, 
                                  confounder.type = 'continuous',conf.cov.cor=0,confounder.eff.mean = 0,
                                  confounder.eff.sd = 0,conf.diff.otu.pct = 0,conf.nondiff.otu.pct = 0,
                                  # Depth setting 
                                  depth.mu = 10000, depth.theta = 5,depth.conf.factor = 0) 
      res2=try({myfit=MYFIT(sim.obj) },silent = F)
      if (inherits(res2,"try-error")==FALSE) break}
  }
  MTR=MytestResult(myfit,sim.obj,cutoff = 0.625)
  fdrvalue=MTR[[1]]
  tprvalue=MTR[[2]]
  fdrvalue.common=MTR[[3]]
  tprvalue.common=MTR[[4]]
  fdrvalue.rare=MTR[[5]]
  tprvalue.rare=MTR[[6]]
  
  method.name=c('Negative Binomial','Poisson','ZicoSeq','bootstrap','sandwich','linda','ancombc2')
  fdr=rbind(fdr,fdrvalue)
  colnames(fdr)=method.name
  tpr=rbind(tpr,tprvalue)
  colnames(tpr)=c(method.name,'commontaxa')
  
  fdr.common=rbind(fdr.common,fdrvalue.common)
  colnames(fdr.common)=method.name
  tpr.common=rbind(tpr.common,tprvalue.common)
  colnames(tpr.common)=c(method.name,'diffcommon')
  
  fdr.rare=rbind(fdr.rare,fdrvalue.rare)
  colnames(fdr.rare)=method.name
  tpr.rare=rbind(tpr.rare,tprvalue.rare)
  colnames(tpr.rare)=c(method.name,'diffrare')
}
})
dim(fdr)
fdr
colMeans(fdr)
colMeans(tpr)
colMeans(fdr.rare)
colMeans(tpr.rare)
colMeans(tpr.rare[complete.cases(tpr.rare),])
fdr.rare <- ifelse(is.nan(fdr.rare), 0, fdr.rare)
colMeans(fdr.common)
colMeans(tpr.common)
tpr.rare <- ifelse(is.nan(tpr.rare), 0, tpr.rare)
tpr.rare
colMeans(tpr.rare)

######################################################
#draw the bar plot

OVERALL = cbind(colMeans(tpr),c(colMeans(fdr),0))
colnames(OVERALL)=c("TPR","FDR")
OVERALL
RARE = cbind(colMeans(tpr.rare),c(colMeans(fdr.rare),0))
colnames(RARE)=c("TPR","FDR")
RARE
COMMON = cbind(colMeans(tpr.common),c(colMeans(fdr.common),0))
colnames(COMMON)=c("TPR","FDR")
COMMON


library(ggplot2)
library(dplyr)
library(tidyr)
# Assuming 'data_overall' is your dataframe for the Overall group
data_overall <- data.frame(
  Method = rownames(OVERALL),
  OVERALL
)
data_overall=data_overall[-10,]
data_overall=data_overall[-c(4,9),]
data_overall

data_common <- data.frame(
  Method = rownames(COMMON),
  COMMON
)
data_common=data_common[-10,]
data_common=data_common[-c(4,9),]
data_common

data_rare <- data.frame(
  Method = rownames(RARE),
  RARE
)
data_rare=data_rare[-10,]
data_rare=data_rare[-c(4,9),]
data_rare

library(ggplot2)
library(dplyr)
library(tidyr)


# Assuming your data frames are already set up and combined
data_combined <- bind_rows(
  mutate(data_overall, Group = "Overall"),
  mutate(data_rare, Group = "Rare"),
  mutate(data_common, Group = "Common")
)

# Convert the 'Group' variable to a factor with the desired order
data_combined$Group <- factor(data_combined$Group, levels = c("Overall", "Rare", "Common"))

# Reshape the data to long format
data_long <- data_combined %>%
  pivot_longer(cols = c(TPR, FDR), names_to = "Metric", values_to = "Value") %>%
  mutate(Method = factor(Method, levels = unique(Method)))  # Order methods as they appear

# Plotting
myplot = ggplot(data_long, aes(x = Method, y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_text(aes(label = round(Value, 3)), vjust = -0.3, position = position_dodge(width = 0.7), size = 2) +  # Smaller size
  facet_grid(Group ~ Metric, switch = "y") +
  geom_hline(data = subset(data_long, Metric == "FDR"), aes(yintercept = 0.05), linetype = "dashed", color = "red") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", colour = "white"), # Set the background color to white
    panel.spacing=unit(1,"lines"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "Comparison of Methods by TPR and FDR across Groups", x = "Method", y = "Value")
myplot
# Save the plot using ggsave
# Plotting
myplot = ggplot(data_long, aes(x = Method, y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_text(aes(label = round(Value, 3)), vjust = -0.3, position = position_dodge(width = 0.7), size = 2) +  # Smaller size
  facet_grid(Group ~ Metric, switch = "y") +
  geom_hline(data = subset(data_long, Metric == "FDR"), aes(yintercept = 0.05), linetype = "dashed", color = "red") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", colour = "white"), # Set the background color to white
    panel.spacing=unit(1,"lines"),
    panel.border = element_rect(colour = "black", fill = NA),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1,size = 8),
    strip.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "Comparison of Methods by TPR and FDR across Groups", x = "Method", y = "Value")
myplot

ggsave(filename = "Barplothighdensity.jpg", plot = myplot, width = 9, height = 8, dpi = 300)




