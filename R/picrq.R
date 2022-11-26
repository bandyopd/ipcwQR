#' @importFrom stats as.formula binomial lm predict sd
NULL
#' Fit the partly interval-censored AFT model with quantile regressions
#'
#' Fit inverse weighted quantile regression with partially interval-censored data
#'
#' @param L left-censoring time, having 0 if left-censored.
#' @param R right-censoring time, having \code{Inf} if right-censored.
#' @param T exactly observed time.
#' @param delta censoring indicator, 1: observed; 0: interval-censored.
#' @param x X matrix of baseline covariates.
#' @param tau quantile level.
#' @param estimation estimating method of partly interval censored, if estimation="dr", doubly robust estimator is estimated.
#' @param wttype weight estimating method, default is "param".
#' @param hlimit bandwidth value, default is 0.5.
#' @param id cluster id. If the data does not have clustered structure, set \code{id=NULL}.
#' @param k index of cluster weight.
#' @param maxit maximum number of iteration for the log-rank estimator, default is 100.
#' @param tol tolerance of iteration for the log-rank estimator, default is 1e-3.
#'
#' @return \code{picrq} returns a data frame containing at least the following components:
#' \itemize{
#'   \item \code{tau}: quantile level.
#'   \item \code{coefficients}: regression estimator.
#'   \item \code{se}: standard error estimates for \code{est}.
#'   \item \code{pvalue}: p-value.
#'   \item \code{lower bd}: lower bound of coefficients under 95% confidence level.
#'   \item \code{upper bd}: upper bound of coefficients under 95% confidence level.
#' }
#'
#' @details
#' see Kim et al., (2022+) for detailed method explanation.
#'
#' @references
#' Chiou, S. H., Kang, S. and Yan, J. (2015). Rank-based estimating equations with general weight for accelerated failure time models: an induced smoothing approach. Statistics in Medicine 34(9): 1495–-1510.
#' 
#' Kim, Y., Choi, T., Park, S., Choi, S. and Bandyopadhyay, D. (2022+). Inverse weighted quantile regression with partially interval-censored data.
#' 
#' Pan, C. (2021). PICBayes: Bayesian Models for Partly Interval-Censored Data. R package. https://CRAN.R-project.org/package=PICBayes.
#'
#'
#' @examples
#' \dontrun{
#' # Simulations
#' set.seed(111)
#' n = 200
#' x1 = runif(n,-1,1)
#' x2 = rbinom(n,1,0.43)
#' x = cbind(x1,x2)
#' T = 2 + x1 + x2 + rnorm(n)
#' U = (1 - 0.25*x1)*runif(n, -6, 5)
#' V = U + (1 - 0.1*x2)*runif(n, 6, 20)
#' U = exp(dplyr::case_when(TRUE ~ T, T>V ~ V, T<U ~ -Inf))
#' V = exp(dplyr::case_when(TRUE ~ T, T>V ~ Inf, T<U ~ U))
#' delta = ifelse(U==V, 1, 0)
#' tau=0.3
#' picrq(L=V,R=U,delta=delta,x=x,tau=tau)
#' picrq(L=V,R=U,delta=delta,x=x,tau=tau,estimation = "dr")
#' 
#' 
#' # Data example
#' library(PICBayes)
#' data("mCRC")
#' d = with(data.frame(mCRC), data.frame(U = ifelse(y==0,R,L),
#'                                       V = ifelse(y==2,L,R),
#'                                       # Cluster weighted data
#'                                       id=(rep(c(table(SITE)),c(table(SITE)))),
#'                                       # Treatment arm: 0 = FOLFIRI alone, 1 = Panitumumab + FOLFIRI.
#'                                       x1= case_when(TRT_C == 0 ~ 0, #Pan et al data
#'                                                     TRT_C == 1 ~ 1),
#'                                       # Tumor KRAS mutation status: 0 = wild-type, 1 = mutant.
#'                                       x2= case_when(KRAS_C == 0 ~ 1,
#'                                                     KRAS_C == 1 ~ 0),
#'                                       delta = case_when(IC == 0 ~ 1,
#'                                                         IC == 1 ~ 0)
#'));
#'L=d$U;R=d$V; delta=d$delta
#'L=(log(d$U));R=log(d$V); delta=d$delta
#'x = cbind(d$x1,d$x2); id=d$id;  tau=0.1;
#'picrq(L,R,delta,x=x,tau=tau)
#'picrq(L,R,delta,x=x,tau=tau,hlimit=0.9)
#'picrq(L,R,delta,x=x,tau=tau,estimation = "dr")
#'picrq(L,R,delta,x=x,tau=tau,id=id)
#' }
#' @export
#'
#'
#'


# library(tidyverse)


picrq=function(L,R,delta,x,tau,estimation=NULL,wttype="param",hlimit=0.5,id=NULL,k=1,maxit=100,tol=1e-3){
  
  
  library(survival)
  library(tidyverse)
  library(extRemes)
  wtfunc2=function(L,R,delta){
    
    Y=pmax(ifelse(delta==0,R,L),1e-8); n=length(Y)
    
    kml = survfit(Surv(L) ~ 1, d)
    kmr = survfit(Surv(R) ~ 1, d)
    ww = rep(0,n)
    
    for (i in 1:n) {
      if(delta[i]==1){
        sl=approx(c(0,kml$time,100),c(1,kml$surv,0), xout=L[i], method = "linear", ties = mean)$y
        sr=approx(c(0,kmr$time,100),c(1,kmr$surv,0), xout=R[i], method = "linear", ties = mean)$y
        ww[i] = 1/pmax(1-(sr-sl),0.001)
      }
    }
    ww
  }
  
  
  
  # tau=0.3; hlimit=0.9
  Berwtfunc = function(L,R,delta,x, hlimit=NULL) {
    library(survival)
    Y=pmax(ifelse(delta==0,R,L),1e-8); n=length(Y)
    ker = dnorm(outer(x[,1],x[,1],"-")/hlimit)
    Wnj = ker / rowSums(ker)
    sr = sl = srl= 0
    denomr = rowSums(outer(Y,Y,">=")*(Wnj))
    denoml = rowSums(outer(Y,Y,"<=")*(Wnj))
    for (i in 1:length(Y)) {
      y0 = Y[i]
      etar = 1*(Y>=y0 & delta==0)
      etal = 1*(Y<=y0 & delta==0)
      nom = Wnj[,i]
      sr = prod((1 - nom/denomr)^etar)
      sl = 1-prod((1 - nom/denoml)^etal)
      srl[i] = 1-(sr-sl)
    }
    del = ifelse(delta==1,1,0)
    del/pmax( srl, 1e-4)
  }
  
  
  Lwtfunc=function(L,R,delta){
    
    Y=pmax(ifelse(delta==0,R,L),1e-8); n=length(Y)
    
    kml = survfit(Surv(L) ~ 1, d)
    kmr = survfit(Surv(R) ~ 1, d)
    ww = rep(0,n)
    
    for (i in 1:n) {
      if(delta[i]==1){
        sl=approx(c(0,kml$time,100),c(1,kml$surv,0), xout=L[i], method = "linear", ties = mean)$y
        sr=approx(c(0,kmr$time,100),c(1,kmr$surv,0), xout=R[i], method = "linear", ties = mean)$y
        ww[i] = 1/pmax((1+sl),0.001)
      }
    }
    ww
  }
  
  
  Rwtfunc=function(L,R,delta){
    
    
    Y=pmax(ifelse(delta==0,R,L),1e-8); n=length(Y)
    
    kml = survfit(Surv(L) ~ 1, d)
    kmr = survfit(Surv(R) ~ 1, d)
    ww = rep(0,n)
    
    for (i in 1:n) {
      if(delta[i]==1){
        sl=approx(c(0,kml$time,100),c(1,kml$surv,0), xout=L[i], method = "linear", ties = mean)$y
        sr=approx(c(0,kmr$time,100),c(1,kmr$surv,0), xout=R[i], method = "linear", ties = mean)$y
        ww[i] = 1/pmax((sr),0.001)
      }
    }
    ww
  }
  
  # eta=1
  PICrq=function(L,R,delta,x,ww){
    
    Y=pmax(ifelse(delta==0,R,L),1e-8); n=length(Y)
    quantreg::rq((Y)~x, weights = eta*ww, tau = tau)$coef #intc, beta1, beta2
  }
  
  # beta=c(1,1,1); Sigma=diag(3)
  Efunc=function(L,R,delta,x,Sigma,beta,tau,ww){
    Y=pmax(ifelse(delta==0,R,L),1e-8); n=length(Y)
    xx = as.matrix(cbind(1,x)); p = ncol(xx)
    ss = sqrt( pmax(1e-3, diag(xx%*%Sigma%*%t(xx))) ) 
    res = as.numeric(Y - xx%*%beta)
    Phi = as.vector( pnorm( -res/ss ) * ww )
    
    U = as.vector( t(xx)%*%(Phi - tau) )
    U/n
  }
  
  # wr=Rwtfunc(L,R,delta);wl=Lwtfunc(L,R,delta);
  DREfunc=function(L,R,delta,x,Sigma,beta,tau,ww,wl,wr){
    Y=pmax(ifelse(delta==0,R,L),1e-8); n=length(Y)
    n=length(Y); 
    ndelta=sum(ifelse(delta==1,1,0))
    ndelta=ifelse(ndelta==0,1,ndelta)
    xx = as.matrix(cbind(1,x)); p = ncol(xx)
    ss = sqrt( pmax(1e-3, diag(xx%*%Sigma%*%t(xx))) ) 
    res = as.numeric(Y - xx%*%beta)
    ind = ifelse(res<=0,1,0)
    wp=ww*pnorm( -(res/ss) )
    U = as.vector( t(xx*eta)%*%(as.numeric(wp-tau)) )
    UR=matrix(0,p,1); UL=matrix(0,p,1);
    
    for(i in 1:n){
      if(delta[i]==1){
        yindr=Y>=Y[i]
        denomr=sum(yindr)
        Rft = as.vector( t(xx*eta*(wr))%*%(as.numeric((yindr*ind-tau)*ind)) )
        UR=UR+((Rft/denomr))
        
        yindl=Y<=Y[i]
        denoml=sum((1-yindl))+1
        Lft = as.vector( t(xx*eta*(wl))%*%(as.numeric(((1-yindl)*ind-tau)*ind)) )
        UL=UL+((Lft/denoml))
      }
    }
    (U-(UR/ndelta)+(UL/ndelta))/n
  }
  
  
  Afunc=function(L,R,delta,x,Sigma,beta,ww){
    Y=pmax(ifelse(delta==0,R,L),1e-8); n=length(Y)
    xx = as.matrix(cbind(1,x)); p = ncol(xx)
    ss = sqrt( pmax(1e-3, diag(xx%*%Sigma%*%t(xx))) ) 
    res = as.numeric(Y - xx%*%beta)
    Phi = as.vector( pnorm( -res/ss ) * ww )
    A = t(Phi * xx) %*% xx + diag(p)*0.05
    A/n
  }
  
  Gfunc=function(L,R,delta,x,Sigma,beta,tau,ww){
    
    Y=pmax(ifelse(delta==0,R,L),1e-8); n=length(Y)
    xx = as.matrix(cbind(1,x)); p = ncol(xx)
    ss = sqrt( pmax(1e-3, diag(xx%*%Sigma%*%t(xx))) ) 
    res = as.numeric(Y - xx%*%beta)
    ind = ifelse(res<=0,1,0)
    resind = ww*ind
    
    Gam=xx*as.numeric(resind-tau)
    Gamma=( t(Gam)%*%(Gam*eta) )
    GammaR=matrix(0,p,p)
    GammaL=matrix(0,p,p)
    
    for(i in 1:n){
      if(delta[i]==0){
        # yind=Y<=Y[i] | Y>=Y[i]
        yind=Y>=Y[i]
        denom=sum(yind)
        num1=as.vector( t(xx*eta)%*% (yind*resind) )
        num2=as.vector( t(xx)%*% (yind*resind) )
        
        Rft1 = num1/denom
        Rft2 = num2/denom
        GammaR=GammaR+(Rft1)%*%t(Rft2)
        
        yind2=Y<=Y[i]
        denom2=sum(yind2)
        num3=as.vector( t(xx*eta)%*% (yind2*resind) )
        num4=as.vector( t(xx)%*% (yind2*resind) )
        
        Lft1 = num3/denom2
        Lft2 = num4/denom2
        GammaL=GammaL+(Lft1)%*%t(Lft2)
      }
    }
    (Gamma-GammaR+GammaL)/n
  }
  
  # update variance estimator
  up_Sigma = function(Y,Afunc,Gfunc){
    n=length(Y)
    invA = solve(Afunc)
    newSigma = (( (invA) %*% Gfunc %*% (invA) ) )
    newSigma/n
  }
  
  Y=pmax(ifelse(delta==0,R,L),1e-8); n=length(Y)
  if(is.null(id)){eta=1}
  else{ci=rep(c(table(id)),c(table(id))); wi=(1/ci); eta=(wi^k)}
  
  if(wttype=="param"){ww=wtfunc2(L,R,delta);}
  if(wttype=="nonparam" && is.null(hlimit)){print("hlimit should be entered.")}
  if(wttype=="nonparam"){ww=Berwtfunc(L,R,delta,x,hlimit);}
  xx = as.matrix(cbind(1,x)); p = ncol(xx)
  # wttype="param"; eta=1; maxit=100; tol=10; estimation=NULL
  # wttype="nonparam"; eta=1; maxit=100; tol=10; estimation="dr"
  old_beta = init = beta = PICrq(L,R,delta,x,ww=ww)
  old_Sigma = Sigma = diag(p)/n
  
  
  i=0; eps=1; max.iter=100; tol = 1e-3; 
  while (i<max.iter & eps >= tol ) {
    Amat = Afunc(L,R,delta,x,beta=c(old_beta),ww=ww,Sigma = old_Sigma)
    if(is.null(estimation)){
      new_beta = c(old_beta) - solve(Amat)%*%Efunc(L,R,delta,x,ww=ww,beta=c(old_beta),Sigma = old_Sigma,tau)/n}
    else if(estimation=="dr"){wr=Rwtfunc(L,R,delta);wl=Lwtfunc(L,R,delta);
    new_beta = c(old_beta) - solve(Amat)%*%DREfunc(L,R,delta,x,ww=ww,wl=wl,wr=wr,beta=c(old_beta),Sigma = old_Sigma,tau)/n}
    # new_beta = BB::dfsane(par=old_beta,fn=U_n,Y=Y,x=x,ww=ww,Sigma=old_Sigma,tau=tau,control=list(trace=FALSE))$par
    Gamma = Gfunc(L,R,delta,x,beta=c(old_beta),Sigma = old_Sigma,ww=ww, tau=tau)
    new_Sigma = up_Sigma(Y,Amat,Gamma)
    
    if (det(new_Sigma) <= 0) {
      new_beta = old_beta; new_Sigma = old_Sigma
    }
    
    eps = max(max(abs(new_beta - old_beta)),
              max(abs(new_Sigma - old_Sigma)))
    old_beta = new_beta; old_Sigma = new_Sigma; 
    i = i+1
    
  }
  se = sqrt(diag(new_Sigma))
  res=data.frame(tau=c(rep(tau,p)),
                 est=new_beta,se=se,
                 pvalue = 1 - pnorm(abs(new_beta/se)),
                 lb = new_beta-1.96*se, ub = new_beta+1.96*se)
  colnames(res)=c("tau","coefficients","se","pvalue","lower bd","upper bd")
  rownames(res)[1]="Intercept"
  round((res), 6)
}

