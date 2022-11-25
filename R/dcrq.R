#' @importFrom stats as.formula binomial lm predict sd
NULL
#' Fit the doubly interval-censored AFT model with quantile regression
#'
#' Fit inverse weighted quantile regression with doubly interval-censored data
#'
#'
#' @param L left-censoring time, having 0 if left-censored.
#' @param R right-censoring time, having \code{Inf} if right-censored.
#' @param T exactly observed time.
#' @param delta censoring indicator, 1: observed; 0: interval-censored.
#' @param x baseline covariates.
#' @param tau quantile level.
#' @param estimation estimating method of partly interval censored, if estimation="dr", doubly robust estimator is estimated.
#' @param wttype weight estimating method, default is "param".
#' @param h bandwidth value, default is 0.5.
#' @param id cluster id. If the data does not have clustered structure, set \code{id=NULL}.
#' @param k index of cluster weight.
#' @param maxit maximum number of iteration for the log-rank estimator, default is 100.
#' @param tol tolerance of iteration for the log-rank estimator, default is 1e-3.
#'
#' @return \code{dcrq} returns a data frame containing at least the following components:
#' \itemize{
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
#' rm(list=ls())
#' library(tidyverse)
#' set.seed(111)
#' cluster=60
#' v=rlnorm(cluster,mean=0,sd=0.5)
#' k = sample(0:9*10,cluster,replace = TRUE)
#' m = ifelse(quantile(v, k/100) <= v & quantile(v, (k+10)/100) <= v,
#'            k/10 + 2, 2)
#' id = as.numeric(unlist(apply(cbind(1:cluster,m),1,function(a) rep(a[1],a[2]))))
#' n=sum(m)
#' x1 = rnorm(n); x2 = rbinom(n,1,0.5); x = cbind(x1,x2)
#' T = exp(2 + x1 + x2 + rnorm(n))
#' L=runif(n,-2.8,1.7); R=L+runif(n,4.2,10); 
#' L=exp(L); R=exp(R)
#' Y=pmin(R,pmax(T,L))
#' delta=case_when(T<L ~ 3, T>R ~ 2, TRUE ~ 1) #observed
#' tau=0.3
#' dcrq(L,R,T,delta,x,tau)
#' }
#' @export
#'
#'
#'


# library(tidyverse)

dcrq=function(L,R,T,delta,x,tau,estimation=NULL,wttype="param",h=0.5,id=NULL,k=1,maxit=100,tol=1e-3){
  
  
  
  wtfunc=function(L,R,T,delta){
    library(survival)
    Y=pmin(R,pmax(T,L))
    n=length(Y);
    kml = survfit(Surv(-Y,delta==3)~1)
    kmr = survfit(Surv(Y,delta==2)~1)
    ww = rep(0,n)
    
    for (i in 1:n) {
      if(delta[i]==1){
        sl = approx( c(0, -kml$time, 100), c(1,1-kml$surv,0), xout=Y[i], method = "linear", ties = mean)$y
        sr = approx( c(0, kmr$time, 100), c(1, kmr$surv, 0), xout=Y[i], method = "linear", ties = mean)$y
        ww[i] = 1/pmax( sr-sl, 1e-3)
      }
    }
    ww
  }
  
  # h=0.9
  Berwtfunc = function(L,R,T,delta,x, h=NULL) {
    library(survival)
    Y=pmin(R,pmax(T,L)); y = Y;  n = length(Y)
    kml = survfit(Surv(-Y,delta==3)~1)
    kmr = survfit(Surv(Y,delta==2)~1)
    ker = dnorm(outer(x[,1],x[,1],"-")/h)
    Wnj = ker / rowSums(ker)
    sr = sl = srl= rep(0,n)
    denomr = rowSums(outer(y,y,"<=")*(Wnj))
    denoml = rowSums(outer(y,y,">=")*(Wnj))
    for (i in 1:length(y)) {
      if(delta[i]==1){
        y0 = y[i]
        etar = 1*(y<=y0 & delta==2)
        etal = 1*(y>=y0 & delta==3)
        nom = Wnj[,i]
        sr = prod((1 - nom/denomr)^etar)
        sl = 1-prod((1 - nom/denoml)^etal)
        srl[i] = 1/pmax( sr-sl, 1e-3)
      }
    }
    srl
  }
  
  
  Lwtfunc=function(L,R,T,delta){
    Y=pmin(R,pmax(T,L))
    n=length(Y);
    kml = survfit(Surv(-Y,delta==3)~1)
    kmr = survfit(Surv(Y,delta==2)~1)
    ww = rep(0,n)
    
    
    for (i in 1:n) {
      if(delta[i]==1){
        sl = approx( c(0, -kml$time, 100), c(1,1-kml$surv,0), xout=Y[i], method = "linear", ties = mean)$y
        sr = approx( c(0, kmr$time, 100), c(1, kmr$surv, 0), xout=Y[i], method = "linear", ties = mean)$y
        ww[i] = 1/pmax( sl, 0.5)
      }
    }
    ww
  }
  
  
  Rwtfunc=function(L,R,T,delta){
    
    Y=pmin(R,pmax(T,L)); Y=log(Y)
    n=length(Y);
    kml = survfit(Surv(-Y,delta==3)~1)
    kmr = survfit(Surv(Y,delta==2)~1)
    ww = rep(0,n)
    
    
    for (i in 1:n) {
      if(delta[i]==1){
        sl = approx( c(0, -kml$time, 100), c(1,1-kml$surv,0), xout=Y[i], method = "linear", ties = mean)$y
        sr = approx( c(0, kmr$time, 100), c(1, kmr$surv, 0), xout=Y[i], method = "linear", ties = mean)$y
        ww[i] = 1/pmax( sr, 1e-3)
      }
    }
    ww
  }
  
  
  # eta=1
  DCrq=function(L,R,T,delta,x,tau,eta){
    Y=pmin(R,pmax(T,L));  Y=log(Y)
    n=length(Y);
    ww=wtfunc(L,R,T,delta)
    quantreg::rq((Y)~x, weights = ww*eta, tau = tau)$coef #intc, beta1, beta2
  }
  
  # Sigma=diag(p)/n; beta=rep(1,3)
  Efunc=function(L,R,T,delta,x,Sigma,beta,tau,ww){
    Y=pmin(R,pmax(T,L)); Y=log(Y)
    n=length(Y);
    xx = as.matrix(cbind(1,x)); p = ncol(xx)
    ss = sqrt( pmax(1e-3, diag(xx%*%Sigma%*%t(xx))) ) 
    res = as.numeric(Y - xx%*%beta)
    Phi = as.vector( pnorm( -res/ss ) * ww )
    U = as.vector( t(xx)%*%(Phi - tau) )
    U/(sqrt(n))
  }
  
  
  DREfunc=function(L,R,T,delta,x,Sigma,beta,tau,ww,wl,wr){
    Y=pmin(R,pmax(T,L));  Y=log(Y)
    n=length(Y); nr=table(delta)[2]; nl=table(delta)[3];
    xx = as.matrix(cbind(1,x)); p = ncol(xx)
    ss = sqrt( pmax(1e-3, diag(xx%*%Sigma%*%t(xx))) ) 
    res = as.numeric(Y - xx%*%beta)
    ind = ifelse(res<=0,1,0)
    wp=ww*pnorm( -(res/ss) )
    U = as.vector( t(xx)%*%(as.numeric(wp-tau)) )
    UR=matrix(0,p,1); UL=matrix(0,p,1);
    
    for(i in 1:n){
      if(delta[i]==2){
        yindr=Y>=Y[i]
        denomr=sum(yindr)
        Rft = as.vector( t(xx*eta*(wr))%*%(as.numeric((yindr*ind-tau)*ind)) )
        UR=UR+((Rft/denomr))
      }
      if(delta[i]==3){  
        yindl=Y<=Y[i]
        denoml=sum((1-yindl))+1
        Lft = as.vector( t(xx*eta*(wl))%*%(as.numeric(((1-yindl)*ind-tau)*ind)) )
        UL=UL+((Lft/denoml))
      }
    }
    (U-(UR)+(UL))/(sqrt(n))
  }
  
  Afunc=function(L,R,T,delta,x,Sigma,beta,ww){
    Y=pmin(R,pmax(T,L));  Y=log(Y)
    n=length(Y);
    xx = as.matrix(cbind(1,x)); p = ncol(xx)
    ss = sqrt( pmax(1e-3, diag(xx%*%Sigma%*%t(xx))) ) 
    res = as.numeric(Y - xx%*%beta)
    phi = as.vector( dnorm( -res/ss )* (ww/ss) )
    A = t(phi * xx) %*% xx + diag(p)*0.05
    A*(sqrt(n))
  }
  
  Gfunc=function(L,R,T,delta,x,Sigma,beta,tau,ww,eta){
    
    Y=pmin(R,pmax(T,L));  Y=log(Y)
    n=length(Y);
    xx = as.matrix(cbind(1,x)); p = ncol(xx)
    ss = sqrt( pmax(1e-3, diag(xx%*%Sigma%*%t(xx))) ) 
    res = as.numeric(Y - xx%*%beta)
    ind = ifelse(res<=0,1,0)
    resind = ww*ind
    
    Gam=xx*as.numeric(resind-tau)
    Gamma=( t(eta*Gam)%*%(Gam) )
    GammaR=matrix(0,p,p)
    GammaL=matrix(0,p,p)
    
    for(i in 1:n){
      if(delta[i]==2){
        yindr=Y>=Y[i]
        denom=sum(yindr)
        nom=as.vector( t(xx)%*% (yindr*resind) )
        Rft = nom/denom
        GammaR=GammaR+(Rft)%*%t(Rft)
      }
      if(delta[i]==3){
        yindl=Y<=Y[i]
        denom=sum((1-yindl))+1
        nom=as.vector( t(xx)%*% ((1-yindl)*resind) )
        Lft = nom/denom
        GammaL=GammaL+(Lft)%*%t(Lft)
      }
    }
    (Gamma-GammaR+GammaL)*(sqrt(n))
  }
  
  # update variance estimator
  up_Sigma = function(Y,Afunc,Gfunc){
    n=length(Y)
    invA = solve(Afunc)
    newSigma = (( (invA) %*% Gfunc %*% (invA) ) )
    newSigma
  }
  
  
  Y=pmin(R,pmax(T,L));  Y=log(Y)
  n=length(Y);
  
  if(is.null(id)){eta=1}
  else{ci=rep(c(table(id)),c(table(id))); eta=(1/(ci^(k)))}
  
  if(wttype=="param"){ww=wtfunc(L,R,T,delta);}
  if(wttype=="nonparam" && is.null(h)){print("h should be entered.")}
  if(wttype=="nonparam"){ww=Berwtfunc(L,R,delta,x,h);}
  # wttype="param"; eta=1; maxit=100; tol=10
  xx = as.matrix(cbind(1,x)); p = ncol(xx)
  old_beta = init = beta = DCrq(L,R,T,delta,x,tau,eta)
  old_Sigma = Sigma = diag(p)/n
  
  
  
  i=0; eps=1; max.iter=maxit; tol = tol;
  while (i<max.iter & eps >= tol ) {
  Amat = Afunc(L,R,T,x=x,beta=c(old_beta),ww=ww,Sigma = old_Sigma)
  
  if(is.null(estimation)){new_beta = c(old_beta) - solve(Amat)%*%Efunc(L,R,T,delta,x=x,ww=ww,beta=c(old_beta),Sigma = old_Sigma,tau=tau)/n}
  else if(estimation=="dr"){wr=Rwtfunc(L,R,T,delta);wl=Lwtfunc(L,R,T,delta);
  new_beta = c(old_beta) - solve(Amat)%*%DREfunc(L,R,T,delta,x=x,ww=ww,wl=wl,wr=wr,beta=c(old_beta),Sigma = old_Sigma,tau=tau)/n}
  Gamma = Gfunc(L,R,T,delta=delta,x=x,beta=c(old_beta),Sigma = old_Sigma,ww=ww, tau=tau,eta=eta)
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
  res=data.frame(tau=c(tau,tau,tau),
                 est=new_beta,se=se,
                 pvalue = 1 - pnorm(abs(new_beta/se)),
                 lb = new_beta-1.96*se, ub = new_beta+1.96*se)
  colnames(res)=c("tau","coefficients","se","pvalue","lower bd","upper bd")
  rownames(res)[1]="Intercept"
  round((res), 6)
}

