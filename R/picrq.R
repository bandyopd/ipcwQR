#' @importFrom stats as.formula binomial lm predict sd
NULL
#' Fit the partly interval-censored AFT model with quantile regressions
#'
#' Fit inverse weighted quantile regression with partially interval-censored data
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
#' @return \code{picrq} returns a data frame containing at least the following components:
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
#' Kim, Y., Choi, T., Pak, S., Choi, S. and Bandyopadhyay, D. (2022+). Inverse weighted quantile regression with partially interval-censored data.
#' 
#' Pan, C. (2021). PICBayes: Bayesian Models for Partly Interval-Censored Data. R package. https://CRAN.R-project.org/package=PICBayes.
#'
#'
#' @examples
#' \dontrun{
#' # Simulations
#' # library(tidyverse)
#' library(ipcwQR)
#' set.seed(111)
#' x1 = rnorm(n)
#' x2 = rbinom(n,1,0.5)
#' x = cbind(x1,x2)
#' T = 2 + x1 + x2 + rnorm(n)
#' U = (1 - 0.25*x1)*runif(n, -6, 3.85)
#' V = U + (1 - 0.5*x2)*runif(n, 6, 19.6)
#' U = exp(dplyr::case_when(TRUE ~ T, T>V ~ V, T<U ~ -Inf))
#' V = exp(dplyr::case_when(TRUE ~ T, T>V ~ Inf, T<U ~ U))
#' delta = ifelse(U==V, 1, 0); tau=0.3
#' picrq(V,U,delta,x,tau)
#' picrq(V,U,delta,x,tau,estimation = "dr")
#' 
#' 
#' # Data example
#' library(PICBayes)
#' library(ipcwQR)
#' data(mCRC)
#' d=data.frame(mCRC)
#' y=d$y; U = ifelse(y==0,R,L); V = ifelse(y==2,L,R)
#' delta = case_when(d$IC == 0 ~ 1, d$IC == 1 ~ 0)
#' x1= case_when(d$TRT_C == 0 ~ 0, d$TRT_C == 1 ~ 1)
#' x2= case_when(d$KRAS_C == 0 ~ 1, d$KRAS_C == 1 ~ 0)
#' x = cbind(x1,x2);  
#' id = rep(c(table(d$SITE)),c(table(d$SITE))); tau=0.3
#' picrq(V,U,delta,x,tau)
#' picrq(V,U,delta,x,tau,wttype = "nonparam",h=0.9)
#' picrq(V,U,delta,x,tau,wttype = "nonparam",id=id,h=0.9)
#' 
#' 
#' d = with(dt0,
#'          data.frame(U = ifelse(is.na(L), 0, L),
#'                     V = ifelse(is.na(R), Inf, R),
#'                     Delta = 1-IC,
#'                     x1 = TRT_C,
#'                     x2 = KRAS_C,
#'                     id = SITE))
#' U = d$U; V = d$V; X = cbind(d$x1,d$x2); Delta = d$Delta; id = d$id
#' aft_rank(U = U, V = V, X = X, Delta = Delta, id = id, 
#'          alpha = 1, type = "gehan", nboot = 10)
#' aft_rank(U = U, V = V, X = X, Delta = Delta, id = id, 
#'          alpha = 1, type = "logrank", nboot = 10)
#' }
#' @export
#'
#'
#'


# library(tidyverse)


picrq=function(L,R,delta,x,tau,estimation=NULL,wttype="param",h=0.5,id=NULL,k=1,maxit=100,tol=1e-3){
  
  
  wtfunc2=function(L,R,delta){
    library(survival)
    Y = pmax(L,1e-8); L = pmax(L,1e-8); R = (R)
    n=length(Y);
    kml = survfit(Surv(L) ~ 1)
    kmr = survfit(Surv(R) ~ 1)
    ww = rep(0,n)
    
    for (i in 1:n) {
      if(delta[i]==1){
        sl=approx(c(0,kml$time,100),c(1,kml$surv,0), xout=L[i])$y
        sr=approx(c(0,kmr$time,100),c(1,kmr$surv,0), xout=R[i])$y
        ww[i] = 1/pmax(1-(sr-sl),0.001)
      }
    }
    ww
  }
  
  
  Berwtfunc = function(L,R,delta,x, h=NULL) {
    library(survival)
    Y = pmax(L,1e-8); delta = delta; y = Y;  n = length(Y)
    ker = dnorm(outer(x[,1],x[,1],"-")/h)
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
  
  
  Lwtfunc=function(L,R,T,delta){
    Y = pmax(L,1e-8); delta = delta; n = length(Y)
    n=length(Y);
    kml = survfit(Surv(-Y,delta==3)~1)
    kmr = survfit(Surv(Y,delta==2)~1)
    ww = rep(0,n)
    
    
    for (i in 1:n) {
      if(delta[i]==1){
        sl = approx( c(0, -kml$time, 100), c(1,1-kml$surv,0), xout=Y[i])$y
        sr = approx( c(0, kmr$time, 100), c(1, kmr$surv, 0), xout=Y[i])$y
        ww[i] = 1/pmax( sl, 0.5)
      }
    }
    ww
  }
  
  
  Rwtfunc=function(L,R,delta){
    
    Y = pmax(L,1e-8); delta = delta; n = length(Y)
    n=length(Y);
    kml = survfit(Surv(-Y,delta==3)~1)
    kmr = survfit(Surv(Y,delta==2)~1)
    ww = rep(0,n)
    
    
    for (i in 1:n) {
      if(delta[i]==1){
        sl = approx( c(0, -kml$time, 100), c(1,1-kml$surv,0), xout=Y[i])$y
        sr = approx( c(0, kmr$time, 100), c(1, kmr$surv, 0), xout=Y[i])$y
        ww[i] = 1/pmax( sr, 1e-3)
      }
    }
    ww
  }
  
  PICrq=function(L,R,delta,x,tau){
    Y = pmax(L,1e-8); L = pmax(L,1e-8); R = (R); 
    n=length(Y);
    if(wttype=="param"){ww=wtfunc2(L,R,delta)}
    if(wttype=="nonparam" && is.null(h)){print("h should be entered.")}
    if(wttype=="nonparam"){ww=Berwtfunc(L,R,delta,x,h)}
    quantreg::rq((Y)~x, weights = ww, tau = tau)$coef #intc, beta1, beta2
  }
  
  Efunc=function(L,R,delta,x,Sigma,beta,tau,ww){
    Y = pmax(L,1e-8); L = pmax(L,1e-8); R = (R); 
    n=length(Y);
    xx = as.matrix(cbind(1,x)); p = ncol(xx)
    ss = sqrt( pmax(1e-3, diag(xx%*%Sigma%*%t(xx))) ) 
    res = as.numeric(Y - xx%*%beta)
    Phi = as.vector( pnorm( -res/ss ) * ww )
    U = as.vector( t(xx*eta)%*%(Phi - tau) )
    U/sqrt((length(Y)))
  }
  
  
  DREfunc=function(L,R,delta,x,Sigma,beta,tau,ww,wl,wr){
    Y = pmax(L,1e-8); L = pmax(L,1e-8); R = (R); 
    n=length(Y); ndelta=sum(ifelse(delta==1,0,1))
    xx = as.matrix(cbind(1,x)); p = ncol(xx)
    ss = sqrt( pmax(1e-3, diag(xx%*%Sigma%*%t(xx))) ) 
    res = as.numeric(Y - xx%*%beta)
    ind = ifelse(res<=0,1,0)
    wp=ww*pnorm( -(res/ss) )
    U = as.vector( t(xx)%*%(as.numeric(wp-tau)) )
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
    (U-(UR/ndelta)+(UL/ndelta))/sqrt((length(Y)))
  }
  
  Afunc=function(L,R,delta,x,Sigma,beta,ww){
    Y = pmax(L,1e-8); L = pmax(L,1e-8); R = (R); 
    n=length(Y);
    xx = as.matrix(cbind(1,x)); p = ncol(xx)
    ss = sqrt( pmax(1e-3, diag(xx%*%Sigma%*%t(xx))) ) 
    res = as.numeric(Y - xx%*%beta)
    phi = as.vector( dnorm( -res/ss )* (ww/ss) )
    A = t(phi * xx) %*% xx + diag(p)*0.05
    A
  }
  
  Gfunc=function(L,R,delta,x,Sigma,beta,tau,ww){
    
    Y = pmax(L,1e-8); L = pmax(L,1e-8); R = (R);  
    n=length(Y);
    xx = as.matrix(cbind(1,x)); p = ncol(xx)
    ss = sqrt( pmax(1e-3, diag(xx%*%Sigma%*%t(xx))) ) 
    res = as.numeric(Y - xx%*%beta)
    ind = ifelse(res<=0,1,0)
    resind = ww*ind
    
    Gam=xx*as.numeric(resind-tau)
    Gamma=( t(Gam)%*%(Gam) )
    GammaR=matrix(0,p,p)
    GammaL=matrix(0,p,p)
    
    for(i in 1:n){
      if(delta[i]==0){
        yindr=Y>=Y[i]
        denom=sum(yindr)
        nom1=as.vector( t(xx*eta)%*% (yindr*resind) )
        nom2=as.vector( t(xx)%*% (yindr*resind) )
        R1 = nom1/denom; R2 = nom2/denom
        GammaR=GammaR+(R1)%*%t(R2)
        
        yindl=Y<=Y[i]
        denom=sum((1-yindl))+1
        nom1=as.vector( t(xx*eta)%*% ((1-yindl)*resind) )
        nom2=as.vector( t(xx)%*% ((1-yindl)*resind) )
        L1 = nom1/denom; L2 = nom2/denom
        GammaL=GammaL+(L1)%*%t(L2)
      }
    }
    (Gamma-GammaR+GammaL)/sqrt((length(Y)))
  }
  
  # update variance estimator
  up_Sigma = function(Y,Afunc,Gfunc){
    n=length(Y)
    invA = solve(Afunc)
    newSigma = (( (invA) %*% Gfunc %*% (invA) ) )
    newSigma
  }
  
  
  Y = pmax(L,1e-8); L = pmax(L,1e-8); R = (R);  
  n=length(Y);
  if(is.null(id)){eta=1}
  else{ci=rep(c(table(id)),c(table(id))); eta=(1/(ci^(k)))}
  
  if(wttype=="param"){ww=wtfunc2(L,R,delta);wr=Rwtfunc(L,R,T,delta)}
  if(wttype=="nonparam" && is.null(h)){print("h should be entered.")}
  if(wttype=="nonparam"){ww=Berwtfunc(L,R,delta,x,h);wr=Rwtfunc(L,R,T,delta)}
  xx = as.matrix(cbind(1,x)); p = ncol(xx)
  # wttype="param"; eta=1; max.iter=100; tol=10
  old_beta = init = beta = PICrq(L,R,delta,x,tau)
  old_Sigma = Sigma = diag(p)/n
  
  i=0; eps=1; max.iter=maxit; tol = tol;
  while (i<max.iter & eps >= tol ) {
    Amat = Afunc(L,R,x=x,beta=c(old_beta),ww=ww,Sigma = old_Sigma)
    if(is.null(estimation)){new_beta = c(old_beta) - solve(Amat)%*%Efunc(L,R,delta,x=x,ww=ww,beta=c(old_beta),Sigma = old_Sigma,tau)/((length(Y)))}
    else if(estimation=="dr"){wr=Rwtfunc(L,R,T,delta);wl=Lwtfunc(L,R,T,delta);
    new_beta = c(old_beta) - solve(Amat)%*%DREfunc(L,R,delta,x=x,ww=ww,wl=wl,wr=wr,beta=c(old_beta),Sigma = old_Sigma,tau)/((length(Y)))}
    Gamma = Gfunc(L,R,delta=delta,x=x,beta=c(old_beta),Sigma = old_Sigma,ww=ww, tau=tau)
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


