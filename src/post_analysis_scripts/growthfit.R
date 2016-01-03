# -------------------------------------------------
# growthfit.R
# Erin Landguth
# Description: Reads in sensitivity runs produced by 
# diagnostics.py GrowthRate.csv that are rows of results
# to compare, y(t) data. Data is then fit using the 
# package grofit and Richards/or logistic growth 
# equations to produce a lambda value. 
# --------------------------------------------------

library(grofit)
library(nls)

library(MASS)
library(mgcv)

# ------------------------
# User data
# ------------------------
file <- 'D:/projects/CDmetaPOP/Seattle/Runs/data373_20140817/Scenario1_MortEggsPatch_1408304994/MatureN.csv'
gen <- 20
batchno <- 5

# -----------------------
# Main work
# -----------------------

# Extract time in grofit format
timepoints <- 1:gen
timepoints <- t(matrix(rep(timepoints, batchno), c(gen, batchno)))

# Read in file
data <- read.table(file, header = FALSE, sep = ",")

# Extract data in grofit format


# ---------------------------------------
# Testing Ellner et al Gradient matching
# ---------------------------------------

x <- 1:gen
y <- as.numeric(as.matrix(data[2,1:gen+1]))

# fit & plot with penalty=2
fit2<-gcv.rss(x,y,hmin=-5,hmax=5,intercept=T,penalty=2,draw=T);

# fit & don't plot with penalty=1
fit1<-gcv.rss(x,y,hmin=-6,hmax=1,intercept=T,penalty=1,draw=T);


library(mgcv); par(mfcol=c(2,2)); 

#load data and estimate the gradient 
x<-as.matrix(read.table("c:\\GRADFuncs\\NichAdults.txt")); 
ad<-x/1000; tvals<-2*(c(1:length(ad))-1); 

dxhat<-bcgrad(tvals,ad,hder=2); 
plot(tvals,1000*ad,type="p",xlab="Time (d)",ylab="Adults"); 
points(tvals,1000*dxhat$smooth,type="l",lty=1); title("(a) Time series & smooth"); 
zhat<-dxhat$unbiased; zraw<-dxhat$raw; 

# plot gradient bias correction curve 
plot(1000*zraw,1000*zhat,type="p",xlab="Raw gradient estimate",ylab="Bias corrected estimate");
ez<-order(zraw); 
points(1000*zraw[ez],1000*zraw[ez],type="l",lty=1);
title(" (b) Gradient bias correction"); 

## Create data vectors: xt1=current density, xt0=lagged density, dadt=current gradient
xt1<-trimr(ad,9,2); xt0<-trimr(ad,2,9); dadt<-trimr(zhat,9,2); 

## Unconstrained regression spline fit using gam() from mgcv library
gamfit<-gam(dadt~s(xt1,20)+s(xt0,20));
newd1<-data.frame(xt1=(0:120)/20,xt0=0*(0:120)/20);
f1<- -predict.gam(gamfit,newd1); f1<-f1-f1[1];
newd2<-data.frame(xt1=0*(0:120)/20,xt0=(0:120)/20);
f0<-predict.gam(gamfit,newd2); f0<-f0-f0[1];
matplot(1000*(0:120)/20, 1000*(f1%~%f0),type="l",lty=1,xlab="Number of adults",
        ylab="Fitted B and D rates"); title("(c) Unconstrained GAM fit");


############## Example of gcv.rssM 
x<-2*runif(250); x<-x[order(x)]; ybar<-(x^20)/(1+x^20); 
y<-ybar+0.5*rnorm(250);
fitM<-gcv.rssM(x,y,hmin=-12,hmax=0,intercept=F,penalty=1,draw=T); 
points(x,ybar,type="l",lty=2,lwd=2);

# fit & plot the same data without monotonicity constraint
fit<-gcv.rss(x,y,hmin=-12,hmax=0,intercept=F,penalty=1,draw=F); 
points(x,fit$fitted.values,type="l",lty=4,lwd=1); 
legend(0,2.2,c("true","monotonic","unconstrained"),lty=c(2,1,4),lwd=c(2,1,1),cex=0.75); 


x<-runif(100)*3-0.1; x<-x[order(x)];
ybar<-5*(exp(-x)-4*exp(-2*x)+3*exp(-3*x)); 
y<-ybar+0.3*rnorm(100);

# fit & plot with penalty=2
fit2<-gcv.rss(x,y,hmin=-5,hmax=5,intercept=T,penalty=2,draw=T);
points(x,ybar,type="l",lty=2,lwd=2);

# fit & don't plot with penalty=1
fit1<-gcv.rss(x,y,hmin=-6,hmax=1,intercept=T,penalty=1,draw=T); 

# predict & plot new values from the model with penalty=1
newx<-range(x)[1]+0.01*c(0:100)*(range(x)[2]-range(x)[1]); 
yhat1<-predict.rss(fit1$model,newx);
points(yhat1$x,yhat1$y,type="l",lty=4,lwd=1); 
legend(0.25,1,c("true","penalty=2","penalty=1"),lty=c(2,1,4),lwd=c(2,1,1),cex=0.75); 

##############

############## Example of gcv.rssM 
x<-2*runif(250); x<-x[order(x)]; ybar<-(x^20)/(1+x^20); 
y<-ybar+0.5*rnorm(250);
fitM<-gcv.rssM(x,y,hmin=-12,hmax=0,intercept=F,penalty=1,draw=T); 
points(x,ybar,type="l",lty=2,lwd=2);

# fit & plot the same data without monotonicity constraint
fit<-gcv.rss(x,y,hmin=-12,hmax=0,intercept=F,penalty=1,draw=F); 
points(x,fit$fitted.values,type="l",lty=4,lwd=1); 
legend(0,2.2,c("true","monotonic","unconstrained"),lty=c(2,1,4),lwd=c(2,1,1),cex=0.75); 

##############


################# Example of bcgrad and gcv.rssadd2 on Nicholson data
library(mgcv); par(mfcol=c(2,2)); 

#load data and estimate the gradient 
x<-as.matrix(read.table("c:\\GRADFuncs\\NichAdults.txt")); 
ad<-x/1000; tvals<-2*(c(1:length(ad))-1); 

dxhat<-bcgrad(tvals,ad,hder=2); 
plot(tvals,1000*ad,type="p",xlab="Time (d)",ylab="Adults"); 
points(tvals,1000*dxhat$smooth,type="l",lty=1); title("(a) Time series & smooth"); 
zhat<-dxhat$unbiased; zraw<-dxhat$raw; 

# plot gradient bias correction curve 
plot(1000*zraw,1000*zhat,type="p",xlab="Raw gradient estimate",ylab="Bias corrected estimate");
ez<-order(zraw); 
points(1000*zraw[ez],1000*zraw[ez],type="l",lty=1);
title(" (b) Gradient bias correction"); 

## Create data vectors: xt1=current density, xt0=lagged density, dadt=current gradient
xt1<-trimr(ad,9,2); xt0<-trimr(ad,2,9); dadt<-trimr(zhat,9,2); 

## Unconstrained regression spline fit using gam() from mgcv library
gamfit<-gam(dadt~s(xt1,20)+s(xt0,20));
newd1<-data.frame(xt1=(0:120)/20,xt0=0*(0:120)/20);
f1<- -predict.gam(gamfit,newd1); f1<-f1-f1[1];
newd2<-data.frame(xt1=0*(0:120)/20,xt0=(0:120)/20);
f0<-predict.gam(gamfit,newd2); f0<-f0-f0[1];
matplot(1000*(0:120)/20, 1000*(f1%~%f0),type="l",lty=1,xlab="Number of adults",
        ylab="Fitted B and D rates"); title("(c) Unconstrained GAM fit");


## Sign-constrained regression spline fit using gcv.rssadd2
## Note: gam() from mgcv is MUCH FASTER because we use a general-purpose
## optimizer (R's optim() with method="Nelder-Mead") for minimizing GCV(alpha). 
##  mgcv is built for speed, this set of functions was built speedily. 
rssfit<-gcv.rssadd2(x=xt1%~%xt0,y=dadt,Boundary.knots1=c(0,range(xt1)[2]),
                    Boundary.knots2=c(0,range(xt0)[2]), signs=c(-1,1),intercept1=F,
                    intercept2=F,penalty=2);

f1c<- -1000*rssfit$fitted.values[,1];
f0c<- 1000*rssfit$fitted.values[,2]; 
e0<-order(xt0); 
plot(1000*c(0,xt0[e0]),c(0,f0c[e0]),type="l",lty=1,xlab="Number of adults x 1000",
     ylab="Fitted B and D rates",xlim=c(0,6000)); title("(d) Sign-constrained fit");
e1<-order(xt1); 
points(1000*c(0,xt1[e1]),c(0,f1c[e1]),type="l",lty=1)


## SIMEX example ########################################################################
## This uses some of the 'utility' routines that are not described in the documentation,
## but the comments in the source code should be sufficient. 
##  rssfit() is similar to gcv.rss except that the smoothing parameter is user-supplied
##	rssbasis() generates the spline basis functions at a set of points, given the knots
##	  and the left endpoint of the interval ("origin" in the call). 

# the "data"; x and ybar are the unobserved truth, xobs and y are the error-corrupted observations.
x<-5*runif(200)-2.5; x<-x[order(x)]; 
ybar<-exp(-(x^2)); y<-ybar+0.01*rnorm(200); 
sigme<-0.35;
xobs<-x+sigme*rnorm(200);    
plot(xobs,y,xlim=c(-3,3), ylim=c(0,1.2),xlab="Observed x values",ylab="Observed y values"); 
points(x,ybar,type="l",lty=1,lwd=2); 

# fit using the observed data
rawfit<-gcv.rss(xobs,y,hmin=-6,hmax=6,intercept=T,penalty=1,draw=F);
exobs<-order(xobs); 
points(xobs[exobs],rawfit$fitted.values[exobs],type="l",lty=2,lwd=2); 

# "freeze" the smoothing parameter at the optimum for the observed data
hopt<-rawfit$hopt;

# extract other parameters of the fit to the observed data
nknots<-length(rawfit$model$knots);
Boundary.knots<-rawfit$model$Boundary.knots;

# refit 1000 times with one added dose of measurment error
fitted2<-0*rawfit$model$coef; fitted3<-fitted2;
for (irep in 1:1000) {
  ei<-sigme*rnorm(200);
  x2<-xobs+ei;
  fit2<-rssfit(x2,y,nknots=nknots,Boundary.knots=Boundary.knots,alpha=hopt,intercept=T); 
  fitted2<-fitted2+fit2$coef;
  x3<-xobs+sqrt(2)*ei; 
  fit3<-rssfit(x3,y,nknots=nknots,Boundary.knots=Boundary.knots,alpha=hopt,intercept=T); 
  fitted3<-fitted3+fit3$coef;
  cat(c(irep, "\n")); 
}
fitted2<-fitted2/1000; fitted3<-fitted3/1000; 

SimexCoef<-6*rawfit$model$coef + 4*fitted3 - 9*fitted2;

xplot<- xobs[order(xobs)]; knots<- rawfit$model$knots; origin<-rawfit$model$origin; 
X<-rssbasis(xplot,knots=knots,origin=origin,intercept=T); 
points(xplot,X%*%SimexCoef,type="l",lty=4,lwd=3); 
legend(-3,1,c("True","Raw fit","SIMEX"),lty=c(1,2,4),lwd=c(2,2,3)); 
#####################################################################



#### ----------------------     Introduction 
## R functions for 
## (0) making R understand a bit of GAUSS. Some of this code was ported from GAUSS.
## (1) Estimating the gradient dx/dt for a time series x(t) with observations
##     at times t0<t1<t2<... 
## (2) Fitting penalized regression splines of several types:
##     univariate, univariate monotonic increasing, and additive with sign
##     constraints. Fitting criterion is GCV with user-specified penalty parameter 
##     as in smooth.spline (the default, penalty=1, gives the usual GCV criterion)

## To run under Splus, uncomment the following line: 
##   library(Matrix); solve<-solve.Matrix;

## Not everything runs under Splus, though adapting them to do so should 
## not be too difficult. The main R-only feature is the use of 
## optim() with method = "Nelder-Mead" to fit the two smoothing
## parameters for the ridge functions in additive models; a call to
## to nlmin() could be substituted, but only if you're willing to trust 
## results from nlmin. 


#### -------------------------   CODE STARTS HERE

library(quadprog); 
library(Mass); 

## Utility routines betraying that this code was ported from GAUSS.

"%~%"<-function(a,b) {cbind(a,b)}
"%|%"<-function(a,b) {rbind(a,b)}
minindc<-function(x) {match(range(x)[1],x) };

stack<-function(xx, yy){
  nx <- length(xx)
  ny <- length(yy)
  z <- matrix(0, nx + ny, 1)
  z[1:nx] <- xx[1:nx]
  z[(nx + 1):(nx + ny)] <- yy[1:ny]
  return(z)
}

trimr<-function (a,n1,n2) {
  da<-dim(a); 
  if(is.null(da)) {a[(n1+1):(length(a)-n2)]}
  else {a[(n1+1):(da[1]-n2),]};
};

zeros<-function(nrow,ncol) {matrix(0,nrow,ncol)};
ones<-function(nrow,ncol) {matrix(1,nrow,ncol)};
rndn<-function(nrow,ncol) {matrix(rnorm(nrow*ncol),nrow,ncol)};
rndu<-function(nrow,ncol) {matrix(runif(nrow*ncol),nrow,ncol)};
rows<-function(a) {
  da<-dim(a)
  if(is.null(da)) {length(a)} else {da[1]}
}
cols<-function(a) {
  da<-dim(a)
  if(is.null(da)) {1} else {da[2]}
}

xy<-function(x,y,...) {
  e<-order(x);
  plot(x[e],y[e],...);
}

xypoints<-function(x,y,...) {
  e<-order(x);
  points(x[e],y[e],...);
}

## Routines for gradient estimation via local polynomial regression, 
## Gaussian kernel truncated at (-10*h,10*h) 

# Utility routine for locregG. Does the fit at a single point.		
locreg1G<-function(xvals,yvals,x,h,degree=2) {
  x0<-(xvals-x);
  e<-(abs(x0)<=10*h);
  yvals<-yvals[e]; x0<-x0[e];
  X<-matrix(0,length(x0),degree); 
  for (i in 1:degree) {
    X[,i]<-x0^i;
  }
  w<-exp(-((x0/h)^2));
  betahat<-lsfit(X,yvals,wt=w,intercept=T)$coef;
  return(as.matrix(betahat));
}

# locregG: given data (x,y), produce a fitted value at a vector
# of locations xnew by local regression smoothing with bandwidth h 
# and selected degree (default 2).
# RETURNS a vector ysmu whose ith column is the ith-degree coefficient
# in the local polynomial fit centered at the data points. In particular
# ysmu[,1] is the fitted value, ysmu[,2] is the estimated derivative.
locregG<-function(x,y,xnew,h,degree=2) {
  ysmu<-matrix(0,length(xnew),degree+1); 
  for (i in 1:length(xnew)) {
    xi<-xnew[i]
    bi<-locreg1G(xvals=x,yvals=y,x=xi,h=h,degree=degree)
    ysmu[i,]<-matrix(bi,1,degree+1)
  }
  return(ysmu)
}


## Routines for fitting penalized cubic regression splines by GCV
## "Penalty" is the cost per effective model df, as in smooth.spline(); 

# The "plus" regression spline basis
# x are the valuation points, knots are the knots, origin is the left
# endpoint of the interval corresponding to Boundary.knots[1] in the
# functions that use the basis, and intercept is logical for inclusion
# of a constant term in the basis
rssbasis<-function(x, knots, origin, intercept) {
  xs<-(x-origin);
  p<-cbind(xs,xs^2,xs^3); 
  if(intercept) p<-cbind( matrix(1,length(x),1),p); 
  nknots<-length(knots); 
  for (i in 1:nknots) {
    si<-pmax(x-knots[i],0);
    p<-cbind(p,si^3);
  }
  return(as.matrix(p))
}

# Fit Penalized Regression Spline with user-supplied smoothing parameter
#	NOTE alpha is the log10 of the smoothing parameter
#	Other parameters are the same as in gcv.rss 
rssfit<-function(x,y,nknots,Boundary.knots,intercept=T,alpha,penalty=1)	{
  origin<-Boundary.knots[1]; 
  fac<-(Boundary.knots[2]-Boundary.knots[1])/(nknots+1)
  knots<- origin+c(1:nknots)*fac; 
  X<-rssbasis(x,knots,origin,intercept);
  XTX<-t(X)%*%X;
  d<-rep(1,dim(XTX)[1]); 
  d[1:3]<-0; if(intercept) d[4]<-0; 
  D<-diag(d); 
  betahat<-solve(XTX+(10^alpha)*D,t(X)%*%y,tol=1.e-16);
  yhat<-X%*%betahat;
  
  msr<-mean((yhat-y)^2); 
  hat<-X%*%solve(XTX+(10^alpha)*D,tol=1.e-16)%*%t(X);
  d1<-sum(diag(hat))/length(y);
  gcv<-msr/((1-penalty*d1)^2);
  if(penalty*d1>1) gcv<-1.e12; 
  return(list(x=x,knots=knots,origin=origin,alpha=alpha,Boundary.knots=Boundary.knots,
              intercept=intercept, coef=betahat, fitted.values=yhat,gcv=gcv,enp=sum(diag(hat))))
}

predict.rss<-function(model,newx){
  X<-rssbasis(newx,model$knots,model$origin,model$intercept);
  yhat<-X%*%(model$coef);
  return(list(x=newx,y=yhat))
}


# Fit Penalized Regression Spline with smoothing parameter chosen by GCV(p)
# Hmin and Hmax are user-supplied brackets on log10(optimal smoothing parameter) 
# Uses golden section search so hmin and hmax MUST bracket a minimum or the
# routine will return garbage. Use draw=T (the default) to make sure you have a bracket. 
gcv.rss<-function(xvals,yvals,hmin,hmax,nknots=NA,Boundary.knots=NA,intercept=T,penalty=1,tol=0.01,draw=T) {
  if(is.na(Boundary.knots)) Boundary.knots<-range(xvals)
  if(is.na(nknots)) nknots<-20;
  C<-0.61803399; R<-1-C;
  dh<-hmax-hmin;
  hvals<-hmin+0.1*c(0:10)*dh;
  Gvals<-matrix(0,11,1);
  for (i in 1:11) {
    Gvals[i]<-rssfit(xvals,yvals,nknots,Boundary.knots,intercept=intercept,
                     hvals[i],penalty=penalty)$gcv;
    print(c(i,hvals[i],Gvals[i]));
  }
  imin<-minindc(Gvals);
  if((imin<=1)||(imin>=11)) {hopt<-hvals[imin]} else
  { 
    ax<-hvals[imin+1]; bx<-hvals[imin]; cx<-hvals[imin-1]; 
    x0<-ax; x3<-cx; 
    if(abs(cx-bx)>abs(bx-ax)) 
    {x1<-bx; x2<-bx+C*(cx-bx)} 
    else {x2<-bx; x1<-bx-C*(bx-ax)}
    f1<-rssfit(xvals,yvals,nknots,Boundary.knots,intercept=intercept,penalty=penalty,alpha=x1)$gcv;	
    f2<-rssfit(xvals,yvals,nknots,Boundary.knots,intercept=intercept,penalty=penalty,alpha=x2)$gcv;
    while(abs(x3-x0)>tol*(abs(x1)+abs(x2))) {
      if(f2<f1) {x0<-x1; x1<-x2; x2<-R*x1+C*x3; f1<-f2; 
                 f2<-rssfit(xvals,yvals,nknots,Boundary.knots,intercept=intercept,penalty=penalty,alpha=x2)$gcv;}
      else {x3<-x2; x2<-x1; x1<-R*x2+C*x0; f2<-f1;
            f1<-rssfit(xvals,yvals,nknots,Boundary.knots,intercept=intercept,penalty=penalty,alpha=x1)$gcv;}
    }		
    if(f1<f2) {hopt<-x1; fmin<-f1} else {hopt<-x2; fmin<-f2}
    print(c(hopt,fmin));
  }
  
  bestfit<-rssfit(xvals,yvals,nknots,Boundary.knots,intercept=intercept,penalty=penalty,alpha=hopt);	
  
  if (draw==T) {
    win.graph(); par(mfrow=c(2,1));
    plot(hvals,log10(Gvals),type="l",xlab="log10(h)",
         ylab="GCV");
    title("Bandwidth selection curve");
    ox<-order(xvals);
    plot(xvals,yvals); points(xvals[ox],bestfit$fitted.values[ox],type="l");
    title("Data and P-spline regression curve");
  }
  
  return(list(hvals=hvals,Gvals=Gvals,hopt=hopt,fitted.values=bestfit$fitted.values,
              penalty=penalty,gcv=bestfit$gcv,model=bestfit))
  
}

# Use Qprog to Fit Penalized Regression Spline with user-supplied 
# smoothing parameter, constrained to be monotonic at the knots
# NOTE x values are sorted into increasing order on return. If x is unsorted
# on input, plotting fitted values versus input x will give garbage. 
rssfitM<-function(x,y,nknots,Boundary.knots,intercept=T,penalty=1,alpha,ncongrid)	{
  ex<-order(x); x<-x[ex]; y<-y[ex]; 
  origin<-Boundary.knots[1]; 
  fac<-(Boundary.knots[2]-Boundary.knots[1])/(nknots+1)
  knots<- origin+c(1:nknots)*fac; 
  X<-rssbasis(x,knots,origin,intercept);
  XTX<-t(X)%*%X;
  
  d<-rep(1,dim(XTX)[1]); 
  d[1:3]<-0; if(intercept) d[4]<-0;
  D<-diag(d); 
  
  # Quadratic programming objective function
  qmat<-XTX+(10^alpha)*D; rmat<-as.vector(t(X)%*%y);  
  
  # Impose monotonicity at all constraining grid points
  cfac<-(Boundary.knots[2]-Boundary.knots[1])/ncongrid;
  aknots<- as.vector(Boundary.knots[1]+c(0:ncongrid)*cfac); 
  cmat<-rssbasis(aknots,knots,origin,intercept); 	
  nrc<-dim(cmat)[1];
  zmat<-cmat[2:nrc,]-cmat[1:(nrc-1),];
  qpfit<-solve.QP(Dmat=qmat,dvec=rmat,Amat=t(zmat),meq=0,factorized=F);
  betahat<-qpfit$solution;
  yhat<-X%*%betahat;
  msr<-mean((yhat-y)^2); 
  iact<-qpfit$iact; nact<-sum(iact>0);
  if(nact<1) {
    hat<-X%*%solve(XTX+(10^alpha)*D,tol=1.e-16)%*%t(X); Z<-1;
  }
  else {		
    iact<-iact[iact>0]; 
    Cact<-zmat[iact,];
    if(nact<2) Cact<-matrix(Cact,1,cols(zmat)) 
    Z<-Null(t(Cact));
    ZT<-t(Z);
    hat<-X%*%Z%*%solve(ZT%*%(XTX+(10^alpha)*D)%*%Z,tol=1.e-16)%*%ZT%*%t(X);
  }
  enp<-sum(diag(hat));
  d1<-sum(diag(hat))/length(y);
  gcv<-msr/((1-penalty*d1)^2);
  if(penalty*d1>1) gcv<-1.e12; 
  
  return(list(x=x,knots=knots,origin=origin,alpha=alpha,
              intercept=intercept, coef=betahat, fitted.values=yhat,gcv=gcv,
              penalty=penalty,Z=Z,enp=enp))
}

# Fit Monotone Penalized Regression Spline with smoothing parameter chosen by GCV(p)
# Hmin and Hmax are user-supplied brackets on the log10(optimal smoothing parameter) 
# Uses golden section search; hmin and hmax must bracket a minimum
gcv.rssM<-function(xvals,yvals,hmin,hmax,nknots=NA,Boundary.knots=NA,intercept=T,penalty=1,ncongrid=NA,tol=0.01,draw=T) {
  if(is.na(Boundary.knots)) Boundary.knots<-range(xvals);
  if(is.na(nknots)) nknots<-20;
  if(is.na(ncongrid)) ncongrid<-50; 
  
  ex<-order(xvals); xvals<-xvals[ex]; yvals<-yvals[ex]; 
  
  C<-0.61803399; R<-1-C;
  hvals<-hmin+0.1*c(0:10)*(hmax-hmin);
  Gvals<-matrix(0,11,1);
  for (i in 1:11) {
    Gvals[i]<-rssfitM(xvals,yvals,nknots,Boundary.knots,intercept=intercept,
                      alpha=hvals[i],penalty=penalty,ncongrid=ncongrid)$gcv;
    print(c(i,hvals[i],Gvals[i]));
  }
  imin<-minindc(Gvals);
  if((imin<=1)||(imin>=11)) {hopt<-hvals[imin]} else
  { 
    ax<-hvals[imin-1]; bx<-hvals[imin]; cx<-hvals[imin+1]; 
    x0<-ax; x3<-cx; 
    if(abs(cx-bx)>abs(bx-ax)) 
    {x1<-bx; x2<-bx+C*(cx-bx)} 
    else {x2<-bx; x1<-bx-C*(bx-ax)}
    f1<-rssfitM(xvals,yvals,nknots,Boundary.knots,intercept=intercept,penalty=penalty,alpha=x1,ncongrid=ncongrid)$gcv;	
    f2<-rssfitM(xvals,yvals,nknots,Boundary.knots,intercept=intercept,penalty=penalty,alpha=x2,ncongrid=ncongrid)$gcv;
    while(abs(x3-x0)>tol*(abs(x1)+abs(x2))) {
      if(f2<f1) {x0<-x1; x1<-x2; x2<-R*x1+C*x3; f1<-f2; 
                 f2<-rssfitM(xvals,yvals,nknots,Boundary.knots,intercept=intercept,penalty=penalty,alpha=x2,ncongrid=ncongrid)$gcv;}
      else {x3<-x2; x2<-x1; x1<-R*x2+C*x0; f2<-f1;
            f1<-rssfitM(xvals,yvals,nknots,Boundary.knots,intercept=intercept,penalty=penalty,alpha=x1,ncongrid=ncongrid)$gcv;}
    }		
    if(f1<f2) {hopt<-x1; fmin<-f1} else {hopt<-x2; fmin<-f2}
    print(c(hopt,fmin));
  }
  
  bestfit<-rssfitM(xvals,yvals,nknots,Boundary.knots,intercept=intercept,penalty=penalty,alpha=hopt,ncongrid=ncongrid);	
  
  if (draw==T) {
    win.graph(); par(mfrow=c(2,1));
    plot(hvals,log10(Gvals),type="l",xlab="log10(alpha)",
         ylab="GCV");
    title("Bandwidth selection curve");
    plot(xvals,yvals); points(xvals,bestfit$fitted.values,type="l");
    title("Data and Monotone P-spline fit");
  }
  
  return(list(hvals=hvals,Gvals=Gvals,hopt=hopt,x=xvals,fitted.values=bestfit$fitted.values,
              penalty=penalty,gcv=bestfit$gcv,model=bestfit))
  
}

bcgrad<-function(tvals,x,hder=NA) {
  xrange<-range(x)[2]-range(x)[1]; xmin<-range(x)[1];
  x<-(x-range(x)[1])/xrange;
  if(is.na(hder)) {hder<-(tvals[2]-tvals[1])}
  
  betahat<-locregG(tvals,x,tvals,hder,4); 
  xsmu<-betahat[,1]; xdot1<-betahat[,2];
  
  xfit<-xsmu;
  nt<-length(tvals); 
  tsim<-0.5*(tvals[2:nt]+tvals[1:(nt-1)]);
  betahat<-locregG(tvals,xfit,tsim,hder,4); 
  dxhat<-betahat[,2];
  
  f<-splinefun(tvals,xfit);
  dx<-0.05*(tvals[2]-tvals[1]);
  xplus<-f(tsim+dx); xminus<-f(tsim-dx);
  dxsim<-(xplus-xminus)/(2*dx);
  
  dxhat<-dxhat[5:(nt-4)];
  dxsim<-dxsim[5:(nt-4)]; 
  
  e<-order(dxhat);
  xvals<-dxhat[e]; yvals<-dxsim[e];
  
  fit<-gcv.rssM(xvals,yvals,hmin=-8,hmax=5,nknots=20,Boundary.knots=range(dxhat),
                intercept=T,penalty=1,ncongrid=50,tol=0.01,draw=F);
  
  f<-splinefun(xvals,fit$fitted.values,method="natural"); 
  yhat<-f(xdot1); 
  return(list(raw=xrange*xdot1,unbiased=xrange*yhat,smooth=xmin+xrange*xfit,hder=hder))
}

rssadd2<-function(x,y,nknots1=20,nknots2=20,Boundary.knots1,Boundary.knots2,intercept1=T,intercept2=T,
                  penalty=1, alpha,signs=c(1,1)) {
  #	Additive model, y=f_1(x[,1])+f_2(x[,2]) + e, cubic penalized regression splines for f1 & f2
  # 	Input variable 'signs' indicates sign constraints: 
  #		f_i >= 0 or <=0 according to whether signs[i] is +ive or -ive. 
  # 	Uses solve.QP to apply sign constraints at the knots
  #
  x1<-x[,1]; x2<-x[,2]; signs<-signs/abs(signs); 
  origin1<-Boundary.knots1[1];
  fac1<-(Boundary.knots1[2]-Boundary.knots1[1])/(nknots1+1);
  knots1<-origin1+c(1:nknots1)*fac1;
  x1mat<-rssbasis(x1,knots1,origin1,intercept1); 
  origin2<-Boundary.knots2[1];
  fac2<-(Boundary.knots2[2]-Boundary.knots2[1])/(nknots2+1);
  knots2<-origin2+c(1:nknots2)*fac2;
  x2mat<-rssbasis(x2,knots2,origin2,intercept2); 
  xmat<-x1mat%~%x2mat;
  XTX<-t(xmat)%*%xmat;
  
  d1<-rep(1,dim(x1mat)[2]); 
  d1[1:3]<-0; if(intercept1) d1[4]<-0;
  d2<-rep(1,dim(x2mat)[2]); 
  d2[1:3]<-0; if(intercept2) d2[4]<-0;
  D<-c((10^alpha[1])*d1,(10^alpha[2])*d2); D<-diag(D);
  
  # qprog matrices
  qmat<-XTX+D; rmat<-as.vector(t(xmat)%*%y);  	
  
  # constrain the sign of each term #
  aknots<-c(Boundary.knots1[1],knots1,Boundary.knots1[2]); 
  cmat1<- signs[1]*rssbasis(aknots,knots1,origin1,intercept1); 
  
  aknots<-c(Boundary.knots2[1],knots2,Boundary.knots2[2]); 
  cmat2<- signs[2]*rssbasis(aknots,knots2,origin2,intercept2); 
  
  cmat<-zeros(rows(cmat1)+rows(cmat2),cols(cmat1)+cols(cmat2));
  cmat[1:rows(cmat1),1:cols(cmat1)]<-cmat1;
  cmat[(1+rows(cmat1)):rows(cmat),(1+cols(cmat1)):cols(cmat)]<-cmat2;
  
  qpfit<-solve.QP(Dmat=qmat,dvec=rmat,Amat=t(cmat),meq=0,factorized=F);
  betahat<-qpfit$solution;
  yhat<-xmat%*%betahat;
  yhat1<-x1mat%*%betahat[1:cols(x1mat)]; 
  yhat2<-x2mat%*%betahat[(1+cols(x1mat)):(cols(x1mat)+cols(x2mat))];
  msr<-mean((yhat-y)^2);
  
  
  iact<-qpfit$iact; nact<-sum(iact>0);
  if(nact<1) {
    hat<-xmat%*%solve(XTX+D,tol=1.e-16)%*%t(xmat); Z<-1;
  }
  else {		
    iact<-iact[iact>0]; 
    Cact<-cmat[iact,];
    if(nact<2) Cact<-matrix(Cact,1,cols(cmat)) 
    Z<-Null(t(Cact)); 
    ZT<-t(Z);
    hat<-xmat%*%Z%*%solve(ZT%*%(XTX+D)%*%Z,tol=1.e-16)%*%ZT%*%t(xmat);
  }
  enp<-sum(diag(hat));
  d1<-sum(diag(hat))/length(y);
  if(penalty*d1>1) {gcv<-1.0e12} else {gcv<-msr/((1-penalty*d1)^2)};
  
  return(list(x=x,knots1=knots1,knots2=knots2,origin1=origin1, origin2=origin2,
              alpha=alpha, 
              intercept1=intercept1, intercept2=intercept2,coef=betahat,fitted.values=yhat1%~%yhat2,gcv=gcv,
              penalty=penalty,enp=enp))
}

gcvscore.rssadd2<-function(par,x,y,nknots1=20,nknots2=20,Boundary.knots1,Boundary.knots2,
                           intercept1=T,intercept2=T,penalty=1,signs=c(1,1)) 
{
  fit<-rssadd2(x,y,nknots1=nknots1,nknots2=nknots2,Boundary.knots1=Boundary.knots1, 
               Boundary.knots2=Boundary.knots2,intercept1=intercept1,intercept2=intercept2,
               penalty=penalty,alpha=par,signs=signs)
  return(log10(fit$gcv))
}

gcv.rssadd2<-function(x,y,nknots1=20,nknots2=20,Boundary.knots1=NA,Boundary.knots2=NA,
                      intercept1=T,intercept2=T,penalty=1,signs=c(1,1),par=c(0,0)) 
{
  if(is.na(Boundary.knots1)) Boundary.knots1<-range(x[,1]);
  if(is.na(Boundary.knots2)) Boundary.knots2<-range(x[,2]);
  bestfit<-optim(par,gcvscore.rssadd2,gr=NULL, method="Nelder-Mead", lower=-Inf, upper=+Inf,
                 control=list(maxit=1000,trace=3), hessian=F,
                 x,y,nknots1,nknots2,Boundary.knots1,
                 Boundary.knots2,intercept1=intercept1,intercept2=intercept2,penalty,signs);
  hopt<-bestfit$par; 
  bestfit<-optim(hopt,gcvscore.rssadd2,gr=NULL, method="Nelder-Mead", lower=-Inf, upper=+Inf,
                 control=list(maxit=1000,trace=3), hessian=F,
                 x,y,nknots1,nknots2,Boundary.knots1,
                 Boundary.knots2,intercept1,intercept2,penalty,signs);
  hopt<-bestfit$par; 
  
  rssfit<-rssadd2(x,y,nknots1,nknots2,Boundary.knots1,Boundary.knots2,intercept1=intercept1,
                  intercept2=intercept2,penalty=penalty,alpha=bestfit$par,signs=signs);
  
  return(list(hopt=bestfit$par,gcv=rssfit$gcv, fitted.values=rssfit$fitted.values,model=rssfit));
}

predict.rssadd2<-function(model, newx) {
  x1<-newx[,1]; x2<-newx[,2];  
  x1mat<-rssbasis(x1,model$knots1,model$origin1,model$intercept1); 
  x2mat<-rssbasis(x2,model$knots2,model$origin2,model$intercept2); 
  betahat<-model$coef; 
  yhat1<-x1mat%*%betahat[1:cols(x1mat)]; 
  yhat2<-x2mat%*%betahat[(1+cols(x1mat)):(cols(x1mat)+cols(x2mat))];
  return(list(x=newx,y=yhat1%~%yhat2));
}



