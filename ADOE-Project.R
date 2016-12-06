##ISYE 7400 Project Code
setwd("/Users/dshartra/Documents/GT-Fall-2016/ADOE")

## Load in data
data = read.table("runprojdata.txt",header=T,sep=" ")
head(data)

##Function to calculate vdot (approximation of you V02Max)
vdotf = function(x){
  permax = 0.8 +0.1894393*exp(-0.012778*x[2])+0.2989558*exp(-0.1932605*x[2]) #calculate percent maximum
  v02 = -4.6 +0.182258*(x[1]/x[2]) + 0.000104*(x[1]/x[2])^2 #calculate your v02
  vdot = round(v02/permax) #calculate vdot
  return(vdot)
  
}

##Analysis of Observed base data
x = as.matrix(data)

##Edit data to prevent ties in time (scaled based on tenths and hundreths of seconds)
p1 = sort(x[,1])
p2 = sort(x[,2])
p2[10] = p2[10] +.017
p2[17] = p2[17] +.017
X = cbind(p1,p2)
X = X[1:18,] #remove cross country times (prediction is much better, makes sense as hard to compare XC to track)
f12 = apply(X,1,vdotf)

##Collect Test Dataset (obtained as logitudinal study, same athletes, same coach, same distances, year following training set)
test = read.csv('runprojtest.txt',sep="")
p3 = sort(test$Distance)
p4 = sort(test$Time)
test = cbind(p3,p4)
test = test[1:25,] #remove cross country times (again for better prediction)

##MARS (regression is easier to predict new time from distance)
library(mda)
a=mars(X,f12, degree=2)
summary(lm(f12~a$x-1))
library(lhs)
y.test=apply(test,1,vdotf)
y.pred=predict(a,test)
sqrt(mean((y.pred-y.test)^2)) ##approximately 0.8777997 
plot(y.pred,y.test)
abline(0,1)

##Kernel Regression
y = apply(X,1,vdotf)
E=(as.matrix(dist(X, diag=T, upper=T)))^2
n=18
one=rep(1,n)
I =diag(n)
mscv=function(theta)
{
  R=exp(-theta*E)+10^(-6)*I
  den=c(R%*%one)
  e=y-diag(1/den)%*%R%*%y
  cv=e/(1-1/den)
  return(mean(cv^2))
}

theta=optimize(mscv,c(.1,100000))$min
theta #35,906.05

r=function(x)
  exp(-theta*(x-X)^2)
fhat=function(x)
  sum(r(x)*y)/sum(r(x))

y.kr = apply(test,1,fhat)
y.test = apply(test,1,vdotf)
sqrt(mean((y.kr-y.test)^2)) ## approximatey 4.734648
plot(y.kr,y.test)
abline(0,1)

#Inverse Distance Weigthing
y = apply(X,1,vdotf)
fhat=function(x)
{
  d=abs(x-X)
  val=sum(y/(d+10^(-6))^2)/sum(1/(d+10^(-6))^2)
  return(val)
}
y.idw = apply(test,1,fhat)
y.test = apply(test,1,vdotf)
sqrt(mean((y.idw-y.test)^2)) ## approximatey 4.734648
plot(y.idw,y.test)
abline(0,1)

##LOESS
## Has issues with the singularties
y = apply(X,1,vdotf)
X1 = as.data.frame(X)
b=loess(y~p1+p2,data=X1,degree=2,span=.5)
y.loess=predict(b,test) ##why is this not working properly but it does in class code
sqrt(mean((y.loess-y.test)^2))
plot(y.loess,y.test)
abline(0,1)

##Radial Basis
y=apply(X,1,vdotf)
E=as.matrix(dist(X, diag=T, upper=T))
n=18
I=diag(n)
MSCV=function(para)
{
  R=exp(-para*E^2)
  Rinv=solve(R+10^(-6)*I)
  Rinvd=diag(Rinv)
  cv=diag(1/Rinvd)%*%Rinv%*%y
  val=mean(cv^2)
  return(val)
}

a.opt=optimize(MSCV,c(1,10000))
theta=a.opt$min
theta #58.24415
R=exp(-theta*E^2)
basis=function(h)
  exp(-theta*sum(h^2))
r=function(x)
{
  A=t(t(X)-x)
  vec=apply(A,1,basis)
  return(vec)
}
coef=solve(R,y)

fhat=function(x)
  t(r(x))%*%coef
y.test=apply(test,1,vdotf)
y.rbf=apply(test,1,fhat)
sqrt(mean((y.rbf-y.test)^2)) #36.41985
plot(y.rbf,y.test)
abline(0,1)

##Radial Basis Function with polynomial precision
y=apply(X,1,vdotf)
E=as.matrix(dist(X, diag=T, upper=T))
n=18
I=diag(n)
one=rep(1,n)
MSCV=function(para)
{
  R=exp(-para*E^2)
  Rinv=solve(R+10^(-10)*I)
  Rinvd=diag(Rinv)
  mu=drop(t(one)%*%Rinv%*%y/(t(one)%*%Rinv%*%one))
  cv=diag(1/Rinvd)%*%Rinv%*%(y-mu)
  val=log(mean(cv^2))
  return(val)
}
a.opt=optimize(MSCV,c(1,10000))
theta=a.opt$min
theta #203.7373

R=exp(-theta*E^2)
Rinv=solve(R+10^(-10)*I)
mu=drop(t(one)%*%Rinv%*%y/(t(one)%*%Rinv%*%one))

coef=Rinv%*%(y-mu)

basis=function(h)
  exp(-theta*sum(h^2))
r=function(x)
{
  A=t(t(X)-x)
  vec=apply(A,1,basis)
  return(vec)
}
fhat=function(x)
  mu+t(r(x))%*%coef

y.test=apply(test,1,vdotf)
y.rbfpp=apply(test,1,fhat)
sqrt(mean((y.rbfpp-y.test)^2)) #3.492863
plot(y.rbfpp,y.test)
abline(0,1)

## Some Kriging Stuff Now
## Ordinary Kriging(should be the same as RBF)
y=apply(X,1,vdotf)
n=18
E=as.matrix(dist(X, diag=T, upper=T))
theta=203.7373
R=exp(-theta*E^2)
one=rep(1,n)
I=diag(n)
Rinv=solve(R+10^(-10)*I)
mu=drop(t(one)%*%Rinv%*%y/(t(one)%*%Rinv%*%one))
coef=Rinv%*%(y-mu)

basis=function(h)
  exp(-theta*sum(h^2))
r=function(x)
{
  A = t(t(X)-x)
  vec=apply(A,1,basis)
  return(vec)
}
fhat=function(x)
  mu+t(r(x))%*%coef
y.ok=apply(test,1,fhat)
y.test=apply(test,1,vdotf)
sqrt(mean((y.ok-y.test)^2)) #3.492863, confirm same as RBF so thats good to see
plot(y.ok,y.test)
abline(0,1)

#Limit Kriging
theta = 203.7373
n=18
E=as.matrix(dist(X, diag=T, upper=T))
R=exp(-theta*E^2)
one=rep(1,n)
I=diag(n)
Rinv=solve(R+10^(-10)*I)
mu=drop(t(one)%*%Rinv%*%y/(t(one)%*%Rinv%*%one))
coef=Rinv%*%y

basis=function(h)
  exp(-theta*sum(h^2))
r=function(x)
{
  A = t(t(X)-x)
  vec=apply(A,1,basis)
  return(vec)
}
den=Rinv%*%one
fhat=function(x){
  t(r(x))%*%coef/(t(r(x))%*%den)
}
y.lk = apply(test,1,fhat)
y.test = apply(test,1,vdotf)
sqrt(mean((y.lk-y.test)^2)) #2.148368
plot(y.lk,y.test)
abline(0,1)

## Universal Kriging
y=apply(X,1,vdotf)
n=18
E=as.matrix(dist(X, diag=T, upper=T))
theta=203.7373
R=exp(-theta*E^2)
one=rep(1,n)
I=diag(n)
Rinv=solve(R+10^(-10)*I)
P=as.matrix(cbind(one,X))
Sinv=solve(t(P)%*%Rinv%*%P)
beta=drop(Sinv%*%t(P)%*%Rinv%*%y)
coef=Rinv%*%(y-P%*%beta)

basis=function(h)
  exp(-theta*sum(h^2))
r=function(x)
{
  A = t(t(X)-x)
  vec=apply(A,1,basis)
  return(vec)
}

fhat=function(x){
  t(c(1,x))%*%beta+t(r(x))%*%coef
}
y.uk = apply(test,1,fhat)
y.test = apply(test,1,vdotf)
sqrt(mean((y.uk-y.test)^2)) #1.857107
plot(y.uk,y.test)
abline(0,1)

par(mfrow=c(3,3))
plot(y.pred,y.test,main="MARS")
abline(0,1)
plot(y.kr,y.test,main="Kernel Regression")
abline(0,1)
plot(y.idw,y.test,main="Inverse Distance Weighting")
abline(0,1)
plot(y.loess,y.test,main="Loess")
abline(0,1)
plot(y.rbf,y.test,main="Radial Basis")
abline(0,1)
plot(y.rbfpp,y.test,main="RBF Polynomial Precision")
abline(0,1)
plot(y.ok,y.test,main="Ordinary Kriging")
abline(0,1)
plot(y.lk,y.test,main="Limit Kriging")
abline(0,1)
plot(y.uk,y.test,main="Universial Kriging")
abline(0,1)
par(mfrow=c(1,1))
## Sensitivity Analysis
min1 = min(test[,1])
max1 = max(test[,1])
min2 = min(test[,2])
max2 = max(test[,2])
z1 = (test[,1]-min1)/(max1-min1)
z2 = (test[,2]-min2)/(max2-min2)
z = cbind(z1,z2)
low = c(1500,3.983)
up = c(5000,17.7)
#Function to scale to [0,1] to do sensitivity analysis (did not do in regular analysis because of various 
# difficulties to be discussed later)
scale01 = function(x){
  x=low+diag(up-low)%*%x
}
dat = apply(z,1,scale01)
datnew = cbind(dat[1,],dat[2,])
f.eval=predict(a,datnew)
f0=mean(f.eval)
x.ev=seq(0,1,length=20)
val=x.ev
n=25
main.effect=function(ind)
{
  D0=z
  for(i in 1:20)
  {
    D0[,ind]=rep(x.ev[i],n)
    D0temp =apply(D0,1,scale01)
    D0new = cbind(D0temp[1,],D0temp[2,])
    val[i]=mean(predict(a,D0new))
  }
  return(val-f0)
}
d=2
M=matrix(0,nrow=20,ncol=d)
for(j in 1:d)
  M[,j]=main.effect(j)
matplot(x.ev,M,type="l",lty=1:d,col=1:d,lwd=2,main="Sensitivity of Variables Distance and Time")
legend(.6,0,c("Distance","Time"),bty="n",lty=1:d,col=1:d,lwd=2)

## Prediction based on "Case Studies"
library(chemCal)
## EXAMPLE: Runner A runs this time in 5k and gets a vdot value from
## one of the prediction models. Run an inverse prediction model
## for all 3 distances (1500,3k,5k) and use the fitted vdot
## to predict the time run for other distances.
w15 =test[1:9,2]
w3 = test[10:18,2]
w5 = test[19:25,2]
pred15 = y.pred[1:9]
pred3 = y.pred[10:18]
pred5 = y.pred[19:25]
inv.pred15 <- lm(pred15 ~ w15)
inv.pred3 <- lm(pred3 ~ w3)
inv.pred5 <- lm(pred5 ~ w5)
##Case studies of team members
##predict kyles 1500/5k time from his 3k vdot from 2015
kyle = cbind(3000,8.483)
kyinv = apply(kyle,1,vdotf)
ky15 = inverse.predict(inv.pred15,kyinv)
ky5 = inverse.predict(inv.pred5,kyinv)
##predict casey 1500/5k time from his 3k vdot from 2015
casey = cbind(3000,10.217)
cainv = apply(casey,1,vdotf)
ca15 = inverse.predict(inv.pred15,cainv)
ca5 = inverse.predict(inv.pred5,cainv)
##predict merlin from 2015 by 3k
merlin = cbind(3000,9.317)
merinv = apply(merlin,1,vdotf)
mer15 = inverse.predict(inv.pred15,merinv)
mer5 = inverse.predict(inv.pred5,merinv)
##predict collin from 2015 by 3k
col = cbind(3000,8.933)
colinv = apply(col,1,vdotf)
col15 = inverse.predict(inv.pred15,colinv)
col5 = inverse.predict(inv.pred5,colinv)
##predict steve from 2014 by 1500
st = cbind(1500,4.083)
stinv = apply(st,1,vdotf)
st3 = inverse.predict(inv.pred3,stinv)
st5 = inverse.predict(inv.pred5,stinv)
##predict chris from 2012 by 5k
ch = cbind(5000,15.883)
chinv = apply(ch,1,vdotf)
ch15 = inverse.predict(inv.pred15,chinv)
ch3 = inverse.predict(inv.pred3,chinv)
##jack bennett from 2012 by 5k
ja = cbind(5000,15.45)
jainv = apply(ja,1,vdotf)
ja15 = inverse.predict(inv.pred15,jainv)
ja3 = inverse.predict(inv.pred3,jainv)

