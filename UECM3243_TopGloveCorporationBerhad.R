library(fBasics)
library(fGarch)
library(forecast)

tg1 = read.csv("UECM3243_TopGlove.csv", header = T)
# Option 1 
# to eliminate null value(s) and check numerical
sapply(X = tg1, FUN = function(x) sum(x=="null"))
tg = tg1[tg1$Open!="null",]
tp_open = as.numeric(levels(tg$Open))[tg$Open]
#tp_open=as.numeric((tg$Open))
str(tp_open)

###### PART I #####
# daily log returns
lrtn = diff(log(tp_open))

######   time series plot       #####
# time stamp
year= c(1:2471)/252 + 2009

par(mfrow=c(1,1))
plot(year, tp_open, main = "Top Glove Corp Bhd", type = "l",
     ylab = "Daily Opening Price")

year1= c(1600:2471)/252 + 2009
plot(year1, tp_open[1600:2471], main = "Top Glove Corp Bhd", type = "l",
     ylab = "Daily Opening Price",xlim=c(2015,2019))

######   ACF PACF tp_open      ######
par(mfrow=c(1,2))
acf(tp_open, main="Top Glove Corp Bhd")
pacf(tp_open, main="Top Glove Corp Bhd")
#acf decays very slowly
#pacf lag-1 is approximately 1.0, first reg diff required

######   log return plot        #######
par(mfrow=c(1,1))
plot(year[-1], lrtn, main = "Top Glove Corp Bhd", type = "l", xlab= "Year", ylab = "Daily Log Returns")
#appears to be weekly stationary in mean and variance
#volatility cluster exists in 2018

year1= c(1600:2471)/252 + 2009
plot(year1[-1], diff(log(tp_open[1600:2471])), main = "Top Glove Corp Bhd", type = "l",
     ylab = "Daily Opening Price",xlim=c(2015,2019))


######   ACF PACF lrtn          ######
par(mfrow=c(1,2))
acf(lrtn, main = "Daily log returns of Top Glove Corp Bhd")
pacf(lrtn, main = "Daily log returns of Top Glove Corp Bhd")

# First 12-lags of acf lie within the error band
# suggesting white noise series and no evidence of non-randomness

######   Histogram & density plot  lrtn   #############
par(mfrow=c(1,1))
hist(lrtn, nclass=40, xlab="Daily Log Returns", main="Top Glove Corp Bhd") 
plot(density(lrtn)$x, density(lrtn)$y, type="l", xlab="Daily Log Returns", ylab="Density", main="Top Glove Corp Bhd") 
range(lrtn)
x = seq(-0.3, 0.12, 0.001)   
y = dnorm(x, mean(lrtn), sd(lrtn))
lines(x, y, lty=2, col="blue") 

######   Statistical Testings (Analysis)   ######
basicStats(lrtn)

# Mean testing 
t.test(lrtn)
# Mean is significant. Need to include mean

# Skewness testing 
t1 = skewness(lrtn)/sqrt(6/length(lrtn))
t1
pv1 = 2*pnorm(t1)
pv1
#  The daily log returns is skewed to the left

# Tail thickness testing 
t2 = kurtosis(lrtn)/sqrt(24/length(lrtn))
t2
pv2 = 2*(1-pnorm(t2))
pv2
# Heavy tails

# JB test 
normalTest(lrtn, method = "jb")
# Normality assumption for lrtn is rejected. 

######   Data Testing           ######

Box.test(lrtn, lag = 12, type = "Ljung")
# Do not rej H_0, serially uncorrelated

par(mfrow=c(1,2))

acf(lrtn,main="Daily Log Return")
# From first 12 acf and pacf, weakly stationary and white noise

acf(abs(lrtn), main = "Absolute Daily Log Return")
Box.test(abs(lrtn), lag = 12, type = "Ljung") 
# Rej H_0, it is dependent series

x1 = lrtn - mean(lrtn) #estimate of epsilon_t 

Box.test(x1^2, lag = 12, type = "Ljung")
# Rej H_0, exist ARCH effect





####   PART II      #####
par(mfrow=c(1,2))
acf(x1^2, main=expression(epsilon[t]^2),ylim=c(-0.1,1))
pacf(x1^2, main=expression(epsilon[t]^2),ylim=c(-0.1,1))

#######   GARCH (1,2)        ######
#based on acf and pacf
a1.1 = garchFit(~garch(1,2), data = lrtn, include.mean = T, trace = F)
summary(a1.1) #mean is insignificant
a1.1 = garchFit(~garch(1,2), data = lrtn, include.mean = F, trace = F)
summary(a1.1)
#Pass SRT,SSRT for first 15 lags.
par(mfrow=c(1,3))
plot(a1.1)
10
11
13 #QQ plot deviates, normality assumption fails
0

a1.2 = garchFit(~garch(1,2), data = lrtn, include.mean = T, trace = F, cond.dist = "std")
summary(a1.2) #mean is insignificant
a1.2 = garchFit(~garch(1,2), data = lrtn, include.mean = F, trace = F, cond.dist = "std")
summary(a1.2)
#Pass SRT,SSRT for first 15 lags.
par(mfrow=c(1,3))
plot(a1.2)
10
11
13 #QQ plot seems to be adequate
0

skewness(lrtn)
a1.3 = garchFit(~garch(1,2), data = lrtn, include.mean = T, trace = F, cond.dist = "sstd")
summary(a1.3) #beta 2 insignificant,reduce to garch (1,1)
a1.3 = garchFit(~garch(1,1), data = lrtn, include.mean = T, trace = F, cond.dist = "sstd")
summary(a1.3) 
#Pass SRT,SSRT for first 15 lags.
par(mfrow=c(1,3))
plot(a1.3)
10
11
13 #QQ plot adequate
0

# all coef sig at 5%, pass SRT , SSRT for Gaussian  mean insig  AIC -5.240867
# all coef sig at 5%, pass SRT , SSRT for std       mean insig  AIC -5.449111
# all coef sig at 5%, pass SRT , SSRT for sstd     beta2 insig  AIC -5.450186

#skew parameter testing
t_a1=(1.066e+00-1)/ 2.679e-02 
t_a1
pv_a1=2*(1-pnorm(t_a1))
pv_a1
#rej, model adequate. 

v1.3 = volatility(a1.3) 
resi1 = residuals(a1.3, standardize = T)
par(mfrow=c(3,1))
plot(year[-1], lrtn, type="l", xlab="year", ylab="log returns", main ="Top Glove Corp Bhd")
plot(year[-1], v1.3, type="l", xlab="year", ylab="volatility", main ="Top Glove Corp Bhd")
plot(year[-1], resi1, type="l", xlab="year", ylab="standardized residuals", main ="Top Glove Corp Bhd")

par(mfrow=c(1,3))
plot(a1.1)
13
0
plot(a1.2)
13
0
plot(a1.3)
13
0

#model a1.3 preferred

#######   GARCH (1,1)        ######

a2.1 = garchFit(~garch(1,1), data = lrtn, include.mean = T, trace = F)
summary(a2.1) #mean is insignificant
a2.1 = garchFit(~garch(1,1), data = lrtn, include.mean = F, trace = F)
summary(a2.1)
#Pass SRT,SSRT for first 15 lags.
par(mfrow=c(1,3))
plot(a2.1)
10
11
13 #QQ plot deviates, normality assumption fails
0

a2.2 = garchFit(~garch(1,1), data = lrtn, include.mean = T, trace = F, cond.dist = "std")
summary(a2.2) #mean is insignificant
a2.2 = garchFit(~garch(1,1), data = lrtn, include.mean = F, trace = F, cond.dist = "std")
summary(a2.2)
#Pass SRT,SSRT for first 15 lags.
par(mfrow=c(1,3))
plot(a2.2)
10
11
13 #QQ plot seems to be adequate
0

a2.3 = garchFit(~garch(1,1), data = lrtn, include.mean = T, trace = F, cond.dist = "sstd")
summary(a2.3)
#Pass SRT,SSRT for first 15 lags.
par(mfrow=c(1,3))
plot(a2.3)
10
11
13 # QQ plot shows adequate
0

# all coef sig at 5%, pass SRT , SSRT for Gaussian  mean insig  AIC -5.239471
# all coef sig at 5%, pass SRT , SSRT for std       mean insig  AIC -5.448725
# all coef sig at 5%, pass SRT , SSRT for sstd                  AIC -5.450186

t_a2=(1.066e+00-1)/2.679e-02
t_a2
pv_a2=2*(1-pnorm(t_a2))
pv_a2
#rej h0, model adequate

v2.3 = volatility(a2.3) 
resi2 = residuals(a2.3, standardize = T)
par(mfrow=c(3,1))
plot(year[-1], lrtn, type="l", xlab="year", ylab="log returns", main ="Top Glove Corp Bhd")
plot(year[-1], v2.3, type="l", xlab="year", ylab="volatility", main ="Top Glove Corp Bhd")
plot(year[-1], resi2, type="l", xlab="year", ylab="standardized residuals", main ="Top Glove Corp Bhd")

par(mfrow=c(1,3))
plot(a2.1)
13
0
plot(a2.2)
13
0
plot(a2.3)
13
0

#model a2.3 preferred

#######   GARCH (2,0)        ######
ar.mle(x1^2)$order
#suggest GARCH(2,0)
auto.arima(x1^2)
#suggest GARCH(0,2), for garch model must at least p =1 , model fails

a3.1 = garchFit(~garch(2,0), data = lrtn, include.mean = T, trace = F)
#cannot run as distribution of yt has heavier tail than normal distribution
summary(a3.1)
#cannot run

a3.2 = garchFit(~garch(2,0), data = lrtn, include.mean = T, trace = F, cond.dist = "std")
summary(a3.2)  #mean insignificant
a3.2 = garchFit(~garch(2,0), data = lrtn, include.mean = F, trace = F, cond.dist = "std")
summary(a3.2)
#Pass SRT,SSRT for first 15 lags.
plot(a3.2)
10
11
13
0
#QQ plot seems to be adequate

a3.3 = garchFit(~garch(2,0), data = lrtn, include.mean = T, trace = F, cond.dist = "sstd")
summary(a3.3)  
#Pass SRT,SSRT for first 15 lags.
plot(a3.3)
10
11
13 #QQ plot seems to be adequate
0

# Fail  for Gaussian  
# all coef sig at 5%, pass SRT first 12 lags, SSRT  for std   mean insig  AIC -5.438258
# all coef sig at 5%, pass SRT first 12 lags, SSRT  for sstd  AIC -5.440216

t_a3=( 1.073e+00-1)/2.731e-02
t_a3
pv_a3 = 2*(1-pnorm(t_a3))
pv_a3
# rej H0, model adequate


v3.3 = volatility(a3.3) 
resi3 = residuals(a3.3, standardize = T)
par(mfrow=c(3,1))
plot(year[-1], lrtn, type="l", xlab="year", ylab="log returns", main ="Top Glove Corp Bhd")
plot(year[-1], v3.3, type="l", xlab="year", ylab="volatility", main ="Top Glove Corp Bhd")
plot(year[-1], resi3, type="l", xlab="year", ylab="standardized residuals", main ="Top Glove Corp Bhd")

par(mfrow=c(1,2))

plot(a3.2)
13
0
plot(a3.3)
13
0
#model a3.3 preferred

#######   Model Preferred          #####
a1.3  #GARCH(1,1) AIC -5.450186
a2.3  #GARCH(1,1) AIC -5.450186
a3.3  #GARCH(2,0) AIC -5.440216 

#Lowest AIC , GARCH(1,1)
# one-step ahead
par(mfrow=c(1,1))
pre=predict(a1.3, n.ahead=1,plot=T)

lcv=qsstd(0.025,nu=3.377e+00,xi=1.066e+00)
lcv
rcv=qsstd(0.975,nu=3.377e+00,xi=1.066e+00)
rcv
#predictive interval
ll=0.0007182676-1.82443*0.02072517
ll
#-0.03710246
ul=0.0007182676+1.989323*0.02072517
ul
#0.04193786

plot(year[-1],lrtn,type="l",xlab="year",main="1 Step Ahead Forecast")
points(year[length(year)],pre$meanForecast , col = "red",pch=4)
points(year[length(year)],ll , col = "blue", cex = 1.5,pch=4)
points(year[length(year)],ul , col = "green", cex = 1.5,pch=4)
legend("bottomleft", legend=c("mean Forecast", 
                        expression("mean Forecast +1.989"+sqrt(MSE)), expression("mean Forecast -1.825"+sqrt(MSE))),
       col=c("red","green","blue"), cex=1.5,box.lty=0,pch=c(4,4,4))

#######   Predictive Interval         ####
a1.3 = garchFit(~garch(1,1), data = lrtn, include.mean = T, trace = F, cond.dist = "sstd")
summary(a1.3)
v1.3 = volatility(a1.3) 
v1.3
mu=7.183e-04 
lcv=qsstd(0.025,nu=3.377e+00,xi=1.066e+00)
lcv
rcv=qsstd(0.975,nu=3.377e+00,xi=1.066e+00)
rcv
UL = mu + rcv*v1.3
LL = mu + lcv*v1.3

year= c(1:2471)/252 + 2009
plot(year[-1],lrtn,type="l",xlab="year",ylab="Daily log returns",
     main="Top Glove Corp")
lines(year[-1],UL,lty=2 , col="navy")
lines(year[-1],LL,lty=2 , col="navy")
abline(h=mu,col='red',lwd=1.5)

####   Part III           ####
am1 = arima(x1^2, order = c(1,0,1), include.mean = T)
am1

mean(lrtn)     #for the purpose of mean equation 

a1.3 = garchFit(~garch(1,1), data = lrtn, include.mean = T, trace = F, cond.dist = "sstd")
ex1 = x1^2 - am1$residuals
v1.3 = volatility(a1.3)
cor(v1.3, sqrt(ex1))

