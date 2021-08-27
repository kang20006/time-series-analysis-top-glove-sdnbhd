library(fBasics)
library(fGarch)
library(forecast)

tg1 = read.csv("UECM3243_TopGlove.csv", header = T)

# Option 1
# to eliminate null value(s) and check numerical
sapply(X = tg1, FUN = function(x) sum(x=="null"))
tg = tg1[tg1$Adj.Close!="null",]
#tp_adjcl = as.numeric(levels(tg$Adj.Close))[tg$Adj.Close]
tp_adjcl=as.numeric((tg$Adj.Close))
str(tp_adjcl)

# Option 2
#run str(tg1)
# Check if adjust.close shows as numeric then do nothing
# if adjust.close shows as factor run
#tg_adjcl=as.numeric(levels(tg$Adj.Close))[tg$Adj.Close]
# if adjust.close shows as char run
#tg_adjcl=as.numeric((tg$Adj.Close))

# daily log returns
lrtn = diff(log(tp_adjcl))

# time stamp
year= c(1:2719)/252 + 2008

######   time series plot       #####
par(mfrow=c(1,1))
plot(year, tp_adjcl, main = "Top Glove Corp Bhd", type = "l",
     ylab = "Daily Adjusted Closing Price")

#plot(c(1:521)/252 + 2001.25, tp_adjcl[1:521], main = "Top Glove Corp Bhd", type = "l", ylab = "Daily Adjusted Closing Price",xlim=c(2001,2004))
#year1=c(1:521)/252 + 2001.25
#plot(year1[-1], diff(log(tp_adjcl[1:521])), main = "Top Glove Corp Bhd", type = "l", ylab = "Daily log returns")

#plot(c(2481:3462)/252 + 2001.25, tp_adjcl[2481:3462], main = "Top Glove Corp Bhd", type = "l", ylab = "Daily Adjusted Closing Price")
#year2=c(2481:3462)/252 + 2001.25
#plot(year2[-1], diff(log(tp_adjcl[2481:3462])), main = "Top Glove Corp Bhd", type = "l", ylab = "Daily log returns")

#plot(c(3463:4454)/252 + 2001.25, tp_adjcl[3463:4454], main = "Top Glove Corp Bhd", type = "l", ylab = "Daily Adjusted Closing Price")
#year3=c(3463:4454)/252 + 2001.25
#plot(year3[-1], diff(log(tp_adjcl[3463:4454])), main = "Top Glove Corp Bhd", type = "l", ylab = "Daily log returns")


######   ACF PACF tp_adjcl      ######
par(mfrow=c(1,2))
acf(tp_adjcl, main="Top Glove Corp Bhd")
pacf(tp_adjcl, main="Top Glove Corp Bhd")
#acf decays very slowly
#pacf lag-1 is approximately 1.0, first reg diff required

######   log return plot        #######
par(mfrow=c(1,1))
plot(year[-1], lrtn, main = "Top Glove Corp Bhd", type = "l", xlab= "Year", ylab = "Daily Log Returns")
#appears to be weekly stationary in mean and variance


######   ACF PACF lrtn          ######
par(mfrow=c(1,2))
acf(lrtn, main = "Daily log returns of Top Glove Corp Bhd")
pacf(lrtn, main = "Daily log returns of Top Glove Corp Bhd")

# First 12-lags of acf lie within the error band
# suggesting white noise series and no evidence of non-randomness

######   Histogram & density plot  lrtn   #############
par(mfrow=c(1,1))
hist(lrtn, nclass=40, xlab="Daily Log Returns", main="Top Glove Corp Bhd") 
plot(density(lrtn)$x, density(lrtn)$y, type="l", xlab="Daily Log Returns", ylab="Density", main="Top Corp Bhd") 

x = seq(-0.3, 0.1, 0.001)   
y = dnorm(x, mean(lrtn), sd(lrtn))
lines(x, y, lty=2, col="blue") 

######   Statistical Testings   ######
basicStats(lrtn)

# Mean testing 
t.test(lrtn)
# Mean is significant. Need to include mean

# Skewness testing 
t1 = skewness(lrtn)/sqrt(6/length(lrtn))
t1
pv1 = 2*pnorm(t1)
pv1
# The daily log returns is not symmetrical

# Tail thickness testing 
t2 = kurtosis(lrtn)/sqrt(24/length(lrtn))
t2
pv2 = 2*(1-pnorm(t2))
pv2
# Heavy tails

# JB test 
t3 = t2^2 + t1^2
normalTest(lrtn, method = "jb")
# Normality assumption for lrtn is rejected. 

######   1(C) Testing           ######

Box.test(lrtn, lag = 12, type = "Ljung")
# Do not rej H_0, serially uncorrelated

par(mfrow=c(1,2))

acf(lrtn,main="Daily Log Return")
pacf(lrtn,main="Daily Log Return")
# From first 12 acf and pacf, weakly stationary and white noise

acf(abs(lrtn), main = "Absolute Daily Log Return")
Box.test(abs(lrtn), lag = 12, type = "Ljung") 
# Rej H_0, it is dependent series

x1 = lrtn - mean(lrtn)

Box.test(x1^2, lag = 12, type = "Ljung")
# Rej H_0, exist ARCH effect

#####                           ######
####   PART II      #####
#######   GARCH (1,2)        #####
par(mfrow=c(1,2))
acf(x1^2, main=expression(epsilon[t]^2))
pacf(x1^2, main=expression(epsilon[t]^2))
a1.1 = garchFit(~garch(1,2), data = lrtn, include.mean = F, trace = F, cond.dist = "norm") 
summary(a1.1)
par(mfrow=c(1,3))
plot(a1.1)
10
11
13
0

a1.2 = garchFit(~garch(1,2), data = lrtn, include.mean = F, trace = F, cond.dist = "std") 
summary(a1.2)
plot(a1.2)
10
11
13
0

a1.3 = garchFit(~garch(1,2), data = lrtn, include.mean = F, trace = F, cond.dist = "sstd") 
summary(a1.3)
plot(a1.3)
10
11
13
0

# all coef sig at 5%, pass SRT , SSRT for Gaussian  mean insig  AIC -5.204402
# all coef sig at 5%, pass SRT , SSRT for std       mean insig  AIC -5.488314
# all coef sig at 5%, pass SRT , SSRT for sstd      mean insig  AIC -5.487817

t_a1 = (9.061e-01-1)/5.095e-02
pv_a1 = 2* pnorm(t_a1)
pv_a1
#inadequate sstd since symmetry assumption do not rej

v1.2 = volatility(a1.2) 
resi1 = residuals(a1.2, standardize = T)
par(mfrow=c(3,1))
plot(year[-1], lrtn, type="l", xlab="year", ylab="log returns", main ="Top Glove Corp Bhd")
plot(year[-1], v1.2, type="l", xlab="year", ylab="volatility", main ="Top Glove Corp Bhd")
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
# Model a1.2 preferred

#######   GARCH (1,1)        ######
a2.1 = garchFit(~garch(1,1), data = lrtn, include.mean = T, trace = F, cond.dist = "norm")
summary(a2.1)
par(mfrow=c(1,3))
plot(a2.1)
10
11
13
0

a2.2 = garchFit(~garch(1,1), data = lrtn, include.mean = F, trace = F, cond.dist = "std")
summary(a2.2)
par(mfrow=c(1,3))
plot(a2.2)
10
11
13
0

a2.3 = garchFit(~garch(1,1), data = lrtn, include.mean = F, trace = F, cond.dist = "sstd")
summary(a2.3)
par(mfrow=c(1,3))
plot(a2.3)
10
11
13
0

# all coef sig at 5%, pass SRT , SSRT for Gaussian              AIC -5.195012
# all coef sig at 5%, pass SRT , SSRT for std       mean insig  AIC -5.487884
# all coef sig at 5%, pass SRT , SSRT for sstd      mean insig  AIC -5.487372

t_a2=(1.016e+00-1)/2.114e-02 
pv_a2=2*(1-pnorm(t_a.2))
pv_a2
#inadequate sstd since symmetry assumption do not rej

v2.2 = volatility(a2.2) 
resi2 = residuals(a2.2, standardize = T)
par(mfrow=c(3,1))
plot(year[-1], lrtn, type="l", xlab="year", ylab="log returns", main ="Top Glove Corp Bhd")
plot(year[-1], v2.2, type="l", xlab="year", ylab="volatility", main ="Top Glove Corp Bhd")
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

#model a2.2 preferred

#######   GARCH (1,0)        ######
ar.mle(x1^2)$order
a3.1 = garchFit(~garch(1,0), data = lrtn, include.mean = T, trace = F)
summary(a3.1)
plot(a3.1)
10
11
13
0

a3.2 = garchFit(~garch(1,0), data = lrtn, include.mean = F, trace = F, cond.dist = "std")
summary(a3.2)
plot(a3.2)
10
11
13
0

a3.3 = garchFit(~garch(1,0), data = lrtn, include.mean = F, trace = F, cond.dist = "sstd")
summary(a3.3)
plot(a3.3)
10
11
13
0

# all coef sig at 5%, pass SRT Q(10), SSRT Q(20) for Gaussian              AIC -5.187530
# all coef sig at 5%, pass SRT Q(10), SSRT Q(20) for std       mean insig  AIC -5.472105
# all coef sig at 5%, pass SRT Q(10), SSRT Q(20) for sstd      mean insig  AIC -5.471632

t_a3=( 1.018e+00-1)/2.149e-02 
pv_a3 = 2*(1-pnorm(t_a3))
pv_a3
#inadequate sstd since symmetry assumption not rej

v3.2 = volatility(a3.2) 
resi3 = residuals(a3.2, standardize = T)
par(mfrow=c(3,1))
plot(year[-1], lrtn, type="l", xlab="year", ylab="log returns", main ="Top Glove Corp Bhd")
plot(year[-1], v3.2, type="l", xlab="year", ylab="volatility", main ="Top Glove Corp Bhd")
plot(year[-1], resi3, type="l", xlab="year", ylab="standardized residuals", main ="Top Glove Corp Bhd")

par(mfrow=c(1,3))
plot(a3.1)
13
0
plot(a3.2)
13
0
plot(a3.3)
13
0

auto.arima(x1^2)
#suggest GARCH(1,0)

#model a3.2 preferred

#######   Part II b          #####
a1.2  #GARCH(1,2) AIC -5.488314  # SMALLEST
a2.2  #GARCH(1,1) AIC -5.487884
a3.2  #GARCH(1,0) AIC -5.472105

# one-step ahead
par(mfrow=c(1,1))
predict(a1.2, n.ahead=1, plot = T)

qstd(0.025, nu = 2.504e+00)
qstd(0.975, nu = 2.504e+00)

#######   Part II c          ####

v_a1.2 = volatility(a1.2)
upp = 0 + 1.96*v_a1.2
low = 0 - 1.96*v_a1.2

plot(year[-1], lrtn, type="l", xlab="year", ylab="log returns", main ="Top Glove Corp Bhd")
lines(year[-1], upp, col="blue")
lines(year[-1], low, col="blue")
abline(h=0, col="red")


