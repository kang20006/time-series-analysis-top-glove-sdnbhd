library(fBasics)
library(fGarch)
library(forecast)

tg1 = read.csv("UECM3243_TopGlove.csv", header = T)

# to eliminate null value(s) and check numerical
sapply(X = tg1, FUN = function(x) sum(x=="null"))
tg = tg1[tg1$Adj.Close!="null",]

str(tg$Adj.Close)
# Check if adjust.close shows as numeric then do nothing
# If adjust.close shows as factor run
tp_adjcl = as.numeric(levels(tg$Adj.Close))[tg$Adj.Close]
# If adjust.close shows as char run
#tp_adjcl=as.numeric((tg$Adj.Close))
str(tp_adjcl)
tp_adjcl

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

range(lrtn)
x = seq(-0.3, 0.2, 0.001)   
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
normalTest(lrtn, method = "jb")
# Normality assumption for lrtn is rejected. 

######   1(C) Testing           ######

Box.test(lrtn, lag = 12, type = "Ljung")
# Do not rej H_0, serially uncorrelated
# X-squared = 12.232, df = 12, p-value = 0.4272

par(mfrow=c(1,2))
acf(lrtn,main="Daily Log Return",ylim=c(-0.04,1))
title(main="(a)",adj=0)
acf(abs(lrtn), main = "Absolute Daily Log Return",ylim=c(-0.04,1))
title(main="(b)",adj=0)
# From first 12 acf and pacf, weakly stationary and white noise

Box.test(abs(lrtn), lag = 12, type = "Ljung") 
# Rej H_0, it is dependent series
# X-squared = 283.85, df = 12, p-value < 2.2e-16

x1 = lrtn - mean(lrtn)

Box.test(x1^2, lag = 12, type = "Ljung")
# Rej H_0, exist ARCH effect
# X-squared = 42.17, df = 12, p-value = 3.118e-05

####   PART II      #####

# maybe we can perform a JB test to check for normality for x1^2
# so we not need to test for so many model
normalTest(x1^2, method = "jb")
# p Value: < 2.2e-16 
# reject h0 normality assumption not true, Gaussian innovation may not be accepted

#######   GARCH (1,2)        #####
a1.1 = garchFit(~garch(1,2), data = lrtn, include.mean = F, trace = F, cond.dist = "norm") 
# mu is not significant so exclude the mean.
summary(a1.1)
#           Estimate  Std. Error  t value Pr(>|t|)    
#   omega  1.871e-05   3.528e-06    5.303 1.14e-07 ***
#   alpha1 1.541e-01   1.857e-02    8.295  < 2e-16 ***
#   beta1  3.728e-01   6.359e-02    5.862 4.56e-09 ***
#   beta2  4.424e-01   5.966e-02    7.414 1.22e-13 ***
par(mfrow=c(1,3))
plot(a1.1)
10 #acf standardized residuals epsilon/sigma
#no spike
11 #acf squared standardized residuals 
#spike at lag 19
13 #qqnorm 
#deviate so much from normal
0

a1.2 = garchFit(~garch(1,2), data = lrtn, include.mean = F, trace = F, cond.dist = "std") 
summary(a1.2)
#         Estimate  Std. Error  t value Pr(>|t|)    
#   omega  1.048e-04   3.455e-05    3.034 0.002417 ** 
#   alpha1 5.674e-01   1.661e-01    3.416 0.000636 ***
#   beta1  4.043e-01   9.631e-02    4.198 2.69e-05 ***
#   beta2  1.412e-01   8.225e-02    1.717 0.085942 .  
#   shape  2.502e+00   1.692e-01   14.785  < 2e-16 ***
plot(a1.2)
10
11
13 
0

a1.3 = garchFit(~garch(1,2), data = lrtn, include.mean = F, trace = F, cond.dist = "sstd") 
summary(a1.3)
#           Estimate  Std. Error  t value Pr(>|t|)    
#   omega  1.044e-04   3.443e-05    3.031 0.002434 ** 
#   alpha1 5.644e-01   1.655e-01    3.410 0.000649 ***
#   beta1  4.045e-01   9.638e-02    4.197 2.71e-05 ***
#   beta2  1.415e-01   8.239e-02    1.718 0.085791 .  
#   skew   1.003e+00   1.467e-02   68.339  < 2e-16 ***
#   shape  2.504e+00   1.700e-01   14.729  < 2e-16 ***
plot(a1.3)
10
11
13
0

# all coef sig at 5%, pass SRT (10)(15)(20), SSRT (10)(15)(20) for Gaussian  mean insig  AIC -5.203805 BIC -5.195111
# all coef sig at 5%, pass SRT (10)(15), SSRT (10)(15)(20)     for std       mean insig  AIC -5.488945 BIC -5.478077
# all coef sig at 5%, pass SRT (10)(15) , SSRT (10)(15)(20)    for sstd      mean, beta 2 insig  AIC -5.488222 BIC -5.475181

t_a1 = (1.017e+00 -1)/2.123e-02
pv_a1 = 2* (1-pnorm(t_a1))
pv_a1
#p-value 0.4232743
#inadequate sstd since symmetry assumption do not rej
#choose std model a1.2


#
v1.2 = volatility(a1.2) 
resi1 = residuals(a1.2, standardize = T)
par(mfrow=c(3,1))
plot(year[-1], lrtn, type="l", xlab="year", ylab="log returns", main ="Log Return")
plot(year[-1], v1.2, type="l", xlab="year", ylab="volatility", main ="Volality")
plot(year[-1], resi1, type="l", xlab="year", ylab="standardized residuals", main ="Standardized residuals")

# compare QQ plot
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
par(mfrow=c(1,2))
acf(x1^2, main=expression(epsilon[t]^2))
pacf(x1^2, main=expression(epsilon[t]^2))
#acf and pacf suggest (1,1)
a2.1 = garchFit(~garch(1,1), data = lrtn, include.mean = F, trace = F, cond.dist = "norm")
summary(a2.1)
#           Estimate  Std. Error  t value Pr(>|t|)    
#   omega  1.404e-05   2.644e-06    5.311 1.09e-07 ***
#   alpha1 1.085e-01   1.390e-02    7.805 6.00e-15 ***
#   beta1  8.682e-01   1.445e-02   60.101  < 2e-16 ***
par(mfrow=c(1,3))
plot(a2.1)
10
11
13
0

a2.2 = garchFit(~garch(1,1), data = lrtn, include.mean = F, trace = F, cond.dist = "std")
summary(a2.2)
#           Estimate  Std. Error  t value Pr(>|t|)    
#   omega  1.098e-04   3.648e-05    3.009 0.002622 ** 
#   alpha1 5.514e-01   1.664e-01    3.314 0.000921 ***
#   beta1  5.534e-01   7.420e-02    7.459 8.73e-14 ***
#   shape  2.491e+00   1.683e-01   14.801  < 2e-16 ***
par(mfrow=c(1,3))
plot(a2.2)
10
11
13
0

a2.3 = garchFit(~garch(1,1), data = lrtn, include.mean = F, trace = F, cond.dist = "sstd")
summary(a2.3)
#           Estimate  Std. Error  t value Pr(>|t|)    
#   omega  1.094e-04   3.636e-05    3.008 0.002630 ** 
#   alpha1 5.487e-01   1.658e-01    3.309 0.000938 ***
#   beta1  5.538e-01   7.430e-02    7.453  9.1e-14 ***
#   skew   1.002e+00   1.460e-02   68.649  < 2e-16 ***
#   shape  2.493e+00   1.691e-01   14.747  < 2e-16 ***
par(mfrow=c(1,3))
plot(a2.3)
10
11
13
0

# all coef sig at 5%, pass SRT (10)(15)(20) , SSRT (10)(15)(20) for Gaussian  mean insig AIC -5.194407 BIC -5.187887
# all coef sig at 5%, pass SRT (10)(15)(20) , SSRT (10)(15)(20) for std       mean insig  AIC -5.488516 BIC -5.479822
# all coef sig at 5%, pass SRT (10)(15)(20), SSRT (10)(15)(20) for sstd       mean insig  AIC -5.487791 BIC -5.476923

t_a2=(1.002e+00-1)/1.460e-02 
pv_a2=2*(1-pnorm(t_a2))
pv_a2
#inadequate sstd since symmetry assumption do not rej
#p-value=0.8910416

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
#           Estimate  Std. Error  t value Pr(>|t|)    
#   mu     8.458e-04   3.236e-04    2.614  0.00896 ** 
#   omega  2.402e-04   9.443e-06   25.432  < 2e-16 ***
#   alpha1 4.367e-01   4.973e-02    8.782  < 2e-16 ***
plot(a3.1)
10
11
13
0

a3.2 = garchFit(~garch(1,0), data = lrtn, include.mean = F, trace = F, cond.dist = "std")
summary(a3.2)
#           Estimate  Std. Error  t value Pr(>|t|)    
#   omega  0.0003965   0.0001010    3.926 8.65e-05 ***
#   alpha1 0.8467813   0.2548581    3.323 0.000892 ***
#   shape  2.4485973   0.1634926   14.977  < 2e-16 ***
plot(a3.2)
10
11
13
0

a3.3 = garchFit(~garch(1,0), data = lrtn, include.mean = F, trace = F, cond.dist = "sstd")
summary(a3.3)
#         Estimate  Std. Error  t value Pr(>|t|)    
# omega  0.0003939   0.0000996    3.955 7.65e-05 ***
# alpha1 0.8386810   0.2516725    3.332 0.000861 ***
# skew   1.0046810   0.0144695   69.434  < 2e-16 ***
# shape  2.4531708   0.1645415   14.909  < 2e-16 ***
plot(a3.3)
10
11
13
0

# all coef sig at 5%, pass SRT Q(10), SSRT Q(20) for Gaussian              AIC -5.187530 BIC -5.181009
# all coef sig at 5%, pass SRT Q(10), SSRT Q(20) for std       mean insig  AIC -5.472105 BIC -5.466276
# all coef sig at 5%, pass SRT Q(10), SSRT Q(20) for sstd      mean insig  AIC -5.471632 BIC -5.463406

t_a3=(1.0046810-1)/0.0144695
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
a1.2  #GARCH(1,2) AIC -5.488945  BIC -5.478077 # SMALLEST AIC
a2.2  #GARCH(1,1) AIC -5.488516  BIC -5.479822 # SMALLEST BIC
a3.2  #GARCH(1,0) AIC -5.472105  BIC -5.466276

# one-step ahead
par(mfrow=c(1,1))
predict(a1.2, n.ahead=1, plot = T)

qstd(0.025, nu = 2.504e+00) #-1.601902
qstd(0.975, nu = 2.504e+00) #1.601902
#predictive interval
#=0+-1.601902*0.02503842
#=0+-0.04006124 

#######   Part II c          ######

v_a1.2 = volatility(a1.2)
upp = 0 + 1.96*v_a1.2
low = 0 - 1.96*v_a1.2

plot(year[-1], lrtn, type="l", xlab="year", ylab="log returns", main ="Top Glove Corp Bhd")
lines(year[-1], upp, col="blue")
lines(year[-1], low, col="blue")
abline(h=0, col="red")
#or
plot(a1.2)
3
0
