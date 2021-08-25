library(fBasics)
library(fGarch)

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
year= c(1:4454)/252 + 2001.25

######   time series plot   #####
par(mfrow=c(1,1))
plot(year, tp_adjcl, main = "Top Glove Corp Bhd", type = "l",
     ylab = "Daily Adjusted Closing Price")

plot(c(1:521)/252 + 2001.25, tp_adjcl[1:521], main = "Top Glove Corp Bhd", type = "l",
     ylab = "Daily Adjusted Closing Price",xlim=c(2001,2004))
year1=c(1:521)/252 + 2001.25
plot(year1[-1], diff(log(tp_adjcl[1:521])), main = "Top Glove Corp Bhd", type = "l",
     ylab = "Daily log returns")

plot(c(2481:3462)/252 + 2001.25, tp_adjcl[2481:3462], main = "Top Glove Corp Bhd", type = "l",
     ylab = "Daily Adjusted Closing Price")
year2=c(2481:3462)/252 + 2001.25
plot(year2[-1], diff(log(tp_adjcl[2481:3462])), main = "Top Glove Corp Bhd", type = "l",
     ylab = "Daily log returns")

plot(c(3463:4454)/252 + 2001.25, tp_adjcl[3463:4454], main = "Top Glove Corp Bhd", type = "l",
     ylab = "Daily Adjusted Closing Price")
year3=c(3463:4454)/252 + 2001.25
plot(year3[-1], diff(log(tp_adjcl[3463:4454])), main = "Top Glove Corp Bhd", type = "l",
     ylab = "Daily log returns")

######   ACF PACF tp_adjcl    ######
par(mfrow=c(1,2))
acf(tp_adjcl, main="Top Glove Corp Bhd")
pacf(tp_adjcl, main="Top Glove Corp Bhd")
#acf decays very slowly
#pacf lag-1 is approximately 1.0, first reg diff required

######   log return plot   #######
par(mfrow=c(1,1))
plot(year[-1], lrtn, main = "Top Glove Corp Bhd", type = "l",
     xlab= "Year", ylab = "Daily Log Returns")
#appears to be weekly stationary in mean and variance


######   ACF PACF lrtn   ######
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
# The daily log returns is symmetrical

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

######   1(C) Testing   ######

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

######   PART II   #####
#test for Garch
m1=garchFit(~garch(5,9),data=lrtn,trace=F,cond.dist = "sstd")
m1
summary(m1)
