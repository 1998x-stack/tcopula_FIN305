#GARCH model
setwd("C:/Users/surface/Desktop/FIN305/CW")

library("zoo")

library("rugarch")

library("xts")
rm(list=ls())
WFC_data=read.csv(file="WFC.csv",header = TRUE)
VOD_data=read.csv(file="VOD.csv",header = TRUE)
Price_WFC=WFC_data[,6]
Price_VOD=VOD_data[,6]
ret_WFC=diff(log(Price_WFC))#ret_EDU:Log return
ret_VOD=diff(log(Price_VOD))#ret_EDU:Log return
#GARCH-type Model estimate
GARCH11=ugarchspec(variance.model = list(model="sGARCH",
                                    garchOrder=c(1,1)),
                       mean.model = list(armaOrder=c(0,0),
                                    include.mean=FALSE),
               distribution.model = "norm")
GARCH12=ugarchspec(variance.model = list(model="sGARCH",
                                    garchOrder=c(2,1)),
                       mean.model = list(armaOrder=c(0,0),
                                    include.mean=FALSE),
               distribution.model = "norm")
GARCH21=ugarchspec(variance.model = list(model="sGARCH",
                                    garchOrder=c(1,2)),
                       mean.model = list(armaOrder=c(0,0),
                                    include.mean=FALSE),
              distribution.model = "norm")

#GARCH_fit_WFC
GARCH11_fit_WFC=ugarchfit(GARCH11,ret_WFC)                    
show(GARCH11_fit_WFC)
coef(GARCH11_fit_WFC)
infocriteria(GARCH11_fit_WFC)

GARCH12_fit_WFC=ugarchfit(GARCH12,ret_WFC)                    
show(GARCH12_fit_WFC)
coef(GARCH12_fit_WFC)
infocriteria(GARCH12_fit_WFC)

GARCH21_fit_WFC=ugarchfit(GARCH21,ret_WFC)                    
show(GARCH21_fit_WFC)
coef(GARCH21_fit_WFC)
infocriteria(GARCH21_fit_WFC)

#GARCH_fit_VOD
GARCH11_fit_VOD=ugarchfit(GARCH11,ret_VOD)                    
show(GARCH11_fit_VOD)
coef(GARCH11_fit_VOD)
infocriteria(GARCH11_fit_VOD)

GARCH12_fit_VOD=ugarchfit(GARCH12,ret_VOD)                    
show(GARCH12_fit_VOD)
coef(GARCH12_fit_VOD)
infocriteria(GARCH12_fit_VOD)

GARCH21_fit_VOD=ugarchfit(GARCH21,ret_VOD)                    
show(GARCH21_fit_VOD)
coef(GARCH21_fit_VOD)
infocriteria(GARCH21_fit_VOD)

#extract the _EDestimated conditional volatility, using GARCH(1,1)

ConVolatility_WFC=sigma(GARCH11_fit_WFC)
return_fitted_WFC=fitted(GARCH11_fit_WFC)
date_WFC=WFC_data[,1]
time_WFC=as.Date(date_WFC[2:length(Price_WFC)])
TSvol_WFC=xts(ConVolatility_WFC,time_WFC)
plot(TSvol_WFC,format.labels="%m-%d-%Y",ylim=c(0,0.13))


ConVolatility_VOD=sigma(GARCH11_fit_VOD)
return_fitted_VOD=fitted(GARCH11_fit_VOD)
date_VOD=VOD_data[,1]
time_VOD=as.Date(date_VOD[2:length(Price_VOD)])
TSvol_VOD=xts(ConVolatility_VOD,time_VOD)
plot(TSvol_VOD,format.labels="%m-%d-%Y",ylim=c(0,0.12))




library("fitdistrplus")
library("copula")
library("readxl")

dt_G <- function(x, mean, sd, nu){#density
  dt((x-mean)/sd,nu)/sd
}
pt_G <- function(q, mean, sd, nu){#CDF
  pt((q-mean)/sd,nu)
}
qt_G <- function(x, mean, sd, nu){#quantile
  qt(x,nu)*sd+mean
}

#standardize residuals xts
u_WFC=residuals(GARCH11_fit_WFC, standardize=T)
u_VOD=residuals(GARCH11_fit_VOD, standardize=T)

#standardize the log returns number
volVals_WFC=GARCH11_fit_WFC@fit$sigma
u_WFC =(ret_WFC-coef(GARCH11_fit_WFC)[1])/volVals_WFC
volVals_VOD=GARCH11_fit_VOD@fit$sigma
u_VOD=(ret_VOD-coef(GARCH11_fit_VOD)[1])/volVals_VOD


#Fitting marginal distribution:
# Fitting Gaussian distribution & t distribution for WFC (three parameters:mean, sd and nu)
# u_WFC is the standardized residuals of log returns
norm_u_WFC <- fitdist(u_WFC,"norm")
#fitted margins
summary(norm_u_WFC)
ft_u_WFC <- fitdist(u_WFC,"t_G", start=list(mean=mean(u_WFC), sd=sd(u_WFC), nu=5))
#fitted margins
summary(ft_u_WFC)

# Fitting Gaussian distribution and t distribution for VOD (three parameters:mean, sd and nu)
# u_VOD is the standardized residuals of log returns
norm_u_VOD=fitdist(u_VOD,"norm", start=list(mean=mean(u_VOD), sd=sd(u_VOD)))
#fitted margins
summary(norm_u_VOD)
ft_u_VOD=fitdist(u_VOD,"t_G", start=list(mean=mean(u_VOD), sd=sd(u_VOD), nu=5))
#fitted margins
summary(ft_u_VOD)



#histogram plot and density of t and Gaussian
histt <- function(data, ft, location = "topright", legend.cex = 1, xlab){
  # parameters of the fitted t distribution
  mean_t <- as.list(ft$estimate)$mean
  sd_t <- as.list(ft$estimate)$sd
  nu_t <- as.list(ft$estimate)$nu
  
  # drawing histogram and assigning to h so that we can get the breakpoints between histogram cells (we will use it!)
  h <- hist(data, breaks=30)
  
  # x sequence for the additional plots on the histogram
  x_seq <- seq(-3,3,length=10000)
  
  # y sequence: density of fitted t distr. at x_seq
  yhistt <- dt_G(x_seq, mean=mean_t, sd=sd_t, nu=nu_t)
  
  # y sequence: density of normal distr. with mean and standard deviation of log-returns at x_seq
  # Note. I did not fit the normal distribution as we would get almost the same mean and standard deviation
  yhistNorm <- dnorm(x_seq, mean=mean(data), sd=sd((data)))
  
  # drawing histogram but this time we draw its density as y axis (freq=FALSE)
  hist(data,freq=FALSE,xlab=xlab, ylab="Density",breaks=h$breaks,
  main=paste(""),ylim=c(0,0.6),xlim=c(-4,4),cex.lab=1.5,las=1,cex.lab=1.5)

  # adding the density of the fitted t distribution at x_seq
lines(x_seq, yhistt, col=4)
# adding the density of the normal distribution at x_seq
lines(x_seq, yhistNorm, lty = "dashed")

# Legend

tmp.text <- c("t", "Gaussian")
legend("topright", legend = tmp.text, cex = legend.cex, lty = c(1,2), col=c(4,1))	
}



# Drawing the histogram using the function that we wrote
histt(u_WFC, ft_u_WFC,xlab="Standardized residuals - WFC")
# Drawing the histogram using the function that we wrote
histt(u_VOD, ft_u_VOD, xlab="Standardized residuals - VOD")


# Fit copula:
u=matrix(nrow=length(u_WFC),ncol=2)
u[,1]=pt_G(u_WFC, mean=as.list(ft_u_WFC$estimate)$mean,
            sd=as.list(ft_u_WFC$estimate)$sd, nu=as.list(ft_u_WFC$estimate)$nu)
u[,2]=pt_G(u_VOD, mean=as.list(ft_u_VOD$estimate)$mean,
            sd=as.list(ft_u_VOD$estimate)$sd, nu=as.list(ft_u_VOD$estimate)$nu)
norm.cop=normalCopula(dim=2,dispstr="un")
n.cop=fitCopula(norm.cop,u,method="ml")
t.copula=tCopula(dim=2,dispstr="un",df=5)
t.cop=fitCopula(t.copula,u,method="ml")
#Comparison between normal copula and t copula:
summary(n.cop)
summary(t.cop)
# the copula correlation between stock log returns
coef(n.cop)
coef(t.cop)
#the density of the normal copula 
persp(normalCopula(dim=2,coef(n.cop)),dCopula)
#the density of the t copula
rho <- coef(t.cop)[1]
df <- coef(t.cop)[2]
persp(tCopula(dim=2,rho,df=df),dCopula)

#É¢µãÍ¼ t
n <- rCopula(8000,tCopula(dim=2,rho,df=df))
plot(n[,1],n[,2],pch='.',col='blue')
#É¢µãÍ¼ normal
rho_t <- coef(n.cop)[1]
b <- rCopula(8000,normalCopula(dim=2,rho_t))
plot(b[,1],b[,2],pch='.',col='blue')



mean_r1=as.list(ft_u_WFC$estimate)$mean
sd_r1=as.list(ft_u_WFC$estimate)$sd
nu_r1=as.list(ft_u_WFC$estimate)$nu

mean_r2=as.list(ft_u_VOD$estimate)$mean
sd_r2=as.list(ft_u_VOD$estimate)$sd
nu_r2=as.list(ft_u_VOD$estimate)$nu

# create the correlation matrix using the fitted normal copula
cor=matrix(data=c(1,0.4149913,0.4149913,1), nrow=2,ncol=2)


d=2 # we have 2 stocks in one portfolio
L=t(chol(cor))
set.seed(1234)
N=10000
Sim_R=matrix(nrow=N,ncol=2)

for (i in 1:N){
  z=rnorm(d)
  z_tilde=L%*%z
  R1=qt_G(pnorm(z_tilde[1]),mean=mean_r1, sd=sd_r1, nu=nu_r1)
  R2=qt_G(pnorm(z_tilde[2]),mean=mean_r2, sd=sd_r2, nu=nu_r2)
  Sim_R[i,]=c(R1,R2)
}
portfolio=rowMeans(Sim_R)
sort_p=sort(portfolio)
alpha=0.05 #confidence level
VaR_p=sort_p[floor(length(portfolio)*alpha)]

#ret_p=0.5*ret_VOD+0.5*ret_WFC
#sum(ret_p[1:2515]<VaR_p)
#sum=64>0.1*586=58.6, therefore we test the following
#1-pbinom(63,586,0.1)
ee=sort_p[1:floor(length(portfolio)*alpha)-1]





