####### Bass Cumulative Adoption Process #####

f <- function(t, m, p, q) {(m - m * exp(-(p + q) * t))/(1 + (q/p) * exp(-(q + p) * t))}


########################### Generation 1 ################################

# data

t0 = 1

t1 <- t0:(t0+23)

a1 <- c(1811, 286, 136, 124, 184, 101, 65, 40, 55, 30, 40, 36, 22, 24, 16, 24, 36, 22, 14, 8, 7, 9, 8, 9)

a_cum_1 <- cumsum(a1)

plot(a1 ~ t1)

plot(a_cum_1 ~ t1)

# Model

AC1<- function(t,m12,p12,q12,d11,d12){
  
    I(t<t0)*f(t,m12,p12,q12)+
    I(t==t0)*(f(t,m12,p12,q12) + d11)+
    I(t0+4>t & t>t0)*(f(t,m12,p12,q12) + d11)+
    I(t==t0+4)*(f(t,m12,p12,q12) + d12 + d11)+
    I(t>t0+4)*(f(t,m12,p12,q12) + d11+ d12)
}

# Fit AC1 with nls

BB1<- nls(a_cum_1 ~ AC1(t1,m12,p12,q12,d11,d12),  
          algorithm = "port", 
          start = list(m12=1367.37, p12=0.221058,q12=0., d11=1508.73,d12=100),
          lower = c(0,0,0,0,0))

summary(BB1)

# Plot the fit

plot(t1,a_cum_1)
lines(t1,predict(BB1),col="blue",lwd=2)

plot(t1, a_cum_1, xlim = c(0,24), ylim = c(0,3500))
lines(t1, predict(BB1), col="blue",lwd=2, xlim = c(0,24), ylim = c(0,3500))

#########################################################################

dm1 <- c(1,rep(0,23))
dm2 <- c(rep(0, 4), 1,rep(0,19))

dat <- cbind(a_cum_1, t1, 1, t1-1, dm1, dm2, (t1-1)^2, (t1-1)^4)

## Moments

moments <- function(par,data) {
  y <- as.numeric(data[,1])
  x <- data.matrix(data[,3:8])
  m <- x * as.vector((y - AC1(data[,2],par[1], par[2], par[3], par[4], par[5])))
  return(cbind(m))
}


#

library(gmm)
library(sandwich)
#  type=c("twoStep","cue","iterative")

GMM_Estimation_G1 <- gmm(moments,
                         x = dat, type = "twoStep",
                         method = "L-BFGS-B",
                         t0=c(m12=1418.43, p12=0.1752, q12=0.00000001, d11=1639, d12=45),
                         lower = c(0,0,0,0,0),
                         upper = c(2000,5,5,2000,200)
                         )

# Methods: “Nelder-Mead”, “BFGS”, “CG”, “L-BFGS-B”, “SANN”, “Brent”

summary(GMM_Estimation_G1)

# Plot the estimated model

TB <- AC1(
  1:24,
  coef(GMM_Estimation_G1)[1], 
  coef(GMM_Estimation_G1)[2], 
  coef(GMM_Estimation_G1)[3], 
  coef(GMM_Estimation_G1)[4], 
  coef(GMM_Estimation_G1)[5]
)

plot(a_cum_1 ~ t1, ylab = "Adoptions", xlab = "Time", ylim = c(0, 4000), xlim = c(0, 25),type="b",lwd=2)
par(new=TRUE)
plot(TB ~ t1, col="red", ylab = "Adoptions", xlab = "Time", ylim = c(0, 4000), xlim = c(0, 25),type="b",lwd=2)



########################### Generation 2 ################################

# data

t0 = 1

t2 <- t0:(t0+23)

a2 <- c(3622, 950, 712, 686, 1644, 910, 487, 390, 234, 358, 
        305, 355, 397, 302, 296, 261, 407, 249, 163, 136, 
        119, 100, 115, 137)

a_cum_2 <- cumsum(a2)

plot(a2 ~ t2)

plot(a_cum_2 ~ t2, ylim = c(0,15000))

# Model

AC2<- function(t,m22,p22,q22,d21,d22){
  
    I(t<t0)*f(t,m22,p22,q22)+
    I(t==t0)*(f(t,m22,p22,q22) + d21)+
    I(t0+4>t & t>t0)*(f(t,m22,p22,q22) + d21)+
    I(t==t0+4)*(f(t,m22,p22,q22) + d22 + d21)+
    I(t>t0+4)*(f(t,m22,p22,q22) + d21+ d22)

}

# Fit AC2 with nls

BB2<- nls(a_cum_2 ~ AC2(t2,m22,p22,q22,d21,d22),  
          algorithm = "port", 
          start = list(m22=1800, p22=0.221058, q22=0., d21=1508.73, d22=100),
          lower = c(0,0,0,0,0))

summary(BB2)

# Plot the fit

plot(t2, a_cum_2)
lines(t2,predict(BB2),col="blue",lwd=2)

########################### G2 GMM ######################################

dm1 <- c(1,rep(0,23))
dm2 <- c(rep(0, 4), 1,rep(0,19))

#Question: generate random variable as IV?
random_instroment <- c(1,20,12, 22, 16, 30, 32,  42,
                       38,  36,  52,  55,  45,  62,  64,
                       72,  70,  80,  82,  85,  90,  96, 100, 102)

dat <- cbind(a_cum_2, t1, 1, t1-1, dm1, dm2, (t1-1)^2, (t1-1)^4, random_instroment^2)

## Moments

moments <- function(par,data) {
  y <- as.numeric(data[,1])
  x <- data.matrix(data[,4:9])
  m <- x * as.vector((y - AC2(data[,2],par[1], par[2], par[3], par[4], par[5])))
  return(cbind(m))
}

#

#library(gmm)
#library(sandwich)
coef(BB2)[5]

GMM_Estimation_G2 <- gmm(moments,
                         x = dat,
                         method = "L-BFGS-B",
                         t0=c(m12=10790, p12=0.07, q12=0.00000001, d11=2961, d12=1252),
                         lower = c(10000,0.01,0,2500,1252),
                         upper = c(12000,5,5,4000,1500)
)

# Methods: “Nelder-Mead”, “BFGS”, “CG”, “L-BFGS-B”, “SANN”, “Brent”

summary(GMM_Estimation_G2)

# Plot the estimated model

TB <- AC1(
  1:24,
  coef(GMM_Estimation_G2)[1], 
  coef(GMM_Estimation_G2)[2], 
  coef(GMM_Estimation_G2)[3], 
  coef(GMM_Estimation_G2)[4], 
  coef(GMM_Estimation_G2)[5]
)

plot(a_cum_2 ~ t1, ylab = "Adoptions", xlab = "Time", ylim = c(0, 15000), xlim = c(0, 25),type="b",lwd=2)
par(new=TRUE)
plot(TB ~ t1, col="green", ylab = "Adoptions", xlab = "Time", ylim = c(0, 15000), xlim = c(0, 25),type="b",lwd=2)

########################### Generation 3 ################################

# data

t0 = 1

t3 <- t0:(t0+23)

a3 <- c(4773, 982, 763, 1133, 
        2428, 840, 432, 432, 
        404, 312, 316, 357, 
        411, 400, 325, 310, 
        296, 369, 245, 144, 
        110, 112, 121, 121)

a_cum_3 <- cumsum(a3)

plot(a3 ~ t3)

plot(a_cum_3 ~ t3)

# The Model

AC3 <- function(t,m32,p32,q32,d31,d32){
  
  I(t<t0)*f(t,m32,p32,q32)+
    I(t==t0)*(f(t,m32,p32,q32) + d31)+
    I(t0+4>t & t>t0)*(f(t,m32,p32,q32) + d31)+
    I(t==t0+4)*(f(t,m32,p32,q32) + d32 + d31)+
    I(t>t0+4)*(f(t,m32,p32,q32) + d31+ d32)
}

# Fit AC2 with nls

BB3 <- nls(a_cum_3 ~ AC3(t3,m32,p32,q32,d31,d32),  
          algorithm = "port", 
          start = list(m32=1800, p32=0.221058, q32=0., d31=1508.73, d32=100),
          lower = c(0,0,0,0,0))

summary(BB3)

# Plot the fit

plot(t3, a_cum_3)
lines(t3,predict(BB3),col="blue",lwd=2)

########################### Gen 3 GMM #################################

dm1 <- c(1,rep(0,23))
dm2 <- c(rep(0, 4), 1,rep(0,19))

random_instroment <- c(1,20,12, 22, 16, 30, 32,  42,
                       38,  36,  52,  55,  45,  62,  64,
                       72,  70,  80,  82,  85,  90,  96, 100, 102)

dat <- cbind(a_cum_3, t1, 1, t1-1, dm1, dm2, (t1-1)^2, (t1-1)^4, random_instroment^2)

## Moments

moments <- function(par,data) {
  y <- as.numeric(data[,1])
  x <- data.matrix(data[,4:9])
  m <- x * as.vector((y - AC3(data[,2],par[1], par[2], par[3], par[4], par[5])))
  return(cbind(m))
}

#

#library(gmm)
#library(sandwich)
coef(BB3)[5]

GMM_Estimation_G3 <- gmm(moments,
                         x = dat,
                         method = "L-BFGS-B",
                         t0=c(m32=11700, p32=0.08, q32=0.00000001, d31=4077, d32=2070),
                         lower = c(10000,0.01,0,2500,1252),
                         upper = c(12000,5,5,4000,1500)
)

# Methods: “Nelder-Mead”, “BFGS”, “CG”, “L-BFGS-B”, “SANN”, “Brent”

summary(GMM_Estimation_G3)

# Plot the estimated model

TB <- AC3(
  1:24,
  coef(GMM_Estimation_G3)[1], 
  coef(GMM_Estimation_G3)[2], 
  coef(GMM_Estimation_G3)[3], 
  coef(GMM_Estimation_G3)[4], 
  coef(GMM_Estimation_G3)[5]
)

plot(a_cum_3 ~ t1, ylab = "Adoptions", xlab = "Time", ylim = c(0, 20000), xlim = c(0, 25),type="b",lwd=2)
par(new=TRUE)
plot(TB ~ t1, col="green", ylab = "Adoptions", xlab = "Time", ylim = c(0, 20000), xlim = c(0, 25),type="b",lwd=2)


########################### Generation 4 ################################

# data

t0 = 1

t4 <- t0:(t0+23)

a4 <- c(5664, 2269, 1106, 1020, 
        2526, 1052, 605, 484, 
        318, 342, 325, 426, 
        422, 366, 312, 347, 
        478, 321, 175, 118, 
        116, 123, 134, 131)

a_cum_4 <- cumsum(a4)

plot(a4 ~ t4)

plot(a_cum_4 ~ t4)

# The Model

AC4 <- function(t,m42,p42,q42,d41,d42){
  
  I(t<t0)*f(t,m42,p42,q42)+
    I(t==t0)*(f(t,m42,p42,q42) + d41)+
    I(t0+4>t & t>t0)*(f(t,m42,p42,q42) + d41)+
    I(t==t0+4)*(f(t,m42,p42,q42) + d42 + d41)+
    I(t>t0+4)*(f(t,m42,p42,q42) + d41+ d42)
}

# Fit AC2 with nls

BB4 <- nls(a_cum_4 ~ AC4(t4,m42,p42,q42,d41,d42),  
           algorithm = "port", 
           start = list(m42=1367.37, p42=0.221058, q42=0., d41=1508.73, d42=100),
           lower = c(0,0,0,0,0))

summary(BB4)

# Plot the fit

plot(t4, a_cum_4)
lines(t4,predict(BB4),col="blue",lwd=2)

########################### Gen 4 GMM ###################################

dm1 <- c(1,rep(0,23))
dm2 <- c(rep(0, 4), 1,rep(0,19))

random_instroment <- c(1,20,12, 22, 16, 30, 32,  42,
                       38,  36,  52,  55,  45,  62,  64,
                       72,  70,  80,  82,  85,  90,  96, 100, 102)

dat <- cbind(a_cum_4, t1, 1, t1-1, dm1, dm2, (t1-1)^2, (t1-1)^4, random_instroment^2)

## Moments

moments <- function(par,data) {
  y <- as.numeric(data[,1])
  x <- data.matrix(data[,4:9])
  m <- x * as.vector((y - AC4(data[,2],par[1], par[2], par[3], par[4], par[5])))
  return(cbind(m))
}

#

#library(gmm)
#library(sandwich)
coef(BB4)[5]

GMM_Estimation_G4 <- gmm(moments,
                         x = dat,
                         method = "L-BFGS-B",
                         t0=c(m42=13200, p42=0.11, q42=0.00000001, d41=4900, d42=1660),
                         lower = c(10000,0.01,0,2500,1252),
                         upper = c(14000,5,5,6000,1500)
)

# Methods: “Nelder-Mead”, “BFGS”, “CG”, “L-BFGS-B”, “SANN”, “Brent”

summary(GMM_Estimation_G4)

# Plot the estimated model

TB <- AC3(
  1:24,
  coef(GMM_Estimation_G4)[1], 
  coef(GMM_Estimation_G4)[2], 
  coef(GMM_Estimation_G4)[3], 
  coef(GMM_Estimation_G4)[4], 
  coef(GMM_Estimation_G4)[5]
)

plot(a_cum_4 ~ t1, ylab = "Adoptions", xlab = "Time", ylim = c(0, 20000), xlim = c(0, 25),type="b",lwd=2)
par(new=TRUE)
plot(TB ~ t1, col="green", ylab = "Adoptions", xlab = "Time", ylim = c(0, 20000), xlim = c(0, 25),type="b",lwd=2)



########################### Generation 5 ################################

# The data

t0 = 1

t5 <- t0:(t0+23)

a5 <- c(
  6290, 1267, 1065, 1645, 
  3458, 1222, 649, 506, 
  387, 428, 399, 693, 
  587, 473, 402, 385, 
  519, 336, 179, 166, 
  136, 196, 159, 191
)

a_cum_5 <- cumsum(a5)

plot(a5 ~ t5)

plot(a_cum_5 ~ t5)

# The Model

AC5 <- function(t,m52,p52,q52,d51,d52){
  
    I(t<t0)*f(t,m52,p52,q52)+
    I(t==t0)*(f(t,m52,p52,q52) + d51)+
    I(t0+4>t & t>t0)*(f(t,m52,p52,q52) + d51)+
    I(t==t0+4)*(f(t,m52,p52,q52) + d52 + d51)+
    I(t>t0+4)*(f(t,m52,p52,q52) + d51+ d52)
  
}

# Fit AC5 with nls

BB5 <- nls(a_cum_5 ~ AC5(t5,m52,p52,q52,d51,d52),  
           algorithm = "port", 
           start = list(m52=12191.7, p52=0.0842829,  q52=0., d51=3745.45,  d52=1694.52),
           lower = c(0,0,0,0,0))

summary(BB5)

# Plot the fit

plot(t5, a_cum_5)
lines(t5,predict(BB5),col="blue",lwd=2)

########################### Gen 5 GMM ###################################

dm1 <- c(1,rep(0,23))
dm2 <- c(rep(0, 4), 1,rep(0,19))

random_instroment <- c(1,20,12, 22, 16, 30, 32,  42,
                       38,  36,  52,  55,  45,  62,  64,
                       72,  70,  80,  82,  85,  90,  96, 100, 102)

dat <- cbind(a_cum_5, t1, 1, t1-1, dm1, dm2, (t1-1)^2, (t1-1)^4, random_instroment^2)

## Moments

moments <- function(par,data) {
  y <- as.numeric(data[,1])
  x <- data.matrix(data[,5:9])
  m <- x * as.vector((y - AC5(data[,2],par[1], par[2], par[3], par[4], par[5])))
  return(cbind(m))
}

#

#library(gmm)
#library(sandwich)
coef(BB5)[5]

GMM_Estimation_G5 <- gmm(moments,
                         x = dat,
                         method = "L-BFGS-B",
                         t0=c(m52=15400, p52=0.088, q52=0.00000001, d51=5190, d52=2990),
                         lower = c(10000,0.01,0,2500,1252),
                         upper = c(16000,5,5,6000,3500)
)

# Methods: “Nelder-Mead”, “BFGS”, “CG”, “L-BFGS-B”, “SANN”, “Brent”

summary(GMM_Estimation_G5)

# Plot the estimated model

TB <- AC5(
  1:24,
  coef(GMM_Estimation_G5)[1], 
  coef(GMM_Estimation_G5)[2], 
  coef(GMM_Estimation_G5)[3], 
  coef(GMM_Estimation_G5)[4], 
  coef(GMM_Estimation_G5)[5]
)

plot(a_cum_5 ~ t1, ylab = "Adoptions", xlab = "Time", ylim = c(0, 25000), xlim = c(0, 25),type="b",lwd=2)
par(new=TRUE)
plot(TB ~ t1, col="green", ylab = "Adoptions", xlab = "Time", ylim = c(0, 25000), xlim = c(0, 25),type="b",lwd=2)


########################### Generation 6 ################################

# Data

t0 = 1

t6 <- t0:(t0+20)

a6 <- c(
  8193, 1306, 1012, 1309,
  3252, 1399, 772, 628, 
  470, 468, 649, 552, 
  470, 355, 412, 428, 
  573, 474, 200, 160, 132
)

a_cum_6 <- cumsum(a6)

plot(a6 ~ t6)

plot(a_cum_6 ~ t6)

# Model

AC6 <- function(t,m62,p62,q62,d61,d62){
  
    I(t<t0)*f(t,m62,p62,q62)+
    I(t==t0)*(f(t,m62,p62,q62) + d61)+
    I(t0+4>t & t>t0)*(f(t,m62,p62,q62) + d61)+
    I(t==t0+4)*(f(t,m62,p62,q62) + d62 + d61)+
    I(t>t0+4)*(f(t,m62,p62,q62) + d61+ d62)
  
}

# Fit AC6 with nls

BB6 <- nls(a_cum_6 ~ AC6(t6,m62,p62,q62,d61,d62),  
           algorithm = "port", 
           start = list(m62=12191.7, p62=0.0842829,  q62=0., d61=3745.45,  d62=1694.52),
           lower = c(0,0,0,0,0))

summary(BB6)

# Plot the fit

plot(t6, a_cum_6)
lines(t6,predict(BB6),col="blue",lwd=2)

########################### Gen 6 GMM ###################################

dm1 <- c(1,rep(0,23))
dm2 <- c(rep(0, 4), 1,rep(0,19))

random_instroment <- c(1,20,12, 22, 16, 30, 32,  42,
                       38,  36,  52,  55,  45,  62,  64,
                       72,  70,  80,  82,  85,  90,  96, 100, 102)

dat <- cbind(a_cum_6, t1, 1, t1-1, dm1, dm2, (t1-1)^2, (t1-1)^4-1, random_instroment^2)

## Moments

moments <- function(par,data) {
  y <- as.numeric(data[,1])
  x <- data.matrix(data[,3:9])
  m <- x * as.vector((y - AC6(data[,2],par[1], par[2], par[3], par[4], par[5])))
  return(cbind(m))
}

#

#library(gmm)
#library(sandwich)
coef(BB6)[5]

GMM_Estimation_G6 <- gmm(moments,
                         x = dat,
                         method = "L-BFGS-B",
                         t0=c(m62=16400, p62=0.088, q62=0.00000001, d61=6900, d62=2750),
                         lower = c(10000,0.01,0,2500,1252),
                         upper = c(17000,5,5,8000,3500)
)

# Methods: “Nelder-Mead”, “BFGS”, “CG”, “L-BFGS-B”, “SANN”, “Brent”

#Question: estimation is not good.
summary(GMM_Estimation_G6)

# Plot the estimated model

TB <- AC5(
  1:24,
  coef(GMM_Estimation_G6)[1], 
  coef(GMM_Estimation_G6)[2], 
  coef(GMM_Estimation_G6)[3], 
  coef(GMM_Estimation_G6)[4], 
  coef(GMM_Estimation_G6)[5]
)

plot(a_cum_6 ~ t6, ylab = "Adoptions", xlab = "Time", ylim = c(0, 25000), xlim = c(0, 25),type="b",lwd=2)
par(new=TRUE)
plot(TB[1:length(t6)] ~ t6, col="green", ylab = "Adoptions", xlab = "Time", ylim = c(0, 25000), xlim = c(0, 25),type="b",lwd=2)



########################### Generation 7 ################################

# Data

t0 = 1

t7 <- t0:(t0+8)

a7 <- c(11081, 1324, 948, 1376, 4118, 1852, 863, 624, 1002)

plot(a7 ~ t7)

a_cum_7 <- cumsum(a7)

plot(a7 ~ t7)

plot(a_cum_7 ~ t7)

# The Model

AC7 <- function(t,m72,p72,q72,d71,d72){
  
  I(t<t0)*f(t,m72,p72,q72)+
    I(t==t0)*(f(t,m72,p72,q72) + d71)+
    I(t0+4>t & t>t0)*(f(t,m72,p72,q72) + d71)+
    I(t==t0+4)*(f(t,m72,p72,q72) + d72 + d71)+
    I(t>t0+4)*(f(t,m72,p72,q72) + d71+ d72)
  
}

# Fit AC6 with nls

BB7 <- nls(a_cum_7 ~ AC7(t7,m72,p72,q72,d71,d72),  
           algorithm = "port", 
           start = list(m72=12191.7, p72=0.0842829,  q72=0., d71=3745.45,  d72=1694.52),
           lower = c(0,0,0,0,0))

summary(BB7)

# Plot the fit

plot(t7, a_cum_7)
lines(t7,predict(BB7),col="blue",lwd=2)

plot(t7, a_cum_7, ylim = c(0, 25000))
lines(t7,predict(BB7),col="blue",lwd=2, ylim = c(0, 25000))

########################### Gen 7 GMM ###################################

dm1 <- c(1,rep(0,23))
dm2 <- c(rep(0, 4), 1,rep(0,19))

random_instroment <- c(1,20,12, 22, 16, 30, 32,  42,
                       38,  36,  52,  55,  45,  62,  64,
                       72,  70,  80,  82,  85,  90,  96, 100, 102)

dat <- cbind(a_cum_7, t1, 1, t1-1, dm1, dm2, (t1-1)^2, (t1-1)^4-1, random_instroment^2)

length(t1)

## Moments

moments <- function(par,data) {
  y <- as.numeric(data[,1])
  x <- data.matrix(data[,3:9])
  m <- x * as.vector((y - AC7(data[,2],par[1], par[2], par[3], par[4], par[5])))
  return(cbind(m))
}

#

#library(gmm)
#library(sandwich)
coef(BB7)[5]

GMM_Estimation_G7 <- gmm(moments,
                         x = dat,
                         method = "L-BFGS-B",
                         t0=c(m72=11000, p72=0.05, q72=0.00000001, d71=10000, d72=2800),
                         lower = c(10000,0.01,0,2500,1252),
                         upper = c(17000,5,5,8000,3500)
)

# Methods: “Nelder-Mead”, “BFGS”, “CG”, “L-BFGS-B”, “SANN”, “Brent”

summary(GMM_Estimation_G6)

# Plot the estimated model

TB <- AC5(
  1:24,
  coef(GMM_Estimation_G7)[1], 
  coef(GMM_Estimation_G7)[2], 
  coef(GMM_Estimation_G7)[3], 
  coef(GMM_Estimation_G7)[4], 
  coef(GMM_Estimation_G7)[5]
)

plot(a_cum_7 ~ t7, ylab = "Adoptions", xlab = "Time", ylim = c(0, 25000), xlim = c(0, 25),type="b",lwd=2)
par(new=TRUE)
plot(TB[1:length(t7)] ~ t7, col="green", ylab = "Adoptions", xlab = "Time", ylim = c(0, 25000), xlim = c(0, 25),type="b",lwd=2)




####### Aggregating the derived trends for the generations and removing the dummy from the initial spike #####

tt <- c(1: 95)

### Aggregaate data with dummy

Aggr_Act_Adp <- c(
  
                  rep(0, 11),
                  
                  AC1(0:11,coef(BB1)[1],coef(BB1)[2],coef(BB1)[3],coef(BB1)[4],coef(BB1)[5]),
                  
                  AC1(12:23,coef(BB1)[1],coef(BB1)[2],coef(BB1)[3],coef(BB1)[4],coef(BB1)[5])+
                    AC2(0:11,coef(BB2)[1],coef(BB2)[2],coef(BB2)[3],coef(BB2)[4],coef(BB2)[5]),
                  
                  AC1(24,coef(BB1)[1],coef(BB1)[2],coef(BB1)[3],coef(BB1)[4],coef(BB1)[5])+
                  AC2(12:23,coef(BB2)[1],coef(BB2)[2],coef(BB2)[3],coef(BB2)[4],coef(BB2)[5])+
                    AC3(0:11,coef(BB3)[1],coef(BB3)[2],coef(BB3)[3],coef(BB3)[4],coef(BB3)[5]),
                  
                  AC1(24,coef(BB1)[1],coef(BB1)[2],coef(BB1)[3],coef(BB1)[4],coef(BB1)[5])+
                    AC2(24,coef(BB2)[1],coef(BB2)[2],coef(BB2)[3],coef(BB2)[4],coef(BB2)[5])+
                  AC3(12:23,coef(BB3)[1],coef(BB3)[2],coef(BB3)[3],coef(BB3)[4],coef(BB3)[5])+
                    AC4(0:11,coef(BB4)[1],coef(BB4)[2],coef(BB4)[3],coef(BB4)[4],coef(BB4)[5]),
                  
                  
                  AC1(24,coef(BB1)[1],coef(BB1)[2],coef(BB1)[3],coef(BB1)[4],coef(BB1)[5])+
                    AC2(24,coef(BB2)[1],coef(BB2)[2],coef(BB2)[3],coef(BB2)[4],coef(BB2)[5])+
                    AC3(24,coef(BB3)[1],coef(BB3)[2],coef(BB3)[3],coef(BB3)[4],coef(BB3)[5])+
                  AC4(12:23,coef(BB4)[1],coef(BB4)[2],coef(BB4)[3],coef(BB4)[4],coef(BB4)[5])+
                    AC5(0:11,coef(BB5)[1],coef(BB5)[2],coef(BB5)[3],coef(BB5)[4],coef(BB5)[5]),
                  
                  AC1(24,coef(BB1)[1],coef(BB1)[2],coef(BB1)[3],coef(BB1)[4],coef(BB1)[5])+
                    AC2(24,coef(BB2)[1],coef(BB2)[2],coef(BB2)[3],coef(BB2)[4],coef(BB2)[5])+
                    AC3(24,coef(BB3)[1],coef(BB3)[2],coef(BB3)[3],coef(BB3)[4],coef(BB3)[5])+
                    AC4(24,coef(BB4)[1],coef(BB4)[2],coef(BB4)[3],coef(BB4)[4],coef(BB4)[5])+
                  AC5(12:23,coef(BB5)[1],coef(BB5)[2],coef(BB5)[3],coef(BB5)[4],coef(BB5)[5])+
                    AC6(0:11,coef(BB6)[1],coef(BB6)[2],coef(BB6)[3],coef(BB6)[4],coef(BB6)[5]),
                  
                  AC1(24,coef(BB1)[1],coef(BB1)[2],coef(BB1)[3],coef(BB1)[4],coef(BB1)[5])+
                    AC2(24,coef(BB2)[1],coef(BB2)[2],coef(BB2)[3],coef(BB2)[4],coef(BB2)[5])+
                    AC3(24,coef(BB3)[1],coef(BB3)[2],coef(BB3)[3],coef(BB3)[4],coef(BB3)[5])+
                    AC4(24,coef(BB4)[1],coef(BB4)[2],coef(BB4)[3],coef(BB4)[4],coef(BB4)[5])+
                    AC5(24,coef(BB5)[1],coef(BB5)[2],coef(BB5)[3],coef(BB5)[4],coef(BB5)[5])+
                  AC6(12:23,coef(BB6)[1],coef(BB6)[2],coef(BB6)[3],coef(BB6)[4],coef(BB6)[5])+
                    AC7(0:11,coef(BB7)[1],coef(BB7)[2],coef(BB7)[3],coef(BB7)[4],coef(BB7)[5])
                  
)


plot(Aggr_Act_Adp ~ tt)

########### Dummies removed


for(i in 13:95) {Aggr_Act_Adp[i] <- Aggr_Act_Adp[i]-coef(BB1)[4]}
for(i in 25:95) {Aggr_Act_Adp[i] <- Aggr_Act_Adp[i]-coef(BB2)[4]}
for(i in 37:95) {Aggr_Act_Adp[i] <- Aggr_Act_Adp[i]-coef(BB3)[4]}
for(i in 49:95) {Aggr_Act_Adp[i] <- Aggr_Act_Adp[i]-coef(BB4)[4]}
for(i in 61:95) {Aggr_Act_Adp[i] <- Aggr_Act_Adp[i]-coef(BB5)[4]}
for(i in 73:95) {Aggr_Act_Adp[i] <- Aggr_Act_Adp[i]-coef(BB6)[4]}
for(i in 85:95) {Aggr_Act_Adp[i] <- Aggr_Act_Adp[i]-coef(BB7)[4]}

plot(Aggr_Act_Adp ~ tt)

############################## Product line level model ####################

dd1 <- coef(BB1)[4]
dd2 <- coef(BB2)[4]
dd3 <- coef(BB3)[4]
dd4 <- coef(BB4)[4]
dd5 <- coef(BB5)[4]
dd6 <- coef(BB6)[4]
dd7 <- coef(BB7)[4]

DUMMIES <- c(dd1, dd2, dd3, dd4, dd5, dd6, dd7)

####### GDMR Cumulative ##########


X <- c(-0.0514718425553177, 0.0514718425553177, -0.1538699136085835, 
       0.1538699136085835, -0.2546369261678899, 0.2546369261678899,
       -0.3527047255308781, 0.3527047255308781, -0.4470337695380892,
       0.4470337695380892, -0.5366241481420199, 0.5366241481420199,
       -0.6205261829892429, 0.6205261829892429, -0.6978504947933158,
       0.6978504947933158, -0.7677774321048262, 0.7677774321048262,
       -0.8295657623827684, 0.8295657623827684, -0.8825605357920527,
       0.8825605357920527, -0.9262000474292743, 0.9262000474292743,
       -0.9600218649683075, 0.9600218649683075, -0.9836681232797472,
       0.9836681232797472, -0.9968934840746495, 0.9968934840746495)

w <- c(0.1028526528935588, 0.1028526528935588, 0.1017623897484055,
       0.1017623897484055, 0.0995934205867953, 0.0995934205867953, 
       0.0963687371746443, 0.0963687371746443, 0.0921225222377861, 
       0.0921225222377861, 0.0868997872010830, 0.0868997872010830,
       0.0807558952294202, 0.0807558952294202, 0.0737559747377052, 
       0.0737559747377052, 0.0659742298821805, 0.0659742298821805, 
       0.0574931562176191, 0.0574931562176191, 0.0484026728305941, 
       0.0484026728305941, 0.0387991925696271, 0.0387991925696271,
       0.0287847078833234, 0.0287847078833234, 0.0184664683110910, 
       0.0184664683110910, 0.0079681924961666, 0.0079681924961666)

##Spouge's approximation for the gamma function

h <- 5

c0 <- 1

C <- function(k) {(1/sqrt(2 * pi)) * (((-1)^(k - 1))/factorial(k - 1)) * ((-k + h)^(k - 1/2)) * exp(-k + h)}

G <- function(x) {
  ((x - 1 + h)^(x - 1/2)) * exp(-(x - 1 + h)) * sqrt(2 * pi) * (c0 + C(1)/(x - 1 + 1) + C(2)/(x - 1 + 2) + C(3)/(x - 1 + 3) + C(4)/(x - 1 + 4))
}

##Bass cumulative adoption process

f <- function(t, m, p, q) {(m - m * exp(-(p + q) * t))/(1 + (q/p) * exp(-(q + p) * t))}

##Derivatives of the Bass process

d_f <- function(t, m, p, q){(exp((p+q)*t)*m*p*(p+q)^2)/(exp((p+q)*t)*p+q)^2}

d2_f <- function(t, m, p, q){-((exp((p + q)*t) *
                                  m * p * (exp((p + q) * t) * p - q)
                                * (p + q)^3)/(exp((p + q) * t) * p + q)^3)}

d3_f <- function(t, m, p, q){(exp((p + q) * t) *
                                m * p * (p + q)^4 
                              * (exp(2 * (p + q) * t) 
                                 * p^2 - 4 * exp((p + q) * t) * p * q + q^2))/(exp((p + q) * t) * p + q)^4}

## Cumulative GDMR

g <- function(s, a, t, m, p, q) {(t - s)^(3 - a) * d3_f(s, m, p, q)}

S <- function(s, a, t, m, p, q) {
  sm = 0
  for(i in 1:30) {
    sm = sm + w[i] * g(((t/2) * X[i] + (t/2)), a, t, m, p, q)
  }
  sm
  #sum(w*g(((t/2) * X + (t/2)), a, t, m, p, q))
}




CGDMR <- function(t, a, p, q, m) {
  (d_f(0, m, p, q) / G(3 - a)) * (t^(2 - a)) +
    (d2_f(0, m, p, q)/G(4 - a)) * (t^(3 - a)) +
    (1  /G(4 - a)) * (t/2) * S(s, a, t, m, p, q)
}


#############################################


PFunc <- function (x,a,p,q,m,q1,q2,q3,q4,q5,q6,q7,q8,d,m8) {
  
    I(0 < x & x < 12) * (CGDMR(x, a, p, q, m) - f(x, d*dd1, p, q1)) +
    
    I(12 <= x & x < 24) * (CGDMR(x, a, p, q, m) - f(x - 12, d*dd2, p, q2) - f(12, d*dd1, p, q1)) +
    
    I(24 <= x & x < 36) * (CGDMR(x, a, p, q, m) - f(x - 24, d*dd3, p, q3) - f(12, d * dd2, p, q2) - f(12, d*dd1, p, q1)) +
    
    I(36 <= x & x < 48) * (CGDMR(x, a, p, q, m) - f(x - 36, d*dd4, p, q4) - f(12, d * dd3, p, q3) - f(12, d * dd2, p, q2) - f(12, d*dd1, p, q1)) +
    
    I(48 <= x & x < 60) * (CGDMR(x, a, p, q, m) - f(x - 48, d*dd5, p, q5) - f(12, d* dd4, p, q4) - f(12, d * dd3, p, q3) - f(12, d * dd2, p, q2) - f(12, d*dd1, p, q1)) +
                             
    I(60 <= x & x < 72) * (CGDMR(x, a, p, q, m) - f(x - 60, d*dd6, p, q6) - f(12, d*dd5, p, q5) - f(12, d* dd4, p, q4) - f(12, d * dd3, p, q3) - f(12, d * dd2, p, q2) - f(12, d*dd1, p, q1)) +
                             
    I(72 <= x & x< 84) * (CGDMR(x, a, p, q, m) - f(x - 72, d*dd7, p, q7) - f(12, d*dd6, p, q6) - f(12, d*dd5, p, q5) - f(12, d* dd4, p, q4) - f(12, d * dd3, p, q3) - f(12, d * dd2, p, q2) - f(12, d*dd1, p, q1)) +
                             
    I(84 <= x ) * (CGDMR(x, a, p, q, m) - f(x - 84, m8, p, q8) - f(12, d*dd7, p, q7) - f(12, d*dd6, p, q6) - f(12, d*dd5, p, q5) - f(12, d* dd4, p, q4) - f(12, d * dd3, p, q3) - f(12, d * dd2, p, q2) - f(12, d*dd1, p, q1))
}

PFUNCITON<-rep(0,96)
for(i in 1:96){
  PFUNCITON[i] <- PFunc (i, 0.79, 0.0016, 0.05, 72000, 0.648, 0.648, 0.648, 0.648, 0.648, 0.648, 0.648, 0.000000001, 1.16, 10000)
}


TTIIMME <- c(1:96)

############## Plot data along 

plot(PFUNCITON[1: 95], ylim = c(0, 100000), xlim = c(0, 96), xlab = "time", ylab = "sales", col="blue")

par(new=TRUE)

plot(Aggr_Act_Adp ~ tt, ylim = c(0, 100000), xlim = c(0, 96), xlab = "time", ylab = "sales",col="red")




############################## Model Fitting ###############################

# NLS model fitting using LM algoritm
library(minpack.lm)
nls_LM_out = nlsLM(Aggr_Act_Adp ~ PFunc(tt,a,p,q,m,q1,q2,q3,q4,q5,q6,q7,q8,d,m8),  
                     start = list(a= 0.79, p= 0.0016, q= 0.05, m= 72000, q1= 0.684, q2= 0.684, q3= 0.684, q4= 0.648, q5= 0.648, q6= 0.648, q7= 0.648, q8= 0.000001, d= 1.16, m8= 10000),
                     lower = c(rep(0, 14)),
                     upper = c(a= 1, p= 0.1, q= 0.1, m=200000, q1= 2, q2= 2, q3= 2, q4= 2, q5= 2, q6= 2, q7= 2, q8= 2, d= 3, m8= 20000),
                     control = nls.control(maxiter = 500, tol=100))

summary(nls_LM_out)
summary(nls_LM_out1)



######### Estimation With Constraints #########

library(nloptr)
# define sum of square error term
obj_fun <-function(param){

  return(sum((Aggr_Act_Adp-PFunc(tt,param[1],param[2],param[3],param[4],param[5],param[6],param[7],param[8],param[9],param[10],param[11],param[12],param[13],param[14]))^2))

}

#define equality constraint
eval_g_eq <- function (param){
  a=param[1]
  p=param[2]
  q=param[3]
  m=param[4]
  q1=param[5]
  q2=param[6]
  q3=param[7]
  q4=param[8]
  q5=param[9]
  q6=param[10]
  q7=param[11]
  q8=param[12]
  d=param[13]
  m8=param[14]
  
  return (c(f(12, d*dd1, p, q1) -dd1 ,
            f(12, d*dd2, p, q2) -dd2 ,
            f(12, d*dd3, p, q3) -dd3 ,
            f(12, d*dd4, p, q4) -dd4 ,
            f(12, d*dd5, p, q5) -dd5 ,
            f(12, d*dd6, p, q6) -dd6 ,
            f(12, d*dd7, p, q7) -dd7 ))
  
}

#define inequality constraint
hin<- function (param){
  a=param[1]
  p=param[2]
  q=param[3]
  m=param[4]
  q1=param[5]
  q2=param[6]
  q3=param[7]
  q4=param[8]
  q5=param[9]
  q6=param[10]
  q7=param[11]
  q8=param[12]
  d=param[13]
  m8=param[14]
  
  return (c(f(12, d*dd1, p, q1) -dd1 ,
            f(12, d*dd2, p, q2) -dd2 ,
            f(12, d*dd3, p, q3) -dd3 ,
            f(12, d*dd4, p, q4) -dd4 ,
            f(12, d*dd5, p, q5) -dd5 ,
            f(12, d*dd6, p, q6) -dd6 ,
            f(12, d*dd7, p, q7) -dd7 ,
          -f(12, d*dd1, p, q1) +dd1 ,
          -f(12, d*dd2, p, q2) +dd2 ,
          -f(12, d*dd3, p, q3) +dd3 ,
          -f(12, d*dd4, p, q4) +dd4 ,
          -f(12, d*dd5, p, q5) +dd5 ,
          -f(12, d*dd6, p, q6) +dd6 ,
          -f(12, d*dd7, p, q7) +dd7,
          dd7-f(12, m8, p, q8),
          dd7 * d -m8)
          )
  
}

opts <- list( "algorithm" = "NLOPT_LN_COBYLA",#"NLOPT_LN_COBYLA", #"NLOPT_GN_ISRES",
              "xtol_rel" = 1.0e-4,
              "maxeval" = 10000
              )

BB_Product_line_Const = nloptr( 
                      x0=c(0.77,0.0016,0.05,69000,0.66,0.66,0.66,0.66,0.66,0.66,0.66,0.5,1.14,18000),#as.numeric(coef(nls_LM_out)),
                      eval_f=obj_fun,
                      lb=c(0,0,0,69000,0,0,0,0,0,0,0,0.5,1,0),
                      ub=c(1,0.008,0.1,Inf, 2,2,2,2,2,2,2,2,5,30000),
                      #eval_g_eq=eval_g_eq,
                      eval_g_ineq=hin,
                      opts = opts)

BB_Product_line_Const$status

#given_param <- c(0.79,0.0016,0.05,72000,0.684,0.684,0.684,0.648,0.648,0.648,0.648,0.000001,1.16,10000)
given_param <- c(0.769997,0.00166661,0.0480419,69000,0.658768,0.658768,0.658768,0.658768,0.658768,0.658768,0.658768,0.5,1.14331,18087.4)

obj_fun(given_param)
eval_g_eq(given_param)
obj_fun(BB_Product_line_Const$solution)
eval_g_eq(BB_Product_line_Const$solution)


PFUNCITON<-rep(0,96)
for(i in 1:96){
  #PFUNCITON[i] <- PFunc (i, 0.769997,0.00166661,0.0480419,69000,0.658768,0.658768,0.658768,0.658768,0.658768,0.658768,0.658768,0.5,1.14331,18087.4)
  PFUNCITON[i] <- PFunc (i,8.290017e-01, 1.580634e-03, 4.639335e-02, 8.757562e+04, 6.611465e-01, 6.611465e-01, 6.611465e-01,6.611465e-01, 6.611465e-01, 6.611465e-01 ,6.611465e-01, 5.331092e-01, 1.147513e+00, 1.811274e+04)
}
plot(PFUNCITON[1: 95], ylim = c(0, 100000), xlim = c(0, 96), xlab = "time", ylab = "sales", col="blue")
par(new=TRUE)
plot(Aggr_Act_Adp ~ tt, ylim = c(0, 100000), xlim = c(0, 96), xlab = "time", ylab = "sales",col="red")

