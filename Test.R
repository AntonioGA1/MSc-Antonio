#Import methods
source("Parameters.R")
source("Algorithms/MCHestonAmerican.R") #LSM method
source("Algorithms/FDHestonAmerican.R") #FD method
source("Algorithms/ExactGBM.R") #GBM method
source("Algorithms/HestonExactGL.R") #Exact European option price

##### Test individual price #####
testLSM = MCHestonAmerican(NS,NT,S0,V0,K,T,r,xi,kappa,theta,rho,opttype,simtype,BStype)
AmericanPriceLSM = testLSM$AmerPrice
EuropeanPriceLSM = testLSM$EuroPrice
AmericanPriceFD  <- FDHestonAmerican(Smax,Vmax,m1,m2,N,c,d,T,r,xi,w,kappa,theta,rho,S0,V0,opttype,method,american)$y
EuropeanPriceFD  <- FDHestonAmerican(Smax,Vmax,m1,m2,N,c,d,T,r,xi,w,kappa,theta,rho,S0,V0,opttype,method,american)$y
ExactEuropeanPrice <- HestonGaussLaguerre(kappa,theta,lambda,rho,xi,S0,K,T,r,q,n,a)$Put
AmericanPriceGBM <- AmerPriceGBM(NB,r,q,S0,sigma,K,T,opttype)
EuropeanPriceGBM <- EuroPriceGBM(NB,S0,K,r,q,sigma,T,opttype)

source("Tables.R")
#### Comparison Tables #####
source()
Svec = c(90,100,110)
Vvec = c(0.04,0.16)
rvec = c(0.05,0.10)

#CV tables
gridlist = list(c(10,5),c(20,10),c(50,20))
Nvec = c(5,20)
FDCVTable <- FDCV(Smax,Vmax,gridlist,Nvec,c,d,T,r,xi,w,kappa,theta,rho,Svec,Vvec,opttype,method,n,a)
NTvec = c(50,200)
NSVec = c(500,2500,20000)
MCCVTable <- MCCV(NTvec,NSVec,T,r,xi,kappa,theta,rho,Svec,Vvec,opttype,simtype,BStype,seeds=NULL,n,a)

#GBM tables
FDGBMTable <- FDGBM(Smax,Vmax,m1,m2,N,c,d,T,rvec,w,Svec,Vvec,opttype,method,NB)
MCGBMTable <- MCGBM(NT,NS,T,rvec,Svec,Vvec,opttype,simtype,BStype,seeds=NULL,NB)

#Comparison table
CompTable <- Comp(Smax,Vmax,m1,m2,N,c,d,T,r,xi,w,kappa,theta,rho,Svec,Vvec,opttype,method,NT,NS,simtype,BStype,seeds=NULL)
