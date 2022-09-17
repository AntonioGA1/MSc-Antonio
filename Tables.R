### Comparison between Heston European and Exact to check convergence using CV method ###
##Tabla 4.1 of the MSc
FDCV <- function(Smax,Vmax,gridlist,Nvec,c,d,T,r,xi,w,kappa,theta,rho,Svec,Vvec,opttype,method,n,a){
  Tabla <- data.frame(S0 = integer(),
                      V0 = double(),
                      N = integer(),
                      m1 = integer(),
                      m2 = integer(),
                      AmerFD = double(),
                      EuroFD = double(),
                      ExactEuro = double()
                      )
  cont = 1
  for(S0 in Svec){
    for(V0 in Vvec){
      for(N in Nvec){
        for(grid in gridlist){
          m1 = grid[1]
          m2 = grid[2]
          AmerFD <- FDHestonAmerican(Smax,Vmax,m1,m2,N,c,d,T,r,xi,w,kappa,theta,rho,S0,V0,opttype,method,TRUE)$y
          EuroFD <- FDHestonAmerican(Smax,Vmax,m1,m2,N,c,d,T,r,xi,w,kappa,theta,rho,S0,V0,opttype,method,FALSE)$y
          if(opttype == 'C'){
            EExact <- HestonGaussLaguerre(kappa,theta,lambda,rho,xi,S0,K,T,r,q,n,a)$Call
          }else{
            EExact <- HestonGaussLaguerre(kappa,theta,lambda,rho,xi,S0,K,T,r,q,n,a)$Put
          }
          Tabla[cont,] <- c(S0,V0,N,m1,m2,AmerFD,EuroFD,EExact)
          cont = cont+1
        }
      }
    }
  }
  Tabla$ErrorEur = abs(Tabla$EuroFD - Tabla$ExactEuro)
  Tabla$RelError = Tabla$ErrorEur/Tabla$ExactEuro
  return(Tabla)
}
##Table 4.2 of the MSc
MCCV <- function(NTvec,NSVec,T,r,xi,kappa,theta,rho,Svec,Vvec,opttype,simtype,BStype,seeds=NULL,n,a){
  Tabla <- data.frame(S0 = integer(),
                      V0 = double(),
                      NS = integer(),
                      NT = integer(),
                      AmerLSM = double(),
                      EuroLSM = double(),
                      ExactEuro = double()
                      )
  cont = 1
  for(S0 in Svec){
    for(V0 in Vvec){
      for(NT in NTvec){
        for(NS in NSVec){
          LSM = MCHestonAmerican(NS,NT,S0,V0,K,T,r,xi,kappa,theta,rho,opttype,simtype,BStype,seeds)
          if(opttype == 'C'){
            EExact <- HestonGaussLaguerre(kappa,theta,lambda,rho,xi,S0,K,T,r,q,n,a)$Call
          }else{
            EExact <- HestonGaussLaguerre(kappa,theta,lambda,rho,xi,S0,K,T,r,q,n,a)$Put
          }
          Tabla[cont,] <- c(S0,V0,NT,NS,LSM$AmerPrice,LSM$EuroPrice,EExact)
          cont = cont+1
        }
      }
    }
  }
  Tabla$ErrorEur = abs(Tabla$EuroLSM - Tabla$ExactEuro)
  Tabla$RelError = Tabla$ErrorEur/Tabla$ExactEuro
  return(Tabla)
}

### Comparison between Heston for American option with no stochastic variance and GBM ###

##Table 4.3
FDGBM <- function(Smax,Vmax,m1,m2,N,c,d,T,rvec,w,Svec,Vvec,opttype,method,NB){
  Tabla <- data.frame(S0 = integer(),
                      r = double(),
                      v0 = double(),
                      EuroFD = double(),
                      EuroGBM = double(),
                      AmerFD = double(),
                      AmerGBM = double()
  )

  cont = 1
  for(s0 in Svec){
    for (r in rvec){
      for(v0 in Vvec){
        AmerFD <- FDHestonAmerican(Smax,Vmax,m1,m2,N,c,d,T,r,0,w,0,0,1,s0,v0,opttype,method,TRUE)$y
        EuroFD <- FDHestonAmerican(Smax,Vmax,m1,m2,N,c,d,T,r,0,w,0,0,1,s0,v0,opttype,method,FALSE)$y
        EuroGBM <- EuroPriceGBM(NB,s0,K,r,q,sqrt(v0),T,opttype)
        AmerGBM <- AmerPriceGBM(NB,r,q,s0,sqrt(v0),K,T,opttype)

        Tabla[cont,] <- c(s0,r,v0,EuroFD,EuroGBM,AmerFD,AmerGBM)
        cont = cont+1
      }
    }
  }

  Tabla$ErrorFDEuro = abs(Tabla$EuroFD - Tabla$EuroGBM)
  Tabla$RErrorFDEuro = Tabla$ErrorFDEuro/Tabla$EuroGBM
  Tabla$ErrorFDAmer = abs(Tabla$AmerFD - Tabla$AmerGBM)
  Tabla$RErrorFDAmer = Tabla$ErrorFDAmer/Tabla$AmerGBM
  return(Tabla)
}

## Table 4.4
MCGBM <- function(NT,NS,T,rvec,Svec,Vvec,opttype,simtype,BStype,seeds=NULL,NB){
  Tabla <- data.frame(S0 = integer(),
                      r = double(),
                      v0 = double(),
                      EuroMC = double(),
                      EuroGBM = double(),
                      AmerMC = double(),
                      AmerGBM = double()
  )

  cont = 1
  for(s0 in Svec){
    for (r in rvec){
      for(v0 in Vvec){
        LSM = MCHestonAmerican(NS,NT,s0,v0,K,T,r,0,0,0,1,opttype,simtype,BStype,seeds)
        EuroGBM <- EuroPriceGBM(NB,s0,K,r,q,sqrt(v0),T,opttype)
        AmerGBM <- AmerPriceGBM(NB,r,q,s0,sqrt(v0),K,T,opttype)

        Tabla[cont,] <- c(s0,r,v0,LSM$EuroPrice,EuroGBM,LSM$AmerPrice,AmerGBM)
        cont = cont+1
      }
    }
  }

  Tabla$ErrorMCEuro = abs(Tabla$EuroMC - Tabla$EuroGBM)
  Tabla$RErrorMCEuro = Tabla$ErrorMCEuro/Tabla$EuroGBM
  Tabla$ErrorMCAmer = abs(Tabla$AmerMC - Tabla$AmerGBM)
  Tabla$RErrorMCAmer = Tabla$ErrorMCAmer/Tabla$AmerGBM
  return(Tabla)
}

#Import methods
source("Parameters.R")
source("MCHestonAmerican.R") #LSM method
source("FDHestonAmerican.R") #FD method
source("ExactGBM.R") #GBM method
source("HestonExactGL.R") #Exact European option price

#Comparison between algorithms
Comp <- function(Smax,Vmax,m1,m2,N,c,d,T,r,xi,w,kappa,theta,rho,Svec,Vvec,opttype,method,NT,NS,simtype,BStype,seeds=NULL){
    Tabla <- data.frame(S0 = integer(),
                      V0 = double(),
                      AmerFD = double(),
                      TimeFD = double(),
                      AmerMC = double(),
                      TimeMC = double()
  )
  cont = 1
  for(S0 in Svec){
    for(V0 in Vvec){
      t1 = proc.time()
      AmerFD <- FDHestonAmerican(Smax,Vmax,m1,m2,N,c,d,T,r,xi,w,kappa,theta,rho,S0,V0,opttype,method,TRUE)$y
      t2 = proc.time()
      FDtime = (t2-t1)[['elapsed']]
      t1 = proc.time()
      AmerMC <- MCHestonAmerican(NS,NT,S0,V0,K,T,r,xi,kappa,theta,rho,opttype,simtype,BStype,seeds)$AmerPrice
      t2 = proc.time()
      MCtime = (t2-t1)[['elapsed']]

      Tabla[cont,] = c(S0,V0,AmerFD,FDtime,AmerMC,MCtime)
      cont = cont+1
    }
  }
  Tabla$AbsDiff = abs(Tabla$AmerFD-Tabla$AmerMC)
  Tabla$RelDiff = Tabla$AbsDiff/max(Tabla$AmerFD,Tabla$AmerMC)
  return(Tabla)
}

CompTable(Smax,Vmax,50,25,20,c,d,T,r,xi,w,kappa,theta,rho,c(90,100,110),c(0.04,0.25),opttype,method,50,5000,simtype,BStype,seeds=NULL)

