#Simulate the Stock price and volatility prices
#Inputs:
## Heston parameters
## NS: Number of paths
## NT: Number of time steps
## simtype: Type of simulation 'E' for Euler, 'M' for Milstein
## BStype: Negative variance corrector
#Output: European and American option price and executed paths at time t

simulation <- function(NS,NT,T,r,xi,kappa,theta,rho,S0,V0,simtype,BStype,seeds){
  dt = T/NT
  set.seed(0)
  #Initialize stock-price and volatility
  V = matrix(0,NS,NT)
  S = matrix(0,NS,NT)
  #Add the initial value
  V[,1] = V0
  S[,1] = S0
  #Select boundary scheme
  if(BStype == 'A'){
    BS = function(x) pmax(0,x)
  }else if(BStype == 'R'){
    BS = function(x) abs(x)
  } else{
    BS = function(x) x
  }

  if(is.null(seeds)){
    Z1 = rnorm(NS*NT,mean=0,sd=1)
    Z2 = rnorm(NS*NT,mean=0,sd=1)
  }else{
    set.seed(seeds[1])
    Z1 = rnorm(NS*NT,mean=0,sd=1)
    set.seed(seeds[2])
    Z2 = rnorm(NS*NT,mean=0,sd=1)
  }

  #Generate correlated random variables
  Zv = matrix( Z1, NS, NT)
  Zs = rho * Zv + sqrt(1-rho^2)* matrix( Z2, NS, NT)
  #Simulate the variance and stock-price
  for (t in 2:NT){
    if (simtype == 'E'){
        V[,t] = BS(V[,t-1] + kappa*dt*(theta-V[,t-1]) + xi*sqrt(V[,t-1]*dt)*Zv[,t-1])
    }
    else if(simtype == 'M'){
        V[,t] = BS(V[,t-1] + kappa*dt*(theta-V[,t-1]) + xi*sqrt(V[,t-1]*dt)*Zv[,t-1] + (1/4)*xi^2*
          dt*(Zv[,t-1]^2-1))
    }else{
      num = 0.5*xi^2*V[,t-1]*(1-exp(-2*kappa*dt))/kappa
      den = (exp(-kappa*dt)*V[,t-1] + (1-exp(-kappa*dt))*theta)^2
      gam = log(1 + num/den)
      V[,t] = (exp(-kappa*dt)*V[,t-1] + (1-exp(-kappa*dt))*theta) * exp(-0.5*gam^2 + gam*Zv[,t])
    }
    #Discretization of the log stock-price
    S[,t] = S[,t-1]*exp( (r-V[,t-1]/2)*dt + sqrt(V[,t-1]*dt)*Zs[,t])
  }
  return(list("S"=S,"V"=V))
}

MCHestonAmerican <- function(NS,NT,S0,V0,K,T,r,xi,kappa,theta,rho,opttype,simtype,BStype,seeds = NULL){
  #Initialization
  dt = T/NT
  C <-E <-ExBool <- matrix(0,NS,NT)  #Cashflow, Exercise price, Exercise boolean
  #First we simulate the price with the function above
  S <- simulation(NS,NT,T,r,xi,kappa,theta,rho,S0,V0,simtype,BStype,seeds)$S

  #Define payoff function. s represents any subset of S
  if(opttype == 'P'){
    fpayoff_ = function(s,K) pmax(0,K-s)
  }else{
    fpayoff_ = function(s,K) pmax(0,s-K)
  }

  #We start by defining the value at maturity
  E[,NT] <- fpayoff_(S[,NT],K)
  C[,NT] <- fpayoff_(S[,NT],K)

  EuroPrice = exp(-r*T)*mean(C[,NT])
  ExBool[,NT] <-1

  #Before maturity t=2 to t=NS-1
  for (t in (NT-1):2){
    E[,t] = fpayoff_(S[,t],K)
    #We look for the ITM stock-paths among the others paths
    #We are storing the prices along the indices in the same df
    df = setNames(data.frame(S[,t],fpayoff_(S[,t],K),E[,t]>0),
                  c("StockPrice","Payoff","ITMBool"))
    I = which(df["ITMBool"] == TRUE)

    if (length(I) != 0){
      #Cash flows at t+1
      B = C[I,t+1] * exp(-r*dt)
      if (length(B)> 2){
        numOP = 2
      }else{
        numOP = length(B)-1
      }
      #Define laguerre polynomials as proposed in the original article
      f_laguerre <- lapply(laguerre.polynomials(numOP,normalized = TRUE),as.function)
      A_laguerre <-sapply(f_laguerre,function(x) sapply(df[I,"StockPrice"],x))
      #Make the fit of the points based on the laguerre base
      fit <- coef(lm(B~0+A_laguerre))
      names(fit) <- NULL
      #Now, we predict the value
      predCF = A_laguerre %*% matrix(fit,ncol = 1,nrow = ifelse(is.null(ncol(A_laguerre)),1,ncol(A_laguerre)))
      index_ = which(df[,"ITMBool"]==TRUE)
      temp_ = rep(NA,NS-length(index_))
      cont = 1
      for(i in index_){
        temp_ = append(temp_,predCF[cont],after = i-1)
        cont = cont+1
      }
      df$predCF = temp_
    }else{
      df$predCF = NA
    }
    #See what path we are exercising or not
    df$Exercise = ( (df[,"ITMBool"] == TRUE) & (df[,"Payoff"] > df[,"predCF"]))
    #Insert the Cashflows
    Ex = which(df[,"Exercise"] == TRUE)
    NEx = which(df[,"Exercise"] == FALSE)
    C[Ex,t] = df[Ex,"Payoff"] #Ejercemos
    ExBool[Ex,t] <- 1
    C[NEx,t] = exp(-r*dt)*C[NEx,t+1] #No Ejercemos y descontamos un periodo
    ExBool[NEx,t] <- 0
  }
  AmerPrice = exp(-r*dt)*mean(C[,2]) #Obtenemos el precio como media de todos los caminos
  return(list('ExBool' = ExBool,'AmerPrice' = AmerPrice,'EuroPrice'= EuroPrice))
}




