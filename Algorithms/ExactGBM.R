#Exact price European option with GBM

#Inputs
## GBM parameters (with sigma = sqrt(V0) )
## NB: number of iteration to compute optimal exercise boundary

#Outputs:
## European price GBM
## American price GBM
European <- function(S,K,r,q,sigma,t,opttype){
  d1 = log(S/K) + (r-q + 0.5 * sigma^2) * (T-t) / (sigma * sqrt(T-t))
  d2 = log(S/K) + (r-q - 0.5 * sigma^2) * (T-t) / (sigma * sqrt(T-t))

  if(opttype == 'C'){
    EurPrice = S * exp(-q * (T-t)) * pnorm(d1) - K *exp(-r*(T-t)) * pnorm(d2)
  }else if(opttype == 'P'){
    EurPrice = K * exp(-r * (T-t)) * pnorm(-d2) - S * exp(-q * (T-t)) * pnorm(-d1)
  }
  return(EurPrice)
}

b_GBM_AmPut <- function (tol = 1e-2, mu = 0, sigma = 1,
                         T = 1, N = 50, t = NULL,t_ini = 0, S = 10,
                         Kernel = K_GBM_Amput) {

  # How many boundaries to compute
  n <- length(sigma)
  # Check if the discretization is given
  if (is.null(t)) {
    t <- log(seq(exp(t_ini), exp(0.5), l = N))
    N <- length(t) - 1
    T <- t[N + 1]                       # Expiration date
    Delta <- t[2:(N + 1)] - t[1:N]
    # Discretization
  } else {
    N <- length(t) - 1                  # Number of subintervals
    T <- t[N + 1]                       # Expiration date
    Delta <- t[2:(N + 1)] - t[1:N]      # Length step matrix
  }
  # Preallocating boundary
  bnd <- matrix(rep(S, (N + 1) * n),  ncol = N + 1)

  # Computing boundary's second last point
  # bnd[, N] <- sigma / 2 * sqrt(pi * (T - t[N]) / 2)  + y

  # Boundary computation
  for ( i in N:1 ) {
    bnd_old <- bnd[, i + 1]
    e <- 1                    # error in the while-loop
    t_u <- t[(i + 1):(N + 1)]       # time from the (i + 1)-th element
    bt_u <- bnd[, (i + 1):(N + 1)]  # boundary from the (i + 1)-th element onwards
    # A term needed that doesn't depend on b(t_i)
    # aux <- 2 * (T - t[i]) / ( (T - t[i]) + (t[N] - t[i]) )
    # Evaluate integral int_{t_{N - 1}}^{T} sqrt((u - t_i) / (T - u)) du
    # H <- ( T - t[i] ) * ( pi / 2 - atan(sqrt((t[N] - t[i]) / (T - t[N]))) ) +
    #  ( T - t[N] ) * sqrt((t[N] - t[i]) / (T - t[N]))
    # Fixed point algorithm
    cont  = 1
    while ( (e > tol) & !(cont > 1000 & e < 1)  ) {
      # Evaluate the kernel
      K <- Kernel(t[i], bnd_old, t_u, bt_u, y, T, mu, sigma)
      # Updating the boundary
      I <- integrate(function(z) pnorm((log((S - z) / bnd_old) - (mu - sigma^2 / 2) * (T - t[i]))
                                       / (sigma * sqrt(T - t[i])),
                                       mean = 0, sd = 1, lower.tail = T), 0, S, rel.tol = 1e-8)
      bnd[, i] <- S - exp(-mu * (T - t[i])) * I$value -
                  mu *  S * rowSums( t(Delta[i:N] * t(K)) )
      # Relative error
      e <- max( abs((bnd[, i] - bnd_old) / bnd[, i]) )
      # Caution: the while will keep running until all errors associated to
      # each of the sigma's are smaller than the tolerange
      # Update
      bnd_old <- bnd[, i]
      cont = cont+1
    }
  }
  return(bnd)
}
# Kernel of the Volterra integral equation using the Mean as optimality criterion
K_GBM_Amput <- function(t, x, t_u, bt_u, y, T, mu, sigma){

  # Number of sigmas and length of the time grid
  n <- length(sigma)
  N <- length(t_u)

  # Compute z
  mu <- matrix(rep(mu, N), ncol = N)
  sigma <- matrix(rep(sigma, N), ncol = N)
  x <- matrix(rep(x, N), ncol = N)
  t_u <-matrix(rep(t_u, each = n), nrow = n)

  # Evaluate Kernel
  z <- (log(bt_u/ x) - (mu - sigma^2/2) * (t_u - t)) / (sigma * sqrt(t_u - t))
  K <- exp(-mu * (t_u - t)) * pnorm(z, mean = 0, sd = 1, lower.tail = T)

  return(K)
}

AmerPriceGBM = function(N,r,q,S0,sigma,K,T,opttype){
  tvec <- log(seq(exp(0), exp(T), l = N))
  b <- b_GBM_AmPut(mu = r, sigma = sigma, N = N,t = tvec, S=S0, T=T)
  f = function(u,t) exp(-r*(u-t))* pnorm( (1/(sigma*sqrt(u-t))) *(log(interp1(tvec,as.vector(b),u)/S0) - (r-sigma^2/2)*(u-t)) )
  AmerGBM = European(S0,K,r,q,sigma,tvec[1],opttype) + r*K * integrate(f,lower = tvec[1] ,upper = T, t = tvec[1])$value
  return(AmerGBM)
}

EuroPriceGBM <- function(N,S,K,r,q,sigma,T,opttype){
  tvec <- log(seq(exp(0), exp(T), l = N))
  return(European(S,K,r,q,sigma,tvec[1],opttype))
}

