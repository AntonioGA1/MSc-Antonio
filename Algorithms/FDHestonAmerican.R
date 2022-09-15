#Inputs:
## Heston model parameters
## m1, m2, N for the stock price, volatility and time grid
## c and d constant for non-uniform grid
## w: from 0 (explicit) to 1 (implicit)
## Method: Choose between 'DO' (Douglas), 'CS' (Craig-Sneyd), 'MCS' (Modiﬁed Craig-Sneyd) or 'HV' (Hundsdorfer-Verwer) schemes
## American: TRUE for pricing an American option or FALSE to price an European option

#Outputs:
##y: American/European option price.
##U: Price and volatility apprx for every point at time N
##Svec and Vvec: Stock price and Volatility grid.

gen_grid <- function(Smax,S0,Vmax,V0,K,m1,m2,c,d){
  Dxi = (1.0 / (m1-1)) * (asinh((Smax - K) / c) - asinh(-K / c))
  s = asinh(-K/c) + 0:(m1-1) * Dxi #Equidistant points
  s_grid =  K+c*sinh(s)
  s_grid = sort(c(s_grid,S0))
  s_grid = s_grid[-length(s_grid)]

  #In this case we use beta instead of eta, that is alredy used for the model
  Deta = (1.0 / (m2-1)) * asinh(Vmax / d)
  v = 0:(m2-1) * Deta
  v_grid = d*sinh(v)
  v_grid = sort(c(v_grid,V0))
  v_grid = v_grid[-length(v_grid)]

  return(list("s_grid" = s_grid,
              "v_grid" = v_grid))
}
delta_coeff <- function(i,vector){
  D1 = vector[i+1] - vector[i]
  D2 = vector[i+2] - vector[i+1]
  L = 2 / (D1 * (D1 + D2))
  C = -2 / (D1 * D2)
  R = 2 / (D2 * (D1 + D2))
  return(c(L,C,R))
}
alpha_coeff <- function(i,vector){
  D0 = vector[i] - vector[i-1]
  D1 = vector[i+1] - vector[i]
  L = D1 / (D0 * (D0 + D1))
  C = (-D0 - D1) / (D0 * D1)
  R = (D0 + 2 * D1) / (D1 * (D0 + D1))
  return(c(L,C,R))
}
beta_coeff <- function(i,vector){
  D1 = vector[i+1] - vector[i]
  D2 = vector[i+2] - vector[i+1]
  L = -D2 / (D1 * (D1 + D2))
  C = (D2 - D1) / (D1 * D2)
  R = D1 / (D2 * (D1 + D2))
  return(c(L,C,R))
}
gamma_coeff <- function(i,vector){
  D2 = vector[i+2] - vector[i+1]
  D3 = vector[i+3] - vector[i+2]

  L = (-2 * D2 - D3) / (D2 * (D2 + D3))
  C = (D2 + D3) / (D2 * D3)
  R = -D2 / (D3 * (D2 + D3))

  return(c(L,C,R))
}
GenMatrix <- function(m1,m2,m,r,rho,xi,theta,kappa,SVector,VVector){
  # Use the definition of derivates of the original article
  A0 = matrix(0,m,m)
  A1 = matrix(0,m,m)
  A2 = matrix(0,m,m)

  `%+=%` = function(e1,e2) eval.parent(substitute(e1 <- e1 + e2))
  #A0 Matrix
  for (j in 1:(m2-2)){
    for (i in 1:(m1-2)){
      c = rho * xi * SVector[i+1]*VVector[j+1]
      for (k in c(-1,0,1)){
        A0[(i+1) + j * m1,seq((i+k+1) + (j-1) * m1,(i+k+1) + (j+1) * m1,by=m1)] %+=% (c * (matrix(beta_coeff(i, SVector),3,1) %*%
          matrix(beta_coeff(j, VVector),1,3))[k+2,])
      }
    }
  }
  #A1 Matrix Definition
  for (j in 0:(m2-1)){
    for(i in 1:(m1-2)){
      a = 0.5*SVector[i+1]^2 * VVector[j+1]
      b = r*SVector[i+1]
      A1[(i+1) + j * m1, (i + j * m1):((i+2) + j * m1)] %+=% (a * delta_coeff(i, SVector) + b * beta_coeff(i, SVector))
      A1[(i+1) + j * m1, (i+1) + j * m1] %+=% (- 0.5 * r)
    }
    A1[m1 + j * m1, m1 + j * m1] %+=% (- 0.5 * r)
  }

  for(j in 0:(m2-3)){
    for (i in 0:(m1-1)){
      temp = kappa*(theta - VVector[j+1])
      temp2 = 0.5*xi^2 * VVector[j+1]
      if (VVector[j+1] > 1){
        A2[(i+1) + (j + 1) * m1, seq((i+1) + m1 * (j-1),(i+1) + m1 * (j + 1),by=m1)] %+=% (temp * alpha_coeff(j+1, VVector))
        A2[(i+1) + (j + 1) * m1, seq((i+1) + m1 * j,(i+1) + m1 * (j + 2),by=m1)] %+=% (temp2 * delta_coeff(j, VVector))
      }
      if (j == 0){
        A2[i+1, seq((i+1),i + 1 + 2*m1,by=m1)] %+=% (temp * gamma_coeff(j+1, VVector))
      } else{
        A2[(i+1) + j * m1, seq( (i+1) + m1 * (j - 1),(i+1) + m1 * (j + 1),by = m1)] %+=% (temp * beta_coeff(j , VVector) + temp2 * delta_coeff(j , VVector))
      }
      A2[(i+1) + j * m1, (i+1) + j * m1] %+=% (- 0.5 * r)
    }
  }
  A = A0+A1+A2
  return(list('A0'=A0,'A1'=A1,'A2'=A2,'A'=A))
}
GenBoundary <- function(m1,m2,m,r,SVector){
  b0 = numeric(m)
  b1 = numeric(m)
  b2 = numeric(m)

  for(j in 1:m2){
    b1[m1*j] = r*SVector[length(SVector)]
  }
  for(i in 2:m1){
    b2[m-m1-1+i] = -0.5*r*SVector[i]
  }

  b = b0+b1+b2

  return(list('b0'=b0,'b1' = b1,'b2' = b2, 'b' = b))
}

FDHestonAmerican <- function(Smax,Vmax,m1,m2,N,c,d,T,r,xi,w,kappa,theta,rho,S0,V0,opttype,method,american){
  #initializa time increment and array
  dt = T/N
  t = 0:T *dt
  m = m1 * m2

  #Make grid
  grid <- gen_grid(Smax,S0,Vmax,V0,K,m1,m2,c,d)
  S <- grid$s_grid
  V <- grid$v_grid

  A_list = GenMatrix(m1,m2,m,r,rho,xi,theta,kappa,S,V)
  A <- A_list$A
  A0 <- A_list$A0
  A1 <- A_list$A1
  A2 <- A_list$A2

  b_list = GenBoundary(m1,m2,m,r,S)
  b <- b_list$b
  b0 <- b_list$b0
  b1 <- b_list$b1
  b2 <- b_list$b2

  listU = list()
  SVector = rep(t(S),m2)
  #Now we want to loop over time
  u = matrix(0,nrow = m,1)
  I = eye(m)
  if (opttype == 'P'){
    U = pmax(0,K-SVector)
  }else{
    U = pmax(0,SVector-K) #Value of U(T)
  }
  listU[[1]] <- reshape(array(U),m1,m2)
  print("Aplicando Metodo")
  for (t in 2:(N+1)){
    u = U
    Y0 = (I + dt*A)%*%u + dt*b
    Y1 = mldivide( (I - w * dt*A1), (Y0 - w*dt*A1%*%u))
    Y2 = mldivide( (I - w * dt*A2), (Y1 - w*dt*A2%*%u))
    if (method == 'DO'){
      U = Y2
    }else if(method == 'CS'){
      Z0 = Y0 + (1/2)*dt*(A0%*%Y2 - A0%*%u)
      Z1 = mldivide( (I-w*dt*A1),(Z0 - w*dt*A1%*%u) )
      Z2 = mldivide( (I-w*dt*A2),(Z1 - w*dt*A2%*%u) )
      U = Z2
    }else if(method == 'MCS'){
      Y0_ = Y0 + w * dt * (A0%*%Y2 - A0%*%u)
      Z0 = Y0_ + (1/2-w)*dt*((A0+A1+A2)%*%Y2 - (A0+A1+A2)%*%u)
      Z1 = mldivide((I - w*dt*A1),(Z0 - w*dt*A1%*%u))
      Z2 = mldivide((I - w*dt*A2),(Z1 - w*dt*A2%*%u))
      U = Z2
    }else if(method == 'HV'){
      Z0 = Y0 + (1/2)*dt*((A0+A1+A2)%*%Y2 - (A0+A1+A2)%*%u)
      Z1 = mldivide((I - w*dt*A1),(Z0 - w*dt*A1%*%Y2))
      Z2 = mldivide((I - w*dt*A2),(Z1 - w*dt*A2%*%Y2))
      U = Z2
    }
    #Include the American early exercise
    if (american == TRUE)
    {temp = reshape(U,m1,m2)
      for (i in 1:m1){
        if (opttype == 'C'){
          temp[i,] = pmax(temp[i,],S[i] - K)
        }else if(opttype == 'P'){
          temp[i,] = pmax(temp[i,],K - S[i])
        }
      }
    U = reshape(temp,m,1)
    }
    listU[[t]] <- reshape(U,m1,m2)
    print(paste0(t-1,"/",N))
  }
  U = reshape(U,m1,m2)
  y = interp2(V,S,U,V0,S0)
  return(list("y"=y,"U"=U,"Svec" = S,"Vvec" = V))
}
