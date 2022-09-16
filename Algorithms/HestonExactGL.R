#Inputs:
## Heston parameters
## Gauss-Laguerre parameters: n number of nodes, a: exponent of x (a>=0)

#Outputs: Put and call price
HestonGaussLaguerre = function(kappa,theta,lambda,rho,xi,S,K,T,r,q,n,a){

  HestonProb <- function(phi,kappa,theta,lambda,rho,xi,tau,K,S,r,q,V0,pnum){
    x = log(S)

    a = kappa * theta
    if(pnum==1){
      u = 0.5
      b = kappa + lambda - rho*xi
    }else{
      u = -0.5
      b = kappa + lambda
    }

    d = sqrt((rho*xi*complex(imaginary = 1)*phi - b)^2 - xi^2*(2*u*complex(imaginary = 1)*phi - phi^2));
    g = (b - rho*xi*complex(imaginary = 1)*phi + d) / (b - rho*xi*complex(imaginary = 1)*phi - d);

    G = (1 - g*exp(d*tau))/(1-g);
    C = (r-q)*complex(imaginary = 1)*phi*tau + a/xi^2*((b - rho*xi*complex(imaginary = 1)*phi + d)*tau - 2*log(G));
    D = (b - rho*xi*complex(imaginary = 1)*phi + d)/xi^2*((1-exp(d*tau))/(1-g*exp(d*tau)));

    # The characteristic function.
    f = exp(C + D*V0 + complex(imaginary = 1)*phi*x);

    # Return the real part of the integrand.
    y = Re(exp(-complex(imaginary = 1)*phi*log(K))*f/complex(imaginary = 1)/phi)
    return(y)
  }

  #Generate Gauss-Laguerre weights and nodes
  GL = gaussLaguerre(n,a)
  W = c()
  X = c()
  for(i in 1:n){
    X = c(X,GL$x[i])
    W = c(W,GL$w[i]*GL$x[i]^(a)*exp(GL$x[i]))
  }

int1 = c()
int2 = c()
for(k in 1:n){
  int1 = c(int1,W[k]*HestonProb(X[k],kappa,theta,lambda,rho,xi,T,K,S0,r,q,V0,1))
  int2 = c(int2,W[k]*HestonProb(X[k],kappa,theta,lambda,rho,xi,T,K,S0,r,q,V0,2))
}

#Define P1 and P2
P1 = 1/2 + (1/pi)*sum(int1);
P2 = 1/2 + (1/pi)*sum(int2);

# The call price
HestonC = S*exp(-q*T)*P1 - K*exp(-r*T)*P2;

# The put price by put-call parity
HestonP = HestonC - S*exp(-q*T) + K*exp(-r*T);
return(list('Call'= HestonC,'Put'=HestonP))
}
