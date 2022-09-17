rm(list = ls())
library(orthopolynom)
library(pracma)
library(Matrix)
library(matlab)

### PARAMETERS ###
#Heston model parameters with q=0 and lambda = 0
S0 = 100 #Real value 2
V0 = 0.25   # initial volatility
K = 100       # strike price
T = 1/2        # time to maturity in years
r = 0.05      # annual risk-free rate
xi = 0.1       # volatility of the volatility
kappa = 3   # rate reversion theta
theta = 0.04    # longterm volatility
rho = -0.1   # Correlation
q=0
lambda = 0
opttype = 'P' # Option Type 'C' or 'P'

### Parameters for LSM algorithm.
NT = 100 #Steps
NS = 10000 #Trials
simtype = 'M' # Simulation Type 'E' for Euler or 'M' for Miller
BStype = 'R' # Boundary Scheme for negative variance 'A' for absolute or 'R' for reflection

### Parameters for FD algorithm
#We assume that the lower limit is Smin = Vmin = t0 = 0, while the upper limit
Smax = 200
Vmax = 0.75

#Number of steps
m1 = 50 #Number of steps for S
m2 = 25 #Number of steps for V
N = 10 #Number of steps for T
w = 0.5 #weight
method = 'CS'
american = TRUE

#Costants for Non-uniform grid
c = K/5
d = Vmax/500

### Parameters for GBM
NB = 100
sigma = sqrt(V0)

### Parameters for Exact Heston
n = 32
a = 0

