# MSc-Antonio
This project implements two different algorithms to price the American and European option for the Heston model, and measure its error using the exact European option price and the exact Geometric Brownian Motion for pricing European and American options.

## Algorithms

### ADI method.
Four different schemes type are implemented for the Alternating Direction Implicit (ADI) method:
 - Douglas
 - Craig–Sneyd
 - Modified Craig–Sneyd
 - Hundsdorfer–Verwer
The implementation is based on In’t Hout and Foulon (2010).

### Monte Carlo Least-Square method.
The algorithm can be divided in two parts. The first one is the simulation of the Heston paths based on the following discretization:
 - Euler
 - Milstein
 - Moment-matched
The details of the schemes can be found in Chapter 7 from Rouah (2013)
 
The second part of the algorithm is the application of the Least-Square Method (LSM) to the generated paths to get the optimal early exercise of the option determined by comparing the immediate exercise with the conditional expectation based on the previous information. It was proposed by Longstaff and Schwartz (2001).

### Exact Heston model.
The implementation of the original Heston (1993) model using Gauss-Laguerre quadrature to approximate the integral.

### Exact Geometric Brownian Motion.
It implements the exact European Geometric Brownian Motion (GBM) and the exact American GBM solve by Peskir(2005). Credits to Abel Guada Azze, who developed the optimal exercise boundary.

## Tables

### Control variate tables
It compares the exact European option price with the ADI and LSM approximation for European option. Then, it approximates the American option price error by the following equation:

![equation](https://latex.codecogs.com/svg.image?\begin{equation}\label{eq:CVEq}\text{Price}_{\text{True}}^{\text{American}}-\text{Price}_{\text{Approx.}}^{\text{American}}=&space;\text{Price}_{\text{True}}^{\text{European}}-\text{Price}_{\text{Approx.}}^{\text{European}}.\end{equation})

This is known as Control Variate (CV) method and can be used as well to minimize the American option price Error.

### GBM tables
It compares the exact GBM for European and American options with the European and American option using the Heston model with constant volatility (we set the volatility parameters to be 0).

### Comparison table
It compares both, ADI and LSM methods and their computational time.

## Bibliography
 - Heston, S. L. (1993). A closed-form solution for options with stochastic volatility with applications to bond and currency options. The review of financial studies, 6(2),327–343.
 - In’t Hout, K., & Foulon, S. (2010). Adi finite difference schemes for option pricing in the heston model with correlation. International Journal of Numerical Analysis & Modeling, 7(2).
 - Longstaff, F. A., & Schwartz, E. S. (2001). Valuing american options by simulation: A simple least-squares approach. The review of financial studies, 14(1), 113–147.
 - Peskir, G. (2005). On the american option problem. Mathematical Finance: An International Journal of Mathematics, Statistics and Financial Economics, 15(1), 169–181.
 - Rouah, F. (2013). The Heston model and its extensions in Matlab and C#. John Wiley & Sons, Inc.
