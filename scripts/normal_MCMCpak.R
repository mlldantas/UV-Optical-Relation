# From Bayesian Models for Astrophysical Data 
# by Hilbe, de Souza & Ishida, 2016, Cambridge Univ. Press

# Chapter 4 - Example of Bayesian linear regression in R using MCMCregress
# 1 response (y) and 1 explanatory variable (x1) 

library(MCMCpack)


nobs= 5000                      # number of obs in model 
x1 <- runif(nobs)               # random uniform variab;le
x2 <- runif(nobs)               # random uniform variable
xb <- 2 + 3*x1       # linear predictor, xb
y <- rnorm(nobs, xb, sd=1)      

posteriors <- MCMCregress(y ~ x1, 
                          thin=1,
                          seed=1056,
                          burnin = 2500,
                          mcmc = 5000,
                          verbose=1)
summary(posteriors)

plot(posteriors)
