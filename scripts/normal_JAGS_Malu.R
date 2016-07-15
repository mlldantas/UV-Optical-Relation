# From Bayesian Models for Astrophysical Data 
require(R2jags)
require(jagstools)
require(ggplot2)

# Read data

UV_dat <- read.table("..//data_trial/output_results.txt",header=T)

y <- UV_dat$fuv_mag-UV_dat$dered_r
x1 <- UV_dat$dered_g-UV_dat$dered_r
nobs<-nrow(UV_dat)

# Prepare data for prediction 
M=500
xx = seq(from =  min(x1), 
         to =  max(x1), 
         length.out = M)


X <- model.matrix(~ 1 + x1)
K <- ncol(X)
jags_data <- list(Y = y,
                 X  = X,
                 K  = K,
                 N  = nobs,
                 M = M,
                 xx= xx)


NORM <-" model{
    # Diffuse normal priors for predictors
    for (i in 1:K) { beta[i] ~ dnorm(0, 0.0001) }
    
   # Uniform prior for standard deviation
     tau <- pow(sigma, -2)       # precision
     sigma ~ dunif(0, 100)       # standard deviation
   
    # Likelihood function 
    for (i in 1:N){
    Y[i]~dnorm(mu[i],tau)
    mu[i]  <- eta[i]
    eta[i] <- inprod(beta[], X[i,])
    }
   # Prediction for new data
   for (j in 1:M){
   etax[j]<-beta[1]+beta[2]*xx[j]
   mux[j]  <- etax[j]
   Yx[j]~dnorm(mux[j],tau)
}
    }"


inits <- function () {
  list(
    beta = rnorm(K, 0, 0.01))

}

params <- c("beta", "sigma","Yx")

jagsfit <- jags(
           data       = jags_data,
           inits      = inits,
           parameters = params,
           model      = textConnection(NORM),
           n.chains   = 3,
           n.iter     = 5000,
           n.thin     = 1,
           n.burnin   = 2500)


print(jagsfit,justify = "left", digits=2)





# Plot
yx <- jagsresults(x=jagsfit, params=c('Yx'))


normdata <- data.frame(x1,y)
gdata <- data.frame(x =xx, mean = yx[,"mean"],lwr1=yx[,"25%"],lwr2=yx[,"2.5%"],upr1=yx[,"75%"],upr2=yx[,"97.5%"])


ggplot(normdata,aes(x=x1,y=y))+ 
 geom_ribbon(data=gdata,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL), alpha=0.95, fill=c("gray60"),show.legend=FALSE) +
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL), alpha=0.35, fill = c("gray80"),show.legend=FALSE) +
  geom_point(size=3,colour="cyan3",alpha=0.4)+
  geom_line(data=gdata,aes(x=xx,y=mean),colour="gray25",linetype="dashed",size=1,show.legend=FALSE)+
  theme_bw()+xlab("g-r")+ylab("FUV-r")+coord_cartesian(xlim=c(-1,1.5),ylim=c(0,6))



