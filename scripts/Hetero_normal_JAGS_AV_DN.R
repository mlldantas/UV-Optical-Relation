# From Bayesian Models for Astrophysical Data 
require(R2jags)
require(jagstools)
require(ggplot2)

# Read data

UV_dat <- read.csv("..//Data/av_dn4000_fuv.csv",header=T)

y <- (UV_dat$fuv_flux_obs-UV_dat$fuv_flux_syn)/UV_dat$fuv_flux_obs
x1 <- UV_dat$extinction
x2 <- UV_dat$dn4000_synth
nobs<-nrow(UV_dat)


X <- model.matrix(~ 1 + x1)
K <- ncol(X)
jags_data <- list(Y = y,
                 X  = X,
                 K  = K,
                 N  = nobs)

NORM <-" model{
    # Diffuse normal priors for predictors
    for (i in 1:K) { beta[i] ~ dnorm(0, 0.0001) }
    for (i in 1:K) { alpha[i] ~ dnorm(0, 0.0001) }

    # Likelihood function 
    for (i in 1:N){
    Y[i]~dnorm(mu[i],pow(sigma[i],-2))
    log(sigma[i]) <- inprod(alpha[], X[i,])
    mu[i]  <- eta[i]
    eta[i] <- inprod(beta[], X[i,])
    }
  
    }"


inits <- function () {
  list(
    beta = rnorm(K, 0, 0.01))

}

params <- c("beta", "sigma","alpha")

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
sigmax <- jagsresults(x=jagsfit, params=c('sigma'))


x1x2 <- expand.grid( x1 = x1, x2 = x2)

y.pred <- matrix(predict(fit, newdata = x1x2), 
                 nrow = grid.lines, ncol = grid.lines)



fit<-sigmax[,"mean"]


# scatter plot with regression plane
scatter3D(x1, x2, y, pch = 19,                    
          cex = 0.5, cex.lab=1.5,  
          theta = 130, phi = 25, ticktype = "detailed",
          col="red2",bty = "b2",t="l",
          xlab="x1",
          ylab="x2",
          zlab="y", 
          surf = list(col="cyan",x = x1, y = x2, z = fit,  
                      facets = NA, fit = fitpoints,lwd=1.5,lty=3),colkey = FALSE)




normdata <- data.frame(x1,y)
gdata <- data.frame(x =xx, mean = yx[,"mean"],lwr1=yx[,"25%"],lwr2=yx[,"2.5%"],upr1=yx[,"75%"],upr2=yx[,"97.5%"])


ggplot(normdata,aes(x=x1,y=y))+ 
 geom_ribbon(data=gdata,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL), alpha=0.95, fill=c("gray60"),show.legend=FALSE) +
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL), alpha=0.35, fill = c("gray80"),show.legend=FALSE) +
  geom_point(size=3,colour="cyan3",alpha=0.4)+
  geom_line(data=gdata,aes(x=xx,y=mean),colour="gray25",linetype="dashed",size=1,show.legend=FALSE)+
  theme_bw()+xlab("g-r")+ylab("FUV-r")+coord_cartesian(xlim=c(-1,1.5),ylim=c(0,6))



