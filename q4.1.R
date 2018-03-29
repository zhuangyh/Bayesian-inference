setwd("~/Dropbox/Biostatics_class/BIOS7717/HW3")
bw <- read.table("bweight.dat", header=TRUE)
lmod = lm(bweight ~ packs + folate, data=bw)
summary(lmod)

y <- bw$bweight
bw$x0 <- 1 
X <- as.matrix(bw[, c("x0", "packs", "folate")])

library(mvtnorm)
library(invgamma)
library(ggplot2)
library(dplyr)
library(tidyr)

# function for blocked gibbs sampler
gibbs_sampler<-function(y, x, iter, burnin, trim, a=0.1, g=1000){
  xprimex_inv<-solve(t(x)%*%x) # vectorization
  tau<-numeric(iter) 
  b<-matrix(nrow=iter, ncol = 3) 
  tau[1]<-1 
  
  # gibbs sampling
  for(i in 2:iter ){
    b[i,]<-rmvnorm(n = 1, 
                   mean = ((xprimex_inv%*%t(x))%*%y), 
                   sigma = tau[i-1]*xprimex_inv )
    
    tau[i]<-rinvgamma(n = 1, 
                      shape = (n/2 + a), 
                      rate = .5*( t((y - x%*%t(t(b[i,])) ))%*%(y - x%*%t(t(b[i,])) ) ) + g)
  }
  
  # burnin and trimming  
  keep_draws<-seq(burnin,iter,trim)
  tau<-tau[keep_draws]
  b<-b[keep_draws,]
  
  joint_post<-data.frame(b=b,tau=tau)
  colnames(joint_post)[1:(ncol(x))]<-paste0('B',0:(ncol(x)-1) )
  
  joint_post<-gather(joint_post,keep_draws) %>%
    rename(param=keep_draws, draw=value) %>%
    mutate(iter=rep(keep_draws,ncol(joint_post)))
  
  return(joint_post)
}

# run gibbs sampler with specified parameters
post_dist<-gibbs_sampler(y = y, x = X, iter = 200000, burnin = 5000, trim= 50)

# calculate posterior summary statistics (stats not used in rest of code)
post_sum_stats<-post_dist %>%
  group_by(param) %>%
  summarise(mean=mean(draw),
            median=median(draw),
            lwr=quantile(draw,.025),
            upr=quantile(draw,.975))
post_sum_stats
write.csv(post_sum_stats, "post_sum_stats.csv")
# merge on summary statistics
post_dist <- post_dist %>%
  left_join(post_sum_stats, by='param')

# plot MCMC Chains
ggplot(post_dist,aes(x=iter,y=draw)) +
  geom_line() +
  facet_grid(param ~ .,scale='free_y',switch = 'y') +
  theme_bw() + 
  xlab('Gibbs Sample Iteration') + ylab('MCMC Chains') + 
  ggtitle('Gibbs Sampler MCMC Chains by Parameter')

# plot Posterior Distributions
ggplot(post_dist,aes(x=draw)) +
  geom_histogram(aes(x=draw),bins=50) +
  facet_grid(. ~ param, scale='free_x',switch = 'y') +
  theme_bw() + 
  xlab('Posterior Distributions') + ylab('Count') + 
  ggtitle('Posterior Distributions of Parameters')


# run gibbs sampler with specified parameters
post_dist_1<-block_gibbs_sampler(y = y, x = X, iter = 200000, burnin = 5000, trim= 50, a=0.001, g=0.001 )
# run gibbs sampler with specified parameters
post_dist_2<-block_gibbs_sampler(y = y, x = X, iter = 200000, burnin = 5000, trim= 50, a=1000, g=1000 )


tau0 <-post_dist[which(post_dist$param=='tau'), "draw"]
tau1 <-post_dist_1[which(post_dist_1$param=='tau'), "draw"]

tau2 <-post_dist_2[which(post_dist_2$param=='tau'), "draw"]
taus <- cbind(tau0, tau1, tau2)
head(taus)
colnames(taus) <- c("referce", "prior1", "prior2")


tau_long = melt(taus)

p_tau <- ggplot(aes(x=value, colour=X2), data=tau_long) 
p_tau 
  geom_density() +
  
  scale_colour_manual(name = "Priors",
                      labels = c("alpha=0.1, gamma=1000", 
                                 "alpha=0.001, gamma=0.001", 
                                 "alpha=1000, g=1000"),
                      values = c("black", "red", "blue")) +
  ggtitle("Impact of tau priors")


head(tau0)
taus <- cbind(post_dist$tau, post_dist_1$tau, post_dist_2$tau)
head(taus)

library("rjags")
# Reference prior
mod1_string = " model {
                  for (i in 1:n) {
                    y[i] ~ dnorm(mu[i], prec)
                    mu[i] = b[1] + b[2]*packs[i] + b[3]*folate[i]
                  }
                    for (i in 1:3) {
                    b[i] ~ dnorm(0.0, 1.0/1.0e6)
                    }
                  prec ~ dgamma(0.001, 0.001)
                  sig2 = 1.0 / prec
                  sig = sqrt(sig2)
                } "

set.seed(72)
bw_jags = list(y=bw$bweight, n=nrow(bw), 
                  packs =bw$packs, 
                  folate = bw$folate)
params1 = c("b", "sig")
mod1 = jags.model(textConnection(mod1_string), data=bw_jags, n.chains=1)
update(mod1, 1000) # burn-in

mod1_sim = coda.samples(model=mod1,
                        variable.names=params1,
                        n.iter=5000)
summary(mod1_sim)


# Informative prior : b[1] ~ dnorm(50, 1.0/1.0e3)
mod2_string = " model {
      for (i in 1:n) {
      y[i] ~ dnorm(mu[i], prec)
      mu[i] = b[1] + b[2]*packs[i] + b[3]*folate[i]
      }
      b[1] ~ dnorm(50, 1.0/1.0e3)
      b[2] ~ dnorm(0.0, 1.0/1.0e6)
      b[3] ~ dnorm(0.0, 1.0/1.0e6)
      prec ~ dgamma(0.001, 0.001)
      sig2 = 1.0 / prec
      sig = sqrt(sig2)
      } "

set.seed(72)
bw_jags = list(y=bw$bweight, n=nrow(bw), 
               packs =bw$packs, 
               folate = bw$folate)
params1 = c("b", "sig")
mod2 = jags.model(textConnection(mod2_string), data=bw_jags, n.chains=1)
update(mod2, 1000) # burn-in

mod2_sim = coda.samples(model=mod2,
                        variable.names=params1,
                        n.iter=5000)
summary(mod2_sim)

# Informative prior 3: b[1] ~ dnorm(50, 1.0/1.0e1)
mod3_string = " model {
    for (i in 1:n) {
    y[i] ~ dnorm(mu[i], prec)
    mu[i] = b[1] + b[2]*packs[i] + b[3]*folate[i]
    }
    b[1] ~ dnorm(50, 1.0/1.0e2)
    b[2] ~ dnorm(0.0, 1.0/1.0e6)
    b[3] ~ dnorm(0.0, 1.0/1.0e6)
    prec ~ dgamma(0.001, 0.001)
    sig2 = 1.0 / prec
    sig = sqrt(sig2)
    } "

set.seed(72)
mod3 = jags.model(textConnection(mod3_string), data=bw_jags, n.chains=1)
update(mod3, 1000) # burn-in

mod3_sim = coda.samples(model=mod3,
                        variable.names=params1,
                        n.iter=5000)
summary(mod3_sim)

# Informative prior 3: b[1] ~ dnorm(25, 1.0/1.0e1)
mod4_string = " model {
for (i in 1:n) {
y[i] ~ dnorm(mu[i], prec)
mu[i] = b[1] + b[2]*packs[i] + b[3]*folate[i]
}
b[1] ~ dnorm(10, 1.0)
b[2] ~ dnorm(0.0, 1.0/1.0e6)
b[3] ~ dnorm(0.0, 1.0/1.0e6)
prec ~ dgamma(0.001, 0.001)
sig2 = 1.0 / prec
sig = sqrt(sig2)
} "

set.seed(72)
mod4 = jags.model(textConnection(mod3_string), data=bw_jags, n.chains=1)
update(mod4, 1000) # burn-in

mod4_sim = coda.samples(model=mod4,
                        variable.names=params1,
                        n.iter=5000)
summary(mod4_sim)


priors_intercept <- cbind(MCMCchains(mod1_sim, params = "b")[, 1],
                MCMCchains(mod2_sim, params = 'b')[, 1],
                MCMCchains(mod3_sim, params = 'b')[, 1],
                MCMCchains(mod4_sim, params = 'b')[, 1])
colnames(priors_intercept) <- c("Prior0","Prior1", "Prior2", "Prior3")

str(mod1_sim)
long = melt(priors_intercept)

p <- ggplot(aes(x=value, colour=X2), data=long)
p + geom_density() +
  geom_vline(aes(xintercept = 1/0.3193,col='red'), show.legend = FALSE) +
  scale_colour_manual(name = "Priors",
                      labels = c("Reference prior (b0=0, tau=0.000001)", 
                                 "b0=50, tau=0.001", 
                                 "b0=50, tau=0.01", 
                                 "b0=10, tau=1.0"),
                      values = c("black", "pink", "red", "blue")) +
  ggtitle("Impact of intercept priors")
