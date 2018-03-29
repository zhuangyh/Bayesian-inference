
setwd("~/Dropbox/Biostatics_class/BIOS7717/HW3")
library(rjags)
library(mcmcse)

hpd = function(x, alpha = 0.05)
{
  n = length(x)
  m = round(n * alpha)
  x = sort(x)
  y = x[(n - m + 1):n] - x[1:m]
  z = min(y)
  k = which(y == z)[1]
  c(x[k], x[n - m + k])
}

n=12
y=c(4.20, 4.36, 4.11, 3.96, 5.63, 4.50, 5.64, 4.38,4.45, 3.67, 5.26, 4.66)
mean(y)
1/var(y)


# lg = function(mu, tau, y, mu0, b, c, d) {
#   n= length(y)
#   lh = -0.5*n*log(tau) - 0.5*tau*sum((y-mu)^2)
#   lh = lh - (mu - mu0)^2/(2*b) - 0.5*log(b)+ (c-1)*log(tau) - d*tau
#   return (lh)
# }



lg = function(mu, tau, y, mu0=0, b=10000, c=0.001, d=0.001) {
  n= length(y)
  lh = (0.5*n + c +1)*log(tau)+ 0.5*log(b) - 0.5*tau*sum((y-mu)^2)
  lh = lh - (mu-mu0)^2/(2*b) - d*tau
  return (lh)
}



mh = function(y, n_iter, mu_init, tau_init, mu0, b, c, d, cand_sd1, cand_sd2) {
  ## Random-Walk Metropolis-Hastings algorithm
  
  ## step 1, initialize
  mu_out = tau_out= numeric(n_iter)
  accpt = 0
  mu_now = mu_init
  tau_now = tau_init
  
  lg_now = lg(mu=mu_now, tau=tau_now, y, mu0, b, c, d)
  
  ## step 2, iterate
  for (i in 1:n_iter) {
    ## step 2a
    mu_cand = rnorm(n=1, mean=mu_now, sd=cand_sd1) # draw a candidate
    tau_cand = rnorm(n=1, mean=tau_now, sd=cand_sd2) # draw a candidate
    
    if (tau_cand > 0){
      ## step 2b
      lg_cand = lg(mu=mu_cand, tau=tau_cand, y, mu0, b, c, d) # evaluate log of g with the candidate
      lalpha = lg_cand - lg_now # log of acceptance ratio
      alpha = exp(lalpha)
      
      ## step 2c
      u = runif(1) # draw a uniform variable which will be less than alpha with probability min(1, alpha)
      if (u < alpha) { # then accept the candidate
        mu_now = mu_cand
        tau_now = tau_cand 
        accpt = accpt + 1 # to keep track of acceptance
        lg_now = lg_cand
      }}
    
    ## collect results
    mu_out[i] = mu_now # save this iteration's value of mu
    tau_out[i] = tau_now
  }
  
  ## return a list of output
  list(mu=mu_out, tau=tau_out, accpt=accpt/n_iter)
  #list(mu=mu_out[-1000,], sigma=sigma_out[-1000,], accpt=accpt/n_iter)
}



n=10000

hist(y, freq=FALSE) # histogram of the data
#curve(dpearson7(x, mu=10, sigma=20), -300, 300, col = "navy", add=TRUE) # prior for mu

post_1 = mh(y, n_iter=1e5, mu_init=0, tau_init=1, mu0=0, b=1e6, c=0.001, d=0.001,
            cand_sd1=0.13, cand_sd2=0.1)
post_1$accpt
mean(post_1$mu)
mean(post_1$tau)


post_2 = mh(y, n_iter=1e5, mu_init=0, tau_init=1, mu0=4.5, b=1e6, c=0.001, d=0.001,
            cand_sd1=0.23, cand_sd2=0.1)
post_2$accpt
mean(post_2$mu)
mean(post_2$tau)


post_3 = mh(y, n_iter=1e5, mu_init=0, tau_init=1, mu0=0, b=1e2, c=0.001, d=0.001,
            cand_sd1=0.23, cand_sd2=0.1)
post_3$accpt
mean(post_3$mu)
mean(post_3$tau)


post_4 = mh(y, n_iter=1e5, mu_init=0, tau_init=1, mu0=0, b=1e2, c=0.1, d=0.1,
            cand_sd1=0.23, cand_sd2=0.1)
post_4$accpt
mean(post_4$mu)
mean(post_4$tau)

post_5 = mh(y, n_iter=1e5, mu_init=0, tau_init=1, mu0=0, b=1e2, c=1, d=1,
            cand_sd1=0.23, cand_sd2=0.1)
post_5$accpt
mean(post_5$mu)
mean(post_5$tau)


mhfit_mu_1 = data.frame(mcse.mat(as.matrix(post_1$mu)))
mhfit_mu_1 <- cbind(mhfit_mu_1, hpd(post_1$mu)[1], hpd(post_1$mu)[2])
mhfit_tau_1 = data.frame(mcse.mat(as.matrix(post_1$tau)))
mhfit_tau_1 <- cbind(mhfit_tau_1, hpd(post_1$tau)[1], hpd(post_1$tau)[2])
mhfit_est_1 <- cbind(mhfit_mu_1, mhfit_tau_1)


mhfit_mu_2 = data.frame(mcse.mat(as.matrix(post_2$mu)))
mhfit_mu_2 <- cbind(mhfit_mu_2, hpd(post_2$mu)[1], hpd(post_2$mu)[2])
mhfit_tau_2 = data.frame(mcse.mat(as.matrix(post_2$tau)))
mhfit_tau_2 <- cbind(mhfit_tau_2, hpd(post_2$tau)[1], hpd(post_2$tau)[2])
mhfit_est_2 <- cbind(mhfit_mu_2, mhfit_tau_2)


mhfit_mu_3 = data.frame(mcse.mat(as.matrix(post_3$mu)))
mhfit_mu_3 <- cbind(mhfit_mu_3, hpd(post_3$mu)[1], hpd(post_3$mu)[2])
mhfit_tau_3 = data.frame(mcse.mat(as.matrix(post_3$tau)))
mhfit_tau_3 <- cbind(mhfit_tau_3, hpd(post_3$tau)[1], hpd(post_3$tau)[2])
mhfit_est_3 <- cbind(mhfit_mu_3, mhfit_tau_3)

mhfit_mu_4 = data.frame(mcse.mat(as.matrix(post_4$mu)))
mhfit_mu_4 <- cbind(mhfit_mu_4, hpd(post_4$mu)[1], hpd(post_4$mu)[2])
mhfit_tau_4 = data.frame(mcse.mat(as.matrix(post_4$tau)))
mhfit_tau_4 <- cbind(mhfit_tau_4, hpd(post_4$tau)[1], hpd(post_4$tau)[2])
mhfit_est_4 <- cbind(mhfit_mu_4, mhfit_tau_4)

mhfit_mu_5 = data.frame(mcse.mat(as.matrix(post_5$mu)))
mhfit_mu_5 <- cbind(mhfit_mu_5, hpd(post_5$mu)[1], hpd(post_5$mu)[2])
mhfit_tau_5 = data.frame(mcse.mat(as.matrix(post_5$tau)))
mhfit_tau_5 <- cbind(mhfit_tau_5, hpd(post_5$tau)[1], hpd(post_5$tau)[2])
mhfit_est_5 <- cbind(mhfit_mu_5, mhfit_tau_5)

names <- c("mu", "mu_se", "mu_HPD_low", "mu_HPD_up", 
  "tau", "tau_se", "tau_HPD_low", "tau_HPD_up")

colnames(mhfit_est_1) = colnames(mhfit_est_2)= colnames(mhfit_est_3) = names
colnames(mhfit_est_4) = colnames(mhfit_est_5)= names

mhfit_est <-rbind( mhfit_est_1, mhfit_est_2, mhfit_est_3, mhfit_est_4, mhfit_est_5)
rownames(mhfit_est) <- c("mu0=0, b=1e6, c=0.001, d=0.001", 
                         "mu0=4.5, b=1e6, c=0.001, d=0.001", 
                         "mu0=0, b=1e2, c=0.001, d=0.001", 
                          "mu0=0, b=1e2, c=0.1, d=0.1",
                         "mu0=0, b=1e2, c=1, d=1")
mhfit_est
write.csv(mhfit_est, "mhfit_est.csv")


library("coda")
length(post_1$mu)
plot(as.mcmc(post_1$mu[100:length(post_1$mu)]))
plot(as.mcmc(post_1$tau[100:length(post_1$tau)]))

data.c <- list(y=c(4.20,4.36,4.11,3.96,5.63,4.50,
                   5.64,4.38,4.45,3.67,5.26,4.66),
               n=12)

model_c0.string <-"
    model{
      for(i in 1:n){ y[i] ~ dnorm(mu, tau) }
      mu ~ dnorm(4.56, 1/1e6)
      tau ~ dgamma(0.001,0.001)
      sigma <- 1/sqrt(tau)
      }
      "
model_c0.spec<-textConnection(model_c0.string)

# Initialize the JAGS model, and do an adaptation run of length 1,000.
c0.jags <- jags.model(model_c0.spec,
                      data = data.c,
                      n.chains=1,
                      n.adapt=1000)
update(c0.jags, 1000)
params <- c("mu", "tau")
# Now draw 10,000 samples from the posterior.
iter = 100000
c0_samples = try(coda.samples(c0.jags, params, iter), TRUE)
summary(c0_samples)


### 95% certainty
model_c1.string <-"
    model{
      for(i in 1:n){ y[i] ~ dnorm(mu, tau) }
      mu ~ dnorm(4.75, 1/0.0163)
      tau ~ dgamma(0.001,0.001)
      sigma <- 1/sqrt(tau)
      }
      "
model_c1.spec<-textConnection(model_c1.string)

# Initialize the JAGS model, and do an adaptation run of length 1,000.
c1.jags <- jags.model(model_c1.spec,
                     data = data.c,
                     n.chains=1,
                     n.adapt=1000)
update(c1.jags, 1000)
params <- c("mu")
c1_samples = try(coda.samples(c1.jags, params, iter), TRUE)
summary(c1_samples)

### 60% certainty
model_c2.string <-"
    model{
        for(i in 1:n){ y[i] ~ dnorm(mu, tau) }
        mu ~ dnorm(4.75, 1/0.0882)
        tau ~ dgamma(0.001,0.001)
        sigma <- 1/sqrt(tau)
      }
      "
model_c2.spec<-textConnection(model_c2.string)

# Initialize the JAGS model, and do an adaptation run of length 1,000.
c2.jags <- jags.model(model_c2.spec,
                      data = data.c,
                      n.chains=1,
                      n.adapt=1000)
update(c2.jags, 1000)
c2_samples = try(coda.samples(c2.jags, params, iter), TRUE)
summary(c2_samples)


### 40% certainty
model_c3.string <-"
    model{
      for(i in 1:n){ y[i] ~ dnorm(mu, tau) }
        mu ~ dnorm(4.75, 1/0.2267)
        tau ~ dgamma(0.001,0.001)
        sigma <- 1/sqrt(tau)
      }
      "
model_c3.spec<-textConnection(model_c3.string)
c3.jags <- jags.model(model_c3.spec,
                      data = data.c,
                      n.chains=1,
                      n.adapt=1000)
update(c3.jags, 1000)
c3_samples = try(coda.samples(c3.jags, c("mu", "tau"), iter), TRUE)
summary(c3_samples)
plot(c3_samples)

table_q3 <- rbind(cbind(summary(c0_samples)$statistics[1], HPDinterval(c0_samples, prob = 0.95)[[1]]),
                  cbind(summary(c1_samples)$statistics[1], HPDinterval(c1_samples, prob = 0.95)[[1]]),
                   cbind(summary(c2_samples)$statistics[1], HPDinterval(c2_samples, prob = 0.95)[[1]]),
                   cbind(summary(c3_samples)$statistics[1], HPDinterval(c3_samples, prob = 0.95)[[1]]))
               
rownames(table_q3) <- c("N=4.56, Var=1e6",
                      "N=4.75, Var=0.0163", 
                     "N=4.75, Var=0.0882", 
                     "N=4.75, Var=0.2267")
colnames(table_q3) <- c("posterior_mean", "HPD_low", "HPD_high")
table_q3 <- as.data.frame(table_q3)
table_q3


write.csv(table_q3, "table_q3.csv")

priors_q3c <- cbind(MCMCchains(c0_samples, params = 'mu'),
                MCMCchains(c1_samples, params = 'mu'),
                MCMCchains(c2_samples, params = 'mu'),
                MCMCchains(c3_samples, params = 'mu'))
colnames(priors_q3c) <- c("Reference Prior","stong_prior","Weak_prior", "Weaker_prior")


long_prior_q3c = melt(priors_q3c)
head(long_prior_q3c)

p_q3c <- ggplot(aes(x=value, colour=X2), data=long_prior_q3c)
p_q3c + geom_density() +
  ggtitle("Impact of mu priors on the posterior mu" ) +
  scale_colour_manual(name = "Priors",
                      labels = c("N=4.56, Var=1e6 (flat prior)", 
                                "N=4.75, Var=0.0163", 
                                 "N=4.75, Var=0.0882", 
                                 "N=4.75, Var=0.2267"),
                      values = c("black", "red", "blue", "orange")) 



model_c0.string <-"
    model{
for(i in 1:n){ y[i] ~ dnorm(mu, tau) }
mu ~ dnorm(4.56, 1/1e6)
tau ~ dgamma(0.001,0.001)
sigma <- 1/sqrt(tau)
}
"
model_c0.spec<-textConnection(model_c0.string)

# Initialize the JAGS model, and do an adaptation run of length 1,000.
c0.jags <- jags.model(model_c0.spec,
                      data = data.c,
                      n.chains=1,
                      n.adapt=1000)
update(c0.jags, 1000)
params <- c("tau")
# Now draw 10,000 samples from the posterior.
iter = 100000
c0_samples = try(coda.samples(c0.jags, params, iter), TRUE)
summary(c0_samples)


### 95% certainty
model_c1.string <-"
model{
for(i in 1:n){ y[i] ~ dnorm(mu, tau) }
mu ~ dnorm(4.75, 1/0.0163)
tau ~ dgamma(0.001,0.001)
sigma <- 1/sqrt(tau)
}
"
model_c1.spec<-textConnection(model_c1.string)

# Initialize the JAGS model, and do an adaptation run of length 1,000.
c1.jags <- jags.model(model_c1.spec,
                      data = data.c,
                      n.chains=1,
                      n.adapt=1000)
update(c1.jags, 1000)

c1_samples = try(coda.samples(c1.jags, params, iter), TRUE)
summary(c1_samples)

### 60% certainty
model_c2.string <-"
model{
for(i in 1:n){ y[i] ~ dnorm(mu, tau) }
mu ~ dnorm(4.75, 1/0.0882)
tau ~ dgamma(0.001,0.001)
sigma <- 1/sqrt(tau)
}
"
model_c2.spec<-textConnection(model_c2.string)

# Initialize the JAGS model, and do an adaptation run of length 1,000.
c2.jags <- jags.model(model_c2.spec,
                      data = data.c,
                      n.chains=1,
                      n.adapt=1000)
update(c2.jags, 1000)
c2_samples = try(coda.samples(c2.jags, params, iter), TRUE)
summary(c2_samples)


### 40% certainty
model_c3.string <-"
model{
for(i in 1:n){ y[i] ~ dnorm(mu, tau) }
mu ~ dnorm(4.75, 1/0.2267)
tau ~ dgamma(0.001,0.001)
sigma <- 1/sqrt(tau)
}
"
model_c3.spec<-textConnection(model_c3.string)
c3.jags <- jags.model(model_c3.spec,
                      data = data.c,
                      n.chains=1,
                      n.adapt=1000)
update(c3.jags, 1000)
c3_samples = try(coda.samples(c3.jags, params, iter), TRUE)
summary(c3_samples)
plot(c3_samples)

table_q3_tau <- rbind(cbind(summary(c0_samples)$statistics[1], HPDinterval(c0_samples, prob = 0.95)[[1]]),
                  cbind(summary(c1_samples)$statistics[1], HPDinterval(c1_samples, prob = 0.95)[[1]]),
                  cbind(summary(c2_samples)$statistics[1], HPDinterval(c2_samples, prob = 0.95)[[1]]),
                  cbind(summary(c3_samples)$statistics[1], HPDinterval(c3_samples, prob = 0.95)[[1]]))

rownames(table_q3_tau) <- c("N=4.56, Var=1e6",
                        "N=4.75, Var=0.0163", 
                        "N=4.75, Var=0.0882", 
                        "N=4.75, Var=0.2267")
colnames(table_q3_tau) <- c("posterior_mean", "HPD_low", "HPD_high")
table_q3_tau <- as.data.frame(table_q3_tau)
table_q3_tau


write.csv(table_q3_tau, "table_q3_tau.csv")

priors_q3c_tau <- cbind(MCMCchains(c0_samples, params = 'tau'),
                    MCMCchains(c1_samples, params = 'tau'),
                    MCMCchains(c2_samples, params = 'tau'),
                    MCMCchains(c3_samples, params = 'tau'))
colnames(priors_q3c_tau) <- c("Reference Prior","stong_prior","Weak_prior", "Weaker_prior")


long_prior_q3c_tau = melt(priors_q3c_tau)
head(long_prior_q3c_tau)

p_q3c_tau <- ggplot(aes(x=value, colour=X2), data=long_prior_q3c_tau)
p_q3c_tau + geom_density() +
  ggtitle("Impact of priors on the posterior tau" ) +
  scale_colour_manual(name = "Priors",
                      labels = c("N=4.56, Var=1e6 (flat prior)", 
                                 "N=4.75, Var=0.0163", 
                                 "N=4.75, Var=0.0882", 
                                 "N=4.75, Var=0.2267"),
                      values = c("black", "red", "blue", "orange")) 


data.g <- list(y=c(4.20,4.36,4.11,3.96,5.63,4.50,
                5.64,4.38,4.45,3.67,5.26,4.66,NA),
                c=0.001, d=0.001,n=13)

model1.string <-"
    model{
      for(i in 1:n){ y[i] ~ dnorm(mu, tau) }
      mu ~ dnorm(4.75, 0.0163)
      tau ~ dgamma(0.001,0.001)
      sigma <- 1/sqrt(tau)
      gamma <- phi((4.4-mu)/sqrt(1/tau))
      prob <- step(4.4 - y[13])
    }
"

model1.spec<-textConnection(model1.string)


# Initialize the JAGS model, and do an adaptation run of length 1,000.
g.jags <- jags.model(model1.spec,
                       data = data.g,
                       n.chains=1,
                       n.adapt=1000)
update(g.jags, 1000)
bayes.model.params <- c("mu", "tau", "sigma", "gamma", "prob")


# Now draw 10,000 samples from the posterior.
iter = 10000
g_samples = try(coda.samples(g.jags, bayes.model.params, iter), TRUE)
head(g_samples)
library(MCMCvis)
prob_13 <- MCMCchains(g_samples, params = 'prob')
length(prob_13)
sum(prob_13)/iter

# summary of modle with logit link 
summary(g_samples)

pnorm((4.4 - 4.56)/sqrt(1/2.5), 0, 1)

y13 <- cbind(MCMCchains(c0_samples, params = 'mu'), 
             MCMCchains(c0_samples, params = 'tau'))
y13 <-as.data.frame(y13)

y13$prob <- pnorm((4.4 - y13$mu)/sqrt(1/y13$tau), 0, 1)
mean(y13$prob)
hpd(y13$prob)



dnorm(0)

model{
  # observation model
  for (t in 1:T){
    y[t] ~ dnorm(mu[t], V.inv)
    mu[t] <- x[t]*beta[t]
  }
  # state model
  for (t in 2:T){
    beta[t] ~ dnorm(beta[t-1], W.inv)
  }
  # settings for t=1
  beta[1] ~ dnorm(10,0.01)
  # priors
  ...
}

ref_summary <-cbind(summary(ref_coda_samples)$statistics[, c(1, 2)], HPDinterval(ref_coda_samples, prob = 0.95)[[1]])
ref_summary <- ref_summary[c("theta1", "theta2", "riskDifference", "relativeRisk", "oddsRatio", "prob_RD"), ] 
ref_summary
write.csv(ref_summary, file = "q1_ref_summary.csv")


### Data Augmentation priors

model2.string <-"
model {
y1 ~ dbin(theta1,n1)
y2 ~ dbin(theta2,n2)
theta1 ~ dbeta(10,10)
theta2 ~ dbeta(10,10)
odds1 <- theta1/(1-theta1)
odds2 <- theta2/(1-theta2)
riskDifference <- theta2-theta1
relativeRisk <- theta2/theta1
oddsRatio <- odds2/odds1
prob_RD <- step(riskDifference)
}
"

model2.spec<-textConnection(model2.string)

# Initialize the JAGS model, and do an adaptation run of length 1,000.
dataAug.jags <- jags.model(model2.spec,
                           data = medfly,
                           n.chains=2,
                           n.adapt=1000)
update(dataAug.jags, 1000)

dataAug_coda_samples = try(coda.samples(dataAug.jags, bayes.model.params, iter), TRUE)

dataAug_summary <-cbind(summary(dataAug_coda_samples)$statistics[, c(1, 2)], HPDinterval(dataAug_coda_samples, prob = 0.95)[[1]])
dataAug_summary <- dataAug_summary[c("theta1", "theta2", "riskDifference", "relativeRisk", "oddsRatio", "prob_RD"), ] 
dataAug_summary
write.csv(dataAug_summary, file = "q1_dataAug_summary.csv")


### Expert priors
### theta1 ~ beta(2.07, 4.83)
qbeta(0.95, 2.07, 4.83)
### theta2 ~ beta(1.53, 1.53)
qbeta(0.95, 1.53, 1.53)


model3.string <-"
model {
y1 ~ dbin(theta1,n1)
y2 ~ dbin(theta2,n2)
theta1 ~ dbeta(2.07, 4.83)
theta2 ~ dbeta(1.53, 1.53)
odds1 <- theta1/(1-theta1)
odds2 <- theta2/(1-theta2)
riskDifference <- theta2-theta1
relativeRisk <- theta2/theta1
oddsRatio <- odds2/odds1
prob_RD <- step(riskDifference)
}
"

model3.spec<-textConnection(model3.string)
# Initialize the JAGS model, and do an adaptation run of length 1,000.
expert.jags <- jags.model(model3.spec,
                          data = medfly,
                          n.chains=2,
                          n.adapt=1000)
update(expert.jags, 1000)
# Now draw 10,000 samples from the posterior.
expert_coda_samples = try(coda.samples(expert.jags, bayes.model.params, iter), TRUE)

expert_summary <-cbind(summary(expert_coda_samples)$statistics[, c(1, 2)], HPDinterval(expert_coda_samples, prob = 0.95)[[1]])
expert_summary <- expert_summary[c("theta1", "theta2", "riskDifference", "relativeRisk", "oddsRatio", "prob_RD"), ] 
expert_summary
ref_summary
write.csv(expert_summary, file = "q1_expert_summary.csv")


### explore prior effect
### Case 0 (reference prior)
## alpha/(alpha + beta) = 0.5, alpha + beta = 2
## alpha = 1, beta = 1
case0.string <-"
model {
y1 ~ dbin(theta1,n1)
y2 ~ dbin(theta2,n2)
theta1 ~ dbeta(0.6, 1.4)
theta2 ~ dbeta(0.6, 1.4)
}
"
case0.spec<-textConnection(case0.string)
# Initialize the JAGS model, and do an adaptation run of length 1,000.
case0.jags <- jags.model(case0.spec,
                         data = medfly,
                         n.chains=2,
                         n.adapt=1000)
update(case0.jags, 1000)
# Now draw 10,000 samples from the posterior.
iter = 1000
case0_samples = try(coda.samples(case0.jags, c("theta1"), iter), TRUE)






### Case 1
## alpha/(alpha + beta) = 0.3, alpha + beta = 2
## alpha = 0.6, beta = 1.4
case1.string <-"
model {
y1 ~ dbin(theta1,n1)
y2 ~ dbin(theta2,n2)
theta1 ~ dbeta(0.6, 1.4)
theta2 ~ dbeta(0.6, 1.4)
}
"
case1.spec<-textConnection(case1.string)
# Initialize the JAGS model, and do an adaptation run of length 1,000.
case1.jags <- jags.model(case1.spec,
                         data = medfly,
                         n.chains=2,
                         n.adapt=1000)
update(case1.jags, 1000)
# Now draw 10,000 samples from the posterior.
iter = 1000
case1_samples = try(coda.samples(case1.jags, c("theta1"), iter), TRUE)


### Case 2
## alpha/(alpha + beta) = 0.3, alpha + beta = 20
## alpha = 6, beta = 14

case2.string <-"
model {
y1 ~ dbin(theta1,n1)
y2 ~ dbin(theta2,n2)
theta1 ~ dbeta(6, 14)
theta2 ~ dbeta(6, 14)
}
"
case2.spec<-textConnection(case2.string)
# Initialize the JAGS model, and do an adaptation run of length 1,000.
case2.jags <- jags.model(case2.spec,
                         data = medfly,
                         n.chains=2,
                         n.adapt=1000)
update(case2.jags, 1000)
# Now draw 1000 samples from the posterior.
iter = 1000
case2_samples = try(coda.samples(case2.jags, c("theta1"), iter), TRUE)



### Case 3
## alpha/(alpha + beta) = 0.3, alpha + beta = 200
## alpha = 60, beta = 140
case3.string <-"
model {
y1 ~ dbin(theta1,n1)
y2 ~ dbin(theta2,n2)
theta1 ~ dbeta(60, 140)
theta2 ~ dbeta(60, 140)
}
"
case3.spec<-textConnection(case3.string)
# Initialize the JAGS model, and do an adaptation run of length 1,000.
case3.jags <- jags.model(case3.spec,
                         data = medfly,
                         n.chains=2,
                         n.adapt=1000)
update(case3.jags, 1000)
case3_samples = try(coda.samples(case3.jags, c("theta1"), iter), TRUE)




### Case 4
## alpha/(alpha + beta) = 0.45, alpha + beta = 2
## alpha = 0.9, beta = 1.1
case4.string <-"
model {
y1 ~ dbin(theta1,n1)
y2 ~ dbin(theta2,n2)
theta1 ~ dbeta(0.9, 1.1)
theta2 ~ dbeta(0.9, 1.1)
}
"
case4.spec<-textConnection(case4.string)
# Initialize the JAGS model, and do an adaptation run of length 1,000.
case4.jags <- jags.model(case4.spec,
                         data = medfly,
                         n.chains=2,
                         n.adapt=1000)
update(case4.jags, 1000)
case4_samples = try(coda.samples(case4.jags, c("theta1"), iter), TRUE)


### Case 5
## alpha/(alpha + beta) = 0.45, alpha + beta = 20
## alpha = 9, beta = 11
case5.string <-"
model {
y1 ~ dbin(theta1,n1)
y2 ~ dbin(theta2,n2)
theta1 ~ dbeta(9, 11)
theta2 ~ dbeta(9, 11)
}
"
case5.spec<-textConnection(case5.string)
# Initialize the JAGS model, and do an adaptation run of length 1,000.
case5.jags <- jags.model(case5.spec,
                         data = medfly,
                         n.chains=2,
                         n.adapt=1000)
update(case5.jags, 1000)
case5_samples = try(coda.samples(case5.jags, c("theta1"), iter), TRUE)



### Case 6
## alpha/(alpha + beta) = 0.45, alpha + beta = 20
## alpha = 90, beta = 110

case6.string <-"
model {
y1 ~ dbin(theta1,n1)
y2 ~ dbin(theta2,n2)
theta1 ~ dbeta(90, 110)
theta2 ~ dbeta(90, 110)
}
"
case6.spec<-textConnection(case6.string)
# Initialize the JAGS model, and do an adaptation run of length 1,000.
case6.jags <- jags.model(case6.spec,
                         data = medfly,
                         n.chains=2,
                         n.adapt=1000)
update(case6.jags, 1000)
case6_samples = try(coda.samples(case6.jags, c("theta1"), iter), TRUE)

library(MCMCvis)
require (reshape)
require (ggplot2)

priors <- cbind(MCMCchains(case0_samples, params = 'theta1'),
                MCMCchains(case1_samples, params = 'theta1'),
                MCMCchains(case2_samples, params = 'theta1'),
                MCMCchains(case3_samples, params = 'theta1'),
                MCMCchains(case6_samples, params = 'theta1'))
colnames(priors) <- c("case0","case1", "case2", "case3", "case6")


long = melt(priors)

p <- ggplot(aes(x=value, colour=X2), data=long)
p + geom_density() +
  scale_colour_manual(name = "Priors",
                      labels = c("alpha=1, beta=1", 
                                 "alpha=0.6, beta=1.4", 
                                 "alpha=6, beta=14", 
                                 "alpha=60, beta=140",
                                 "alpha=90, beta=110"),
                      values = c("black", "pink", "orange", "red", "blue")) 


table <- rbind(cbind(summary(case0_samples)$statistics[1], HPDinterval(case0_samples, prob = 0.95)[[1]]),
               cbind(summary(case1_samples)$statistics[1], HPDinterval(case1_samples, prob = 0.95)[[1]]),
               cbind(summary(case2_samples)$statistics[1], HPDinterval(case2_samples, prob = 0.95)[[1]]),
               cbind(summary(case3_samples)$statistics[1], HPDinterval(case3_samples, prob = 0.95)[[1]]),
               cbind(summary(case6_samples)$statistics[1], HPDinterval(case6_samples, prob = 0.95)[[1]]))

rownames(table) <- c("alpha=1, beta=1", 
                     "alpha=0.6, beta=1.4", 
                     "alpha=6, beta=14", 
                     "alpha=60, beta=140",
                     "alpha=90, beta=110")
colnames(table) <- c("posterior_mean", "HPD_low", "HPD_high")

write.csv(table, "q1_table.csv")

save.image("q1.Rdata")
