setwd("~/Dropbox/Biostatics_class/BIOS7717/HW3")

asthma <- read.table("fevasthma.dat")
colnames(asthma) <- c("id", "FEV1", "height", "age", "male", "asthma", "no2high" )
head(asthma)
tail(asthma)

library(tableone)
#Create a variable list which we want in Table 1
listVars <- c("FEV1", "male", "age","height", "no2high")

#Define categorical variables
catVars <- c("asthma","male", "no2high")

table5.1 <- CreateTableOne(vars = listVars, data = asthma, factorVars = catVars, strata = "asthma")
table5.1

table5.1 <- print(table5.1, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
## Save to a CSV file
write.csv(table5.1, "table5.1.csv")


asthma_sub <- asthma[ !duplicated(asthma$id), ]

table5.2 <- CreateTableOne(vars = listVars, data = asthma_sub, factorVars = catVars, strata = "asthma")
table5.2

table5.2 <- print(table5.2, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
## Save to a CSV file
write.csv(table5.2, "table5.2.csv")


head(asthma_sub)
### Refernce: MLE estimation

asthma_sub$sc_FEV1 <- scale(asthma_sub$FEV1)

freq.logit_mod_sc <- glm(asthma ~ sc_FEV1 + male ,
                         family=binomial(link='logit'),data=asthma_sub)
summary(freq.logit_mod_sc)

freq.logit_mod <- glm(asthma ~ FEV1 + male ,
                         family=binomial(link='logit'),data=asthma_sub)
summary(freq.logit_mod)



head(asthma)
X<-as.matrix(cbind(1,asthma_sub[,c(8, 5)])) # model matrix
Y<-matrix(asthma_sub$asthma, ncol=1) # outcome vector
p<-ncol(X)
head(X)

#library(LaplacesDemon)
library(invgamma)
library(MASS)

set.seed(10)

p<-ncol(X)


### sample from conditional posterior of tau
new_post_tau<-function(beta, alpha, gamma, lambda, p){
  post_alpha<-alpha + p/2
  post_gamma<-gamma + .5*t(beta - lambda)%*%(beta - lambda)
  draw_tau<-invgamma::rinvgamma(n = 1, shape = post_alpha, rate = post_gamma)
  return(draw_tau)
}

new_post_beta<-function(beta, tau, lambda, X, Y){
  # calculate likelihood
  lik<-0
  for(i in 1:length(Y)){
    xb <- X[i,] %*% beta
    xb<-ifelse(xb>10, 10, ifelse(xb< (-10),-10, xb))
    p_i<-invlogit(xb)
    lik<-lik + Y[i]*log(p_i) + (1 - Y[i])*log(1 - p_i)
  }
  
  # calculate prior 
  pr <- -.5 * (1/tau)*( t(beta - lambda)%*%(beta - lambda) )
  
  log_cond_post <- lik + pr
  return(log_cond_post)
}


new_post_beta(beta=c(0.001, 0.001, 0.001), tau=2, lambda=c(0, 0, 0), X, Y)

# Metropolis Hastings algorithm
beta_mh<-function(beta_0, tau, lambda, X, Y, mh_trials,j_sig){
  
  for(i in 1:mh_trials){
    # draw from proposal distribution
    beta_c <- mvrnorm(1,beta_0,Sigma = diag(j_sig))
    
    # calculate ratio of conditional posterior densities
    r_num <- new_post_beta(beta_c, tau, lambda, X, Y )
    r_denom <- new_post_beta(beta_0, tau, lambda, X, Y )
    
    print (r_num - r_denom)
    a <- exp(r_num - r_denom)
  
    a <-min(a,1)
    
    # accept or reject proposal
    accept<-0
    if(a >= 1){ 
      beta_0<-beta_c 
      accept<-1
    }else{ 
      if(rbinom(1,1,a)==1){ 
        beta_0<-beta_c
        accept<-1
      }
    }
  }
  
  return(c(new_beta=beta_0, accept=accept)  )
}

### Gibbs Sampler
# hyperparameter values for tau
alpha<-5
gamma<-2

# hyperparameter values for betas
lambda<-c(0,0,0)
tau<-0.01 # initialize 

# shell for storing results
gibbs_iter<-2000 + 1
gibbs_res<-matrix(nrow=gibbs_iter, ncol=p+2)

# initialize 
gibbs_res[1,1:p]<-c(0.1, 0.01, 0.15)
head(gibbs_res)

for(i in 2:gibbs_iter){
  # sample from posterior of tau
  gibbs_res[i,p+1] <- new_post_tau(gibbs_res[i-1,1:p], 
                                     alpha, gamma, lambda, p)
  # sample from posterior of beta vector ( using MH )
  mh_draw <- beta_mh(gibbs_res[i-1,1:p], gibbs_res[i,p+1], 
                                lambda, X, Y, mh_trials=5, j_sig=c(0.001,0.001,0.001))
  # store results
  gibbs_res[i,1:p] <- mh_draw[1:p]
  gibbs_res[i,p+2] <- mh_draw[p+1]
}

tail(gibbs_res)
head(gibbs_res)
length(gibbs_res)
sum(gibbs_res[2:2001, 5])


################################################################################
### 3 - Plot Results
################################################################################

par(mfrow=c(2,2))
plot(gibbs_res[,1][1:2001],type='l',xlab='MCMC Iterations',ylab=c('Coefficient Draw'),
     main='Intercept')
plot(gibbs_res[,2][1:2001],type='l',xlab='MCMC Iterations',ylab=c('Coefficient Draw'),
     main='FEV1')
plot(gibbs_res[,3][1:2001],type='l',xlab='MCMC Iterations',ylab=c('Coefficient Draw'),
     main='male')

# calculate posterior means and credible intervals
post_burn_trim<-gibbs_res[seq(1000, gibbs_iter,100),]
colMeans(post_burn_trim)
apply(post_burn_trim, 2, quantile, p=c(.025,.975))


ddexp = function(x, mu, tau) {
  0.5*tau*exp(-tau*abs(x-mu)) 
}

mod1_string = "model {
    for (i in 1:length(y)) {
      y[i] ~ dbern(p[i])
      logit(p[i]) = b[1] + b[2]*sc_FEV1[i] + b[3]*male[i] 
    }
      for (j in 1:3) {
      b[j] ~ dnorm(0, tau[j]) 
      tau[j] ~ dgamma(0.001, 0.001)
      }
      } "

set.seed(92)

X.df <-as.data.frame(X)
head(X.df)
data_jags = list(y=asthma_sub$asthma, sc_FEV1=X.df$sc_FEV1, male=X.df[, "male"])
str(data_jags)
str(X)
params = c("b")

mod1 = jags.model(textConnection(mod1_string), data=data_jags, n.chains=3)
update(mod1, 1e3)

mod1_sim = coda.samples(model=mod1,
                        variable.names=params,
                        n.iter=5e3)

summary(mod1_sim)
mod1_csim = as.mcmc(do.call(rbind, mod1_sim))

summary(mod1_csim)


beta_summary <- cbind(summary(mod1_csim)$statistics, interv)
write.csv(beta_summary, "beta_summary.csv")


## convergence diagnostics
plot(mod1_sim, ask=TRUE)

gelman.diag(mod1_sim)
autocorr.diag(mod1_sim)
autocorr.plot(mod1_sim)
effectiveSize(mod1_sim)

## calculate DIC
dic1 = dic.samples(mod1, n.iter=1e3)






