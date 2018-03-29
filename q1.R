
setwd("~/Dropbox/Biostatics_class/BIOS7717/HW3")
library(rjags)
library(mcmcse)

medfly <- list(n1=134,y1=54,n2=400,y2=224)

### Reference prior
model1.string <-"
  model {
    y1 ~ dbin(theta1,n1)
    y2 ~ dbin(theta2,n2)
    theta1 ~ dbeta(1,1)
    theta2 ~ dbeta(1,1)
    odds1 <- theta1/(1-theta1)
    odds2 <- theta2/(1-theta2)
    riskDifference <- theta2-theta1
    relativeRisk <- theta2/theta1
    oddsRatio <- odds2/odds1
    prob_RD <- step(riskDifference)
  }
  "

model1.spec<-textConnection(model1.string)

# Initialize the JAGS model, and do an adaptation run of length 1,000.
ref.jags <- jags.model(model1.spec,
                           data = medfly,
                           n.chains=2,
                           n.adapt=1000)
update(ref.jags, 1000)
bayes.model.params <- c("theta1", "theta2", "riskDifference", "relativeRisk", "oddsRatio", "prob_RD")

# Now draw 10,000 samples from the posterior.
iter = 10000
ref_coda_samples = try(coda.samples(ref.jags, bayes.model.params, iter), TRUE)
# summary of modle with logit link 
summary(ref_coda_samples)
# HPD of betas
HPDinterval(ref_coda_samples, prob = 0.95)

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

### Data augment priors 2: Beta(100, 100)
model3.string <-"
  model {
    y1 ~ dbin(theta1,n1)
    y2 ~ dbin(theta2,n2)
    theta1 ~ dbeta(100,100)
    theta2 ~ dbeta(100,100)
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
dataAug2.jags <- jags.model(model3.spec,
                           data = medfly,
                           n.chains=2,
                           n.adapt=1000)
update(dataAug2.jags, 1000)

dataAug2_coda_samples = try(coda.samples(dataAug2.jags, bayes.model.params, iter), TRUE)

dataAug2_summary <-cbind(summary(dataAug_coda_samples)$statistics[, c(1, 2)], HPDinterval(dataAug_coda_samples, prob = 0.95)[[1]])
dataAug2_summary <- dataAug_summary[c("theta1", "theta2", "riskDifference", "relativeRisk", "oddsRatio", "prob_RD"), ] 
dataAug2_summary
write.csv(dataAug2_summary, file = "q1_dataAug2_summary.csv")





### Expert priors
### theta1 ~ beta(2.07, 4.83)
qbeta(0.95, 2.07, 4.83)
### theta2 ~ beta(1.53, 1.53)
qbeta(0.95, 1.53, 1.53)


model4.string <-"
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

model4.spec<-textConnection(model3.string)
# Initialize the JAGS model, and do an adaptation run of length 1,000.
expert.jags <- jags.model(model4.spec,
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
  ggtitle("Impact of priors on theta1 ")+
  scale_colour_manual(name = "Priors",
                      labels = c("reference prior: a=1, b=1", 
                                 "a=0.6, =1.4", 
                                 "a=6, b=14", 
                                 "a=60, b=140",
                                 "a=90, b=110"),
                      values = c("black", "pink", "orange", "red", "blue")) 


table <- rbind(cbind(summary(case0_samples)$statistics[1], HPDinterval(case0_samples, prob = 0.95)[[1]]),
               cbind(summary(case1_samples)$statistics[1], HPDinterval(case1_samples, prob = 0.95)[[1]]),
               cbind(summary(case2_samples)$statistics[1], HPDinterval(case2_samples, prob = 0.95)[[1]]),
               cbind(summary(case3_samples)$statistics[1], HPDinterval(case3_samples, prob = 0.95)[[1]]),
               cbind(summary(case6_samples)$statistics[1], HPDinterval(case6_samples, prob = 0.95)[[1]]))

rownames(table) <- c("a=1, b=1", 
                     "a=0.6, b=1.4", 
                     "a=6, b=14", 
                     "a=60, b=140",
                     "a=90, b=110")
colnames(table) <- c("posterior_mean", "HPD_low", "HPD_high")

write.csv(table, "q1_table.csv")

save.image("q1.Rdata")
