### BIOS 7717
### HW4 Question 1

library(foreign)
library(rjags)

pol <- read.dta ("nes5200_processed_voters_realideo.dta")
pol <- pol[pol$year==2000, ]

pol_sub <- pol[, c("partyid7", "ideo7", "income", "state")]
pol_sub <- na.omit(pol_sub)
head(pol_sub)
dim(pol_sub)

# pol_sub1 <- pol[, c("partyid7", "ideo7", "income",  "state")]
# pol_sub1 <- na.omit(pol_sub1)
# dim(pol_sub1)

pol_sub$partyid7 <- as.numeric(substr(pol_sub$partyid7, 0, 1))
pol_sub$ideo7 <- as.numeric(substr(pol_sub$ideo7, 0, 1))
# pol_sub$race.adj <- ifelse (pol_sub[,"white"]== 1, 0, 0.5)
# pol_sub$race.adj <- ifelse (pol_sub[,"black"]== 1, 1, pol_sub[,"race.adj"])
# pol_sub$educ1 <- as.numeric(substr(pol_sub$educ1, 0, 1))
pol_sub$income <- as.factor(substr(pol_sub$income, 0, 1))
#pol_sub$age.discrete <- as.numeric (cut (pol_sub$age, c(0,29.5, 44.5, 64.5, 200)))
#pol_sub$race = pol_sub$white = pol_sub$black = pol_sub$age = NULL
head(pol_sub)

pol_sub$state <- as.factor(pol_sub$state)
pol_sub$income <- as.numeric(pol_sub$income)
unique(pol_sub$income)
print(levels(pol_sub$income))

library(ggplot2)
# Change point shapes and colors
dim(pol_sub)
set.seed(10)

ggplot(pol_sub, aes(x=state, y=partyid7)) +
  geom_jitter(width = 0.40, height = 0.40)+
  geom_boxplot(outlier.colour="black", outlier.shape=16,
              outlier.size=2, notch=FALSE)+
  ggtitle("Box plot of party identification among states")


ggplot(pol_sub, aes(x=ideo7, y=partyid7, color= income)) +
  geom_jitter(width = 0.40, height = 0.40)+
  geom_smooth(method = "lm", se = FALSE)+
  ggtitle("Scatter plot of party identification vs. ideology")+
  scale_color_manual(values=c("red","orange", "black", "blue", "green"))


ggplot(pol_sub, aes(x=ideo7, y=partyid7, color= income)) +
  geom_jitter(width = 0.40, height = 0.40)+
  #geom_smooth(method = "lm", se = FALSE)+
  ggtitle("Associations between party-ID and ideology in different income categroies")+
  scale_color_manual(values=c("red","orange", "black", "blue", "green"))


ggplot(pol_sub, aes(x=ideo7, y=partyid7, color= income)) +
  geom_jitter(width = 0.40, height = 0.40)+
  geom_smooth(method = "lm", se = FALSE)+
  ggtitle("Associations between party-ID and ideology in different income categroies")+
  scale_color_manual(values=c("red","orange", "black", "blue", "green"))



#Create a variable list which we want in Table 1
listVars <- c("partyid7", "ideo7", "income")

#Define categorical variables
catVars <- c("state", "income")

#Total Population
library(tableone)
table1 <- CreateTableOne(vars = listVars, data = pol_sub, factorVars = catVars)

table1 <- print(table1,  exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
## Save to a CSV file
table1
write.csv(table1, file = "HW4_table1.csv")



library(lme4)
MLeModel <- lmer(partyid7 ~ ideo7 + (1 + ideo7 | state), 
                  data = pol_sub)
summary(MLeModel)

MLeModel_inc <- lmer(partyid7 ~ ideo7 + income+ (1 + ideo7 + income | state), 
                 data = pol_sub)
summary(MLeModel_inc)

y <- pol_sub$partyid7
n <- length(y)
ideo7 <- pol_sub$ideo7
income <- pol_sub$income

# get state index variable
state.name <- as.vector(pol_sub$state)
uniq <- unique(state.name)
uniq
J <- length(uniq)
state <- rep (NA, J)
for (i in 1:J){
  state[state.name==uniq[i]] <- i
}

# no predictors
ybarbar = mean(y)

sample.size <- as.vector (table (state))
sample.size.jittered <- sample.size*exp (runif (J, -.1, .1))
cty.mns = tapply(y,state,mean)
cty.vars = tapply(y,state,var)
cty.sds = mean(sqrt(cty.vars[!is.na(cty.vars)]))/sqrt(sample.size)
cty.sds.sep = sqrt(tapply(y,state,var)/sample.size)


# varying-intercept and slope model no correlation in Random effects
## Code for calling the "modeling the correlation" Bugs code

pol.data <- list (ideo7=as.numeric(ideo7),
                  income=as.numeric(income),
                  y=as.numeric(y), 
                  state=as.numeric(state), 
                  n=as.numeric(n), 
                  J=as.numeric(J))
str(pol.data)


## Varying-intercept, varying-slope models ##
# Gamma priors on tau

gamma_varyingint_string <-"
   model {
      for (i in 1:n){
          y[i] ~ dnorm (y.hat[i], tau.y)
          y.hat[i] <- a[state[i]] + 
                      b_ideo7[state[i]]*ideo7[i]
          }
    tau.y ~ dgamma(0.001, 0.001)
    sigma.y  = tau.y^(-0.5)
    
    for (j in 1:J){
        a[j] ~ dnorm (a.hat[j], tau.a)
        b_ideo7[j] ~ dnorm (b_ideo7.hat[j], tau.b_ideo7)

        a.hat[j] <- mu.a
        b_ideo7.hat[j] <- mu.b_ideo7
        }
    mu.a ~ dnorm (0, .0001)
    mu.b_ideo7 ~ dnorm (0, .0001)

    tau.a ~ dgamma(0.001, 0.001)
    tau.b_ideo7 ~ dgamma(0.001, 0.001)
    
    sigma.a <- tau.a^(-0.5)
    sigma.b_ideo7 <- tau.b_ideo7^(-0.5)
} 
"
gamma_varyingint<-textConnection(gamma_varyingint_string)

gamma.fit <- jags.model(gamma_varyingint, data=pol.data, n.chains=1, n.adapt=5000)
update(gamma.fit, 5000)

params <- c("a", "b_ideo7", 
            "mu.a", "mu.b_ideo7", 
            "sigma.a","sigma.b_ideo7") #set parameters you want to extract from model
#samps1 <- coda.samples(jags1.fit, params, n.iter=2000)
gamma.samps <- coda.samples(gamma.fit, params, n.iter=100000)

#summary(window(samps1, start=1001)) #see summary of params

gamma_csim = as.mcmc(do.call(rbind, gamma.samps))
gamma_csim <- gamma_csim[,c( "mu.a", "mu.b_ideo7", 
                           "sigma.a","sigma.b_ideo7")]

gamma_csim1 <- gamma_csim[,c( "mu.a", "mu.b_ideo7")]
gamma_csim2 <- gamma_csim[,c( "sigma.a","sigma.b_ideo7")]

## convergence diagnostics
plot(gamma_csim1, ask=TRUE)
plot(gamma_csim2, ask=TRUE)


autocorr.diag(gamma_csim)
effectiveSize(gamma_csim)

interv <- HPDinterval(gamma_csim)
gamma_summary <- cbind(summary(gamma_csim)$statistics[, 1:2], interv)
gamma_summary
write.csv(gamma_summary, "hw4q1_gamma_summary.csv")



### Uniform prior for tau
unif_varyingint_string <-"
  model {
    for (i in 1:n){
      y[i] ~ dnorm (y.hat[i], tau.y)
      y.hat[i] <- a[state[i]] + b_ideo7[state[i]]*ideo7[i]
    }
    tau.y <- pow(sigma.y, -2)
    sigma.y ~ dunif (0, 100)
    
    for (j in 1:J){
      a[j] ~ dnorm (a.hat[j], tau.a)
      b_ideo7[j] ~ dnorm (b_ideo7.hat[j], tau.b_ideo7)
      
      a.hat[j] <- mu.a
      b_ideo7.hat[j] <- mu.b_ideo7
    }
    mu.a ~ dnorm (0, .0001)
    mu.b_ideo7 ~ dnorm (0, .0001)
    
    tau.a <- pow(sigma.a, -2)
    tau.b_ideo7 <- pow(sigma.b_ideo7, -2)
    
    sigma.a ~ dunif (0, 100)
    sigma.b_ideo7 ~ dunif (0, 100)
  } 
"

unif_varyingint<-textConnection(unif_varyingint_string)
unif.fit <- jags.model(unif_varyingint, data=pol.data, n.chains=1, n.adapt=5000)
update(unif.fit , 5000)

params <- c("a", "b_ideo7", 
            "mu.a", "mu.b_ideo7", 
            "sigma.a","sigma.b_ideo7") #set parameters you want to extract from model
unif.samps <- coda.samples(unif.fit, params, n.iter=100000)
unif_csim = as.mcmc(do.call(rbind, unif.samps))
unif_csim <- unif_csim[,c( "mu.a", "mu.b_ideo7", 
                           "sigma.a","sigma.b_ideo7")]

unif_csim1 <- unif_csim[,c("mu.a", "mu.b_ideo7")]
unif_csim2 <- unif_csim[,c("sigma.a","sigma.b_ideo7")]

summary(unif_csim)
## convergence diagnostics
plot(unif_csim1, ask=TRUE)
plot(unif_csim2, ask=TRUE)

autocorr.diag(unif_csim)
effectiveSize(unif_csim)

interv <- HPDinterval(unif_csim)
unif_summary <- cbind(summary(unif_csim)$statistics[, 1:2], interv)
unif_summary
write.csv(unif_summary, "hw4q1_unif_summary.csv")


######## Add income covariate
gamma_inc_string <-"
   model {
      for (i in 1:n){
        y[i] ~ dnorm (y.hat[i], tau.y)
        y.hat[i] <- a[state[i]] + b_ideo7[state[i]]*ideo7[i] + b_income[state[i]]*income[i]        
        }
      tau.y <- pow(sigma.y, -2)
      sigma.y ~ dunif (0, 100)
      
      for (j in 1:J){
      a[j] ~ dnorm (a.hat[j], tau.a)
      b_ideo7[j] ~ dnorm (b_ideo7.hat[j], tau.b_ideo7)
      b_income[j] ~ dnorm (b_income.hat[j], tau.b_income)
      
      a.hat[j] <- mu.a
      b_ideo7.hat[j] <- mu.b_ideo7
      b_income.hat[j] <- mu.b_income
      }
      mu.a ~ dnorm (0, .0001)
      mu.b_ideo7 ~ dnorm (0, .0001)
      mu.b_income ~ dnorm (0, .0001)
      
      tau.a <- pow(sigma.a, -2)
      tau.b_ideo7 <- pow(sigma.b_ideo7, -2)
      tau.b_income <- pow(sigma.b_income, -2)
      
      sigma.a ~ dunif (0, 100)
      sigma.b_ideo7 ~ dunif (0, 100)
      sigma.b_income ~ dunif (0, 100)
      } 
      "
gamma_inc <-textConnection(gamma_inc_string)

gamma_inc.fit <- jags.model(gamma_inc, data=pol.data, n.chains=1, n.adapt=5000)
update(gamma_inc.fit, 5000)

params_inc <- c("a", "b_ideo7", 
            "mu.a", "mu.b_ideo7", "mu.b_income", 
            "sigma.a","sigma.b_ideo7", "sigma.b_income") #set parameters you want to extract from model

gamma_inc.samps <- coda.samples(gamma_inc.fit, params_inc, n.iter=100000)

gamma_inc_csim = as.mcmc(do.call(rbind, gamma_inc.samps))
gamma_inc_csim <- gamma_inc_csim[,c( "mu.a", "mu.b_ideo7", "mu.b_income", 
                           "sigma.a","sigma.b_ideo7", "sigma.b_income")]

summary(gamma_inc_csim)
## convergence diagnostics
plot(gamma_inc_csim, ask=TRUE)

autocorr.diag(gamma_inc_csim)
effectiveSize(gamma_inc_csim)

interv <- HPDinterval(gamma_inc_csim)
gamma_inc_summary <- cbind(summary(gamma_inc_csim)$statistics[, 1:2], interv)
gamma_inc_summary
write.csv(gamma_inc_summary, "hw4q1_gamma_inc_summary.csv")

### Unif priors for tau
### Adjust for income

unif_inc_string <-"
  model {
  for (i in 1:n){
    y[i] ~ dnorm (y.hat[i], tau.y)
    y.hat[i] <- a[state[i]] + b_ideo7[state[i]]*ideo7[i] + b_income[state[i]]*income[i]
    }
  tau.y <- pow(sigma.y, -2)
  sigma.y ~ dunif (0, 100)
  
  for (j in 1:J){
    a[j] ~ dnorm (a.hat[j], tau.a)
    b_ideo7[j] ~ dnorm (b_ideo7.hat[j], tau.b_ideo7)
    b_income[j] ~ dnorm (b_income.hat[j], tau.b_income)
    
    a.hat[j] <- mu.a
    b_ideo7.hat[j] <- mu.b_ideo7
    b_income.hat[j] <- mu.b_income
    }
  mu.a ~ dnorm (0, .0001)
  mu.b_ideo7 ~ dnorm (0, .0001)
  mu.b_income ~ dnorm (0, .0001)
  
  tau.a <- pow(sigma.a, -2)
  tau.b_ideo7 <- pow(sigma.b_ideo7, -2)
  tau.b_income <- pow(sigma.b_income, -2)
  
  sigma.a ~ dunif (0, 100)
  sigma.b_ideo7 ~ dunif (0, 100)
  sigma.b_income ~ dunif (0, 100)
  } 
  "

unif_inc<-textConnection(unif_inc_string)
unif_inc.fit <- jags.model(unif_inc, data=pol.data, n.chains=1, n.adapt=5000)
update(unif_inc.fit, 5000)


unif_inc.samps <- coda.samples(unif_inc.fit, params_inc, n.iter=100000)
unif_inc_csim = as.mcmc(do.call(rbind, unif_inc.samps))
unif_inc_csim <- unif_inc_csim[,c( "mu.a", "mu.b_ideo7", "mu.b_income",
                           "sigma.a","sigma.b_ideo7", "sigma.b_income")]

summary(unif_inc_csim)
## convergence diagnostics
plot(unif_inc_csim, ask=TRUE)

autocorr.diag(unif_inc_csim)
effectiveSize(unif_inc_csim)

interv <- HPDinterval(unif_inc_csim)
unif_inc_summary <- cbind(summary(unif_inc_csim)$statistics[, 1:2], interv)
unif_inc_summary
write.csv(unif_inc_summary, "hw4q1_unif_inc_summary.csv")

save.image("HW4.1.RData")
load("HW4.1.RData")
