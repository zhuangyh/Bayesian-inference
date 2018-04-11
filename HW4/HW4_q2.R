### Bayesian Inference: HW4 -2

library(foreign)
library(rjags)
library(ggplot2)
library(tableone)
library(tidyverse)
library(magrittr)
library(cowplot)
library(ggthemes)
library(haven)

# http://www.stat.columbia.edu/~gelman/arm/examples/arsenic/
all <-read.dta("all.dta", convert.factors=F)
wells <- read.table("wells.dat")

### Village was missing
### Village information was retrevied (Relevant coded shared by Cuining Liu)

wellsdat_gelmanARM <-
  read_csv("fulldata1.csv") %>%
  select(UnicefWellId, VillageCode) %>%
  mutate(UnicefWellId = as.numeric(UnicefWellId))

# assume "ID" = "UnicefWellId"
village_data <-
  all %>%
  mutate(id = as.numeric(id)) %>%
  left_join(., wellsdat_gelmanARM, 
            by = c("id" = "UnicefWellId")) %>%
  filter(!duplicated(id))

village_data <-
  village_data %>%
  select(id, distnearest, as, VillageCode) %>%
  filter(complete.cases(.))

village_data$distnearest=village_data$as=NULL

# Retreive well id from "all data set"
missing <- is.na (all[,"func"] + all[,"as"] + all[,"distnearest"] + all[,"assn"] + all[,"ed"])
# Include only the wells that are functioning (func==1) and "unsafe" (as>50)
keep <- all[,"func"]==1 & all[,"as"]>50
wells_sub <- (all[!missing & keep,])

# Give convenient names to the variables
id <- wells_sub$id
switch <- wells_sub$switch
arsenic <- wells_sub$as/100
dist <- wells_sub$distnearest
logdist <- log(dist)
assoc <- ifelse (wells_sub$assn>0,1,0)
educ <- wells_sub$ed

newwells <- cbind (id, switch, arsenic, dist, logdist, assoc, educ)
newwells <- merge(newwells,village_data,by="id")
head(newwells)
wells_village <- newwells
wells_village$logdist=wells_village$distnearest= wells_village$as <- NULL
write.table(wells_village, "wells_village.dat")


#Create a variable list which we want in Table 1
listVars <- c("arsenic", "dist", "assoc","educ")

#Define categorical variables
catVars <- c("assoc")

table7.1 <- CreateTableOne(vars = listVars, data = newwells, factorVars = catVars, strata = "switch")
table7.1

table7.1 <- print(table7.1, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
## Save to a CSV file
write.csv(table7.1, "table7.1.csv")

par(mfrow=c(1,2))
boxplot(arsenic~switch,data=newwells, main="Arsenic level of wells", 
        xlab="switch", ylab="arsenic", col=c("gold","darkgreen"))

boxplot(logdist~switch,data=newwells, main="Log of well distance", 
        xlab="switch", ylab="Distance", col=c("gold","darkgreen"))
dev.off()

ggplot(newwells, aes(x=arsenic, y=logdist, color= as.factor(switch))) +
  geom_jitter(width = 0.40, height = 0.40)+
  #geom_smooth(method = "lm", se = FALSE)+
  ggtitle("")



head(newwells)
library(lme4)

# random intercept
logModel <- glmer(switch ~ arsenic + logdist + (1 | VillageCode), data = newwells, family = binomial, control = glmerControl(optimizer = "bobyqa"),
                       nAGQ = 1)
summary(logModel)

# random intercept and slope
logModel_slope <- glmer(switch ~ arsenic + logdist + (1+ logdist +arsenic| VillageCode), data = newwells, family = binomial, control = glmerControl(optimizer = "bobyqa"),
                  nAGQ = 1)
summary(logModel_slope)

y <- newwells$switch
n <- length(y)
arsenic <- newwells$arsenic
dist <- newwells$dist

# get village index variable
village.name <- as.vector(newwells$VillageCode)
uniq <- unique(village.name)
uniq
J <- length(uniq)
village <- rep (NA, J)
for (i in 1:J){
  village[village.name==uniq[i]] <- i
}


# no predictors
ybarbar = mean(y)
length(y)
length(village)
sample.size <- as.vector (table (village))
sample.size.jittered <- sample.size*exp (runif (J, -.1, .1))
cty.mns = tapply(y,village,mean)
cty.vars = tapply(y,village,var)
cty.sds = mean(sqrt(cty.vars[!is.na(cty.vars)]))/sqrt(sample.size)
cty.sds.sep = sqrt(tapply(y,village,var)/sample.size)


# varying-intercept and slope model no correlation in Random effects
## Code for calling the "modeling the correlation" Bugs code

wells.data <- list(arsenic=as.numeric(arsenic),
                  logdist=as.numeric(logdist),
                  y=as.numeric(y), 
                  village=as.numeric(village), 
                  n=as.numeric(n), 
                  J=as.numeric(J))

str(wells.data)

### Random intercept 
lg_varyingint_string <-"
    model {
    for (i in 1:n){
      y[i] ~ dbern(p[i])
      logit(p[i]) = a[village[i]] + b_arsenic*arsenic[i] + b_logdist*logdist[i]
      }
    
    for (j in 1:J){
      a[j] ~ dnorm (a.hat[j], tau.a)
      a.hat[j] <- mu.a
       }
    mu.a ~ dnorm (0, .0001)
    b_arsenic ~ dnorm (0, tau.b_arsenic)
    b_logdist ~ dnorm (0, tau.b_logdist)
    
    sigma.a <- tau.a^(-0.5)
    sigma.b_arsenic <- tau.b_arsenic^(-0.5)
    sigma.b_logdist <- tau.b_logdist^(-0.5)

    tau.a ~ dgamma(0.001, 0.001)
    tau.b_arsenic ~ dgamma(0.001, 0.001)
    tau.b_logdist ~ dgamma(0.001, 0.001)
    } 
    "


lg_varyingint<-textConnection(lg_varyingint_string)

lg_varyingint.fit <- jags.model(lg_varyingint, data=wells.data, n.chains=1, n.adapt=5000)
update(lg_varyingint.fit, 5000)

params <- c("a", "b_arsenic", 
            "mu.a", "b_arsenic", "b_logdist", 
            "sigma.a","sigma.b_arsenic","sigma.b_logdist" ) #set parameters you want to extract from model
#samps1 <- coda.samples(jags1.fit, params, n.iter=2000)
lg_varyingint.samps <- coda.samples(lg_varyingint.fit, params, n.iter=100000)

#summary(window(lg_varyingintslo.samps, start=1001)) #see summary of params

lg_varyingint_csim = as.mcmc(do.call(rbind, lg_varyingint.samps))
lg_varyingint_csim <- lg_varyingint_csim[,c( "mu.a", "b_arsenic", "b_logdist", 
                                                   "sigma.a","sigma.b_arsenic", "sigma.b_logdist")]

summary(lg_varyingint_csim)
## convergence diagnostics

lg_varyingint_csim <- lg_varyingint_csim[,c( "mu.a", "b_arsenic", "b_logdist", 
                                             "sigma.a","sigma.b_arsenic", "sigma.b_logdist")]


lg_varyingint_csim1 <- lg_varyingint_csim[,c( "mu.a", "b_arsenic", "b_logdist")]

lg_varyingint_csim2 <- lg_varyingint_csim[,c("sigma.a","sigma.b_arsenic", "sigma.b_logdist")]

plot(lg_varyingint_csim1, ask=TRUE)
plot(lg_varyingint_csim2, ask=TRUE)

interv <- HPDinterval(lg_varyingint_csim)
lg_int_summary <- cbind(summary(lg_varyingint_csim)$statistics[, 1:2], interv)
lg_int_summary
write.csv(lg_int_summary, "lg_int_summary.csv")


### Random intercept with uniform prior
unif_varyingint_string <-"
    model {
    for (i in 1:n){
      y[i] ~ dbern(p[i])
      logit(p[i]) = a[village[i]] + b_arsenic*arsenic[i] + b_logdist*logdist[i]
      }
      
    for (j in 1:J){
      a[j] ~ dnorm (a.hat[j], tau.a)
      a.hat[j] <- mu.a
      }
    mu.a ~ dnorm (0, .0001)
    b_arsenic ~ dnorm (0, tau.b_arsenic)
    b_logdist ~ dnorm (0, tau.b_logdist)
    
    tau.a <- pow(sigma.a, -2)
    tau.b_arsenic <- pow(sigma.b_arsenic, -2)
    tau.b_logdist <- pow(sigma.b_logdist, -2)

    sigma.a ~ dunif (0, 100)
    sigma.b_arsenic ~ dunif (0, 100)
    sigma.b_logdist ~ dunif (0, 100)
    } 
    "


unif_varyingint<-textConnection(unif_varyingint_string)

unif_varyingint.fit <- jags.model(unif_varyingint, data=wells.data, n.chains=1, n.adapt=5000)
update(unif_varyingint.fit, 5000)

params <- c("a", "b_arsenic", 
            "mu.a", "b_arsenic", "b_logdist", 
            "sigma.a","sigma.b_arsenic","sigma.b_logdist" ) #set parameters you want to extract from model
#samps1 <- coda.samples(jags1.fit, params, n.iter=2000)
unif_varyingint.samps <- coda.samples(unif_varyingint.fit, params, n.iter=100000)

#summary(window(lg_varyingintslo.samps, start=1001)) #see summary of params

unif_varyingint_csim = as.mcmc(do.call(rbind, unif_varyingint.samps))
unif_varyingint_csim <- unif_varyingint_csim[,c( "mu.a", "b_arsenic", "b_logdist", 
                                             "sigma.a","sigma.b_arsenic", "sigma.b_logdist")]

unif_varyingint_csim1 <- unif_varyingint_csim[,c( "mu.a", "b_arsenic", "b_logdist")]

unif_varyingint_csim2 <- unif_varyingint_csim[,c( "sigma.a","sigma.b_arsenic", "sigma.b_logdist")]


summary(unif_varyingint_csim)
## convergence diagnostics
plot(unif_varyingint_csim1, ask=TRUE)
plot(unif_varyingint_csim2, ask=TRUE)

interv <- HPDinterval(unif_varyingint_csim)
unif_int_summary <- cbind(summary(unif_varyingint_csim)$statistics[, 1:2], interv)
unif_int_summary
write.csv(unif_int_summary, "unif_int_summary.csv")



### Random intercept and slope
lg_varyingintslo_string <-"
      model {
      for (i in 1:n){
          y[i] ~ dbern(p[i])
          logit(p[i]) = a[village[i]] + b_arsenic[village[i]]*arsenic[i] + b_logdist[village[i]]*logdist[i]
      }

      for (j in 1:J){
          a[j] ~ dnorm (a.hat[j], tau.a)
          b_arsenic[j] ~ dnorm (b_arsenic.hat[j], tau.b_arsenic)
          b_logdist[j] ~ dnorm (b_logdist.hat[j], tau.b_logdist)
        
          a.hat[j] <- mu.a
          b_arsenic.hat[j] <- mu.b_arsenic
          b_logdist.hat[j] <- mu.b_logdist
        }
      mu.a ~ dnorm (0, .0001)
      mu.b_arsenic ~ dnorm (0, .0001)
      mu.b_logdist ~ dnorm (0, .0001)

      sigma.a <- tau.a^(-0.5)
      sigma.b_arsenic <- tau.b_arsenic^(-0.5)
      sigma.b_logdist <- tau.b_logdist^(-0.5)
      
      tau.a ~ dgamma(0.001, 0.001)
      tau.b_arsenic ~ dgamma(0.001, 0.001)
      tau.b_logdist ~ dgamma(0.001, 0.001)
      } 
      "


lg_varyingintslo<-textConnection(lg_varyingintslo_string)

lg_varyingintslo.fit <- jags.model(lg_varyingintslo, data=wells.data, n.chains=1, n.adapt=5000)
update(lg_varyingintslo.fit, 5000)

params <- c("a", "b_arsenic", 
            "mu.a", "mu.b_arsenic", "mu.b_logdist", 
            "sigma.a","sigma.b_arsenic","sigma.b_logdist" ) #set parameters you want to extract from model
#samps1 <- coda.samples(jags1.fit, params, n.iter=2000)
lg_varyingintslo.samps <- coda.samples(lg_varyingintslo.fit, params, n.iter=100000)

#summary(window(lg_varyingintslo.samps, start=1001)) #see summary of params

lg_varyingintslo_csim = as.mcmc(do.call(rbind, lg_varyingintslo.samps))
lg_varyingintslo_csim <- lg_varyingintslo_csim[,c( "mu.a", "mu.b_arsenic", "mu.b_logdist", 
                           "sigma.a","sigma.b_arsenic", "sigma.b_logdist")]

summary(lg_varyingintslo_csim)
## convergence diagnostics
lg_varyingintslo_csim1 <- lg_varyingintslo_csim[,c( "mu.a", "mu.b_arsenic", "mu.b_logdist")]

lg_varyingintslo_csim2 <- lg_varyingintslo_csim[,c("sigma.a","sigma.b_arsenic", "sigma.b_logdist")]

plot(lg_varyingintslo_csim1, ask=TRUE)
plot(lg_varyingintslo_csim2, ask=TRUE)


autocorr.diag(lg_varyingintslo_csim)
effectiveSize(lg_varyingintslo_csim)

interv <- HPDinterval(lg_varyingintslo_csim)
lg_intslo_summary <- cbind(summary(lg_varyingintslo_csim)$statistics[, 1:2], interv)
lg_intslo_summary
write.csv(lg_intslo_summary, "lg_intslo_summary.csv")



### Random intercept/slope with uniform prior
unif_varyingintslo_string <-"
    model {
    for (i in 1:n){
      y[i] ~ dbern(p[i])
      logit(p[i]) = a[village[i]] + b_arsenic[village[i]]*arsenic[i] + b_logdist[village[i]]*logdist[i]
      }
    
    for (j in 1:J){
      a[j] ~ dnorm (a.hat[j], tau.a)
      b_arsenic[j] ~ dnorm (b_arsenic.hat[j], tau.b_arsenic)
      b_logdist[j] ~ dnorm (b_logdist.hat[j], tau.b_logdist)
      
      a.hat[j] <- mu.a
      b_arsenic.hat[j] <- mu.b_arsenic
      b_logdist.hat[j] <- mu.b_logdist
      }
    mu.a ~ dnorm (0, .0001)
    mu.b_arsenic ~ dnorm (0, .0001)
    mu.b_logdist ~ dnorm (0, .0001)
    
    tau.a <- pow(sigma.a, -2)
    tau.b_arsenic <- pow(sigma.b_arsenic, -2)
    tau.b_logdist <- pow(sigma.b_logdist, -2)

    sigma.a ~ dunif (0, 100)
    sigma.b_arsenic ~ dunif (0, 100)
    sigma.b_logdist ~ dunif (0, 100)
    } 
    "


unif_varyingintslo<-textConnection(unif_varyingintslo_string)

unif_varyingintslo.fit <- jags.model(unif_varyingintslo, data=wells.data, n.chains=1, n.adapt=5000)
update(unif_varyingintslo.fit, 5000)

params <- c("a", "b_arsenic", 
            "mu.a", "mu.b_arsenic", "mu.b_logdist", 
            "sigma.a","sigma.b_arsenic","sigma.b_logdist" ) #set parameters you want to extract from model
#samps1 <- coda.samples(jags1.fit, params, n.iter=2000)
unif_varyingintslo.samps <- coda.samples(unif_varyingintslo.fit, params, n.iter=100000)

#summary(window(lg_varyingintslo.samps, start=1001)) #see summary of params

unif_varyingintslo_csim = as.mcmc(do.call(rbind, unif_varyingintslo.samps))
unif_varyingintslo_csim <- unif_varyingintslo_csim[,c( "mu.a", "mu.b_arsenic", "mu.b_logdist", 
                                                   "sigma.a","sigma.b_arsenic", "sigma.b_logdist")]

unif_varyingintslo_csim1 <- unif_varyingintslo_csim[,c( "mu.a", "mu.b_arsenic", "mu.b_logdist")]
unif_varyingintslo_csim2 <- unif_varyingintslo_csim[,c( "sigma.a","sigma.b_arsenic", "sigma.b_logdist")]



summary(unif_varyingintslo_csim)
## convergence diagnostics
plot(unif_varyingintslo_csim1, ask=TRUE)
plot(unif_varyingintslo_csim2, ask=TRUE)

autocorr.diag(unif_varyingintslo_csim)
effectiveSize(unif_varyingintslo_csim)

interv <- HPDinterval(unif_varyingintslo_csim)
unif_intslo_summary <- cbind(summary(unif_varyingintslo_csim)$statistics[, 1:2], interv)
unif_intslo_summary
write.csv(unif_intslo_summary, "unif_intslo_summary.csv")

save.image("HW4-2.RData")

load("HW4-2.RData")













