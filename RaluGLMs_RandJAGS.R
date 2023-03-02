library(tidyverse)

#, import data from csv file
spp.sites <- read.csv("spp.sites.csv")

## pull out just the frog data and covariates we might want
ralu.sites <- spp.sites %>%
  filter(species=="ralu") %>%
  select(site, date, community, julianday, species, tads, sals, tads, bubos, 
         ralus, PC1, PC2,beav, m30max, m30mean, positive, total, prev, permanent, 
         canopy_cover, elev, depth_mean_1m, log.area, ph_mean_1m, m7max, 
         t30precip, t7precip, log_chap_bubo, log_chap_ralu, log_chap_total, 
         mean_bod_cond, total_density, ralu_density, bubo_density)%>%
  mutate("prevalence" = prev/100)

### 2nd BEST RALU MODEL: ralu density, PC1, PC2, canopy cover, and hydroperiod (best doesn't include hydro)
ralu.prev.best2 <- glm(prevalence ~ ralu_density + permanent + PC2 + PC1 + canopy_cover, 
                      family = "binomial", data = ralu.sites, weights = total)


summary(ralu.prev.best2)

##R^2
cor.test(predict(ralu.prev.best2), ralu.sites$prevalence)


# pull out coefficients
coef <- ralu.prev.best2$coefficients
#(Intercept) ralu_density    permanent          PC2          PC1 canopy_cover 
#-0.47022868   3.51033958   0.43709159   0.45448826   0.25806481   0.01326909 



# predict values (untransformed):
pred.u <- coef['(Intercept)'] + 
  coef['ralu_density'] * ralu.sites$ralu_density +
  coef['permanent'] * ralu.sites$permanent + ### this is binary 1 or 0 # need example of non-binary - I code 3 levels as 1,2,3 - then use that to reference which value to use in a vector of values for that covariate
  coef['PC2'] * ralu.sites$PC2 +
  coef['PC1'] * ralu.sites$PC1 +
  coef['canopy_cover'] * ralu.sites$canopy_cover

#logit <- function(x){ log(x/(1-x)) }
inv.logit <- function(x){
  exp(x)/(1+exp(x))
}

# now do the inverse logit transform to get correct scale
pred <- inv.logit(pred.u)

plot(pred, predict(ralu.prev.best2))
plot(pred.u, predict(ralu.prev.best2)) ## Weird! The predict function doesn't back transform

# need to say predict(ralu.prev.best2, type="response") to get it to transform 


# now that I know predict() doesn't back transform, redo R^2 calculation from before
mod.cor <- cor.test(inv.logit(predict(ralu.prev.best2)), ralu.sites$prevalence)
mod.cor$estimate^2  # R squared = 0.375, slightly better



###############################################################################
# 
# Plot marginal effects and explore effect of beavers
#
###############################################################################


# Explore marginal effects of ralu density and hydroperiod  -----------------#
# keep other variables constant. PCs=0 (mean values), canopy cover = 50%
new.dat.ralu <- data.frame(ralu_density = rep(seq(0, 0.4, length=30),2) , 
                      permanent = c(rep(0, 30),rep(1,30)),
                      PC2 = rep(0, 60), PC1 = rep(0, 60),
                      canopy_cover = rep(50, 60), total = rep(10, 60) )

pred.ralu.hydro <- predict.glm(ralu.prev.best2, newdata = new.dat.ralu, 
                               type = "response", se.fit = TRUE)


# add predictions and se to dataframe
new.dat.ralu$permanent <- factor(new.dat.ralu$permanent)
new.dat.ralu$fit <- pred.ralu.hydro$fit
new.dat.ralu$se.fit <- pred.ralu.hydro$se.fit
new.dat.ralu$CI.lwr <- pred.ralu.hydro$fit - pred.ralu.hydro$se.fit*1.96
new.dat.ralu$CI.upr <- pred.ralu.hydro$fit + pred.ralu.hydro$se.fit*1.96


ggplot(new.dat.ralu, aes(x = ralu_density, y = fit, color = permanent )) +
  #geom_point() +
  geom_ribbon( aes(ymin = CI.lwr, ymax = CI.upr, fill = permanent, color = NULL), alpha = .15) +
  geom_line( aes(y = fit)) +
  ylab("Predicted Prevalence") 




# Explore marginal effects of canopy cover and hydroperiod  -----------------#
# keep other variables constant. PCs=0 (mean values), ralu density = 0.2
new.dat.canopy <- data.frame(ralu_density = rep(0.2, 60) , 
                           permanent = c(rep(0, 30),rep(1,30)),
                           PC2 = rep(0, 60), PC1 = rep(0, 60),
                           canopy_cover = rep(seq(0,100,length=30),2), total = rep(10, 60) )
pred.canopy.hydro <- predict.glm(ralu.prev.best2, newdata = new.dat.canopy, 
                               type = "response", se.fit = TRUE)


# add predictions and se to dataframe
new.dat.canopy$permanent <- factor(new.dat.canopy$permanent)
new.dat.canopy$fit <- pred.canopy.hydro$fit
new.dat.canopy$se.fit <- pred.canopy.hydro$se.fit
new.dat.canopy$CI.lwr <- pred.canopy.hydro$fit - pred.canopy.hydro$se.fit*1.96
new.dat.canopy$CI.upr <- pred.canopy.hydro$fit + pred.canopy.hydro$se.fit*1.96


ggplot(new.dat.canopy, aes(x = canopy_cover, y = fit, color = permanent )) +
  #geom_point() +
  geom_ribbon( aes(ymin = CI.lwr, ymax = CI.upr, fill = permanent, color = NULL), alpha = .15) +
  geom_line( aes(y = fit)) +
  ylab("Predicted Prevalence") 



####################### Explore effect of beavers
# Assume Beaver ponds are permanent and have canopy cover of 21 
# & non-beaver ponds are not permanent and have canopy cover of 42
# probably should make beaver ponds permanent or non-permanent in proportion to data
# assume PC1 and PC2 values of 0 (average)
# explore the effect of ralu_density

# How do beavers affect prevalence?  ------------------------------------------#
# keep other variables constant. PCs=0 (mean values), ralu density = 0.2
new.dat.beaver <- data.frame(ralu_density = rep(seq(0, 0.4, length=30),2) , 
                           permanent = rep(1, 60),
                           PC2 = rep(0, 60), PC1 = rep(0, 60),
                           canopy_cover = c(rep(21, 30), rep(42, 30)), total = rep(10, 60) )
pred.beaver <- predict.glm(ralu.prev.best2, newdata = new.dat.beaver, 
                                 type = "response", se.fit = TRUE)

new.dat.beaver$beaver <- c(rep("beaver pond",30), rep("non-beaver pond",30))
# add predictions and se to dataframe
new.dat.beaver$permanent <- factor(new.dat.beaver$permanent)
new.dat.beaver$fit <- pred.beaver$fit
new.dat.beaver$se.fit <- pred.beaver$se.fit
new.dat.beaver$CI.lwr <- pred.beaver$fit - pred.beaver$se.fit*1.96
new.dat.beaver$CI.upr <- pred.beaver$fit + pred.beaver$se.fit*1.96


ggplot(new.dat.beaver, aes(x = ralu_density, y = fit, color = beaver )) +
  geom_ribbon( aes(ymin = CI.lwr, ymax = CI.upr, fill = beaver, color = NULL), alpha = .15) +
  geom_line( aes(y = fit)) +
  ylab("Predicted Prevalence") 




###############################################################################
### Now Bayesian!!
###############################################################################

library(R2jags)
library(MCMCvis)
library(mcmcplots)

############################## Specify Model in JAGS code to a text file
sink("GLM_Raluprev.jags")
cat("
model {

# Priors
beta.0 ~ dunif(-5, 5)
beta.ralu ~ dunif(-10, 10)
beta.permanent ~ dunif(-5, 5)
beta.PC1 ~ dunif(-5, 5)
beta.PC2 ~ dunif(-5, 5)
beta.canopy ~ dunif(-5, 5)

# Likelihood: 
for (i in 1:n){

   logit(lambda[i]) <- beta.0 + 
      beta.ralu * ralu[i] + 
      beta.permanent * permanent[i] + 
      beta.PC1 * PC1[i] +
      beta.PC2 * PC2[i] +
      beta.canopy * canopy_cover[i]

   # Compare to data with random binomial errors  
   positive[i] ~ dbin(lambda[i], total[i])    #dbin wants prob first, then numb trails (unlike R's dbinom())
   
   } #i
}
",fill = TRUE)
sink()

############################## End specify model



##### Now R code to run the model

# Bundle data
jags.data <- list(positive = ralu.sites$positive, 
                 total = ralu.sites$total,
                 n = length(ralu.sites$positive), 
                 ralu = ralu.sites$ralu_density,
                 permanent = ralu.sites$permanent,
                 PC1 = ralu.sites$PC1,
                 PC2 = ralu.sites$PC2,
                 canopy_cover = ralu.sites$canopy_cover/100 # got an error (node inconsistent with parents without the division)
                 )

# Initial values
inits <- function() list(beta.0 = runif(1, -5, 5),
                         beta.ralu = runif(1, -5, 5),
                         beta.permanent = runif(1,-5, 5),
                         beta.PC1 = runif(1, -5, 5),
                         beta.PC2 = runif(1, -5, 5),
                         beta.canopy = runif(1, -5, 5))
                         


# Parameters monitored
params <- c("beta.0","beta.ralu","beta.permanent","beta.PC1","beta.PC2","beta.canopy","lambda")

# MCMC settings
ni <- 20000
nt <- 2
nb <- 5000
nc <- 3

# Call JAGS from R
out <- jags(data = jags.data, inits = inits, parameters.to.save = params, 
            model.file = "GLM_Raluprev.jags", n.chains = nc, n.thin = nt, 
            n.iter = ni, n.burnin = nb, working.directory = getwd())


####### Make some plots
# look at distributions, etc
mcmcplot(out,parms=params[1:6])

# plot of estimates
MCMCplot(out,params=params[1:6])

## add the GLM points estimates in red
points(coef[c(1:3,5,4)],6:2,col="white",bg="red",pch=22)
points(coef[6]*100,1,col="white",bg="red",pch=22) # need to mult canopy cover back to 100

## Add GLM Confidence intervals which are 1.96*SE
SE <- summary(ralu.prev.best2)$coefficients[,2]
lines(c(coef[1]+1.96*SE[1],coef[1]-1.96*SE[1]),c(6,6),col="red",lwd=3)
lines(c(coef[2]+1.96*SE[2],coef[2]-1.96*SE[2]),c(5,5),col="red",lwd=3)
lines(c(coef[3]+1.96*SE[3],coef[3]-1.96*SE[3]),c(4,4),col="red",lwd=3)
lines(c(coef[5]+1.96*SE[5],coef[5]-1.96*SE[5]),c(3,3),col="red",lwd=3)
lines(c(coef[4]+1.96*SE[4],coef[4]-1.96*SE[4]),c(2,2),col="red",lwd=3)
lines(c(coef[6]*100+1.96*SE[6]*100,coef[6]*100-1.96*SE[6]*100),c(1,1),col="red",lwd=3)






##############################################################################
## Now ask the Beaver question but based on the whole
## posterior parameter distributions (instead of point estimates)
##############################################################################

# first just look at one mean density of frogs
ralu1 <- mean(ralu.sites$ralu_density) #0.133461

# full posteriors of all parameters (so 12000 values for each parameter)
beta.0.post = out$BUGSoutput$sims.list$beta.0
beta.ralu.post = out$BUGSoutput$sims.list$beta.ralu
beta.permanent.post = out$BUGSoutput$sims.list$beta.permanent
beta.PC1.post = out$BUGSoutput$sims.list$beta.PC1
beta.PC2.post = out$BUGSoutput$sims.list$beta.PC2
beta.canopy.post = out$BUGSoutput$sims.list$beta.canopy


### Predict prevalence in beaver ponds at mean frog density of 0.133
beav.pred.post <- beta.0.post + 
  beta.ralu.post * ralu1 +
  beta.permanent.post + 
  #coef['PC2'] * ralu.sites$PC2 +  # assume average value of 0 
  #coef['PC1'] * ralu.sites$PC1 +  # assume average value of 0 
  beta.canopy.post * 0.21
beav.pred.post <- inv.logit(beav.pred.post)
mean(beav.pred.post)
# 0.6707449

### Predict prevalence in non-beaver ponds (assuming they are not permanent)
nonbeav.pred.post <- beta.0.post + 
  beta.ralu.post * ralu1 +
  #coef['permanent']  +  # assume not permanent, so 0
  #coef['PC2'] * ralu.sites$PC2 +  # assume average value of 0 
  #coef['PC1'] * ralu.sites$PC1 +  # assume average value of 0 
  beta.canopy.post * 0.42
nonbeav.pred.post <- inv.logit(nonbeav.pred.post)
mean(nonbeav.pred.post)
# 0.6444863

#### It would probably be a good idea to do 2 separate non-beaver distributions
#### One for permanent and another for not permanent.

#################### Violin Plot for distribution comparison

# put in one data frame for violin plot function
ralu1comparison <- data.frame(prevalence=rbind(beav.pred.post, nonbeav.pred.post),
                              type=c(rep("beaver", length(beav.pred.post)), 
                              rep("nonbeaver", length(beav.pred.post))))
#violin plot
ggplot(ralu1comparison, aes(x=type, y= prevalence)) +
 geom_violin() +
  ggtitle("Comparison at mean frog density of 0.13") +
  xlab("") +
  ylab("Posterior Predicted Prevalence")
  
  
  
  