library(tidyverse)

spp.sites <- read.csv("spp.sites.csv")

## RALU DATA
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



# predict values:
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

pred <- inv.logit(pred.u)

plot(pred, predict(ralu.prev.best2))
plot(pred.u, predict(ralu.prev.best2)) ## Weird! The predict function doesn't back transform

# Assume Beaver ponds are permanent and have canopy cover of 21 
# & non-beaver ponds are not permanent and have canopy cover of 42

# assume PC1 and PC2 values of 0 (sort of average)
# explore the effect of ralu_density

frog.density <- seq(0, 0.4, length=30)

beav.pred1 <- coef['(Intercept)'] + 
  coef['ralu_density'] * frog.density +
  coef['permanent']  + 
  #coef['PC2'] * ralu.sites$PC2 +  # assume average value of 0 
  #coef['PC1'] * ralu.sites$PC1 +  # assume average value of 0 
  coef['canopy_cover'] * 21

beav.pred <- inv.logit(beav.pred1)

nonbeav.pred1 <- coef['(Intercept)'] + 
  coef['ralu_density'] * frog.density +
  #coef['permanent']  +  # assume not permanent, so 0
  #coef['PC2'] * ralu.sites$PC2 +  # assume average value of 0 
  #coef['PC1'] * ralu.sites$PC1 +  # assume average value of 0 
  coef['canopy_cover'] * 42

nonbeav.pred <- inv.logit(nonbeav.pred1)

plot(frog.density, beav.pred,type="l", col="steelblue",lwd=2, ylab="predicted prevalence")
lines(frog.density, nonbeav.pred, col="olivedrab", lwd=2)
legend("bottomright", c("beaver", "non-beaver"), lwd=2, col=c("steelblue","olivedrab"),bty="n")







###############################################################################
### Now Bayesian!!
###############################################################################

library(R2jags)
library(MCMCvis)
library(mcmcplots)

sink("GLM_Raluprev.jags")
cat("
model {

# Priors
beta.0 ~ dunif(-5, 5)
beta.ralu ~ dunif(-5, 5)
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
ni <- 10000
nt <- 2
nb <- 2000
nc <- 3

# Call JAGS from R
out <- jags(data = jags.data, inits = inits, parameters.to.save = params, 
            model.file = "GLM_Raluprev.jags", n.chains = nc, n.thin = nt, 
            n.iter = ni, n.burnin = nb, working.directory = getwd())



mcmcplot(out,parms=params[1:6])

MCMCplot(out,params=params[1:6])

## add the GLM points estimates
points(coef[c(1:3,5,4)],6:2,col="red")
points(coef[6]*100,1,col="red") # need to mult canopy cover back to 100

