###############################################################################
## Code Run robust design mark-recapture models on formatted data
## code to format data and specify model are in separate files
###############################################################################

library(R2jags)
library(MCMCvis)

## load in formatted data
load("formattedCHdata.RData")

## source functions
source("RDMS_Functions.R")

## source model code
source("MSinfRDarray_model_specification2.R") # create a separate file with model specification 

n.months <- max(months.trapped)

# make a version not robust design for the known state and initial value functions
CH.primary <- primary.MSch.array.fun(CH.secondary)  

# create a vector of first marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH.primary, 1, get.first)


#delete individuals that were caught for the first time at the last primary occasion
# they don't provide data and mess with code
rms <- which(f==dim(CH.primary)[2])
if(length(rms)>0){
  CH.secondary <- CH.secondary[-rms, , ]
  CH.primary <- CH.primary[-rms, ]
  f <- f[-rms]
  individual.covariates <- individual.covariates[-rms,]
}

# create a separate f in terms of months captured (not full months)
f.prim <- temporal.covariates$Prim[f]


n.sec.occasions <- apply(CH.secondary, 3, function(x){
   sum(length(which(is.finite(x[1,]))))
})

# These models need a 3 for not seen instead of 0
CH.secondary3 <- replace(CH.secondary, CH.secondary==0, 3)


##### Bundle data
bugs.data=list(y = CH.secondary3,
               f = f,
               f.prim = f.prim,
               f.state = individual.covariates$f.state,
               nind = dim(CH.secondary)[1],
               n.sec.occ = n.sec.occasions, # contains 0's for months not trapped
               n.months = n.months, # total span of months (incl not trapped)
               prim.occasions = months.trapped, # months out of n.months actually trapped (if all trapped then same as previous line)
               sex = as.numeric(individual.covariates$sex), 
               season  = temporal.covariates$season, # 1:4  
               z = known.state.SImsInf(CH.primary))



## Initial values ----------------------------------------------------------- ##


# supply initial values
inits <- function() {
  list(
    z                  = MSinf.init.z(CH.primary,rep(n.months,bugs.data$nind)),
    alpha.0            = runif(1, 0, 1), 
    alpha.season       = runif(3, 0, 1), 
    alpha.inf          = runif(1, 0, 1),
    alpha.inf.male     = runif(1, 0, 1),
    sigma.0            = runif(1, 0, 1),
    sigma.inf          = runif(1, 0, 1),
    beta.0             = runif(1, 0, 1),
    beta.male          = runif(1, 0, 1)
  )
}


## Parameters monitored ----------------------------------------------------- ##


# parameters monitored - we can add things that aren't listed
# won't monitor initial values or parameters that aren't specified in the code
# see Zuni model run for discussion on monitoring phi
parameters <- c("alpha.0", 
                "alpha.season", 
                "alpha.inf",
                "alpha.inf.male" ,
                "sigma.0", 
                "sigma.inf",
                "beta.0",
                "beta.male"
)


## Run the Model --------------------------------------------------------- ##


#MCMCsettings  # set to low numbers for de-bugging - set higher for model run
ni=10000
nt=6
nb=3000
nc=3



date()
## Call JAGS from R
robust.MSinf=jags(bugs.data,
                inits,
                parameters,
                "MSinfRD_specification2.bug", # name from your sourced model specification file
                n.chains=nc,
                n.thin=nt,
                n.iter=ni,
                n.burnin=nb)
date() # to tell how long it ran


save(robust.MSinf, file="ModelOutput.RData")

#summarize posteriors
print(robust.MSinf,digits=3)  



library(mcmcplots)
library(MCMCvis)
# basic plots to check convergence and credible intervals,
mcmcplot(robust.MSinf,
         parms =  parameters
)


# MCMC vis package   
MCMCplot(robust.MSinf,
         #xlim = c(-1, 1),
         ref_ovl = TRUE,
         params = parameters
         
)

