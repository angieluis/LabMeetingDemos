###############################################################################
## Code Run robust design mark-recapture models on formatted data
## code to format data and specify model are in separate files
###############################################################################

library(R2jags)
library(MCMCvis)

## load in formatted data
load("formattedCHdata.RData")

## source functions
source("RDCJS_Functions.R")

## source model code
source("RDCJS_BasicModel.R") # create a separate file with model specification 

n.months <- max(months.trapped)

# make a version not robust design for the known state and intital value functions
CH.primary <- apply(CH.secondary,c(1,3),sum,na.rm=TRUE)
CH.primary[,which(1 - (1:n.months %in% months.trapped)==1)] <- NA # make sure keep NAs when whole month not trapped 
CH.primary <- replace(CH.primary,CH.primary>1,1)

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

##### Bundle data
bugs.data=list(y = CH.secondary,
               f = f,
               f.prim = f.prim,
               nind = dim(CH.secondary)[1],
               n.secondary.occasions = n.sec.occasions, # contains 0's for months not trapped
               n.primary.occasions = n.months, # total span of months (incl not trapped)
               prim.occasions = months.trapped, # months out of n.months actually trapped (if all trapped then same as previous line)
               z = known.state.cjs(CH.primary))



#initial values
inits=function(){list(
  z=cjs.init.z(CH.primary,f),
  mean.phi=runif(1,0,1),
  mean.p=runif(1,0,1),
  mean.c=runif(1,0,1))}

#parameters monitored
parameters=c("mean.phi","mean.p","mean.c")

#MCMCsettings
ni=10000
nt=6
nb=5000
nc=3



date()
## Call JAGS from R
robust.cjs=jags(bugs.data,
                inits,
                parameters,
                "robust_cjs.bug", # name from your sourced model specification file
                n.chains=nc,
                n.thin=nt,
                n.iter=ni,
                n.burnin=nb)
date() # to tell how long it ran
# Basic model on all (46 months) of GrandCanyon T data took 13 minutes

save(robust.cjs, file="ModelOutput.RData")

#summarize posteriors
print(robust.cjs,digits=3)  



library(mcmcplots)
library(MCMCvis)
# basic plots to check convergence and credible intervals,
mcmcplot(robust.cjs,
         parms =  c("mean.c", 
                    "mean.p",
                    "mean.phi")
)


# MCMC vis package   
MCMCplot(robust.cjs,
         #xlim = c(-1, 1),
         ref_ovl = TRUE,
         params = c("mean.c", 
                    "mean.p",
                    "mean.phi")
         
)

