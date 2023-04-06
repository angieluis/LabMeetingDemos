########################################################################################
### Code to run Robust Design multistate infection model 
### This is on simulated data in an array format
########################################################################################

## This contains the covariates from the  Grand Canyon Max Model  ------------------- ##
 # where ndvi12, prcp6, tmin0, and tmax6 were selected as variables to include 


 # packages
   library(R2jags)
   library(MCMCvis)
   library(mcmcplots)



 # load RData specific to model run
   load("MSinfRDarray_SimulatedData.RData")

 # source functions and model
   #source("01_sw_data_functions_more.R")
   source("MSinfRDarray_model_specification.R")

   
  
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
   }
   f.state <- numeric()
   for(i in 1:dim(CH.primary)[1]){
     f.state[i] <- CH.primary[i,f[i]]
   }

   # create a separate f in terms of months captured (not full months)
   Prim <- match(1:n.months, prim.occasions)
   f.prim <- Prim[f]
   

## Bundle data -------------------------------------------------------------- ##


 # list of bugs data
   bugs.data <- list(
      nind                  = dim(CH.secondary)[1],
      y                     = CH.secondary, # 1=observed as S; 2=observed as I; 3=not observed
      f                     = f, 
      f.state               = f.state,
      f.prim                = f.prim,
      sex                   = as.numeric(indiv.data$sex), 
      n.sec.occ             = n.sec.occ, 
      monthlyCH             = CH.primary,
      z                     = known.state.SImsInf(CH.primary),
      MSinf.init.z          = MSinf.init.z,
      n.months              = dim(CH.secondary)[3],
      season                = tempcovdata$season, # 1:4  
      prim.occasions        = prim.occasions,
      ndvi                  = tempcovdata$ndvi12
       )


## Initial values ----------------------------------------------------------- ##


 # supply initial values
   inits <- function() {
      list(
         z                  = bugs.data$MSinf.init.z(bugs.data$monthlyCH,rep(bugs.data$n.months, bugs.data$nind)), 
         alpha.0            = runif(1, 0, 1), 
         alpha.season       = runif(3, 0, 1), 
         #alpha.male         = runif(1, 0, 1), 
         alpha.ndvi         = runif(1, 0, 1),
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
                   #"alpha.male", 
                   "alpha.season", 
                   "alpha.ndvi",
                   "alpha.inf",
                   "alpha.inf.male" ,
                   
                   "sigma.0", 
                   "sigma.inf",
                   "beta.0",
                   "beta.male"
                   
                    )


## Run the Mdodel --------------------------------------------------------- ##


 # date before run  
   date()


   MSinfRD.Model <- jags(data     = bugs.data,
                                            inits, 
                                            parameters, 
                                            "MSinfRD_specification.bug", 
                                            n.chains = 3, 
                                            n.thin   = 6, 
                                            n.iter   = 10000, 
                                            n.burnin = 3000)

   
    
   
 # date after run
   date() # 1 hr 20min for 10000,3000


 # save session data
   save.image("MSinfRDModelOutput.RData")

   


   
## -------------------------------------------------------------------------- ##
   
   
 # basic review
   
 
 # basic plots to check convergence and credible intervals,
   mcmcplot(MSinfRD.Model,
            parms =  c("alpha.0", 
                       #"alpha.male", 
                       "alpha.season", 
                       "alpha.ndvi",
                       "alpha.inf ",
                       "alpha.inf.male" ,
                       
                       "sigma.0", 
                       "sigma.inf",
                       "beta.0",
                       "beta.male")
                       
            )
 
   
 # MCMC vis package   
   MCMCplot(MSinfRD.Model,
            xlim = c(-4, 4),
            ref_ovl = TRUE,
            params = c("alpha.0", 
                       #"alpha.male", 
                       "alpha.season", 
                       "alpha.ndvi",
                       "alpha.inf ",
                       "alpha.inf.male" ,
                       
                       "sigma.0", 
                       "sigma.inf",
                       "beta.0",
                       "beta.male")
            
              )
   