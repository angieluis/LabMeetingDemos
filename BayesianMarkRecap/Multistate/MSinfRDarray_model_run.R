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
   CH.secondary <- CH.secondary[-rms, , ]
   CH.primary <- CH.primary[-rms, ]
   
   f <- f[-rms]
   f.state <- numeric()
   for(i in 1:dim(CH.primary)[1]){
     f.state[i] <- CH.primary[i,f[i]]
   }

   # create a separate f in terms of months captured (not full months)
   Prim <- match(1:n.months, prim.occasions)
   f.prim <- Prim[f]
   
## functions for known state and initial values for z ----------------------- ##
   MSinf.init.z <-  function(ch, #primary CH
            n.months){ # n.months is the length of months in the dataset by individual because can differ by web [i,m] . I think assumes that all sites start at the same time?
     kn.state <- known.state.SImsInf(ms = ch)
     f <- apply(ch,1,function(x){min(which(x > 0))})
     state <- matrix(NA, nrow = dim(ch)[1], ncol = dim(ch)[2]) 
     # fill in with first state caught
     for(i in 1:dim(ch)[1]){
       f.state <- ch[i,f[i]]
       state[i,] <- rep(f.state,dim(ch)[2]) 
     }
     # remove those that are in the known state
     state <- replace(state,!is.na(kn.state),NA)
     
     for(i in 1:(dim(state)[1])){
       state[i,1:f[i]] <- NA # put NA for when first caught (in likelihood)
       
       if(length(which(kn.state[i,] == 2)) > 0){ 
         maxI <- max(which(kn.state[i,] == 2))
         if(maxI < dim(state)[2] ){
           state[i, (maxI + 1):dim(state)[2]] <- 2 # all after caught as I are I (2)
         }
       }
       if(n.months[i]!=max(n.months)){
         state[i,(n.months[i]+1):dim(ch)[2]] <- NA # replace all after last month in the dataset with NA
       }
     }
     return(state)
   }
     
   known.state.SImsInf <- function(ms){ # ms is multistate primary capture history
     # notseen: label for 'not seen' #here is 3
     state <- ms
     state[state == 0] <- NA
     for(i in 1:dim(ms)[1]){
       n1 <- min(which(ms[i,]>0))
       if(length(which(ms[i, ] == 2)) > 0){ #filling in I's where can
         minI <- min(which(ms[i, ] == 2)) #I's are observation 2
         maxI <- max(which(ms[i, ] == 2))
         state[i, minI:maxI] <- 2}         # I's are state 3
       if(length(which(ms[i, ] == 1)) > 0){  #filling in S's where can
         minS <- min(which(ms[i, ] == 1))  # S's are observation 1
         maxS <- max(which(ms[i, ] == 1))
         state[i, minS:maxS] <- 1}         # S's are state 2
       state[i,n1] <- NA
     }
     
     return(state)
   }
   
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


## Model specification ------------------------------------------------------ ##


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
   