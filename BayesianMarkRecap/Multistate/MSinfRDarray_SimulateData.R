#######################################################################
## Simulate Data for Robust design multistate infection models
## to make sure all the code is working
#######################################################################

library(tidyverse)

# Use sites grandcanyon.e as a template
# so use covariate data for those sites for the first nmonths

#source("01_sw_data_functions_more.R")
sw.temp <- read.csv("southwest_covariates_norm.csv")
sw.temp$date2 <- lubridate::dmy(paste("1", sw.temp$date))
sw.temp$site.web <- paste(sw.temp$site,sw.temp$web,sep=".")
# add season covariate 0=winter
sw.temp$season <- ifelse(sw.temp$month == 12 | 
                              sw.temp$month == 1 | 
                              sw.temp$month == 2 | 
                              sw.temp$month == 3, 1, 
                            ifelse(sw.temp$month == 4 | 
                                     sw.temp$month == 5 | 
                                     sw.temp$month == 6, 2, 
                                   ifelse(sw.temp$month == 7 | 
                                            sw.temp$month == 8 | 
                                            sw.temp$month == 9, 3, 4)))

webs <- "grandcanyon.e"

ninds <- 300
ntot <- sum(ninds)
n.months <- 40   
prim.occasions <- c(1:4,6:40) # months trapped. So not trapped in month 5
n.sec.occ <- c(rep(3,4), 0, rep(3,5), 4, rep(3,29)) # all but month 11 have 3 days, it has 4 days of trapping
start.date <- sw.temp$date2[20] # Aug 1994
end.date <- start.date + months(n.months-1)
n.webs <- length(ninds)
maxI <- 20 # to standardize I.dat so beta not so small
#prim.dates <- seq(start.date,end.date,by="month")

sessions <- character()
for(i in 1:length(prim.dates)){
  sessions[i] <- ifelse(month(prim.dates[i])<10 ,paste(year(prim.dates[i]),"0",month(prim.dates[i]),sep=""),paste(year(prim.dates[i]),month(prim.dates[i]),sep=""))
}


### set parameter values
alpha.0            = 1 
alpha.season       = c(1,-0.5,0) 
alpha.ndvi         = 0.2
alpha.inf          = -1
alpha.inf.male     = -0.5
sigma.0            = -0.1 #-1.5
sigma.inf          = 0.3
beta.0             = -2
beta.male          = 0.5
#beta.I             = 2 # should really be a function of I but simplify here

alpha.season.use       = c(0,alpha.season)

rev.logit <- function(x) {
  exp(x) / (1 + exp(x))
}

## Simulate data for each site separately, then paste them together----------------------#
tempcovdata <- list()
indiv.data <- list()
z <- list()

for(w in 1:n.webs){
  web <- webs[w]
  tempcovdata[[w]] <- sw.temp[which(sw.temp$site.web==web),]
  s.ind <- which(tempcovdata[[w]]$date2==start.date) 
  tempcovdata[[w]] <- tempcovdata[[w]][s.ind:(s.ind+(n.months-1)),]

  n <- ninds[w]
  
  # first half female (0), second half male (1)
  indiv.data[[w]] <- data.frame(ID = 1:n, site.web = web, 
                                sex = c(rep(0,n/2), rep(1,n/2)))
  sex <- indiv.data[[w]]$sex
  
  # Model params as functions of covariates ------------------------#

  ############ Survival
  phiS <- matrix(NA,n,n.months-1)
  phiI <- matrix(NA,n,n.months-1)
  for(i in 1:n) {
    for(m in 1:(n.months - 1)) {            
      
      ### Phi for uninfected
      phiS[i, m] <-  rev.logit(          
        
        alpha.0 +
        alpha.season.use[tempcovdata[[w]]$season[m]] +       

        # environmental covariates
        alpha.ndvi * tempcovdata[[w]]$ndvi[m] 
          )
      
      ### Phi for infected
      phiI[i, m] <-   rev.logit(           
        
        alpha.0 +
          alpha.season.use[tempcovdata[[w]]$season[m]] +       

          # environmental covariates
          alpha.ndvi * tempcovdata[[w]]$ndvi[m] +

        # infection terms
        alpha.inf +
        alpha.inf.male * sex[i])
      
    } # m for months
  } # i for individual
  
  ## Need to calculate psi during simulation
   
   
  # Simulate Z -----------------------------------------------------#
  z[[w]] <- matrix(0,n,n.months)
  # simulate when each indiv first enters pop
  f <- sample(1:n.months,size=n,replace=TRUE)
  # simulate state at first capture. 98% uninfected, 2% infected
  f.st <- sample(1:2,size=n,replace=TRUE,prob=c(0.98,0.02)) 
  
  # add those alive the first month
  alive.first.month <- which(f==1)
  z[[w]][alive.first.month,1] <- f.st[alive.first.month]
  
  for(m in 2: n.months){
    I.lasttime <- sum(length(which(z[[w]][,m-1]==2)) )
    for(i in 1:n){  
      if(m==f[i]){ 
        z[[w]][i,m] <- f.st[i]
      } 
      if(m>f[i]){
        last.state <- z[[w]][i,m-1]
      
        # if last.state uninfected
        if(last.state==1){
          # calculate prob of becoming infected
           
          psiSI <- rev.logit(beta.0 + beta.male * sex[i]) # + beta.I * I.lasttime/maxI)
          
          # stay uninfected (1) with prob phiS*(1-psiSI)
          # become infected (2) with prob phiS*psiSI
          # die (0) with prob 1-phiS
          z[[w]][i,m] <- sample(c(1,2,0),size=1,prob=c(
            phiS[i,m-1] * (1-psiSI),
            phiS[i,m-1] * psiSI,
            1-phiS[i,m-1])
          )
        }
        
        # if last.state infected
        if(last.state==2){
         # survive (stay 2) with prob phiI
         # die (0) with prob 1-phiI
          z[[w]][i,m] <- rbinom(1,1,phiI[i,m-1])*2
        }
        #if last.state dead
        if(last.state==0){
          z[[w]][i,m] <- 0
        }
      } # if m>f 
    } #i
  } #m

} #webs

z <- z[[1]] 
tempcovdata <- tempcovdata[[1]]
tempcovdata$long.month <- 1:n.months
indiv.data <- indiv.data[[1]]




# Simulate observations in a array format -----------------------------------------#
# 1 site

y <- array(NA, dim=c(ninds, max(n.sec.occ), n.months))

pS <- matrix(NA,ninds,n.months)
pI <- matrix(NA,ninds,n.months)

  for(i in 1:ninds){
    f <- min(which(z[i,]>0))
    
    for(m in f:n.months){
        pS[i, m] <- rev.logit(
          sigma.0 )                       # intercept
            # could add time varying or individual covariates                

     
        pI[i, m] <- rev.logit(
          sigma.0 + # intercept
          sigma.inf # infection term
          )
    } #m
    
    for(m in prim.occasions){
        
      # simulate secondary captures and paste onto obs
        z.state <- z[i,m]
        
        p <- ifelse(z.state==1, pS[i,m], ifelse(z.state==2, pI[i,m], 0))
        
        obs.state <- rbinom(n.sec.occ[m], 1, p) * z.state
        
        
        y[i, 1:n.sec.occ[m], m] <- obs.state
          
   } #m
 } #i


# remove individuals never caught (because of incomplete sampling)
never.caught <- which(apply(y,1,sum,na.rm=TRUE)==0)

# remove from secondary history
CH.secondary <- y[-never.caught , , ]

# remove from individual covariate data
indiv.data <- indiv.data[-never.caught, ]

# create primary history
CH.primary <- apply(CH.secondary,c(1,3), function(x){
  ifelse( 2 %in% x , 2,
          ifelse(1 %in% x, 1, 0))
  })
CH.primary[,which(1 - (1:n.months %in% prim.occasions )==1)] <- NA # make sure keep NAs when whole month not trapped 

## Change so now 3 represents not caught
CH.secondary <- replace(CH.secondary,CH.secondary==0, 3) # 3 for not caught 





save.image(file="MSinfRDarray_SimulatedData.RData")


