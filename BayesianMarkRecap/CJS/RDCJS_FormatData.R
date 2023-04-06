###############################################################################
## Code to format data for robust design mark-recapture models, including
## creating capture histories
## formatting individual and temporal covariates
## & other necessary variables to run the model (run is code in a another file)
###############################################################################


## Source file with functions for formatting data
source("RDCJS_Functions.R")

## Load in  your cleaned data (and dirty if there are dates not in cleaned dataset)
load("AllCaptureData.RData") # or your cleaned capture data 

## Create a list of monthly matrices 
grandcanyonT.CH.list <- CJS.RDcapture.history.list.fun(
  cleaned.data = southwest.final.clean, # 'cleaned' e.g. captures with no tag numbers removed - formatted
  site.webs = "GrandCanyon.T", # your data need a site column and web column and here it is pasted together
  species = "PM",
  dirty.data = NULL # might need the 'dirty' data because it has all the 
  # dates including those with no captures.
  # if your data are already 'clean' but there 
  # weren't any months where no deer mice were caught, 
  # then leave NULL.
) 

Ch.list <- grandcanyonT.CH.list$Ch.list
# list of length 46 (1 for each month)
# each element is a matrix where rows represent individuals
# and columns represent secondary occasion (days)

##### These do not include months not trapped - but later functions want that so need to update
sec.occ.list <- grandcanyonT.CH.list$sec.occ.list
# list of secondary occasions (dates) per primary occasion (session)
n.sec.occ <- grandcanyonT.CH.list$n.sec.occ
# number of secondary occasions for each primary occasion

sessions <- grandcanyonT.CH.list$all.sessions

# this just has basic covariates like month, season, year, but function can 
# also put in temporal covariates like ndvi 
covariate.data <- monthly.covariate.fun(cleaned.data = southwest.final.clean,
                      CH.secondary = Ch.list,
                      tags = rownames(Ch.list[[1]]),
                      sessions = sessions,
                      site.webs = "Grandcanyon.T")

individual.covariates <- covariate.data$individual.covariates
temporal.covariates <- covariate.data$temporal.covariates

months.trapped <- temporal.covariates$long.month[which(is.finite(temporal.covariates$session))]
n.months <- max(months.trapped) #months spanned

# reformat CH from list to array
CH.secondary <- list.to.array.fun(Ch.list, temporal.covariates)

save(CH.secondary, individual.covariates, temporal.covariates, months.trapped, file="formattedCHdata.RData")
