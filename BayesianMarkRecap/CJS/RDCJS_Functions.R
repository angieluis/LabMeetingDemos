
library(tidyverse)

CJS.RDcapture.history.list.fun <- function(
  # assuming this dirty data has all the dates 
  # including those trapped but no animals/pema were caught
  dirty.data = NULL, # or if not all dates are in the clean data : southwest.dirty, 
  cleaned.data = southwest.final.clean, # 'cleaned'
  site.webs = "GrandCanyon.E", # chracter string with site and web pasted together with a period 
  species = "PM") {
  
  ## function returns a list of monthly matrices
  ## first matrix represents the first month of trapping
  ## where rows represent individuals
  ## and columns represent secondary occasion (days)
  
  # If there are NOT dirty data, make some:
  if(length(dirty.data) == 0){
    dirty.data <- cleaned.data
    dirty.data$date <- paste(as.character(month(cleaned.data$date)),
                             as.character(day(cleaned.data$date)),
                             as.character(year(cleaned.data$date)),sep="/")

  }
  # make everything lowercase
  site.webs <- tolower(site.webs)
  species <- tolower(species)
  names(dirty.data) <- tolower(names(dirty.data))
  dirty.data$site <- tolower(dirty.data$site)
  dirty.data$web <- tolower(dirty.data$web)
  names(cleaned.data) <- tolower(names(cleaned.data))
  cleaned.data$site <- tolower(cleaned.data$site)
  cleaned.data$letter_2 <- tolower(cleaned.data$letter_2)
  cleaned.data$web <- tolower(cleaned.data$web)
  
  # create site.web column in data
  dirty.data$site.web <- paste(dirty.data$site, 
                               dirty.data$web, 
                               sep = ".")
  cleaned.data$site.web <- paste(cleaned.data$site, 
                                 cleaned.data$web, 
                                 sep = ".")
  
  
  
  # cut data to just the site.webs wanted
  ind.dirty <- numeric()
  ind.clean <- numeric()
  for (i in 1:length(site.webs)) {
    ind.dirty <- c(ind.dirty, which(dirty.data$site.web == site.webs[i]))
    ind.clean <- c(ind.clean, which(cleaned.data$site.web == site.webs[i]))
  }
  dirty.data <- dirty.data[ind.dirty, ]
  cleaned.data <- cleaned.data[ind.clean, ]
  
  # new date column in right format
  dirty.data$date1 <- as.Date(gsub(" ", 
                                   "", 
                                   dirty.data$date), 
                              format = "%m/%d/%Y")
  cleaned.data$date1 <- as.Date(gsub(" ", 
                                     "", 
                                     cleaned.data$date), 
                                format = "%m/%d/%Y")
  
  # new site.tag column of data
  cleaned.data$site.tag <- paste(cleaned.data$site, 
                                 cleaned.data$tag, 
                                 sep = ".")
  all.sessions <- sort(unique(dirty.data$session))
  web.dates <- list()
  
  for (i in 1:length(site.webs)) {
    web.dates[[i]] <- sort(unique(dirty.data$date1[which(dirty.data$site.web == site.webs[i])]))
  }
  names(web.dates) <- paste("web", site.webs, sep = ".")
  
  # calculate the number of secondary occasions for each sessions at each site
  sec.occ.list <- list()
  for (i in 1:length(site.webs)) {
    
    # unique sessions for that site.web
    s <- sort(unique(dirty.data$session[dirty.data$site.web == site.webs[[i]]])) 
    sec.occ.list[[i]] <- list()
    for (m in 1:length(s)) {
      sec.occ.list[[i]][[m]] <- sort(unique(dirty.data$date1[dirty.data$site.web == site.webs[[i]] & dirty.data$session == s[m]]))
    }
    names(sec.occ.list[[i]]) <- s
  }
  names(sec.occ.list) <- site.webs
  n.sec.occ.list <- list()
  for (i in 1:length(site.webs)) {
    n.sec.occ.list[[i]] <- unlist(lapply(sec.occ.list[[i]], length))
  }
  names(n.sec.occ.list) <- site.webs
  
  # separate out deermice. These are the site.tag numbers I want to keep (are pema).
  sp.tags <- sort(unique(cleaned.data$site.tag[which(cleaned.data$letter_2 == species)]))
  if (length(which(sp.tags == -9)) > 0) {
    sp.tags <- sp.tags[-which(sp.tags == -9)]
  }
  
  # what are the rows of data that have those tags?
  sp.ind <- numeric()
  for (i in 1:length(sp.tags)) {
    sp.ind <- c(sp.ind, which(cleaned.data$site.tag == sp.tags[i]))
  }
  
  # cleaned data for sites/webs with species in which i'm interested
  sp.data <- cleaned.data[sp.ind, ] 
  
  # max number of secondary occasions per primary (among all sites)
  max.sec.occ <- numeric() 
  for (i in 1:length(all.sessions)) {
    max.sec.occ[i] <- max(unlist(lapply(n.sec.occ.list, function(x) {
      x[which(names(x) == all.sessions[i])]
    })))
  }
  
  
  # first set up blank capture history list to fill in below
  all.IDs <- sort(unique(sp.data$site.tag))
  Ch.list <- list()
  for (i in 1:length(all.sessions)) {
    mat <- matrix(NA, nrow = length(all.IDs), ncol = max.sec.occ[i])
    rownames(mat) <- all.IDs
    colnames(mat) <- 1:(dim(mat)[2])
    Ch.list[[i]] <- mat
  }
  
  # fill in by web
  for (w in 1:length(site.webs)) {
    web.IDs <- unique(sp.data$site.tag[which(sp.data$site.web == site.webs[w])])
    Session.days <- sec.occ.list[[w]]
    web.dat <- sp.data[which(sp.data$site.web == site.webs[w]), ]
    for (m in 1:length(Session.days)) {
      session <- which(all.sessions == names(Session.days)[m])
      days <- Session.days[[m]]
      for (d in 1:length(days)) {
        for (i in 1:length(web.IDs)) {
          indiv <- which(rownames(Ch.list[[1]]) == web.IDs[i])
          id.dat <- web.dat[which(web.dat$date == days[d] & web.dat$site.tag == web.IDs[i]), ]
          Ch.list[[session]][indiv, d] <- ifelse(dim(id.dat)[1] > 0, 1, 0)
        } # i
      } # d
      cat("web = ", w, "; session = ", m, "\n")
    } # m
  } # w
  
  big.list <- list(Ch.list, sec.occ.list[[1]], n.sec.occ.list[[1]], all.sessions)
  names(big.list) <- c("Ch.list", "sec.occ.list", "n.sec.occ", "all.sessions")
  return(big.list)
}


## Create a primary monthly ch from a secondary list ------------------------ ##

# CJS
primary.ch.fun <- function(CH.secondary) { # as list of monthly matrices
  CH.primary <- matrix(NA, 
                       nrow = dim(CH.secondary[[1]])[1], 
                       ncol = length(CH.secondary))
  for (i in 1:dim(CH.primary)[2]) {
    CH.primary[, i] <- apply(CH.secondary[[i]], 1, function(x) {
      ifelse(length(which(is.na(x))) == length(x), NA, sum(x, na.rm = TRUE))
    })
  }
  CH.primary <- replace(CH.primary, CH.primary > 1, 1)
  return(CH.primary)
}


# Multistate 
primary.MSch.fun <- function(CH.secondary){ # as list of monthly matrices 
  CH.primary <- matrix(NA,nrow=dim(CH.secondary[[1]])[1], ncol=length(CH.secondary))
  for(m in 1:length(CH.secondary)){
    for(i in 1:dim(CH.secondary[[1]])[1]){
      chs <- CH.secondary[[m]][i,]
      if(length(which(is.na(chs))) == length(chs)){
        CH.primary[i,m] <- NA
      } else{
        chs <- chs[is.finite(chs)]
        CH.primary[i,m] <- ifelse(sum(chs) == 0, 0, unique(chs[which(chs > 0)]))
      } 
    }
  }                     
  
  return(CH.primary)
}
 








## Monthly covariate function ------------------------------------- ##

# format individual covariates and temporal covariate (season, ndvi, etc)
monthly.covariate.fun <- function(
  cleaned.data = southwest.final.clean,
  
  # as monthly list
  CH.secondary = Ch.list,
  
  # tag names that line up to CH.secondary
  tags = rownames(Ch.list[[1]]), 
  
  # if TRUE, then tags are specified by site e.g., "Zuni.1101"
  by.sitetag = TRUE, 
  
  # all sessions to include e.g. 199806 (sessions trapped even no pm caught
  sessions = session.list$all.sessions,
  
  # data frame of monthly temporal data, with either a column
  # called date or yearmon (must include months not trapped),
  # must be in long format like sw.temp.data and have all the
  # sites and all temporal data with no time lags
  temporal.data = NULL, # sw.temp.data
  
  multistate = FALSE, # if multistate model, 
  # need state at first capture
  
  
  # long data frame of species MNAs and diversities 
  # species 2 letter codes in lowercase
  # cov.list below must match column headers
  # input should be scaled (max of any variable=1 across all sites)
  diversity.data = NULL, #scaled.MNAs.diversity.longdata, 
  
  #  e.g. c("Zuni.1","Zuni.2", "Navajo.1","Navajo.2"
  site.webs = c("Grandcanyon.T"),
  
  # list of temporal covariates and their time lags,
  # e.g, list(ndvi=0,ndvi=1,tmax=3) means use ndvi with no lag
  # and with a lag 1 and tmax with lag 3
  cov.list = NULL, #list(ndvi = 8, prcp = 0, tmax = 0), 
  remove.na = FALSE #
) {
  
  names(cleaned.data) <- tolower(names(cleaned.data))
  tags <- tolower(tags)
  site.webs <- tolower(site.webs)
  cleaned.data$site.tag <- paste(tolower(cleaned.data$site), 
                                 tolower(cleaned.data$tag), sep = ".")
  cleaned.data$site.web <- paste(tolower(cleaned.data$site), 
                                 tolower(cleaned.data$web), sep = ".")
  if(length(temporal.data) > 0){
  temporal.data$site.web <- paste(tolower(temporal.data$site), 
                                  tolower(temporal.data$web), sep = ".")
  }
  nind <- length(tags)
  ind.clean <- numeric()
  for (i in 1:length(site.webs)) {
    ind.clean <- c(ind.clean, which(cleaned.data$site.web == site.webs[i]))
  }
  cleaned.data <- cleaned.data[ind.clean, ]
  
  ## individual covariate data frame to keep track of sex 
  ic <- data.frame(ID = 1:nind, tag = tags)
  site.webi <- character()
  sex <- numeric() # 1 male, 0 female
  for (i in 1:nind) {
    if (by.sitetag == FALSE) {
      ind <- which(cleaned.data$tag == tags[i])
    } else {
      ind <- which(cleaned.data$site.tag == tags[i])
    }
    x <- cleaned.data[ind, ]
    x <- x[order(x$session), ]
    site.webi[i] <- as.character(x$site.web[1])
    sex[i] <- max(x$sex) # they aren't NAs but -9
  }
  ic$web <- factor(site.webi)
  ic$sex <- replace(sex, sex == -9, 0.5) # split the difference for unknown sexes
  
  
  ## create monthly data frame & number sessions and find missing sessions
  sessions.trapped <- sort(unique(sessions))
  session.dates <- ym(sessions.trapped)
  prim.dates <- seq(session.dates[1], session.dates[length(session.dates)], by="month")
  month.data <- data.frame(long.month = 1:length(prim.dates),
                           prim.date=prim.dates, 
                           session = sessions.trapped[match(prim.dates,session.dates)],
                           year = year(prim.dates),
                           month = month(prim.dates),
                           # covariate.prim? not sure
                           Prim = match(prim.dates,session.dates))
 
  month.data$season <- ifelse(month.data$month == 12 | 
                                month.data$month == 1 | 
                                month.data$month == 2 | 
                                month.data$month == 3, 1, 
                              ifelse(month.data$month == 4 | 
                                       month.data$month == 5 | 
                                       month.data$month == 6, 2, 
                                     ifelse(month.data$month == 7 | 
                                              month.data$month == 8 | 
                                              month.data$month == 9, 3, 4)))
  
  
  if (length(temporal.data) > 0) { # |length(diversity.data)>0
    ls <- length(site.webs)
    datas <- temporal.data[which(temporal.data$site.web == site.webs[1]), ]
    if (ls > 1) {
      for (i in 2:ls) {
        datas <- rbind(datas, temporal.data[which(temporal.data$site.web == site.webs[i]), ])
      }
    }
    
    # make a wide data frame with date/yearmon going from first trapped 
    # session to last trapped session
    wdate <- lubridate::dmy(paste("1", 
                                  month.data$month, 
                                  month.data$year, 
                                  sep = "-"))
    data.w <- data.frame(date = sort(unique(wdate)))
    data.w$year <- lubridate::year(data.w$date)
    data.w$month <- lubridate::month(data.w$date)
    datas$date <- lubridate::dmy(paste("1", datas$date))
    
    ##### Now paste in diversity data if present
    if (length(diversity.data) > 0) {
      
      diversity.data$date <- lubridate::ymd(paste(diversity.data$year,diversity.data$month,"1"))
      diversity.data$site.web <- paste(diversity.data$site,diversity.data$web,sep=".")
      diversity.data$site <- as.character(diversity.data$site)
      diversity.data$web <- as.character(diversity.data$web)
      
      # reduce diversity data to just these site.webs
      ls <- length(site.webs)
      div.data <- diversity.data[which(diversity.data$site.web == site.webs[1]), ]
      if (ls > 1) {
        for (i in 2:ls) {
          div.data <- rbind(div.data, diversity.data[which(diversity.data$site.web == site.webs[i]), ])
        }
      }
      
      
      datas <- dplyr::left_join(datas, div.data) 
    }
    
    
    cl <- length(cov.list)
    for (c in 1:cl) {
      nam <- paste(names(cov.list)[c], cov.list[[c]], sep = "_")
      col <- which(names(datas) == names(cov.list)[c])
      fd <- which(datas$date == data.w$date[1])
      ld <- which(datas$date == data.w$date[length(data.w$date)])
      #if (length(fd) == 1) {
      #  data.w <- cbind(data.w, 
      #                  datas[(fd - cov.list[[c]]):(ld - cov.list[[c]]), 
      #                         col])
      #  names(data.w)[dim(data.w)[2]] <- nam
      #}
      #if (length(fd) > 1) { ## if not >1 then prob with '.web' below. I think can do this no matter the length(fd), don't need the bit above
      for (i in 1:length(fd)) {
        data.w <- cbind(data.w, 
                        datas[(fd[i] - cov.list[[c]]):(ld[i] - cov.list[[c]]), 
                              col])
        names(data.w)[dim(data.w)[2]] <- paste(nam, 
                                               paste("web", 
                                                     datas$site.web[fd[i]], 
                                                     sep = ""), 
                                               sep = ".")
      }
      #}
    }
    month.data <- dplyr::left_join(month.data, data.w)
    ind.na <- which(is.na(month.data), arr.ind = TRUE)
    
    # these columns have NAs that need to be filled for CJS models to run, 
    # they are likely the first time points when there is a lag, so fill in 
    # the next value (this is after removing the first 7 column because that 
    # is just the session numbers etc. not the covariate data)
    # May need to consider interpolating instead?
    if(remove.na==TRUE){
      ind.na <- ind.na[-which(ind.na[, 2] < 9), ]
      if (dim(ind.na)[1] > 0) {
        for (i in 1:dim(ind.na)[1]) {
          # find the next finite value after this one and plug it in
          vi <- which(is.finite(month.data[, ind.na[i, 2]]))
          vin <- vi[min(which(vi > ind.na[i, 1]))]
          month.data[ind.na[i, 1], ind.na[i, 2]] <- month.data[vin, ind.na[i, 2]]
        }
      }
    }
    #-----------------------------------------------------------
    
  } #if covariate data
  
  ## add first capture to individual covariates data
  
  CH.primary <- primary.ch.fun(CH.secondary)
  first.caught <- apply(CH.primary, 1, function(x) {
    min(which(x > 0))
  }) # gives primary occasion first caught (not long.month)
  ic$f.longmonth <- month.data$long.month[match(first.caught, month.data$Prim)]
  if(multistate==TRUE){
    # for multistate models, need to state at first caught:
    ms.CH.primary <-  primary.MSch.fun(CH.secondary)
    f.state <- numeric()
    for(i in 1:dim(ic)[1]){
      f.state[i] <- ms.CH.primary[i,first.caught[i]]
    }
    ic$f.state <- f.state
  }
  
  individual.covariates <- ic
  covariate.data <- list(individual.covariates = individual.covariates, 
                         temporal.covariates = month.data)

  return(covariate.data)
}





## Function to change secondary capture histories from list to array format-- ##

list.to.array.fun <- function(Ch.list, #as list
                              temporal.covariates){ # dataframe from the list output from covariate function
  
  nind <- dim(Ch.list[[1]])[1]
  months.trapped <- temporal.covariates$long.month[which(is.finite(temporal.covariates$session))]

  nprim <- max(months.trapped)
  n.sec.occ <- unlist(lapply(Ch.list, function(x){dim(x)[2]}))
  y <- array(NA, dim = c(nind, max(n.sec.occ), nprim))
  
  #need a months.trapped variable e.g, c(1,2,5,6)
  for(m in months.trapped){
   y[ , 1:n.sec.occ[temporal.covariates$long.month[which(months.trapped==m)]], m] <- Ch.list[[temporal.covariates$long.month[which(months.trapped==m)]]] 
   
  }
  return(y)
}





## Function to create matrix with info about known latent state z ----------- ##
known.state.cjs=function(ch){
  state=ch
  for(i in 1:dim(ch)[1]){
    n1=min(which(ch[i,]==1))
    n2=max(which(ch[i,]==1))
    state[i,n1:n2]=1
    state[i,n1]=NA			#only filling in those that were 0s but we know were alive because caught before and after
  }
  state[state==0]=NA
  return(state)
}



## Function to create matrix of initial values for latent state z ----------- ##
# we shouldn't give initial values for those elements of z whose value is specified in the data.
# they get an NA
cjs.init.z <- function(ch,f,prim.occasions){
  ch <- replace(ch,is.na(ch),0)
  for(i in 1:dim(ch)[1]){
    if(sum(ch[i,],na.rm=TRUE)==1) next
    n2=max(which(ch[i,]==1))
    ch[i,f[i]:n2]=NA
  }
  for(i in 1:dim(ch)[1]){
    ch[i,1:f[i]]=NA
  }
  return(ch)	
}

