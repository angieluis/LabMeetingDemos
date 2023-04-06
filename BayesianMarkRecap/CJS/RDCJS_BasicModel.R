
######################################################################################
# Specify model in BUGS/JAGS language
######################################################################################


sink("robust_cjs.bug")
cat("					######<--------------------- uncomment 
model{
	
###############Priors and constraints
mean.phi ~ dunif(0, 1) # prior for phi
mean.p ~ dunif(0, 1)   # prior for p
mean.c ~ dunif(0, 1)   # prior for c
    
for(i in 1:nind){
	for(m in f[i]:n.primary.occasions){  ### for p need every time for phi need -1

		# phi has only 2 dimensions [indiv, and primary occasions]
    phi[i,m] <- mean.phi   # could specify covariates here

    # p and c have 3 dimensions [indiv, secondary, primary]		
    for(d in 1:max(n.secondary.occasions)){
      p[i, d, m] <- mean.p  # could specify temporal or individual covariates here
      c[i, d, m] <- mean.c  # ditto
	   } #d for days
    } #m for months
	} #i for individual


#############Likelihood 		
for(i in 1:nind){
	#define latent state at first capture 
  # dimensions [individual, primary session (month)]
	z[i,f[i]] <- 1		# z is true (latent) state alive or dead, know alive at first capture
	
	for(m in  (f[i]+1):n.primary.occasions){  
		
	  #state process				# alive or dead
	  z[i, m] ~ dbern(mu1[i, m]) 		#mu1 is probability alive
		mu1[i, m] <- phi[i, m] * z[i, m-1] # this assures that animals stay dead

	} # m
	
	for(m in prim.occasions[f.prim[i]:length(prim.occasions)] ) { 
	  #observation process			# caught or not
	  #first secondary occasion within a primary occasion:
	  y[i, 1, m] ~ dbern(p.eff[i, 1, m])
	  p.eff[i, 1, m] <- z[i, m] * p[i, 1, m]   
	  
	  #loop over rest of secondary occasions per primary
	  for(d in 2:n.secondary.occasions[m]){ ## add index for number of secondary occasions if not consistent
		  y[i, d, m] ~ dbern(p.eff[i, d, m]) 		# p.eff is prob of capture
		  p.eff[i, d, m] <- z[i, m] * ifelse(sum(y[i, 1:(d-1), m])==0, p[i, d, m], c[i, d, m])	
      # capture prob= (p if not caught previously that session or c if was caught that session) 
      # multiply by if it was alive (so can't capture animal that's not alive)
		 
	   } #d
		} #m
	} #i
}
",fill=TRUE)  #####<----------------uncomment this
sink()


