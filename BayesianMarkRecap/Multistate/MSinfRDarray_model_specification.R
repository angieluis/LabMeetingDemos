## Multistate Infection Robust Design Model specification 


## To do: ----------------------------------------------------------##
## -----------------------------------------------------------------##


 # specify model
   sink("MSinfRD_specification.bug")
   
   cat("
     
     model {
       
    # -------------------------------------------------
       # States (S):
       # 1 alive as S
       # 2 alive as I
       # 3 dead

       # Observations (O):  
       # 1 seen as S 
       # 2 seen as I
       # 3 not seen
       # -------------------------------------------------



       ##### PRIORS AND CONSTRAINTS #####
      
      
      ##### PRIORS FOR PHI #####
      alpha.0          ~ dnorm(0, 0.4)T(-10, 10)    # prior for survival intercept
      alpha.ndvi       ~ dnorm(0, 0.4)T(-10, 10)    # prior for survival coef on ndvi12
      alpha.inf        ~ dnorm(0, 0.4)T(-10, 10)    # prior for infection on survival
      alpha.inf.male   ~ dnorm(0, 0.4)T(-10, 10)    # prior for infection males on survival


      ## seasonal covariates
      # don't estimate a beta coefficient for winter (assume 0),
      # then estimate a coefficent for all other seasons
      for(m in 1:3) {
        alpha.season[m] ~ dnorm(0, 0.4)T(-10, 10)   # prior for coef on season
       }
      alpha.season.use <- c(0, alpha.season) 

       
      ##### PRIORS FOR PSI #####
      
      beta.0    ~ dnorm(0, 0.4)T(-10, 10)
      beta.male ~ dnorm(0, 0.4)T(-10, 10)
      # beta.I  ~ dnorm(0, 1)T(-10, 10)

	  
      ##### PRIORS FOR RECAPTURE #####
      sigma.0     ~ dnorm(0, 0.4)T(-10, 10)  # prior for intercept on p
      sigma.inf   ~ dnorm(0, 0.4)T(-10, 10)  # prior for infection on p
 

        
        ##### MODEL FOR PHI #####
        for(i in 1:nind) {
          for(m in 1:n.months) {              #I'm defining more Phis and Psis than I will use, since specifying last time step
            # phi has  2 dimensions [indiv, and months]
            
            ### Phi for uninfected
            logit(phiS[i, m]) <-            
              
              alpha.0 +
              alpha.season.use[season[m]] +       

              # environmental covariates
              alpha.ndvi * ndvi[m] 

            ### Phi for infected
            logit(phiI[i, m]) <-              
              
              alpha.0 +
              alpha.season.use[season[m]] +       

              # environmental covariates
              alpha.ndvi * ndvi[m] +

              # infection terms
              alpha.inf +
              alpha.inf.male * sex[i]
       

        
      ##### MODEL FOR PSI #####
        ##### 2 dimensions [indiv, month]

            logit(psiSI[i, m]) <- beta.0 + beta.male * sex[i]  
                                  # + beta.I * I.dat[m]  
  
  
  
        ##### MODEL FOR P #####
        ##### here doesn't vary by day so has 2 dimensions [indiv, month, day] (but coul have day as well) 
 
              logit(pS[i, m]) <-
                sigma.0 


              logit(pI[i, m]) <-
                sigma.0 +                       
                sigma.inf 
       


             
            } # m for month
          } # i for individual
        
      
      ################## Define state-transition and observation matrices
      for (i in 1:nind){  
        # Define probabilities of State(m+1) given State(m)
        for (m in 1:n.months){ 
        ps[1,i,m,1] <- phiS[i,m] * (1-psiSI[i,m]) # survival of S to S
        ps[1,i,m,2] <- phiS[i,m] * psiSI[i,m]     # survival and transition from S to I
        ps[1,i,m,3] <- 1-phiS[i,m]                # S to dead
        ps[2,i,m,1] <- 0                          # I to S (can't happen, tho check data)
        ps[2,i,m,2] <- phiI[i,m]                  # survival of I to I
        ps[2,i,m,3] <- 1-phiI[i,m]                # I to dead
        ps[3,i,m,1] <- 0                          # dead to S
        ps[3,i,m,2] <- 0                          # dead to I
        ps[3,i,m,3] <- 1                          # dead stay dead


        for(d in 1:n.sec.occ[m]) { 
          # Define probabilities of Observation(m) given State(m)
          # first index is state, last index is observation
          # could potentially include observation as wrong state (false neg or pos)

          po[1,i,m,d,1] <- pS[i,m]          # in S and observed as S
          po[1,i,m,d,2] <- 0                # in S and observed as I 
          po[1,i,m,d,3] <- 1-pS[i,m]        # in S and not observed 
          po[2,i,m,d,1] <- 0                # in I and observed as S
          po[2,i,m,d,2] <- pI[i,m]          # in I and observed as I
          po[2,i,m,d,3] <- 1-pI[i,m]        # in I and not observed
          po[3,i,m,d,1] <- 0                # dead and observed as S
          po[3,i,m,d,2] <- 0                # dead and observed as I
          po[3,i,m,d,3] <- 1                # dead and not observed
         } #d
        } #m
       } #i
       

        ##### LIKELIHOOD #####
        
        ##### STATE PROCESS
        for(i in 1:nind) {
          
          # define latent state at first capture 
          # dimensions [individual, month] spans study period
          # z is true (latent) state, know state at first capture
          z[i,f[i]] <- f.state[i]    
          
          for(m in (f[i] + 1):n.months) {  ### can go backwards if only caught last day?
            z[i, m] ~ dcat(ps[z[i, m-1], i, m-1, ])
          } # m 
        
        
        
          ##### OBSERVATION PROCESS
          for(m in prim.occasions[f.prim[i]:length(prim.occasions)] ) {  
            for(d in 1:n.sec.occ[m]){
              
              y[i, d, m] ~ dcat(po[z[i, m],  i, m, d, ])
                   
            } #d
          } #m
        } #i
  
        
        
        
        } # model
    ", fill = TRUE)
    
    sink()
     
  
## -------------------------------------------------------------------------- ##