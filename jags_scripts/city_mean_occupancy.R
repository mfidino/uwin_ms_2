model{
  # Latent state model
  for(site in 1:nsite){
    logit(psi[site]) <- B0[city_vec[site]] + B[city_vec[site]]* bx[site] + 
      b_2016[city_vec[site]] * in_2016[site] + 
      b_2018[city_vec[site]] * in_2018[site] + B_diff[city_vec[site]]
    z[site] ~ dbern(psi[site]* has_species[site])
    
  }
  # Observation model
  for(site in 1:nsite){
    logit(rho[site]) <- D0[city_vec[site]]
    y[site] ~ dbin(rho[site] * z[site], J[site])
  }
  # priors
  for(city in 1:ncity){
    B0[city] <- inprod(b1,U[city,])
    B[city] <- inprod(b2, U[city,])
    D0[city] ~ dnorm(Dmu, rho_tau)
    B_diff[city] ~ dnorm(0, psi_tau)
  }
  #Bmu ~ dlogis(0,1)
  Dmu ~ dlogis(0,1)
  psi_tau ~ dgamma(0.001,0.001)
  psi_sd <- 1 / sqrt(psi_tau )
  rho_tau ~ dgamma(0.001,0.001)
  rho_sd <- 1 / sqrt(rho_tau)
  
  for(pars in 1:npatch_covs){
    b1[pars] ~ dlogis(0,1)
    b2[pars] ~ dlogis(0,1) 
  }
  # to accomodate repeat sampling in cities. 
  #  Chicago and Manhattan in 2016
  b_2016[1] <- 0
  b_2016[2] ~ dlogis(0,1) # Chicago
  b_2016[3] <- 0
  b_2016[4] <- 0
  b_2016[5] <- 0
  b_2016[6] <- 0
  b_2016[7] <- 0
  b_2016[8] ~ dlogis(0,1) # Manhattan
  b_2016[9] <- 0
  b_2016[10] <- 0
  # Repeat sampling in 2018
  b_2018[1] <- 0
  b_2018[2] ~ dlogis(0,1) # Chicago
  b_2018[3] <- 0
  b_2018[4] <- 0
  b_2018[5] ~ dlogis(0,1) # Iowa city
  b_2018[6] ~ dlogis(0,1) # Indianapolis
  b_2018[7] <- 0
  b_2018[8] <- 0
  b_2018[9] ~ dlogis(0,1) # Madison
  b_2018[10] <- 0
  # code for likelihood, used to calculate CPO.
  
  for(site in 1:nsite){
    # Trick to calculate the binomial coefficient.
  BinCo[site] <- exp(logfact(J[site]) - (logfact(y[site]) + 
                                       logfact(J[site] - y[site])))
  # the likelihood of our data given the model
  lik[site] <- ifelse(equals(y[site],0), 
                psi[site]*((1-rho[site])^J[site]) + (1-psi[site]),
                BinCo[site]*psi[site]*(rho[site]^y[site]) * 
                  (1-rho[site])^(J[site]-y[site]))
  }
  
  
}
