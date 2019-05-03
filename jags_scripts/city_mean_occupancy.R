model{
  # Latent state model
  for(site in 1:nsite){
    logit(psi[site]) <- Bmu + inprod(B, bx[site,]) + 
      b_2016[city_vec[site]] * in_2016[site] + 
      b_2018[city_vec[site]] * in_2018[site] + B0[city_vec[site]]
    z[site] ~ dbern(psi[site]* has_species[site])
    
  }
  # Observation model
  for(site in 1:nsite){
    logit(rho[site]) <- Dmu + D0[city_vec[site]]
    y[site] ~ dbin(rho[site] * z[site], J[site])
  }
  # priors
  for(city in 1:ncity){
    B0[city] ~ dnorm(0, psi_tau)
    D0[city] ~ dnorm(0, rho_tau)
  }
  Bmu ~ dlogis(0,1)
  Dmu ~ dlogis(0,1)
  psi_tau ~ dgamma(1,1)
  psi_sd <- 1 / sqrt(psi_tau)
  rho_tau ~ dgamma(1,1)
  rho_sd <- 1 / sqrt(rho_tau)
  
  for(pars in 1:npatch_covs){
    B[pars] ~ dlogis(0,1)
  }
  # to accomodate repeat sampling in locations. 
  #  Chicago in 2016 (no other cities sampled in 2016)
  b_2016[1] <- 0
  b_2016[2] ~ dlogis(0,1)
  for(more_zed in 3:10){
    b_2016[more_zed] <- 0
  }
  # Chicago in 2018
  b_2018[1] <- 0
  b_2018[2] ~ dlogis(0,1)
  b_2018[3] <- 0
  b_2018[4] <- 0
  b_2018[5] ~ dlogis(0,1)
  b_2018[6] <- 0
  b_2018[7] <- 0
  b_2018[8] <- 0
  b_2018[9] ~ dlogis(0,1)
  b_2018[10] <- 0
  
}
