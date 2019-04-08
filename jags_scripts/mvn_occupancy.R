model{
  # Latent state model
  for(site in 1:nsite){
    logit(psi[site]) <- inprod(B[city_vec[site], ], bx[site, ]) + b_2016[city_vec[site]] + b_2018[city_vec[site]]
    z[site] ~ dbern(psi[site]* has_species[site])

  }
  # Observation model
  for(site in 1:nsite){
    logit(rho[site]) <- inprod(D[city_vec[site],], dx[site, ])
    y[site] ~ dbin(rho[site] * z[site], J[site])
  }
  #
  # Setting up the correlation structure for latent state.
  #  Species intercept and URB slope depend on city specific variables.

for(city in 1:ncity){
  B[city, 1:npatch_covs] ~ dmnorm(B.hat[city, ], Tau.B[, ])
  for(patch_covs in 1:npatch_covs){
    B.hat[city, patch_covs] <- inprod(G[patch_covs, ], U[city, ])
  }
}

for(patch_covs in 1:npatch_covs){
  for(city_covs in 1:ncity_covs){
    G[patch_covs, city_covs]  ~ dlogis(0, 1)
  }
}
  
Tau.B[1:npatch_covs, 1:npatch_covs] <- inverse(Sigma.B[,])

Sigma.B[1,1] <- pow(sigma.int, 2)
sigma.int ~ dgamma(1,1)
Sigma.B[2,2] <- pow(sigma.urb, 2)
sigma.urb ~ dgamma(1,1)
Sigma.B[1,2] <- rho_mat * sigma.int * sigma.urb
Sigma.B[2,1] <- Sigma.B[1,2]
rho_mat ~ dunif(-1, 1)

# Priors for detection model
for(city in 1:ncity){
  for(det_covs in 1:ndet_covs){
  D[city,det_covs] ~ dnorm(rho_mu[det_covs], rho_tau[det_covs])
  }
}
for(det_covs in 1:ndet_covs){
  rho_mu[det_covs] ~ dlogis(0, 1)
  rho_tau[det_covs] ~ dgamma(1, 1)
  rho_sigma[det_covs] <- 1 / sqrt(rho_tau[det_covs])
}

# to accomodate repeat sampling in locations. 
#  Chicago in 2016 (no other cities sampled in 2016)
b_2016[1] <- 0
b_2016[2] ~ dlogis(0,1)
for(more_zed in 3:10){
b_2016[more_zed] <- 0
b_2018[more_zed] <- 0
}
# Chicago in 2018
b_2018[1] <- 0
b_2018[2] ~ dlogis(0,1)


}
