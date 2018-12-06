model{
  # Latent state model
  for(site in 1:nsite){
    logit(psi[site]) <- inprod(B[city_vec[site], ], bx[site, ])
    z[site] ~ dbern(psi[site])
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
# this is the inverse of the variance covariance matrix
Tau.B[1:npatch_covs, 1:npatch_covs] <- inverse(Sigma.B[,])
#
# standard deviation of the three city wide variables
sigma.int ~ dgamma(1,1)
sigma.urb ~ dgamma(1,1)
sigma.lat ~ dgamma(1,1)
# Fill the diagonal of the matrix with variance terms
Sigma.B[1,1] <- pow(sigma.int, 2)
Sigma.B[2,2] <- pow(sigma.urb, 2)
Sigma.B[3,3] <- pow(sigma.lat, 2)
# For a 3 x 3 matrix we have 3 correlation terms
rho12 ~ dunir(-1, 1)
rho13 ~ dunir(-1, 1)
rho23 ~ dunir(-1, 1)
# fill up the off diagonal of the variance covariance matrix
Sigma.B[1,2] <- rho12 * sigma.int * sigma.urb
Sigma.B[1,3] <- rho13 * sigma.int * sigma.lat
Sigma.B[2,3] <- rho23 * sigma.urb * sigma.lat
Sigma.B[2,1] <- rho12 * sigma.int * sigma.urb
Sigma.B[3,1] <- rho13 * sigma.int * sigma.lat
Sigma.B[3,2] <- rho23 * sigma.urb * sigma.lat
#
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
}
