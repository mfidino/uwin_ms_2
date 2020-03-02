package_load<-function(packages = NULL, quiet=TRUE, 
                       verbose=FALSE, warn.conflicts=FALSE){
  
  # download required packages if they're not already
  
  pkgsToDownload<- packages[!(packages  %in% installed.packages()[,"Package"])]
  if(length(pkgsToDownload)>0)
    install.packages(pkgsToDownload, repos="http://cran.us.r-project.org", 
                     quiet=quiet, verbose=verbose)
  
  # then load them
  for(i in 1:length(packages))
    require(packages[i], character.only=T, quietly=quiet, 
            warn.conflicts=warn.conflicts)
}

HDIofMCMC = function( sampleVec , credMass=0.95 ) {
  # Computes highest density interval from a sample of representative values,
  #   estimated as shortest credible interval.
  # Arguments:
  #   sampleVec
  #     is a vector of representative values from a probability distribution.
  #   credMass
  #     is a scalar between 0 and 1, indicating the mass within the credible
  #     interval that is to be estimated.
  # Value:
  #   HDIlim is a vector containing the limits of the HDI
  sortedPts = sort( sampleVec )
  ciIdxInc = ceiling( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin ,median(sampleVec), HDImax )
  return( HDIlim )
}


inits <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = z,
      B = array(rnorm(ncity * npatch_covs), 
                dim = c(ncity, npatch_covs)),
      G = array(rnorm(ncity_covs * npatch_covs), 
                dim = c(npatch_covs, ncity_covs)),
      sigma.int = rgamma(1, 1, 1),
      sigma.urb = rgamma(1, 1, 1),
      rho_cor = runif(1, -1, 1),
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Marsaglia-Multicarry",
                         "3" = "base::Super-Duper",
                         "4" = "base::Mersenne-Twister",
                         "5" = "base::Wichmann-Hill",
                         "6" = "base::Marsaglia-Multicarry",
                         "7" = "base::Super-Duper",
                         "8" = "base::Mersenne-Twister"),
      .RNG.seed = sample(1:1e+06, 1)
    )
  }
  return(switch(chain,           
                "1" = gen_list(chain),
                "2" = gen_list(chain),
                "3" = gen_list(chain),
                "4" = gen_list(chain),
                "5" = gen_list(chain),
                "6" = gen_list(chain),
                "7" = gen_list(chain),
                "8" = gen_list(chain)
  )
  )
}


initsimp <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = z,
     # B = rnorm(data_list$ncity,-1),
      Bmu = rnorm(1, -1),
      B_diff = rnorm(data_list$ncity, -1),
      D0 = rnorm(data_list$ncity, -2),
      Dmu = rnorm(1, -2),
      psi_tau = rgamma(1,1),
      rho_tau = rgamma(1,1),
      alpha = runif(1,.2,.8),
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Marsaglia-Multicarry",
                         "3" = "base::Super-Duper",
                         "4" = "base::Mersenne-Twister",
                         "5" = "base::Wichmann-Hill",
                         "6" = "base::Marsaglia-Multicarry",
                         "7" = "base::Super-Duper",
                         "8" = "base::Mersenne-Twister"),
      .RNG.seed = sample(1:1e+06, 1)
    )
  }
  return(switch(chain,           
                "1" = gen_list(chain),
                "2" = gen_list(chain),
                "3" = gen_list(chain),
                "4" = gen_list(chain),
                "5" = gen_list(chain),
                "6" = gen_list(chain),
                "7" = gen_list(chain),
                "8" = gen_list(chain)
  )
  )
}


calc_waic <- function(posterior = NULL, 
                      data = NULL ){
  cat("Calculating WAIC\n")
  with(data,{
    chains <- length(posterior$mcmc)
    store <- dim(posterior$mcmc[[1]])[1]
    Nobs <- length(y)
    L <- array(dim=c(Nobs, chains, store))
    L_bar <- rep(NA, Nobs)
    var_LL <- rep(NA, Nobs)
    pb <- txtProgressBar(min = 0, max = Nobs, style = 3)
    for (i in 1:Nobs){
      for (j in 1:chains){
        post_sims <- posterior$mcmc[[j]]
        indx_z <- which(dimnames(post_sims)[[2]] == 
                          paste0("z[", i, "]"))
        zvals <- post_sims[, indx_z]
        which_city <- city_vec[i]
        
        indx_p <- which(dimnames(post_sims)[[2]] == 
                          paste0("D[", city_vec[i],",1", "]"))
        pvals <- plogis(post_sims[, indx_p])
        L[i, j, ] <- dbinom(rep(y[i], store), 
                            size=J[i], 
                            prob = zvals * pvals, 
                            log=TRUE)
        if(is.na(has_species[i])){
          L[i, j, ] <- NA
        }
        setTxtProgressBar(pb, i)
      }
      L_bar[i] <- mean(exp(c(L[i, , ])), na.rm = TRUE)
      var_LL[i] <- var(c(L[i, , ]), na.rm = TRUE)
    }
    
    lppd <- sum(log(L_bar), na.rm = TRUE)
    p_WAIC <- sum(var_LL, na.rm = TRUE)
    WAIC <- -2 * (lppd - p_WAIC)
    return(list(lppd=lppd, p_WAIC=p_WAIC, WAIC=WAIC, L_bar))
  })
}


calc_cpo <- function(posterior = NULL ,data = NULL){
  chains <- length(posterior$mcmc)
  store <- dim(posterior$mcmc[[1]])[1]
  Nobs <- length(data$y)
  lik <- as.matrix(posterior$mcmc, chains = TRUE)
  lik <- lik[,grep("lik", colnames(lik))]
  cpo_vec <- apply(1/lik, 2, sum)
  cpo_vec <- nrow(lik) / cpo_vec
  cpo_vec[which(data$has_species == 0)] <- NA
  cpo <- -sum(log(cpo_vec), na.rm = TRUE)
  return(cpo)
}
