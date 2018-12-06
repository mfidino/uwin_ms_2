package_load<-function(packages = NULL, quiet=TRUE, verbose=FALSE, warn.conflicts=FALSE){
	
	# download required packages if they're not already
	
	pkgsToDownload<- packages[!(packages  %in% installed.packages()[,"Package"])]
	if(length(pkgsToDownload)>0)
		install.packages(pkgsToDownload, repos="http://cran.us.r-project.org", quiet=quiet, verbose=verbose)
	
	# then load them
	for(i in 1:length(packages))
		require(packages[i], character.only=T, quietly=quiet, warn.conflicts=warn.conflicts)
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
      sigma.lat = rgamma(1, 1, 1),
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