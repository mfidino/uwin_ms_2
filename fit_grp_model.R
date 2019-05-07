
################################
# Load and fit multi-city models
#
# Code written by Mason Fidino
################################

source("sourcer.R")


# Settings for the jags models
to_monitor <- c("B0", "B", "Bmu", "b_2016", "b_2018",
                "D0", "Dmu", "psi_sd", "rho_sd", "z")

my_model <- "./jags_scripts/city_mean_occupancy.R"
nchain <- 6
nadapt <- 10000
nburnin <- 50000
nsample <- 30000
nthin <- 2
my_method <- "parallel"

# get the names of the species we are fitting
my_species <- colnames(has_species)



# Fit global model
model <- "rerun"
for(species in 1:8) {
  print(my_species[species])
  
  data_list <- list(
    y = det_data[,my_species[species]],
    has_species = has_species[,my_species[species]],
    J = det_data$J,
    ncity = ncity,
    nsite = nsite,
    npatch_covs = npatch_covs,
    bx = bx,
    city_vec = city_vec,
    in_2016 = my_2016,
    in_2018 = my_2018
  )
  
  z <- data_list$y
  z[z>1] <- 1
  #if(sum(det_events[,2 + species]) == ncity){
  #  print("model already exists")
  #  model_output <- readRDS(paste0("./results/global/",my_species[species],
  #                                 ".RDS"))
  #  model_waic <- calc_waic(model_output, data_list)
  #  saveRDS(model_waic, paste0("./results/",model,"/", my_species[species], "_waic.RDS"))
  #  rm(model_output)
  #  rm(model_waic)
  #  next
  #}
  
  # make y NA if species is not present in a given city
  if(any(data_list$has_species == 0)){
    data_list$y[which(data_list$has_species == 0)] <- NA
  }
  
  
  
  
  model_output <- run.jags(
    model = my_model,
    data = data_list,
    n.chains = nchain,
    monitor = to_monitor,
    adapt = nadapt,
    burnin = nburnin,  sample = nsample, thin = nthin, method = my_method,
    inits = initsimp, summarise = FALSE, modules = "glm")
  
  saveRDS(model_output, paste0("./results/",model,"/", my_species[species], "_ranefshift.RDS"))
  #model_waic <- calc_waic(model_output, data_list)
  #saveRDS(model_waic, paste0("./results/",model,"/", my_species[species], "_waic.RDS"))
  rm(model_output)
  #rm(model_waic)
  
}
