
################################
# Load and fit multi-city models
#
# Code written by Mason Fidino
################################

source("sourcer.R")


# Settings for the jags models


my_model <- "./jags_scripts/city_mean_occupancy2.R"
nchain <- 5
nadapt <- 60000
nburnin <- 60000
nsample <- 1e5
nthin <- 10
my_method <- "parallel"

# get the names of the species we are fitting
my_species <- colnames(has_species)

to_monitor <- c("B0", "B", "b1", 'b2', "b_2016", "b_2018",
                "D0", "Dmu", "psi_sd", "rho_sd", "B_diff", 'z')


# Fit global model
model <- "reparam_longrun"
for(species in 2:8) {
  print(my_species[species])
  
  data_list <- list(
    y = det_data[,my_species[species]],
    has_species = has_species[,my_species[species]],
    J = det_data$J,
    ncity = ncity,
    nsite = nsite,
    npatch_covs = 3,
    bx = bx[,c(1)],
    city_vec = city_vec,
    in_2016 = my_2016,
    in_2018 = my_2018,
    U=U
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



tmp_mat <- as.matrix(as.mcmc.list(model_output), chains = TRUE)


tmp <- det_data
bxtmp <- bx[,1]

data_list <- list(
  y = tmp[,my_species[species]],
  has_species = has_species[,my_species[species]],
  J = tmp$J,
  ncity = n_distinct(tmp$city),
  nsite = nrow(tmp),
  npatch_covs = 3,
  bx = bxtmp,
  city_vec = as.numeric(factor(tmp$city)),
  in_2016 = my_2016,
  in_2018 = my_2018,
  U = U
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
#if(any(data_list$has_species == 0)){
#  data_list$y[which(data_list$has_species == 0)] <- NA
#}





model_output2 <- run.jags(
  model = my_model,
  data = data_list,
  n.chains = nchain,
  monitor = to_monitor,
  adapt = nadapt,
  burnin = nburnin,  sample = nsample, thin = nthin, method = my_method,
  inits = initsimp, summarise = FALSE, modules = "glm")
