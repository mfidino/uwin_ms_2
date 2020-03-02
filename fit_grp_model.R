
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
nsample <- 2.5e4
nthin <- 10
my_method <- "parallel"


# get the names of the species we are fitting
my_species <- colnames(has_species)

to_monitor <- c("B0", "B", "b1", 'b2', "b_2016", "b_2018",
                "D0", "Dmu", "psi_sd", "rho_sd", "B_diff", 'z', 'lik')


# Fit habitat model
model <- "habitat"
to_drop <- 3
my_cpo <- data.frame(species = my_species, 
                     cpo = NA)
for(species in 1:8) {
  print(my_species[species])
  if(!is.na(to_drop)){
  my_U <- U[,-to_drop]
  } else {
    my_U <- U
  }
  
  data_list <- list(
    y = det_data[,my_species[species]],
    has_species = has_species[,my_species[species]],
    J = det_data$J,
    ncity = ncity,
    nsite = nsite,
    npatch_covs = ncol(my_U),
    bx = bx[,c(1)],
    city_vec = city_vec,
    in_2016 = my_2016,
    in_2018 = my_2018,
    U=my_U
  )
  
  z <- data_list$y
  z[z>1] <- 1
  #if(sum(det_events[,2 + species]) == ncity){
  #  print("model already exists")
  #  model_output <- readRDS(paste0("./results/global/",my_species[species],
  #                                 ".RDS"))

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
  
  my_cpo$cpo[species] <- calc_cpo(model_output, data_list)
  #model_cpo <- calc_cpo(model_output, data_list)
  #saveRDS(model_cpo, paste0("./results/",model,"/", my_species[species], "_cpo.RDS"))
  write.csv(my_cpo, paste0("./results/",model,"/","cpo.csv"))
  saveRDS(model_output, paste0("./results/",model,"/", my_species[species], ".RDS"))
  #model_waic <- calc_waic(model_output, data_list)
  #saveRDS(model_waic, paste0("./results/",model,"/", my_species[species], "_waic.RDS"))
  rm(model_output)
  #rm(model_waic)
  
}



# Fit global model
model <- "global"
to_drop <- NA
my_cpo <- data.frame(species = my_species, 
                     cpo = NA)
for(species in 1:8) {
  print(my_species[species])
  if(!is.na(to_drop)){
    my_U <- U[,-to_drop]
  } else {
    my_U <- U
  }
  
  data_list <- list(
    y = det_data[,my_species[species]],
    has_species = has_species[,my_species[species]],
    J = det_data$J,
    ncity = ncity,
    nsite = nsite,
    npatch_covs = ncol(my_U),
    bx = bx[,c(1)],
    city_vec = city_vec,
    in_2016 = my_2016,
    in_2018 = my_2018,
    U=my_U
  )
  
  z <- data_list$y
  z[z>1] <- 1
  #if(sum(det_events[,2 + species]) == ncity){
  #  print("model already exists")
  #  model_output <- readRDS(paste0("./results/global/",my_species[species],
  #                                 ".RDS"))
  
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
  
  my_cpo$cpo[species] <- calc_cpo(model_output, data_list)
  #model_cpo <- calc_cpo(model_output, data_list)
  #saveRDS(model_cpo, paste0("./results/",model,"/", my_species[species], "_cpo.RDS"))
  write.csv(my_cpo, paste0("./results/",model,"/","cpo.csv"))
  saveRDS(model_output, paste0("./results/",model,"/", my_species[species], ".RDS"))
  #model_waic <- calc_waic(model_output, data_list)
  #saveRDS(model_waic, paste0("./results/",model,"/", my_species[species], "_waic.RDS"))
  rm(model_output)
  #rm(model_waic)
  
}

# Fit global model
model <- "housing_density"
to_drop <- 2
my_cpo <- data.frame(species = my_species, 
                     cpo = NA)
for(species in 1:8) {
  print(my_species[species])
  if(!is.na(to_drop)){
    my_U <- U[,-to_drop]
  } else {
    my_U <- U
  }
  
  data_list <- list(
    y = det_data[,my_species[species]],
    has_species = has_species[,my_species[species]],
    J = det_data$J,
    ncity = ncity,
    nsite = nsite,
    npatch_covs = ncol(my_U),
    bx = bx[,c(1)],
    city_vec = city_vec,
    in_2016 = my_2016,
    in_2018 = my_2018,
    U=my_U
  )
  
  z <- data_list$y
  z[z>1] <- 1
  #if(sum(det_events[,2 + species]) == ncity){
  #  print("model already exists")
  #  model_output <- readRDS(paste0("./results/global/",my_species[species],
  #                                 ".RDS"))
  
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
  
  my_cpo$cpo[species] <- calc_cpo(model_output, data_list)
  #model_cpo <- calc_cpo(model_output, data_list)
  #saveRDS(model_cpo, paste0("./results/",model,"/", my_species[species], "_cpo.RDS"))
  write.csv(my_cpo, paste0("./results/",model,"/","cpo.csv"))
  saveRDS(model_output, paste0("./results/",model,"/", my_species[species], ".RDS"))
  #model_waic <- calc_waic(model_output, data_list)
  #saveRDS(model_waic, paste0("./results/",model,"/", my_species[species], "_waic.RDS"))
  rm(model_output)
  #rm(model_waic)
  
}

# Fit global model
model <- "null"
to_drop <- c(2,3)
my_cpo <- data.frame(species = my_species, 
                     cpo = NA)
for(species in 1:8) {
  print(my_species[species])
  if(!any(is.na(to_drop))){
    my_U <- U[,-to_drop]
    if(is.vector(my_U)){
      my_U <- matrix(my_U, ncol = 1, nrow = length(my_U))
    }
  } else {
    my_U <- U
  }
  
  data_list <- list(
    y = det_data[,my_species[species]],
    has_species = has_species[,my_species[species]],
    J = det_data$J,
    ncity = ncity,
    nsite = nsite,
    npatch_covs = ncol(my_U),
    bx = bx[,c(1)],
    city_vec = city_vec,
    in_2016 = my_2016,
    in_2018 = my_2018,
    U=my_U
  )
  
  z <- data_list$y
  z[z>1] <- 1
  #if(sum(det_events[,2 + species]) == ncity){
  #  print("model already exists")
  #  model_output <- readRDS(paste0("./results/global/",my_species[species],
  #                                 ".RDS"))
  
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
  
  my_cpo$cpo[species] <- calc_cpo(model_output, data_list)
  #model_cpo <- calc_cpo(model_output, data_list)
  #saveRDS(model_cpo, paste0("./results/",model,"/", my_species[species], "_cpo.RDS"))
  write.csv(my_cpo, paste0("./results/",model,"/","cpo.csv"))
  saveRDS(model_output, paste0("./results/",model,"/", my_species[species], ".RDS"))
  #model_waic <- calc_waic(model_output, data_list)
  #saveRDS(model_waic, paste0("./results/",model,"/", my_species[species], "_waic.RDS"))
  rm(model_output)
  #rm(model_waic)
  
}




my_rds <- list.files("./results/housing_density/", "RDS", full.names = TRUE)

model <- "habitat"
to_drop <- 2
for(species in 1:length(my_rds)){
  print(my_species[species])
  if(!is.na(to_drop)){
    my_U <- U[,-to_drop]
  } else {
    my_U <- U
  }
  
  data_list <- list(
    y = det_data[,my_species[species]],
    has_species = has_species[,my_species[species]],
    J = det_data$J,
    ncity = ncity,
    nsite = nsite,
    npatch_covs = ncol(my_U),
    bx = bx[,c(1)],
    city_vec = city_vec,
    in_2016 = my_2016,
    in_2018 = my_2018,
    U=my_U
  )
  
  tmp_rds <- my_rds[grep(my_species[species], my_rds)]
  tmp <- readRDS(tmp_rds)
  my_cpo$cpo[species] <- calc_cpo(model_output, data_list)
  #model_cpo <- calc_cpo(model_output, data_list)
  #saveRDS(model_cpo, paste0("./results/",model,"/", my_species[species], "_cpo.RDS"))
  write.csv(my_cpo, paste0("./results/",model,"/", my_species[species], "_cpo.csv"))
}
habitat_results <- my_cpo

global_waic <- list.files("./results/global/", "waic.RDS", full.names = TRUE )
habitat_waic <- list.files("./results/habitat/","waic.RDS", full.names = TRUE)

globals <- sapply(global_waic, function(x) readRDS(x)$WAIC)
habitat <- sapply(habitat_waic, function(x) readRDS(x)$WAIC)
cbind(globals, habitat)

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
