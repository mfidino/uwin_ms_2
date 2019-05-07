source("sourcer.R")

which_folder <- "global"

my_files <- list.files(paste0("./results/", which_folder),
                       pattern = "matrix", 
                       full.names = TRUE)

sp_name <- gsub(".*/(\\w+)_matrix.csv", "\\1", my_files)

# check if a sub-folder exists
f_exist <- dir.exists(paste0("./plots/", sp_name))
for(ex in 1:length(f_exist)){
  if(!f_exist[ex]){
    dir.create(paste0("./plots/", sp_name[ex]))
  }
}

for(species in 1:length(sp_name)){
  # read in the mcmc matrix
  results_matrix <- data.table::fread(my_files[species],
                                      data.table = FALSE) %>% 
                    as.matrix(.)
  
  # psi_mu(population density)
  new_data <- data.frame(hden = 0,
                         prophab_btwn = 0,
                         hden_btwn = seq(300, 1400, 10))
  
  new_city_data <- data.frame(hden = 0,
                              prophab_btwn = 0,
                              hden_btwn = city_data$hden)
  row.names(new_city_data) <- city_data$city
  
  
  test <- predict.city(mmat = results_matrix,
                       new_data = new_data,
                       city_data = cdat,
                       new_city_data = new_city_data,
                       species_there = det_events[,sp_name[species]])
  
  plot(test$mu[,2] ~ new_data$hden_btwn, ylim = c(0,1), type = 'l',
       main = sp_name[species], xlab = 'average housing density', ylab = 'occupancy',
       bty = 'l')
  lines(test$mu[,1] ~ new_data$hden_btwn, lty = 2)
  lines(test$mu[,3] ~ new_data$hden_btwn, lty = 2)
  points(test$cmu[,2] ~ new_city_data$hden_btwn[test$cities], pch = 19)
  where_sp <- test$cities
  for(i in 1:length(where_sp)){
    lines(x = rep(new_city_data$hden_btwn[where_sp[i]], 2), y = test$cmu[i,-2])
  }
}

for(species in 1:length(sp_name)){
  # read in the mcmc matrix
  results_matrix <- data.table::fread(my_files[species],
                                      data.table = FALSE) %>% 
    as.matrix(.)
  
  # psi_mu(population density)
  new_data <- data.frame(hden = 0,
                         prophab_btwn = seq(0.15, 0.65, 0.005),
                         hden_btwn = 0)
  
  new_city_data <- data.frame(hden = 0,
                              prophab_btwn = city_data$habitat,
                              hden_btwn = 0)
  row.names(new_city_data) <- city_data$city
  
  
  test <- predict.city(mmat = results_matrix,
                       new_data = new_data,
                       city_data = cdat,
                       new_city_data = new_city_data,
                       species_there = det_events[,sp_name[species]])
  
  plot(test$mu[,2] ~ new_data$prophab_btwn, ylim = c(0,1), type = 'l',
       main = sp_name[species], xlab = 'Proportion green space', ylab = 'occupancy',
       bty = 'l')
  lines(test$mu[,1] ~ new_data$prophab_btwn, lty = 2)
  lines(test$mu[,3] ~ new_data$prophab_btwn, lty = 2)
  points(test$cmu[,2] ~ new_city_data$prophab_btwn[test$cities], pch = 19)
  where_sp <- test$cities
  for(i in 1:length(where_sp)){
    lines(x = rep(new_city_data$prophab_btwn[where_sp[i]], 2), y = test$cmu[i,-2])
  }
}

  
  results_matrix %>% 
    gen_preds.popdens(., occ_data = det_data, 
                      city_covs = scale_cdat$pop10_dens,
                      species = sp_name[species]) %>% 
    make_plot.popdens(., species = sp_name[species])
  
  # psi_mu(proportion habitat)
  results_matrix %>% 
    gen_preds.habitat(., occ_data = det_data, 
                      city_covs = scale_cdat$habitat,
                      species = sp_name[species]) %>% 
    make_plot.habitat(., species = sp_name[species])
    
  # psi_mu(latitude)
  results_matrix %>% 
    gen_preds.latitude(., occ_data = det_data, 
                       city_covs = scale_cdat$latitude,
                       species = sp_name[species]) %>% 
    make_plot.latitude(., species = sp_name[species])
  
  # psi_mu(intercept)
results_matrix %>%
    gen_preds.nocorrelates(., occ_data = det_data,
                           species = sp_name[species] ) %>% 
  make_plot.nocorrelates(., species = sp_name[species],
                         city_names = cdat$city)
    
  # response to urbanization(population density)
  results_matrix %>% 
    gen_preds.popdens(., occ_data = det_data, 
                      city_covs = scale_cdat$pop10_dens,
                      species = sp_name[species], intercept = FALSE) %>% 
    make_plot.popdens(., species = sp_name[species], intercept = FALSE)
  
  results_matrix %>% 
    gen_preds.habitat(., occ_data = det_data, 
                      city_covs = scale_cdat$habitat,
                      species = sp_name[species], intercept = FALSE) %>% 
    make_plot.habitat(., species = sp_name[species], intercept = FALSE)
  
  results_matrix %>% 
    gen_preds.latitude(., occ_data = det_data, 
                       city_covs = scale_cdat$latitude,
                       species = sp_name[species], intercept = FALSE) %>% 
    make_plot.latitude(., species = sp_name[species], intercept = FALSE)
  
  results_matrix %>%
    gen_preds.nocorrelates(., occ_data = det_data,
                           species = sp_name[species], intercept = FALSE ) %>% 
    make_plot.nocorrelates(., species = sp_name[species], intercept = FALSE,
                           city_names = cdat$city)
  
}

  
  
  m1 <- data.table::fread("./results/raccoon_matrix.csv",
                          data.table = FALSE) %>%
    as.matrix(.) 
  
  
  mguess <- array(m1[,1:16], dim = c(nrow(m1), 8, 2))
  
  dim(mguess)
  
  t0 <- c(1,0)
  t1 <- c(1,1)
  my_res <- matrix(0, ncol = 3, nrow = 8)
  
  for( i in 1:8){
    b0 <- plogis(mguess[,i,] %*% c(1,0))
    b1 <- plogis(mguess[,i,] %*% c(1,1))
    my_res[i, ] <- HDIofMCMC( b1 - b0)
  }
  
  m2 <- data.table::fread("./results/f_squ_matrix.csv",
                          data.table = FALSE) %>%
    as.matrix(.) 
  
  mguess <- array(m2[,1:14], dim = c(nrow(m1), 7, 2))
  
  my_res <- matrix(0, ncol = 3, nrow = 7)
  for( i in 1:7){
    b0 <- plogis(mguess[,i,] %*% c(1,0))
    b1 <- plogis(mguess[,i,] %*% c(1,1))
    my_res[i, ] <- HDIofMCMC( b1 - b0)
  }
  plot(my_res[,2] ~ cdat$pop10_dens[1:7], 
       xlab = "city population density", 
       ylab = "percent change in patch occupancy from 1 unit change in URB", 
       bty = 'l', pch = 16, cex = 1.5, ylim = c(-.20, .20))
  for(i in 1:7){
    lines(x = rep(cdat$pop10_dens[i], 2),
         y = my_res[i,-2])
  }
  