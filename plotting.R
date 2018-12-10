
which_folder <- "full_model"

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
  