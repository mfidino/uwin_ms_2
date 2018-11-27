
rac_res <- as.matrix(as.mcmc.list(raccoon), chains = TRUE)
r1 <- as.mcmc.list(raccoon)
diagMCMC(r1,"rho_mat")

t1 <- apply(rac_res, 2, HDIofMCMC)

plot(rac_res[,3] ~ rac_res[,10])



##########
# RACCOON
##########

# Raccoon intercept population density
data.table::fread("./results/raccoon_matrix.csv",
                  data.table = FALSE) %>%
  as.matrix(.) %>% 
 gen_preds.popdens(., occ_data = dat, 
                    city_covs = scale_cdat$pop10_dens,
                   species = "raccoon") %>% 
  make_plot.popdens(., species = "raccoon")
  
# Raccoon intercept habitat
data.table::fread("./results/raccoon_matrix.csv",
                                   data.table = FALSE) %>% 
  as.matrix(.) %>% 
  gen_preds.habitat(., occ_data = dat, 
                    city_covs = scale_cdat$habitat,
                    species = "raccoon") %>% 
make_plot.habitat(., species = "raccoon")

# Raccoon slope popdens
data.table::fread("./results/raccoon_matrix.csv",
                  data.table = FALSE) %>% 
  as.matrix(.) %>% 
  gen_preds.popdens(., occ_data = dat, 
                    city_covs = scale_cdat$pop10_dens,
                    species = "raccoon", intercept = FALSE) %>% 
make_plot.popdens(., species = "raccoon", intercept = FALSE)

# Racoon slope habitat
data.table::fread("./results/raccoon_matrix.csv",
                  data.table = FALSE) %>% 
  as.matrix(.) %>% 
  gen_preds.habitat(., occ_data = dat, 
                    city_covs = scale_cdat$habitat,
                    species = "raccoon", intercept = FALSE) %>% 
make_plot.habitat(., species = "raccoon", intercept = FALSE)

#########
# COYOTE
#########
  
# Coyote intercept popdens
data.table::fread("./results/coyote_matrix.csv",
                  data.table = FALSE) %>% 
  as.matrix(.) %>% 
    gen_preds.popdens(., occ_data = dat, 
                      city_covs = scale_cdat$pop10_dens,
                      species = "coyote") %>% 
  make_plot.popdens(., species = "coyote")
  
# Coyote intercept habitat
data.table::fread("./results/coyote_matrix.csv",
                  data.table = FALSE) %>% 
  as.matrix(.) %>% 
    gen_preds.habitat(., occ_data = dat, 
                      city_covs = scale_cdat$habitat,
                      species = "coyote") %>% 
  make_plot.habitat(., species = "coyote")
  
# Coyote slope popdens
data.table::fread("./results/coyote_matrix.csv",
                  data.table = FALSE) %>% 
  as.matrix(.) %>% 
  gen_preds.popdens(., occ_data = dat, 
                    city_covs = scale_cdat$pop10_dens,
                    species = "coyote", intercept = FALSE) %>% 
make_plot.popdens(., species = "coyote", intercept = FALSE)

# Coyote slope habitat
data.table::fread("./results/coyote_matrix.csv",
                  data.table = FALSE) %>% 
  as.matrix(.) %>% 
  gen_preds.habitat(., occ_data = dat, 
                    city_covs = scale_cdat$habitat,
                    species = "coyote", intercept = FALSE) %>% 
make_plot.habitat(., species = "coyote", intercept = FALSE)

#######
# SKUNK
#######

# skunk intercept popdens
data.table::fread("./results/skunk_matrix.csv",
                  data.table = FALSE) %>% 
  as.matrix(.) %>% 
    gen_preds.popdens(., occ_data = dat, 
                      city_covs = scale_cdat$pop10_dens,
                      species = "skunk") %>% 
  make_plot.popdens(., species = "skunk")

# skunk intercept habitat
data.table::fread("./results/skunk_matrix.csv",
                  data.table = FALSE) %>% 
  as.matrix(.) %>% 
    gen_preds.habitat(., occ_data = dat, 
                      city_covs = scale_cdat$habitat,
                      species = "skunk") %>% 
  make_plot.habitat(., species = "skunk")

# skunk slope popdens
data.table::fread("./results/skunk_matrix.csv",
                  data.table = FALSE) %>% 
  as.matrix(.) %>% 
  gen_preds.popdens(., occ_data = dat, 
                    city_covs = scale_cdat$pop10_dens,
                    species = "skunk", intercept = FALSE) %>% 
make_plot.popdens(., species = "skunk", intercept = FALSE)

# Skunk slope habitat
data.table::fread("./results/skunk_matrix.csv",
                  data.table = FALSE) %>% 
  as.matrix(.) %>% 
  gen_preds.habitat(., occ_data = dat, 
                    city_covs = scale_cdat$habitat,
                    species = "skunk", intercept = FALSE) %>% 
make_plot.habitat(., species = "skunk", intercept = FALSE)
  
#########
# Red fox
#########

# fox intercept popdens
data.table::fread("./results/fox_matrix.csv",
                  data.table = FALSE) %>% 
  as.matrix(.) %>% 
    gen_preds.popdens(., occ_data = dat, 
                      city_covs = scale_cdat$pop10_dens,
                      species = "redfox") %>% 
  make_plot.popdens(., species = "redfox")

# fox intercept habitat
  data.table::fread("./results/fox_matrix.csv",
                    data.table = FALSE) %>% 
    as.matrix(.) %>% 
    gen_preds.habitat(., occ_data = dat, 
                      city_covs = scale_cdat$habitat,
                      species = "redfox") %>% 
  make_plot.habitat(., species = "redfox")

# fox slope popdens
  data.table::fread("./results/fox_matrix.csv",
                    data.table = FALSE) %>% 
    as.matrix(.) %>% 
  gen_preds.popdens(., occ_data = dat, 
                    city_covs = scale_cdat$pop10_dens,
                    species = "redfox", intercept = FALSE) %>% 
make_plot.popdens(., species = "redfox", intercept = FALSE)

# fox slope habitat
data.table::fread("./results/fox_matrix.csv",
                  data.table = FALSE) %>% 
  as.matrix(.) %>% 
  gen_preds.habitat(., occ_data = dat, 
                    city_covs = scale_cdat$habitat,
                    species = "redfox", intercept = FALSE) %>% 
make_plot.habitat(., species = "redfox", intercept = FALSE)

####################
# eastern cottontail
####################

# rabbit intercept popdens
data.table::fread("./results/rabbit_matrix.csv",
                  data.table = FALSE) %>% 
  as.matrix(.) %>% 
    gen_preds.popdens(., occ_data = dat, 
                      city_covs = scale_cdat$pop10_dens,
                      species = "rabbit") %>% 
  make_plot.popdens(., species = "rabbit")
  
# rabbit intercept habitat
data.table::fread("./results/rabbit_matrix.csv",
                  data.table = FALSE) %>% 
  as.matrix(.) %>% 
    gen_preds.habitat(., occ_data = dat, 
                      city_covs = scale_cdat$habitat,
                      species = "rabbit") %>% 
  make_plot.habitat(., species = "rabbit")
  
# rabbit slope popdens
data.table::fread("./results/rabbit_matrix.csv",
                  data.table = FALSE) %>% 
  as.matrix(.) %>% 
  gen_preds.popdens(., occ_data = dat, 
                    city_covs = scale_cdat$pop10_dens,
                    species = "rabbit", intercept = FALSE) %>% 
make_plot.popdens(., species = "rabbit", intercept = FALSE)

# rabbit slope habitat
data.table::fread("./results/rabbit_matrix.csv",
                  data.table = FALSE) %>% 
  as.matrix(.) %>% 
  gen_preds.habitat(., occ_data = dat, 
                    city_covs = scale_cdat$habitat,
                    species = "rabbit", intercept = FALSE) %>% 
make_plot.habitat(., species = "rabbit", intercept = FALSE)
  
#########
# opossum
#########

# opossum intercept popdens
aa <- scale_cdat$pop10_dens[-c(3,4)]
attributes(aa) <- list(dim = c(6,1), 
                        `scaled:center` = as.numeric(attributes(scale_cdat$pop10_dens)[2]),
                        `scaled:scale`=   as.numeric(attributes(scale_cdat$pop10_dens)[3]))
 
data.table::fread("./results/opossum_5city_matrix.csv",
                  data.table = FALSE) %>% 
  as.matrix(.) %>% 
    gen_preds.popdens(., occ_data = dat[-grep("foco|deco", dat$city),], 
                      city_covs = aa,
                      species = "opossum") %>% 
  make_plot.popdens(., species = "opossum5")
  
# opossum intercept habitat
  bb <- scale_cdat$habitat[-c(3:4)]
  attributes(bb) <- list(dim = c(6,1), 
                         `scaled:center` = as.numeric(attributes(scale_cdat$habitat)[2]),
                         `scaled:scale`=   as.numeric(attributes(scale_cdat$habitat)[3]))
  data.table::fread("./results/opossum_5city_matrix.csv",
                    data.table = FALSE) %>% 
    as.matrix(.) %>% 
    gen_preds.habitat(., occ_data = dat[-grep("foco|deco", dat$city),], 
                      city_covs = bb,
                      species = "opossum") %>% 
  make_plot.habitat(., species = "opossum5")
  
  # opossum slope popdens
  data.table::fread("./results/opossum_5city_matrix.csv",
                    data.table = FALSE) %>% 
    as.matrix(.) %>% 
    gen_preds.popdens(., occ_data = dat[-grep("foco|deco", dat$city),], 
                      city_covs = aa,
                      species = "opossum", intercept = FALSE) %>% 
  make_plot.popdens(., species = "opossum5", intercept = FALSE)
  
  # opossum slope habitat
  data.table::fread("./results/opossum_5city_matrix.csv",
                    data.table = FALSE) %>% 
    as.matrix(.) %>% 
    gen_preds.habitat(., occ_data = dat[-grep("foco|deco", dat$city),], 
                      city_covs = bb,
                      species = "opossum", intercept = FALSE) %>% 
  make_plot.habitat(., species = "opossum5", intercept = FALSE)
  
  ##############
  # fox squirrel
  ##############
  aa <- scale_cdat$pop10_dens[-c(8)]
  attributes(aa) <- list(dim = c(7,1), 
                         `scaled:center` = as.numeric(attributes(scale_cdat$pop10_dens)[2]),
                         `scaled:scale`=   as.numeric(attributes(scale_cdat$pop10_dens)[3]))
  
  # fox squirrel intercept popdens
  data.table::fread("./results/f_squ_matrix.csv",
                    data.table = FALSE) %>% 
    as.matrix(.) %>% 
    gen_preds.popdens(., occ_data = dat, 
                      city_covs = aa,
                      species = "foxsquirrel") %>% 
  make_plot.popdens(., species = "foxsquirrel")
  
  # fox squirrel intercept habitat
  
  bb <- scale_cdat$habitat[-8]
  attributes(bb) <- list(dim = c(7,1), 
                         `scaled:center` = as.numeric(attributes(scale_cdat$habitat)[2]),
                         `scaled:scale`=   as.numeric(attributes(scale_cdat$habitat)[3]))
  
  data.table::fread("./results/f_squ_matrix.csv",
                    data.table = FALSE) %>% 
    as.matrix(.) %>% 
    gen_preds.habitat(., occ_data = dat, 
                      city_covs = bb,
                      species = "foxsquirrel") %>% 
  make_plot.habitat(., species = "foxsquirrel")
  
  # fox squirrel slope popdens
  data.table::fread("./results/f_squ_matrix.csv",
                    data.table = FALSE) %>% 
    as.matrix(.) %>% 
    gen_preds.popdens(., occ_data = dat, 
                      city_covs = aa,
                      species = "foxsquirrel", intercept = FALSE) %>% 
  make_plot.popdens(., species = "foxsquirrel", intercept = FALSE)
  
  # fox squirrel slope habitat
  data.table::fread("./results/f_squ_matrix.csv",
                    data.table = FALSE) %>% 
    as.matrix(.) %>% 
    gen_preds.habitat(., occ_data = dat, 
                      city_covs = bb,
                      species = "foxsquirrel", intercept = FALSE) %>% 
  make_plot.habitat(., species = "foxsquirrel", intercept = FALSE)
  
  
  
  
  
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
  