
################################
# Load and fit multi-city models
#
# Code written by Mason Fidino
################################

source("sourcer.R")


# Read in the data
det_data <- read.csv("./data/detection_data.csv", stringsAsFactors = FALSE)
det_data <- det_data[order(det_data$city),]

# number of sites
nsite <- nrow(det_data) - 2 # lost two  austin sites

# number of cities
ncity <- length(unique(det_data$city))

# number of patch level parameters (including intercept)
npatch_covs <- 2

# number of city level parameters (including intercept)
ncity_covs <- 4

# number of detection parameters (including intercept)
ndet_covs <- 1

# number of species
nspecies <- 8

# make the patch level occupancy covariates
#  this is for each sampling location per city
bx <- matrix(1, ncol = npatch_covs, nrow = nsite)

# locate file path for each cities covariates
#  they all end with "covs.csv" and are nested in the data folder
patch_cov_path <- list.files(path = "./data/",
                             pattern = "covs.csv", full.names = TRUE)

# read in the patch covs and turn list of dataframes into a single dataframe
patch_covs <- lapply(patch_cov_path, read.csv, header= TRUE, stringsAsFactors = FALSE)
patch_covs <- bind_rows(patch_covs)

#### continue from here ###########

# read in all the sites as we need the locationID
all_sites <- read.csv("./data/uwin_all_sites.csv", stringsAsFactors = FALSE)[,-2]

patch_covs <- left_join(patch_covs, all_sites[,c(1,2)], by = c("site" = "LocationName" ))

# work on making the site code in patch_covs
sc_tmp <- stringr::str_pad(patch_covs$LocationID, width = 3, pad = "0")
sc_tmp <- paste(tolower(patch_covs$city), sc_tmp, sep = "-")
patch_covs$site_code <- sc_tmp


# compile down to just the sites that we have data for
patch_covs <- patch_covs[which(patch_covs$site_code %in% as.character(det_data$site_code)),]

ds <-as.character(det_data$site_code)

det_data <- det_data[-which(det_data$site_code %in% ds[-which(ds %in% patch_covs$site_code)]),]

# bring in the number
det_events <- read.csv("./data/species_in_cities.csv")


has_species <- matrix(0, ncol = nspecies, nrow = nrow(det_data))

for(species in 1:nspecies){
  has_species[,species] <- det_events[, species +2] %>% 
    unlist %>% 
    as.numeric %>% 
    rep(., times = as.numeric(table(det_data$city)) )

}

has_species <- data.frame(has_species)
colnames(has_species) <- colnames(det_data)[4:11]
# make a numeric vector for which city it is
city_vec <- as.numeric(factor(det_data$city))


# make the urbanization covariate


urb500 <- prcomp(patch_covs[,c(7,10)], scale. = TRUE)
urb1000 <- prcomp(patch_covs[,c(8,11)], scale. = TRUE)
urb4000 <- prcomp(patch_covs[,c(9,12)], scale. = TRUE)

urb <- data.frame(urb = urb500$x[,1], urb1 = urb1000$x[,1], urb4 = urb4000$x[,1])

my_meds <- data.frame(city = patch_covs$city, urb = urb$urb1) %>% 
  group_by(city) %>% 
  summarise(um = median(urb))

to_plot <- data.frame(city = factor(patch_covs$city, levels = my_meds$city),
                      urb = urb$urb1)


bx[,2] <- urb$urb1
# make the patch level detection covariates
dx <- matrix(1, ncol = ndet_covs, nsite)

# bring in the city covs
cdat <- read.csv("data/city_level_data.csv")
cdat <- cdat[order(cdat$city),]

# make the city level covariates
U <- matrix(1, ncol = ncity_covs, nrow = ncity)

# scale cdat
scale_cdat <- cdat %>% mutate_if(is.numeric, scale)

# get only the covars we want
to_keep_city <- c("habitat", "pop10_dens", "latitude")
U[,-1] <- scale_cdat[,to_keep_city] %>% as.matrix

# do raccoon analysis
data_list <- list(
  y = det_data$raccoon,
  has_species = has_species$raccoon,
  J = det_data$J,
  ncity = ncity,
  nsite = nsite,
  npatch_covs = npatch_covs,
  ncity_covs = ncity_covs,
  ndet_covs = ndet_covs,
  bx = bx,
  dx = dx,
  U = U,
  city_vec = city_vec
)

z <- data_list$y
z[z>1] <- 1


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
# Settings for the jags models
to_monitor <- c("B", "G", "D",
              "rho_mu", "rho_sigma", "rho_cor",
              "sigma.int", "sigma.urb", "z")
my_model <- "./jags_scripts/mvn_occupancy_3preds.R"
nchain <- 6
nadapt <- 2e6
nburnin <- 2e6
nsample <- 83333
nthin <- 2
my_method <- "parallel"


raccoon <- run.jags(
  model = my_model,
  data = data_list,
  n.chains = nchain,
  monitor = to_monitor,
  adapt = nadapt,
  burnin = nburnin,  sample = nsample, thin = nthin, method = my_method,
  inits = inits, summarise = FALSE, modules = "glm")

saveRDS(raccoon, "./results/raccoon.RDS")
rm(raccoon)

# do rabbit
data_list <- list(
  y = det_data$rabbit,
  J = det_data$J,
  ncity = ncity,
  nsite = nsite,
  npatch_covs = npatch_covs,
  ncity_covs = ncity_covs,
  ndet_covs = ndet_covs,
  bx = bx,
  dx = dx,
  U = U,
  city_vec = city_vec
)

z <- data_list$y
z[z>1] <- 1


rabbit <- run.jags(
  model = my_model,
  data = data_list,
  n.chains = nchain,
  monitor = to_monitor,
  adapt = nadapt,
  burnin = nburnin,  sample = nsample, thin = nthin, method = my_method,
  inits = inits, summarise = FALSE, modules = "glm")

saveRDS(rabbit, "./results/rabbit.RDS")
rm(rabbit)

# do coyote
data_list <- list(
 y = det_data$coyote,
 J = det_data$J,
 ncity = ncity,
 nsite = nsite,
 npatch_covs = npatch_covs,
 ncity_covs = ncity_covs,
 ndet_covs = ndet_covs,
 bx = bx,
 dx = dx,
 U = U,
 city_vec = city_vec
)


z <- data_list$y
z[z>1] <- 1

coyote <- run.jags(
  model = my_model,
  data = data_list,
  n.chains = nchain,
  monitor = to_monitor,
  adapt = nadapt,
  burnin = nburnin,  sample = nsample, thin = nthin, method = my_method,
  inits = inits, summarise = FALSE, modules = "glm")

saveRDS(coyote, "./results/coyote.RDS")
rm(coyote)

# do opossum
# oy <- det_data[-grep("deco|foco", as.character(det_data$city)),]
# data_list <- list(
#   y = oy$opossum,
#   J = oy$J,
#   ncity = ncity-2,
#   nsite = nrow(oy),
#   npatch_covs = npatch_covs,
#   ncity_covs = ncity_covs,
#   ndet_covs = ndet_covs,
#   bx = bx[-grep("deco|foco", as.character(det_data$city)),],
#   dx = dx[-grep("deco|foco", as.character(det_data$city)),],
#   U = U[-c(3,4),],
#   city_vec = as.numeric(factor(as.character(oy$city)))
# )
# 
# 
# z <- data_list$y
# z[z>1] <- 1
# inits <- function(chain){
#   gen_list <- function(chain = chain){
#     list( 
#       z = z,
#       B = array(rnorm(ncity-2 * npatch_covs), 
#                 dim = c(ncity-2, npatch_covs)),
#       G = array(rnorm(ncity_covs * npatch_covs), 
#                 dim = c(npatch_covs, ncity_covs)),
#       sigma.int = rgamma(1, 1, 1),
#       sigma.urb = rgamma(1, 1, 1),
#       rho_mat = runif(1, -1, 1),
#       .RNG.name = switch(chain,
#                          "1" = "base::Wichmann-Hill",
#                          "2" = "base::Marsaglia-Multicarry",
#                          "3" = "base::Super-Duper",
#                          "4" = "base::Mersenne-Twister",
#                          "5" = "base::Wichmann-Hill",
#                          "6" = "base::Marsaglia-Multicarry",
#                          "7" = "base::Super-Duper",
#                          "8" = "base::Mersenne-Twister"),
#       .RNG.seed = sample(1:1e+06, 1)
#     )
#   }
#   return(switch(chain,           
#                 "1" = gen_list(chain),
#                 "2" = gen_list(chain),
#                 "3" = gen_list(chain),
#                 "4" = gen_list(chain),
#                 "5" = gen_list(chain),
#                 "6" = gen_list(chain),
#                 "7" = gen_list(chain),
#                 "8" = gen_list(chain)
#   )
#   )
# }
# hm <-table(z, det_data$city)
# 
# hm <- hm[2,] / colSums(hm)
# 
# opossum <- run.jags(
#   model = "./jags_scripts/mvn_int_slope2.R",
#   data = data_list,
#   n.chains = 6,
#   monitor = c("B", "G", "D",
#               "rho_mu", "rho_sigma", "rho_mat","sigma.int", "sigma.urb"),
#   adapt = 2e6,
#   burnin = 2e6,  sample = 166667, method = 'parallel',
#   inits = inits, summarise = FALSE, modules = "glm")
# saveRDS(opossum, "./results/opossum_5city.RDS")
# rm(opossum)

# do skunk
data_list <- list(
  y = det_data$skunk,
  J = det_data$J,
  ncity = ncity,
  nsite = nsite,
  npatch_covs = npatch_covs,
  ncity_covs = ncity_covs,
  ndet_covs = ndet_covs,
  bx = bx,
  dx = dx,
  U = U,
  city_vec = city_vec
)


z <- data_list$y
z[z>1] <- 1

skunk <- run.jags(
  model = my_model,
  data = data_list,
  n.chains = nchain,
  monitor = to_monitor,
  adapt = nadapt,
  burnin = nburnin,  sample = nsample, thin = nthin, method = my_method,
  inits = inits, summarise = FALSE, modules = "glm")

saveRDS(skunk, "./results/skunk.RDS")
rm(skunk)


# do redfox
data_list <- list(
  y = det_data$redfox,
  J = det_data$J,
  ncity = ncity,
  nsite = nsite,
  npatch_covs = npatch_covs,
  ncity_covs = ncity_covs,
  ndet_covs = ndet_covs,
  bx = bx,
  dx = dx,
  U = U,
  city_vec = city_vec
)


z <- data_list$y
z[z>1] <- 1

hm <-table(z, det_data$city)

hm <- hm[2,] / colSums(hm)

fox <- run.jags(
  model = my_model,
  data = data_list,
  n.chains = nchain,
  monitor = to_monitor,
  adapt = nadapt,
  burnin = nburnin,  sample = nsample, thin = nthin, method = my_method,
  inits = inits, summarise = FALSE, modules = "glm")

saveRDS(fox, "./results/fox.RDS")
rm(fox)

# do g squ
data_list <- list(
  y = det_data$graysquirrel,
  J = det_data$J,
  ncity = ncity,
  nsite = nsite,
  npatch_covs = npatch_covs,
  ncity_covs = ncity_covs,
  ndet_covs = ndet_covs,
  bx = bx,
  dx = dx,
  U = U,
  city_vec = city_vec
)


z <- data_list$y
z[z>1] <- 1
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
      rho_mat = runif(1, -1, 1),
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
hm <-table(z, det_data$city)

hm <- hm[2,] / colSums(hm)

g_squ <- run.jags(
  model = "./jags_scripts/mvn_int_slope2.R",
  data = data_list,
  n.chains = 6,
  monitor = c("B", "G", "D",
              "rho_mu", "rho_sigma", "rho_mat","sigma.int", "sigma.urb","z"),
  adapt = 2e6,
  burnin = 2e6,  sample = 166667, method = 'parallel',
  inits = inits, summarise = FALSE, modules = "glm")


saveRDS(g_squ, "./results/g_squ.RDS")
rm(g_squ)

# do f squ
oy <- det_data[-grep("mawi", as.character(det_data$city)),]
data_list <- list(
  y = oy$foxsquirrel,
  J = oy$J,
  ncity = ncity-1,
  nsite = nrow(oy),
  npatch_covs = npatch_covs,
  ncity_covs = ncity_covs,
  ndet_covs = ndet_covs,
  bx = bx[-grep("mawi", as.character(det_data$city)),],
  dx = dx[-grep("mawi", as.character(det_data$city)),],
  U = U[-8,],
  city_vec = as.numeric(factor(as.character(oy$city)))
)


z <- data_list$y
z[z>1] <- 1
inits <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = z,
      B = array(rnorm(ncity-1 * npatch_covs), 
                dim = c(ncity-1, npatch_covs)),
      G = array(rnorm(ncity_covs * npatch_covs), 
                dim = c(npatch_covs, ncity_covs)),
      sigma.int = rgamma(1, 1, 1),
      sigma.urb = rgamma(1, 1, 1),
      rho_mat = runif(1, -1, 1),
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
hm <-table(z, det_data$city)

hm <- hm[2,] / colSums(hm)

f_squ <- run.jags(
  model = "./jags_scripts/mvn_int_slope2.R",
  data = data_list,
  n.chains = 6,
  monitor = c("B", "G", "D",
              "rho_mu", "rho_sigma", "rho_mat","sigma.int", "sigma.urb","z"),
  adapt = 2e6,
  burnin = 2e6,  sample = 166667, method = 'parallel',
  inits = inits, summarise = FALSE, modules = "glm")

saveRDS(f_squ, "./results/f_squ.RDS")
rm(f_squ)