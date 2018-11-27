
source("sourcer.R")

# required packages
packs <- c("dplyr", "reshape2", "runjags", "MCMCpack", "mcmcplots",
  'parallel', "coda")

# load the packages
package_load(packs)


# read in the data
dat <- read.csv("./data/detection_data.csv")


# number of sites
nsite <- nrow(dat) - 2 # lost two  austin sites

# number of cities
ncity <- length(unique(dat$city))

# number of patch level parameters (including intercept)
npatch_covs <- 2

# number of city level parameters (including intercept)
ncity_covs <- 4

# number of detection parameters (including intercept)
ndet_covs <- 1

# number of species
nspecies <- 8

# make the patch level occupancy covariates
bx <- matrix(1, ncol = npatch_covs, nrow = nsite)

# read in the patch covariates

my_covs <- list.files(path = "./data/", pattern = "covs.csv", full.names = TRUE)

my_dfs <- lapply(my_covs, read.csv, header= TRUE, stringsAsFactors = FALSE)
my_dfs <- bind_rows(my_dfs)
# read in all the sites as we need the locationID
all_sites <- read.csv("./data/uwin_all_sites.csv", stringsAsFactors = FALSE)[,-2]

my_dfs <- left_join(my_dfs, all_sites[,c(1,2)], by = c("site" = "LocationName" ))

# work on making the site code in my_dfs
sc_tmp <- stringr::str_pad(my_dfs$LocationID, width = 3, pad = "0")
sc_tmp <- paste(tolower(my_dfs$city), sc_tmp, sep = "-")
my_dfs$site_code <- sc_tmp


# compile down to just the sites that we have data for
my_dfs <- my_dfs[which(my_dfs$site_code %in% as.character(dat$site_code)),]

ds <-as.character(dat$site_code)

dat <- dat[-which(dat$site_code %in% ds[-which(ds %in% my_dfs$site_code)]),]
# make a numeric vector for which city it is
city_vec <- as.numeric(dat$city)


# make the urbanization covariate


urb500 <- prcomp(my_dfs[,c(4,7,10)], scale. = TRUE)
urb1000 <- prcomp(my_dfs[,c(5,8,11)], scale. = TRUE)
urb4000 <- prcomp(my_dfs[,c(6,9,12)], scale. = TRUE)

urb <- data.frame(urb = urb500$x[,1], urb1 = urb1000$x[,1], urb4 = urb4000$x[,1])

my_meds <- data.frame(city = my_dfs$city, urb = urb$urb1) %>% 
  group_by(city) %>% 
  summarise(um = median(urb))

to_plot <- data.frame(city = factor(my_dfs$city, levels = my_meds$city),
                      urb = urb$urb1)

boxplot(urb ~ city, data = to_plot)

bx[,2] <- urb$urb1
# make the patch level detection covariates
dx <- matrix(1, ncol = ndet_covs, nsite)
dx[,2] <- as.numeric(scale(log(dat$photos)))
dx[,3] <- dx[,2] ^ 2

# bring in the city covs
cdat <- read.csv("data/city_level_data.csv")
cdat <- cdat[order(cdat$city),]

# make the city level covariates
U <- matrix(1, ncol = ncity_covs, nrow = ncity)

# scale cdat
scale_cdat <- cdat %>% mutate_if(is.numeric, scale)

# get only the covars we want
to_keep_city <- c("habitat", "pop10_dens")
U[,-1] <- scale_cdat[,to_keep_city] %>% as.matrix

# do raccoon analysis
data_list <- list(
  y = dat$raccoon,
  J = dat$J,
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
hm <-table(z, dat$city)

hm <- hm[2,] / colSums(hm)

raccoon <- run.jags(
  model = "./jags_scripts/mvn_int_slope2.R",
  data = data_list,
  n.chains = 6,
  monitor = c("B", "G", "D",
              "rho_mu", "rho_sigma", "rho_mat","sigma.int", "sigma.urb","z"),
  adapt = 2e6,
  burnin = 2e6,  sample = 166667, method = 'parallel',
  inits = inits, summarise = FALSE, modules = "glm")

saveRDS(raccoon, "./results/raccoon8.RDS")
rm(raccoon)


# do rabbit
data_list <- list(
  y = dat$rabbit,
  J = dat$J,
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
hm <-table(z, dat$city)

hm <- hm[2,] / colSums(hm)

rabbit <- run.jags(
  model = "./jags_scripts/mvn_int_slope2.R",
  data = data_list,
  n.chains = 6,
  monitor = c("B", "G", "D",
              "rho_mu", "rho_sigma", "rho_mat","sigma.int", "sigma.urb","z"),
  adapt = 2e6,
  burnin = 2e6,  sample = 166667, method = 'parallel',
  inits = inits, summarise = FALSE, modules = "glm")

saveRDS(rabbit, "./results/rabbit.RDS")
rm(rabbit)

# do coyote
#data_list <- list(
#  y = dat$coyote,
#  J = dat$J,
#  ncity = ncity,
#  nsite = nsite,
#  npatch_covs = npatch_covs,
#  ncity_covs = ncity_covs,
#  ndet_covs = ndet_covs,
#  bx = bx,
#  dx = dx,
#  U = U,
#  city_vec = city_vec
#)
#
#
#z <- data_list$y
#z[z>1] <- 1
#inits <- function(chain){
#  gen_list <- function(chain = chain){
#    list( 
#      z = z,
#      B = array(rnorm(ncity * npatch_covs), 
#                dim = c(ncity, npatch_covs)),
#      G = array(rnorm(ncity_covs * npatch_covs), 
#                dim = c(npatch_covs, ncity_covs)),
#      sigma.int = rgamma(1, 1, 1),
#      sigma.urb = rgamma(1, 1, 1),
#      rho_mat = runif(1, -1, 1),
#      .RNG.name = switch(chain,
#                         "1" = "base::Wichmann-Hill",
#                         "2" = "base::Marsaglia-Multicarry",
#                         "3" = "base::Super-Duper",
#                         "4" = "base::Mersenne-Twister",
#                         "5" = "base::Wichmann-Hill",
#                         "6" = "base::Marsaglia-Multicarry",
#                         "7" = "base::Super-Duper",
#                         "8" = "base::Mersenne-Twister"),
#      .RNG.seed = sample(1:1e+06, 1)
#    )
#  }
#  return(switch(chain,           
#                "1" = gen_list(chain),
#                "2" = gen_list(chain),
#                "3" = gen_list(chain),
#                "4" = gen_list(chain),
#                "5" = gen_list(chain),
#                "6" = gen_list(chain),
#                "7" = gen_list(chain),
#                "8" = gen_list(chain)
#  )
#  )
#}
#
#
#coyote <- run.jags(
#  model = "./jags_scripts/mvn_int_slope2.R",
#  data = data_list,
#  n.chains = 6,
#  monitor = c("B", "G", "D",
#              "rho_mu", "rho_sigma", "rho_mat","sigma.int", "sigma.urb","z"),
#  adapt = 2e6,
#  burnin = 2e6,  sample = 166667, method = 'parallel',
#  inits = inits, summarise = FALSE, modules = "glm")
#
#saveRDS(coyote, "./results/coyote.RDS")
#rm(coyote)

# do opossum
oy <- dat[-grep("deco|foco", as.character(dat$city)),]
data_list <- list(
  y = oy$opossum,
  J = oy$J,
  ncity = ncity-2,
  nsite = nrow(oy),
  npatch_covs = npatch_covs,
  ncity_covs = ncity_covs,
  ndet_covs = ndet_covs,
  bx = bx[-grep("deco|foco", as.character(dat$city)),],
  dx = dx[-grep("deco|foco", as.character(dat$city)),],
  U = U[-c(3,4),],
  city_vec = as.numeric(factor(as.character(oy$city)))
)


z <- data_list$y
z[z>1] <- 1
inits <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = z,
      B = array(rnorm(ncity-2 * npatch_covs), 
                dim = c(ncity-2, npatch_covs)),
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
hm <-table(z, dat$city)

hm <- hm[2,] / colSums(hm)

opossum <- run.jags(
  model = "./jags_scripts/mvn_int_slope2.R",
  data = data_list,
  n.chains = 6,
  monitor = c("B", "G", "D",
              "rho_mu", "rho_sigma", "rho_mat","sigma.int", "sigma.urb"),
  adapt = 2e6,
  burnin = 2e6,  sample = 166667, method = 'parallel',
  inits = inits, summarise = FALSE, modules = "glm")
saveRDS(opossum, "./results/opossum_5city.RDS")
rm(opossum)

# do skunk
data_list <- list(
  y = dat$skunk,
  J = dat$J,
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
hm <-table(z, dat$city)

hm <- hm[2,] / colSums(hm)

skunk <- run.jags(
  model = "./jags_scripts/mvn_int_slope2.R",
  data = data_list,
  n.chains = 6,
  monitor = c("B", "G", "D",
              "rho_mu", "rho_sigma", "rho_mat","sigma.int", "sigma.urb","z"),
  adapt = 2e6,
  burnin = 2e6,  sample = 166667, method = 'parallel',
  inits = inits, summarise = FALSE, modules = "glm")

saveRDS(skunk, "./results/skunk.RDS")
rm(skunk)


# do redfox
data_list <- list(
  y = dat$redfox,
  J = dat$J,
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
hm <-table(z, dat$city)

hm <- hm[2,] / colSums(hm)

fox <- run.jags(
  model = "./jags_scripts/mvn_int_slope2.R",
  data = data_list,
  n.chains = 6,
  monitor = c("B", "G", "D",
              "rho_mu", "rho_sigma", "rho_mat","sigma.int", "sigma.urb","z"),
  adapt = 2e6,
  burnin = 2e6,  sample = 166667, method = 'parallel',
  inits = inits, summarise = FALSE, modules = "glm")

saveRDS(fox, "./results/fox.RDS")
rm(fox)

# do g squ
data_list <- list(
  y = dat$graysquirrel,
  J = dat$J,
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
hm <-table(z, dat$city)

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
oy <- dat[-grep("mawi", as.character(dat$city)),]
data_list <- list(
  y = oy$foxsquirrel,
  J = oy$J,
  ncity = ncity-1,
  nsite = nrow(oy),
  npatch_covs = npatch_covs,
  ncity_covs = ncity_covs,
  ndet_covs = ndet_covs,
  bx = bx[-grep("mawi", as.character(dat$city)),],
  dx = dx[-grep("mawi", as.character(dat$city)),],
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
hm <-table(z, dat$city)

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