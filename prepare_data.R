# Read in the data
det_data <- read.csv(
  "./data/detection_data.csv",
  stringsAsFactors = FALSE
)
det_data <- det_data[order(det_data$city,det_data$year),]

# number of sites
nsite <- nrow(det_data)

# number of cities
ncity <- length(unique(det_data$city))

# number of years
nyear <- length(unique(det_data$year))

# number of patch level parameters (excluding intercept)
npatch_covs <- 5

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

# read in all the sites as we need the locationID
all_sites <- read.csv("./data/uwin_all_sites.csv", stringsAsFactors = FALSE)[,-2]


patch_covs <- left_join(patch_covs, all_sites[,c(1,2)], by = c("site" = "LocationName" ))

# to join these we grab the unique name of the city and then the site code info.
sc_temp <- paste(tolower(patch_covs$city), patch_covs$LocationID, 2, sep="-")

det_temp <- strsplit(det_data$site, "-") %>% sapply(., function(x)paste(x[1], x[2], sep="-") )
det_temp <- paste(det_data$city, det_temp, sep = "-")

det_data$site_code <- det_temp

patch_covs$site_code <- sc_temp


# compile down to just the sites that we have data for
patch_covs <- patch_covs[which(patch_covs$site_code %in% as.character(det_data$site_code)),]

# replicate patch covs for varying years

patch_covs <- left_join(det_data[,c("site_code", "year")], patch_covs[,c("site_code", "city", "hd_1000")],
                        by = "site_code") %>% na.omit %>% arrange(year, site_code)

ds <-as.character(det_data$site_code)

det_data <- det_data %>%  arrange(year,site_code)

# bring in the number
#det_events <- read.csv("./data/species_in_cities.csv")

det_events <- det_data %>% group_by(city) %>% 
  summarise_if(.predicate = is.numeric, sum)
det_events <- det_events[,c(3:10)]


det_events<- apply(det_events, 2, function(x) ifelse(x>0, 1, 0))

has_species <- matrix(0, ncol = nspecies, nrow = nrow(det_data))

# make a numeric vector for which city it is
city_vec <- as.numeric(factor(det_data$city))

for(species in 1:nspecies){
  has_species[,species] <- det_events[city_vec, species]
}

has_species <- data.frame(has_species)
colnames(has_species) <- colnames(det_data)[4:11]


# group center housing density

# need to do a better job figuring out what goes where now.
city_mu <- patch_covs %>% group_by(city) %>% 
  summarise(mu = mean(hd_1000))


hd_cwc <- patch_covs %>% group_by(city) %>% 
  mutate(hd_1000 = (hd_1000 - mean(hd_1000)) / 1000) %>% 
  dplyr::select(., site_code,hd_1000, year )

bx[,1] <- hd_cwc$hd_1000
# make the patch level detection covariates
dx <- matrix(1, ncol = ndet_covs, nsite)

# bring in the city covs
cdat <- read.csv("data/city_level_data.csv")
cdat <- cdat[order(cdat$city),]

#cdat$habitat <- (cdat$habitat - 0.3) 

test <- patch_covs %>% group_by(city) %>% 
  summarise(hden = mean(hd_1000))

cdat$hden <- test$hden
rm(test)

#cdat$hden <- cdat$hden / 1000

# make the city level covariates
U <- matrix(1, ncol = 3, nrow = ncity)

# scale cdat
scale_cdat <- cdat %>% mutate_if(is.numeric, scale)

# get only the covars we want
to_keep_city <- c("habitat", "hden")
U[,-1] <- scale_cdat[,to_keep_city] %>% as.matrix


bx[,2] <- U[city_vec, 2]
bx[,3] <- U[city_vec, 3]
bx[,4] <- bx[,1] * bx[,2]
bx[,5] <- bx[,1] * bx[,3]

# count number of repeat years at a city
city_has_data <- (table(hd_cwc$city, hd_cwc$year) > 0) %>% apply(., 2, as.numeric)

my_2016 <- my_2018<-  rep(0, nrow(bx))
my_2016[det_data$year == 2016] <- 1
my_2018[det_data$year == 2018] <- 1
my_2018[det_data$city == "lbca"] <- 0
my_2018[det_data$city == "wide"] <- 0

rep_count <- rowSums(city_has_data)
