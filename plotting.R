source("sourcer.R")



which_folder <- "reparam_longrun"

my_files <- list.files(paste0("./results/", which_folder),
                       pattern = "matrix", 
                       full.names = TRUE)

sp_name <- gsub(".*/(\\w+)_matrix.csv", "\\1", my_files)

pretty_sp <- c('coyote', 'fox squirrel', 'gray squirrel',
               'opossum', 'e. cottontail', 'raccoon',
               'red fox', 'skunk')




# check if a sub-folder exists
f_exist <- dir.exists(paste0("./plots/", sp_name))
for(ex in 1:length(f_exist)){
  if(!f_exist[ex]){
    dir.create(paste0("./plots/", sp_name[ex]))
  }
}

# plots for intercept

for(species in 1:length(sp_name)){
  
  # as a function of housing density
  
  # read in the mcmc matrix
  results_matrix <- data.table::fread(my_files[species],
                                      data.table = FALSE) %>% 
                    as.matrix(.)
  
  # the covariate for prediction
  new_data <- data.frame(hden = 0,
                         prophab_btwn = 0,
                         hden_btwn = seq(200, 1400, 10))
  
  # the city data
  new_city_data <- data.frame(hden = 0,
                              prophab_btwn = 0,
                              hden_btwn = cdat$hden)
  row.names(new_city_data) <- cdat$city
  
  preds <- predict.intercept(mmat = results_matrix,
                       new_data = new_data,
                       city_data = cdat,
                       new_city_data = new_city_data,
                       species_there = det_events[,sp_name[species]])
  makeplot.hdens(preds = preds, species = sp_name[species],x = new_data$hden_btwn,
                   cityx = cdat$hden, species_there = det_events[,sp_name[species]],
                 pp = FALSE)
  
  # for habitat intercept
  new_data <- data.frame(hden = 0,
                         prophab_btwn = seq(0.1, 0.7, 0.005),
                         hden_btwn = 0)
  
  new_city_data <- data.frame(hden = 0,
                              prophab_btwn = cdat$habitat,
                              hden_btwn = 0)
  row.names(new_city_data) <- cdat$city
  
  
  preds <- predict.intercept(mmat = results_matrix,
                                  new_data = new_data,
                                  city_data = cdat,
                                  new_city_data = new_city_data,
                                  species_there = det_events[,sp_name[species]])
  
  makeplot.habitat(preds = preds, species = sp_name[species],x = new_data$prophab_btwn,
                   cityx = cdat$habitat, species_there = det_events[,sp_name[species]],
                   pp = FALSE)
  
}
# plots for slopes

for(species in 1:length(sp_name)){
  
  # as a function of housing density
  
  # read in the mcmc matrix
  results_matrix <- data.table::fread(my_files[species],
                                      data.table = FALSE) %>% 
    as.matrix(.)
  
  # the covariate for prediction
  
  new_data <- data.frame(hden = 0,
                         prophab_btwn = 0,
                         hden_btwn = seq(200, 1400, 10))
  
  
  # the city data
  new_city_data <- data.frame(hden = 0,
                              prophab_btwn = 0,
                              hden_btwn = cdat$hden)
  row.names(new_city_data) <- cdat$city
  
  preds <- predict.slope(mmat = results_matrix,
                             new_data = new_data,
                             city_data = cdat,
                             new_city_data = new_city_data,
                             species_there = det_events[,sp_name[species]])
  
  makeplot.hdens(preds = preds, species = sp_name[species],x = new_data$hden_btwn,
                 cityx = cdat$hden, species_there = det_events[,sp_name[species]],
                 intercept = FALSE, window = FALSE, pp = FALSE)
}

for(species in 1:length(sp_name)){
  
  # as a function of housing density
  
  # read in the mcmc matrix
  results_matrix <- data.table::fread(my_files[species],
                                      data.table = FALSE) %>% 
    as.matrix(.)
  
  # for habitat intercept
  new_data <- data.frame(hden = 0,
                         prophab_btwn = seq(0.1, 0.7, 0.005),
                         hden_btwn = 0)
  
  new_city_data <- data.frame(hden = 0,
                              prophab_btwn = cdat$habitat,
                              hden_btwn = 0)
  row.names(new_city_data) <- cdat$city
  
  
  preds <- predict.slope(mmat = results_matrix,
                             new_data = new_data,
                             city_data = cdat,
                             new_city_data = new_city_data,
                             species_there = det_events[,sp_name[species]])
  
  makeplot.habitat(preds = preds, species = sp_name[species],x = new_data$prophab_btwn,
                   cityx = cdat$habitat,
                   intercept = FALSE, window = FALSE, pp = FALSE)
  
  
}




# alpha diveristy plot

my_z_files <- list.files(paste0("./results/", which_folder),
                       pattern = "zed", 
                       full.names = TRUE)

my_medians <- data.frame(city = det_data$city,
                         site = det_data$site,
                         matrix(0, ncol = length(sp_name),
                                nrow = nrow(det_data)))
colnames(my_medians)[-c(1:2)] <- sp_name

z_matrix <- array(0, dim  =c(500000, 691))

for(species in 1:length(sp_name)){
  zed <- data.table::fread(my_z_files[species],
                           data.table = FALSE, nrows = 500000) %>% 
                             as.matrix(.) %>% 
    sweep(., 2, has_species[,sp_name[species]], "*")
  
  #my_medians[,sp_name[species]] <- apply(zed, 2, median)
  z_matrix[,,species] <-z_matrix + zed
  rm(zed)
}

# do it again but paste together each species
to_skip <- seq(0, 450000, by = 50000)
to_skip[-1] <- to_skip[-1] + 1



for(iter in 2:length(to_skip)){
  z_matrix <- array(0, dim  =c(50000, 691,8))
for(species in 1:length(sp_name)){
  zed <- data.table::fread(my_z_files[species],
                           data.table = FALSE, nrows = 50000,
                           skip=to_skip[iter]) %>% 
    as.matrix(.) %>% 
    sweep(., 2, has_species[,sp_name[species]], "*")
  
  #my_medians[,sp_name[species]] <- apply(zed, 2, median)
  z_matrix[,,species] <-zed
  rm(zed)
}
  species_paste <- apply(z_matrix, c(1,2), paste, collapse = '-')
  readr::write_csv(data.frame(species_paste, stringsAsFactors = FALSE),
                   "tmp_paste.csv", append = TRUE)
  rm(species_paste)
  rm(z_matrix)
  
}



test <- apply(z_matrix, c(1,2), paste, collapse = "-")

# do average species richness in patch per city
unq_city <- sort(unique(det_data$city))
city_rich_low <- matrix(0, ncol = 3, nrow = 10)
city_rich_high <- matrix(0, ncol = 3, nrow = 10)
city_low_gradient <- city_high_gradient <- matrix(0, ncol = 6, nrow = 10)

hd_cwc %>% group_by(city) %>% 
  summarise(min = min(hd_1000),
            max = max(hd_1000))

# pretty city names, ordered by average housing density

cplot <- cdat
cplot$ord <- order(cplot$city)

cplot$pretty <- c("Austin,\nTexas",
                  "Chicago,\nIllinois",
                  "Denver,\nColorado",
                  "Fort Collins,\nColorado",
                  "Iowa City,\nIowa",
                  "Indianapolis,\nIndiana",
                  "Long Beach,\nCalifornia",
                  "Manhattan,\nKansas",
                  "Madison,\nWisconson",
                  "Wilmington,\nDeleware")
cplot <- cplot[order(cplot$hden),]
for(city in 1:ncity){

  city_rich_low[city,] <- quantile(z_matrix[,which(det_data$city == unq_city[city] &
                                                 hd_cwc$hd_1000 <= 0)], 
                         probs = c(0.025,0.5,0.975))  
  city_rich_high[city,] <- quantile(z_matrix[,which(det_data$city == unq_city[city] &
                                                     hd_cwc$hd_1000 > 0)], 
                                   probs = c(0.025,0.5,0.975))
  city_low_gradient[city,] <-  quantile(z_matrix[,which(det_data$city == unq_city[city] &
                                                            hd_cwc$hd_1000 <= 0)], 
                                          probs = c(0.025,0.05,.1,0.9, 0.95, 0.975))
  city_high_gradient[city,] <- quantile(z_matrix[,which(det_data$city == unq_city[city] &
                                                          hd_cwc$hd_1000 > 0)], 
                                        probs = c(0.025,0.05,.1,0.9, 0.95, 0.975))
  
}
#city_rich_low <- sweep(city_rich_low,1, sp_rich, "/")
#city_rich_high <- sweep(city_rich_high,1, sp_rich, "/")
windows(8,5, xpos = 10)
par(mar = c(5,7,0.5,0.5))
m <- matrix(c(1,1,2), ncol = 1, nrow = 3)
layout(m)
plot(1~1, type ="n", bty = 'n', xlab = "", ylab = "", xaxt = "n",
     yaxt = "n", ylim = c(0,7), xlim = c(1,10))

#axis(1, at = seq(1, 10, 1), labels = F, tck = -0.025)


mtext(text = cplot$pretty, 1, at = 1:10, cex = 0.65, line = 0.7)

axis(2, at = 0:7, labels =F, tck = -0.025)
mtext(text = 0:7,2, las = 1, at = 0:7, line = 1)

cgrad_low <- colorRampPalette(c('white', '#ffb226ff'))
cgrad_high<- colorRampPalette(c('white', '#24d5f7ff'))
my_col_low <- cgrad_low(10)[c(2,4,6)]
my_col_high <- cgrad_high(10)[c(2,4,6)]

my_loc <- matrix(c(1:3, rev(4:6)), ncol = 2, nrow = 3)
my_wid <- c(0.025, 0.025*1.5, 0.025*2)
# do a color strip 

for(city in 1:ncity){
  for(color in 1:3){
    rect(city-0.15-my_wid[color], city_low_gradient[cplot$ord[city], my_loc[color,1]],
         city-0.15+my_wid[color], city_low_gradient[cplot$ord[city], my_loc[color,2]],
         border = FALSE, col = my_col_low[color])
    rect(city+0.15-my_wid[color], city_high_gradient[cplot$ord[city], my_loc[color,1]],
         city+0.15+my_wid[color], city_high_gradient[cplot$ord[city], my_loc[color,2]],
         border = FALSE, col = my_col_high[color])
    
  }
}

for(city in 1:ncity){
  #lines(city_rich_low[cplot$ord[city],-2] ~ rep(city-0.15,2),
  #col = "#ffb226ff", lwd = 2)
  #lines(city_rich_high[cplot$ord[city],-2] ~ rep(city+0.15,2),
  #     col = "#24d5f7ff", lwd = 2)
  #lines(c(city_rich_low[cplot$ord[city],2], city_rich_high[cplot$ord[city],2]) ~ 
  #       c(city - 0.15, city + 0.15),
  #    col = 'black', lwd = 5)
  lines(c(city_rich_low[cplot$ord[city],2], city_rich_high[cplot$ord[city],2]) ~ 
          c(city - 0.15, city + 0.15),
        col = '#99adbfff', lwd = 3)
}

# add lines between

points(city_rich_low[cplot$ord,2] ~ c(1:10 - 0.15), pch = 24, bg = "#ffb226ff", cex = 2)
points(city_rich_high[cplot$ord,2] ~ c(1:10 + 0.15), pch = 22, bg = "#24d5f7ff", cex = 2)


# do plot of common species
par(mar = c(0.5,7,0,0.5), xpd = NA)
plot(1~1, type ="n", bty = 'n', xlab = "", ylab = "", xaxt = "n",
     yaxt = "n", ylim = c(1,8), xlim = c(1,10))
mtext(text = sp_name, 2, at = 1:8, cex = 0.65, line = 0.95, las = 1)

# add thick line along each row
for(sp in 1:8){
  lines(c(0.85,10.15), rep(sp,2), col = '#cad3dd', lwd = 4)
  # add points at each spot
  points(c(1:10 - 0.15), y = rep(sp,10), pch = 16, col = '#95a9bcff', cex = 3.5 )
  points(c(1:10 + 0.15), y = rep(sp,10), pch = 16, col = '#95a9bcff', cex = 3.5 )
  #points(c(1:10 - 0.15), y = rep(sp,10), pch = 16, col = 'white', cex = 2.5 )
  #points(c(1:10 + 0.15), y = rep(sp,10), pch = 16, col = 'white', cex = 2.5 )
}

rect(-0.5, 8.6,10.35,10.25, col = "#e9d7b8ff", border = FALSE)
text(x = -0.5, mean(c(8.6,10.25)), labels = 'Avg. housing\n     density\n(units per 1km)', pos = 4, cex = 0.85)
pretty_hd <- round(cplot$hden, -1)

#tester <- seq(9, 10, length.out = 10)
tester <- rep(mean(c(8.6,10.25)), 10)

text(x = 1:10, y = tester, labels = pretty_hd, pos = NULL, cex = 1.2 )

# figure out most common community in low and high

# step 1. select a city
my_city <- 'chil'
city_results <- vector("list", length = ncity)
n_sims <- 50
max_at_species <- function(x, y){
  x <- sort(x)
  species_per <- strsplit(names(x),"-") %>% sapply(.,function(x) sum(as.numeric(x)))
  to_return <- names(x[max(which(species_per == y))])
  return(to_return)
}

filter_rich <- function(x, y){
  species_per <- strsplit(names(x),"-") %>% sapply(.,function(x) sum(as.numeric(x)))
  return(x[which(species_per == y)])
}

to_skip <- seq(0, 490000, by = 10000)
n_sim <- 9999
for(iter in 1:length(to_skip)){
for(city in 1:ncity){
  my_city <- tolower(cplot$city[city])
  tmp_city <- array(0, dim = c(n_sims,length(which(det_data$city == my_city)),8))
  tmp_cwc <- hd_cwc[which(hd_cwc$city == cplot$city[city]),]
  lower_sites <- which(tmp_cwc$hd_1000<=0)
  upper_sites <- which(tmp_cwc$hd_1000>0)
for(species in 1:length(sp_name)){
  
  zed <- data.table::fread(my_z_files[species],
                           data.table = FALSE, select=which(det_data$city == my_city),
                           nrows = n_sims, skip = to_skip[iter]) %>% 
    as.matrix(.) %>% 
    sweep(., 2, has_species[,sp_name[species]][which(det_data$city == my_city)], "*") 
  tmp_city[,,species] <- zed[1:n_sims,]
  rm(zed)
}
  lower_rich <- city_rich_low[cplot$ord[city],2]
  upper_rich <- city_rich_high[cplot$ord[city],2]
  
  lower_comp <- apply(tmp_city[,lower_sites,],c(1,2),paste,collapse='-')
  upper_comp <- apply(tmp_city[,upper_sites,],c(1,2),paste,collapse='-')
  lower_comp <- filter_rich(table(lower_comp),lower_rich)
  upper_comp <- filter_rich(table(upper_comp),upper_rich)

  
  if(iter == 1){
  city_results[[city]] <- list(city = my_city,
    lower = lower_comp,
    upper = upper_comp)
  } else {
    city_results[[city]]$lower_comp <- tapply(c(city_results[[city]]$lower, lower_comp), 
                                    names(c(city_results[[city]]$lower, lower_comp)), sum)
    city_results[[city]]$upper_comp <- tapply(c(city_results[[city]]$upper, upper_comp), 
                                              names(c(city_results[[city]]$upper, upper_comp)), sum)
    
  }
}
}
  
tapply(c(city_results[[city]]$lower, lower_comp), 
       names(c(city_results[[city]]$lower, lower_comp)), sum)



      #max_at_species(lower_comp, lower_rich),
       #                      upper = max_at_species(upper_comp, upper_rich))
   
}
}







my_medians$low <- hd_cwc$hd_1000 < 0

my_medians$comp <- apply(my_medians[,3:10], 1, paste, collapse = '-')
my_medians$rich <- rowSums(my_medians[,3:10])


chicago <- my_medians[my_medians$city =='lbca',]
clow <- chicago[chicago$rich == 1 & chicago$low == TRUE,]

chigh <- chicago[chicago$rich == 2 & chicago$low == FALSE,]
cat('low'); sort(prop.table(table(clow$comp)));cat('high');sort(prop.table(table(chigh$comp)));colnames(my_medians)[3:10]

for(i in 1)

the_winners <- my_medians %>% group_by(city, low) %>% 
  summarise_if(is.numeric, function(x) sum(x)/length(x))
# get ever detected
lower <- round(the_winners[the_winners$low == TRUE,3:10], 2) * 100
lower[lower == 0] <- NA
lower[!is.na(lower)] <- 75
lower$city <- cdat$city
upper <- round(the_winners[the_winners$low == FALSE,3:10],2) * 100
upper[upper == 0] <- NA
upper[!is.na(upper)] <- 75
upper$city <- cdat$city


lower <- lower[cplot$ord,]
upper <- upper[cplot$ord,]
detplot <- det_events
detplot <- detplot[cplot$ord,order(colnames(detplot))]

my_col_low <- cgrad_low(100)
my_col_high <- cgrad_high(100)
for(i in 1:8){
  for(city in 1:10){
  points((city - 0.15),y=i, pch = 16, col = my_col_low[lower[city,i]], cex = 2.5 )
  points((city + 0.15),y=i, pch = 16, col = my_col_high[upper[city,i]], cex = 2.5 )
  if(detplot[city,i] == 0){
    points(city - 0.15, y=i, cex = 2.5, col = '#95a9bcff', pch = 16)
    points(city + 0.15, y=i, cex = 2.5, col = '#95a9bcff', pch = 16)
  }

  }
  
}
abline(v = 6.5, lty = 2)


for(i in 1:nrow(lower)){
  # get biggest winners

  
  }
  

  
  summarise(most_common = names(sort(table(comp))[1]))

my_medians$sp_rich <- rowSums(my_medians[,3:10])


plot(city_rich_low[,2] ~ c(cdat$hden - 10) , ylim = c(0,1), bty = 'l',
     ylab = "Proportion of species seen", xlab = "average housing density in city",
     pch = 19, col = "blue")
for(city in 1:ncity){
  lines(c(city_rich_low[city,2], city_rich_high[city,2]) ~ c(cdat$hden[city] - 10, cdat$hden[city]+10))
}

points(city_rich_high[cplot$ord,2] ~ c(cdat$hden + 10), pch = 15, col = 'red')

abline(v = 800, lty = 2)

legend("topright", legend = c("Less urban patches in city",
                              "More urban patches in city"),col = c("blue", "red"),pch = c(19,15))

# plots for intercept, 
for(species in 1:length(sp_name)){
  # read in the mcmc matrix
  results_matrix <- data.table::fread(my_files[species],
                                      data.table = FALSE) %>% 
    as.matrix(.)
  
  # psi_mu(population density)

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
  