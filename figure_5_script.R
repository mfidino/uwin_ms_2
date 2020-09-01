

# alpha diveristy plot (i.e., Figure 5)

# Note: This is just a rough svg plot, I ended up tweaking the final
#  a bit in inkscape.


source("sourcer.R")

cpo <- read.csv("cpo_scores.csv", 
                stringsAsFactors = FALSE
)
cpo <- cpo[order(cpo$species),]

which_folder <- "best"

my_files <- list.files(paste0("./results/", which_folder),
                       pattern = "matrix", 
                       full.names = TRUE)

sp_name <- gsub(".*/(\\w+)_matrix.csv", "\\1", my_files)

pretty_sp <- c('coyote', 'fox squirrel', 'gray squirrel',
               'opossum', 'e. cottontail', 'raccoon',
               'red fox', 'skunk')

my_z_files <- list.files(paste0("./results/", which_folder),
                         pattern = "zed", 
                         full.names = TRUE)

my_medians <- data.frame(city = det_data$city,
                         site = det_data$site,
                         matrix(0, ncol = length(sp_name),
                                nrow = nrow(det_data)))
colnames(my_medians)[-c(1:2)] <- sp_name

z_matrix <- array(0, dim  =c(125000, 808))

for(species in 1:length(sp_name)){
  zed <- data.table::fread(my_z_files[species],
                           data.table = FALSE, nrows = 125000) %>% 
    as.matrix(.) %>% 
    sweep(., 2, has_species[,sp_name[species]], "*")
  
  #my_medians[,sp_name[species]] <- apply(zed, 2, median)
  z_matrix <-z_matrix + zed
  rm(zed)
}

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
                  "Orange County,\nCalifornia",
                  "Manhattan,\nKansas",
                  "Madison,\nWisconsin",
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

# make lists for 95% CI
low_list <- high_list <-  vector("list", length = ncity)

for(city in 1:ncity){
  
  low_list[[city]] <- z_matrix[,which(det_data$city == unq_city[city] &
                                        hd_cwc$hd_1000 <= 0)]  
  high_list[[city]] <- z_matrix[,which(det_data$city == unq_city[city] &
                                         hd_cwc$hd_1000 > 0)]
  
}

all_lows <- unlist(low_list)
all_highs <- unlist(high_list)


# step 1. select a city
city_results <- vector("list", length = ncity)

filter_rich <- function(x, y){
  species_per <- strsplit(names(x),"-") %>% sapply(.,function(x) sum(as.numeric(x)))
  return(x[which(species_per == y)])
}

to_skip <- seq(0, 115000, by = 10000)
n_sims <- 9999
for(iter in 1:length(to_skip)){
  print(iter)
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


best_low <- best_high <- rep(NA, ncity)

for(i in 1:ncity){
  best_low[[i]] <- names(city_results[[i]]$lower_comp)[which.max(city_results[[i]]$lower_comp)]
  best_high[[i]] <- names(city_results[[i]]$upper_comp)[which.max(city_results[[i]]$upper_comp)]
  
}



# Do the plot now that everything is summarized

svg('./plots/sp_rich_figure_best.svg', height = 6, width = 9)
par(mar = c(5,7,0.5,7))
m <- matrix(c(1,1,2), ncol = 1, nrow = 3)
layout(m)
plot(1~1, type ="n", bty = 'n', xlab = "", ylab = "", xaxt = "n",
     yaxt = "n", ylim = c(0,7), xlim = c(1,10))

#axis(1, at = seq(1, 10, 1), labels = F, tck = -0.025)


mtext(text = cplot$pretty, 1, at = 1:10, cex = 0.7, line = 0.7)

axis(2, at = 0:7, labels =F, tck = -0.025)
mtext(text = 0:7,2, las = 1, at = 0:7, line = 1)
mtext(text = "Species richness", 2, at = 3.5, line = 4, cex = 1.25)

cgrad_low <- colorRampPalette(c('white', '#ffb226ff'))
cgrad_high<- colorRampPalette(c('white', '#24d5f7ff'))
my_col_low <- cgrad_low(15)[c(6,8,10)]
my_col_high <- cgrad_high(15)[c(6,8,10)]

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
  lines(c(city_rich_low[cplot$ord[city],2], city_rich_high[cplot$ord[city],2]) ~ 
          c(city - 0.15, city + 0.15),
        col = '#99adbfff', lwd = 3)
}

# add lines between

points(city_rich_low[cplot$ord,2] ~ c(1:10 - 0.15), pch = 24, bg = "#ffb226ff", cex = 2)
points(city_rich_high[cplot$ord,2] ~ c(1:10 + 0.15), pch = 22, bg = "#24d5f7ff", cex = 2)


# do plot of common species
par(mar = c(3,7,0,7), xpd = NA)
plot(1~1, type ="n", bty = 'n', xlab = "", ylab = "", xaxt = "n",
     yaxt = "n", ylim = c(1,8), xlim = c(1,10))


lowmat <- matrix(as.numeric(unlist(str_split(best_low, "-"))), ncol = 8, nrow = 10, byrow = TRUE)
richmat <- matrix(as.numeric(unlist(str_split(best_high, "-"))), ncol = 8, nrow = 10, byrow = TRUE)

test <- rbind(lowmat, richmat)
new_sp_ord <- order(colSums(rbind(lowmat, richmat)), decreasing = TRUE)

mtext(text = pretty_sp[new_sp_ord], 2, at = 1:8, cex = 0.65, line = 0.95, las = 1)
mtext(text = pretty_sp[new_sp_ord], 4, at = 1:8, cex = 0.65, line = 0.9, las = 1)

# add a box for every other
#for(sp in seq(1,10,1)){
#  rect(sp - 0.25, 1, sp + 0.25, 8, col = "#cad3dd", border = FALSE)
#}

# add thick line along each row

for(sp in 1:8){
  lines(c(0.75,10.25), rep(sp,2), col = '#99adbfff', lwd = 4)
  # add points at each spot
  #points(c(1:10 - 0.15), y = rep(sp,10), pch = 16, col = '#99adbfff', cex = 3.5 )
  #points(c(1:10 + 0.15), y = rep(sp,10), pch = 16, col = '#99adbfff', cex = 3.5 )
  #points(c(1:10 - 0.15), y = rep(sp,10), pch = 16, col = 'white', cex = 2.5 )
  #points(c(1:10 + 0.15), y = rep(sp,10), pch = 16, col = 'white', cex = 2.5 )
}



rect(-0.5, 8.4,10.35,10.25, col = "#e9d7b8ff", border = FALSE)
text(x = -0.5, mean(c(8,10.25)), labels = 'Avg. housing\n     density\n(units per 1km)', pos = 4, cex = 0.9)
pretty_hd <- round(cplot$hden, -1)

#tester <- seq(9, 10, length.out = 10)
tester <- rep(mean(c(8.4,10.25)), 10)

text(x = 1:10, y = tester, labels = pretty_hd, pos = NULL, cex = 1.2 )
mtext("Most likely wildlife community at median species richness", 1, cex = 1.25, line = 1.5)




my_col_low <- cgrad_low(100)[100]
my_col_high <- cgrad_high(100)[100]

nlow <- city_rich_low[cplot$ord, 2]
nhigh <- city_rich_high[cplot$ord,2]  
#best_low <- best_low[cplot$ord]
#best_high <- best_high[cplot$ord]

#test <- matrix(sapply(strsplit(best_low,"-"), ))

# reorder based off of how common the species is


for(city in 1:10){
  tmp_low <- as.numeric(unlist(strsplit(best_low[city],"-")))[new_sp_ord]
  the_sp_low <- which(tmp_low ==1)
  points(rep(city - 0.15, nlow[city]),y=the_sp_low, pch = 24, bg = my_col_low, cex = 2 )
  tmp_high <- as.numeric(unlist(strsplit(best_high[city],"-")))[new_sp_ord]
  the_sp_high <- which(tmp_high ==1)
  points(rep(city + 0.15, nhigh[city]),y=the_sp_high, pch = 22, bg = my_col_high, cex = 2 )
  #if(detplot[city,i] == 0){
  # points(city - 0.15, y=i, cex = 2.5, col = '#95a9bcff', pch = 16)
  #points(city + 0.15, y=i, cex = 2.5, col = '#95a9bcff', pch = 16)
}

fright <- 11.1
ylm = 20

lines(x = c(6.5,6.5), y = c(1,26), lty = 2)
#abline(v = 6.5, lty = 2)

legend(x = 4, y = 26.6, legend = c(expression(underline('<') *' avg. housing density'),
                                   '> avg. housing density'),
       pch = c(24,22), pt.bg = c(my_col_low, my_col_high),
       title = 'Habitat patches surrounded by:', horiz = TRUE, cex = 1.2, pt.cex = 1.8,
       bg = 'white', box.col = 'white')


text(x = fright - 0.17, y = ylm + 2.4, 'Credible\ninterval', cex = 1.2)
my_col_low <- cgrad_low(15)[c(6,8,10)]
my_col_high <- cgrad_high(15)[c(6,8,10)]

cgrad_null<- colorRampPalette(c('white', '#95a9bcff'))
my_col_null <- cgrad_null(15)[c(6,8,10)]
my_new_wid <- my_wid #- 0.004
my_new_wid <- c(0.02,0.04,0.06)
hm1 <- 10.7
rect(hm1 - my_new_wid[1], ylm+1.2,hm1 + my_new_wid[1], ylm - 4.2 ,
     col = my_col_low[1], border = FALSE)

rect(hm1 - my_new_wid[2], ylm+0.1, hm1 + my_new_wid[2], ylm -3.1  ,
     col = my_col_low[2], border = FALSE)

rect(hm1 - my_new_wid[3], ylm-0.9, hm1 + my_new_wid[3], ylm -2.2  ,
     col = my_col_low[3], border = FALSE)

text(x = hm1+0.025, ylm - 1.5, "80%",pos = 4)
text(x = hm1+0.025, ylm - 0.5, "90%",pos = 4)
text(x = hm1+0.025, ylm +0.5, "95%",pos = 4)
text(x = hm1+0.025, ylm - 2.7, "90%",pos = 4)
text(x = hm1+0.025, ylm - 3.7, "95%",pos = 4)

hm2 <- 11.2
rect(hm2 - my_new_wid[1], ylm+1.2,hm2 + my_new_wid[1], ylm - 4.2 ,
     col = my_col_high[1], border = FALSE)

rect(hm2 - my_new_wid[2], ylm+0.1, hm2 + my_new_wid[2], ylm -3.1  ,
     col = my_col_high[2], border = FALSE)

rect(hm2 - my_new_wid[3], ylm-0.9, hm2 + my_new_wid[3], ylm -2.2  ,
     col = my_col_high[3], border = FALSE)



dev.off()


