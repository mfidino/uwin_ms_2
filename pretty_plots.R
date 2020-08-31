
library(USAboundaries)
library(sf)

cities <- read.csv("./data/uscitiesv1.5.csv", stringsAsFactors = FALSE)

# remove cities without population

cities <- cities[!is.na(cities$population),]

# remove cities with less than 2500
cities <- cities[which(cities$population>=2500),]

# drop alaska and hawaii
cities <- cities[-which(cities$state_name %in% c('Alaska', 'Hawaii')),]

cities$size <- NA

cities$area <- cities$population / cities$density


# calculate proportion of population in cities at different sizes

my_ranges <- c(0,5e5, 1e6, 5e6, 1e7,1e200)

city_sum <- matrix(0, ncol = 4, nrow = 5)


for(ran in 2:length(my_ranges)){
  # sum of populations
  city_sum[ran-1,1] <- sum(cities$population[between(cities$population,
                                                     my_ranges[ran-1],
                                                     my_ranges[ran])])
  # count of cities
  city_sum[ran-1,2] <- sum(between(cities$population,
                               my_ranges[ran-1],
                               my_ranges[ran]))
  # area of cities
  city_sum[ran-1,3] <- sum(cities$area[between(cities$population,
                                                     my_ranges[ran-1],
                                                     my_ranges[ran])])
}

# make in millions
city_sum[,1] <- city_sum[,1]/1e6
windows()
plot(cities$lat ~ cities$lng)
identify(cities$lng,cities$lat)
# categorize cities
cities$size[between(cities$population,
                    my_ranges[1],
                    my_ranges[2])] <- "Fewer than 500,000"
cities$size[between(cities$population,
                    my_ranges[2],
                    my_ranges[3])] <- "500,000 to 1 million"
cities$size[between(cities$population,
                    my_ranges[3],
                    my_ranges[4])] <- "1 to 5 million"
cities$size[between(cities$population,
                    my_ranges[4],
                    my_ranges[5])] <- "5 to 10 million"
cities$size[which(cities$population >1e7)] <- "10 million or more"


# color switch

cols <- c("#24d5f7ff", "#5aa80eff","#ffb226ff", '#ff4d4eff', "#8e69fcff" )

col_switch <- function(x){
  switch(x,
         "Fewer than 500,000" = "#24d5f7ff",
         "500,000 to 1 million" = "#5ee38bff",
         "1 to 5 million" = "#ffb226ff",
         '10 million or more' = '#9eb1c4ff') %>% 
    return
}
unq_cat <- unique(cities$size)[c(1,4,2,3,5)]


ms <- state_codes[-which(state_codes$jurisdiction_type == "territory"),]

ms <- ms[-which(ms$state_name %in% c('Hawaii', 'Alaska')),]
us_map <- us_states(resolution = 'high',
                    states = ms$state_abbr)


windows(10,8)
tiff('./plots/pretty_us.tiff', height = 8, width = 10, units = 'in', 
     res = 300, compression = 'lzw')

plot(cities$lat[cities$size == unq_cat[1]] ~ cities$lng[cities$size == unq_cat[1]],
     pch = 19, col = scales::alpha(cols[1], 0.5 ),
     bty = 'n', xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex = 1.4)
plot(st_geometry(us_map), add = TRUE, lwd = 3)

points(cities$lat[cities$size == unq_cat[4]] ~ cities$lng[cities$size == unq_cat[4]],
       pch = 21, bg = cols[4], cex = 2.4)
points(cities$lat[cities$size == unq_cat[3]] ~ cities$lng[cities$size == unq_cat[3]],
       pch = 21, bg= cols[3], cex = 2)
points(cities$lat[cities$size == unq_cat[2]] ~ cities$lng[cities$size == unq_cat[2]],
       pch = 21, bg = cols[2], cex = 1.6)
points(cities$lat[cities$size == unq_cat[5]] ~ cities$lng[cities$size == unq_cat[5]],
       pch = 21, bg = cols[5], cex = 2.8)
legend("bottomleft", pch = c(19,21,21,21,21), col = c(cols[1], rep("black",4)),
       pt.bg = cols, pt.cex = c(1.4, 1.6, 2, 2.4, 2.8), legend = unq_cat,
       bty = 'n', cex = 1.55)
dev.off()


windows(4,4)
tiff('./plots/pretty_citycount.tiff', height = 4, width = 5, units = 'in',
     res = 300, compression = "lzw")
par(mar = c(3.5, 9, 0, 2), xpd = NA)

plot(1~1, type ="n", bty = 'n', xlab = "", ylab = "", xaxt = "n",
     yaxt = "n", ylim = c(0.5,5.5), xlim = c(0,4))

u <- par('usr')
for(i in 1:5){
  rect(0, 6-i-0.4, log10(city_sum[i,2]), 6-i + 0.4, col = cols[i])
  text(x = log10(city_sum[i,2])-0.1, y = 6-i, pos = 4, labels = city_sum[i,2], cex = 1.15)
}
axis(1, at = seq(0,4, 1), labels = F, tck = -0.025)
axis(1, at = seq(0,4, 0.5), labels = F, tck = -0.025/2)

mtext(unq_cat, 2, at = rev(1:5), line = 0.1, las = 2, cex = 1.15)
mtext(expression(log[10] * ' number of cities'), 1, at = 2, line = 2.2, cex = 1.15)
mtext(0:4, 1, at = 0:4, line = 0.5)
dev.off()
windows(5,4)
tiff('./plots/pretty_popcount.tiff', height = 4, width = 5, units = 'in',
     res = 300, compression = "lzw")
svg('./plots/pretty_popcount.svg', height = 4, width = 5)
par(mar = c(3.5, 9, 0, 2), xpd = NA)

plot(1~1, type ="n", bty = 'n', xlab = "", ylab = "", xaxt = "n",
     yaxt = "n", ylim = c(0.5,5.5), xlim = c(0,200))

u <- par('usr')
for(i in 1:5){
  rect(0, 6-i-0.4, city_sum[i,1], 6-i + 0.4, col = cols[i])
  text(x = city_sum[i,1], y = 6-i, pos = 4, labels = round(city_sum[i,1]), cex = 1.15)
}
axis(1, at = seq(0,200, 50), labels = F, tck = -0.025)
axis(1, at = seq(0,200, 25), labels = F, tck = -0.025/2)

mtext(unq_cat, 2, at = rev(1:5), line = 0.1, las = 2, cex = 1.15)
mtext('Number of people (millions)', 1, at = 100, line = 2.2, cex = 1.15)
mtext(seq(0,200,50), 1, at = seq(0,200,50), line = 0.5)
dev.off()

library(waffle)
windows(8,5)

tmp <- round(city_sum[,3]/1000)
names(tmp) <- unq_cat
svg("./plots/area_waffle.svg", width = 7, height = 5)
waffle(tmp, rows = 14, 
       title = '',
       colors = cols, size = 1)
dev.off()
tiff('./plots/pretty_areacount.tiff', height = 4, width = 5, units = 'in',
     res = 300, compression = "lzw", )
par(mar = c(3.5, 9, 0, 2), xpd = NA)

plot(1~1, type ="n", bty = 'n', xlab = "", ylab = "", xaxt = "n",
     yaxt = "n", ylim = c(0.5,5.5), xlim = c(0,250000))

u <- par('usr')
for(i in 1:5){
  rect(0, 6-i-0.4, city_sum[i,3], 6-i + 0.4, col = cols[i])
  text(x = city_sum[i,3], y = 6-i, pos = 4, labels = round(city_sum[i,3]/1e3), cex = 1.15)
}
axis(1, at = seq(0,250000,50000), labels = F, tck = -0.025)
axis(1, at = seq(0,250000,25000), labels = F, tck = -0.025/2)

mtext(unq_cat, 2, at = rev(1:5), line = 0.1, las = 2, cex = 1.15)
mtext(expression('Cumulative area (thousands of '*km^2*')'), 1, at = 125000, line = 2.2, cex = 1.15)
mtext(seq(0,250,50), 1, at = seq(0,250000,50000), line = 0.5)
dev.off()


city_density <- cities %>% group_by(size) %>% 
  summarise(mden = median(density), se_den = sd(density) / sqrt(length(density)))
city_density$size <- factor(city_density$size, levels = unq_cat)
city_density <- city_density[order(city_density$size),]
city_density$high <- city_density$mden + city_density$se_den 
city_density$low <- city_density$mden - city_density$se_den 

windows(8,4)
tiff('./plots/pretty_density.tiff', height = 4, width = 8, units = 'in',
     res = 300, compression = "lzw")
par(mar = c(3.5, 9, 0, 2), xpd = NA)

plot(1~1, type ="n", bty = 'n', xlab = "", ylab = "", xaxt = "n",
     yaxt = "n", ylim = c(0.5,5.5), xlim = c(0,15000))

u <- par('usr')
for(i in 1:5){
  tmp_city <- cities$density[which(cities$size == unq_cat[i])]
  

  #points(tmp_city, y = jitter(rep(6-i, length(tmp_city)), 5), pch = 19, 
   #      col = scales::alpha(cols[i], 0.4 ))
  arrows(city_density$low[i],  6-i, city_density$high[i], 6-i, length = 0.1, angle = 90,code =3,
         col = cols[i], lwd = 4)
  #points(city_density$mden[i], 6-i, pch = 19, 
   #      col = 'black', cex = 4)
  points(city_density$mden[i], 6-i, pch = 21, 
         bg = cols[i], cex = 3.5)
  #points(city_density$mden[i], 6-i, pch = 3, 
   #       cex = 2.5)
}

axis(1, at = seq(0,15000,2500), labels = F, tck = -0.025)
axis(1, at = seq(0,15000,1250), labels = F, tck = -0.025/2)



mtext(unq_cat, 2, at = rev(1:5), line = -1, las = 2, cex = 1.15)
mtext(expression('Avg. population density (thousands per '*km^2*')'), 1, at = 7500, line = 2.2, cex = 1.15)
mtext(seq(0,15,2.5), 1, at = seq(0,15000,2500), line = 0.5)
dev.off()


barplot(log10(city_sum[,2]), horiz = TRUE, xlim = c(0,4), col = cols,
        names.arg = unq_cat, las = 2)
mtext(unq_cat, 2, at = c(1:5)-0.5, line = 2, las = 2)


