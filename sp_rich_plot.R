
svg('./plots/sp_rich_figure.svg', height = 6, width = 9)
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



# figure out most common community in low and high

# step 1. select a city





#lower <- lower[cplot$ord,]
#upper <- upper[cplot$ord,]
#detplot <- det_events
#detplot <- detplot[cplot$ord,order(colnames(detplot))]

my_col_low <- cgrad_low(100)[100]
my_col_high <- cgrad_high(100)[100]

#nlow <- city_rich_low[cplot$ord, 2]
#nhigh <- city_rich_high[cplot$ord,2]  
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

  # rect(fright - 0.4-0.15-my_new_wid[1], ylm ,
  #      fright-0.4 -0.15+my_new_wid[1], ylm + 1,
  #      border = FALSE, col = my_col_high[1])
  # 
  # rect(fright - 0.4-0.15-my_new_wid[2], ylm - 1.25 ,
  #      fright-0.4 -0.15+my_new_wid[2], ylm -.25 ,
  #      border = FALSE, col = my_col_high[2])
  # 
  # rect(fright - 0.4-0.15-my_new_wid[3], ylm - 2.5 ,
  #      fright-0.4 -0.15+my_new_wid[3], ylm -1.5 ,
  #      border = FALSE, col = my_col_high[3])
  # 
  


# text(x = fright - 0.17, y = ylm + 1.8, 'Habitat patches\nsurrounded by:', cex = 1.2)
# lines(x = c(fright - 0.8, fright + 0.4), y = c(ylm + 0.8, ylm + 0.8))
# 
# text(x = fright, y = ylm,'avg. housing\ndensity')
# text(x = fright, y = ylm-1.5,'avg. housing\ndensity')
# text(x = fright - 0.55, y = ylm, expression(underline('<')), cex = 1.8)
# text(x = fright - 0.55, y = ylm-1.5, '>', cex = 1.8)
# points(x = fright - 0.75, y =ylm, pch = 24, bg = my_col_low, cex = 2)
# points(x = fright - 0.75, y = ylm-1.5, pch = 22, bg = my_col_high, cex = 2)
# 
# text(x = fright - 0.17, y = ylm - 4, 'Confidence', cex = 1.2)
# lines(x = c(fright - 0.8, fright + 0.4), y = c(ylm - 4.5, ylm - 4.5))

