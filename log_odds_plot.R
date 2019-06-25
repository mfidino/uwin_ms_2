# source('sourcer.R')
# 
# # bring in the species results, grab just the columns we need
# #  and store the quantiles in a data.frame
# 
# 
# which_folder <- "reparam_longrun"
# 
# my_files <- list.files(paste0("./results/", which_folder),
#                        pattern = "matrix", 
#                        full.names = TRUE)
# 
# sp_name <- gsub(".*/(\\w+)_matrix.csv", "\\1", my_files)
# 
# sp_results <- vector('list', length= nspecies)
# names(sp_results) <- sp_name
# for(species in seq_len(nspecies)){
#   sp_results[[species]] <- data.table::fread(my_files[species],
#                                       data.table = FALSE,
#                                       select = c('b1[1]', 'b1[2]', 'b1[3]',
#                                                  'b2[1]', 'b2[2]', 'b2[3]')) %>% 
#     as.matrix(.) %>% 
#     apply(., 2, quantile, probs = c(0.025,0.5,0.975)) %>% 
#     t
# }
# 
# windows(7,5)

#library(jpeg)
my_images <- c("./images/coyote.JPG",
               "./images/squirrel.JPG",
               "./images/squirrel.JPG",
               "./images/opossum.JPG",
               "./images/rabbit.JPG",
               './images/raccoon.JPG',
               './images/red_fox.jpg',
               './images/skunk.JPG')

tiff("./plots/log_odds.tiff", height = 5, width = 7, units = "in",
     res = 300, compression = 'lzw')

 m <- matrix(1:2, ncol = 2, nrow = 1)
 layout(m)
par(mar = c(5,7,0.5,0), xpd = NA)

plot(1~1, type ="n", bty = 'n', xlab = "", ylab = "", xaxt = "n",
     yaxt = "n", ylim = c(0,nspecies), xlim = c(-3,3),
     usr = c(-3.24,3.24,0.5,8.5))

u <- par("usr")
to_drop <- (abs(u[3])*2)/nspecies
for(box in seq(2,8,2)){
  rect(xleft = u[1], ybottom = box - 1, xright = u[2], ytop = box,
       col = '#C9C6C3', border = NA)
}
#text(x = u[1], y = ci-0.125, labels = pretty_sp[ci], cex = 0.8,
 #    pos = 4)
lines(x = c(0,0), y = c(u[3],8), lty = 2)
for(ci in seq_len(8)){
  my_img <- readJPEG(my_images[ci])
  rasterImage(my_img, -4.8, ci-0.83, -3.5,ci-0.165 )
  }

mtext(text = pretty_sp, 2, las = 2, at = (1:8)-0.5, line = 3.1, cex = 0.8)

for(ci in seq_len(8)){

  tmp_res <- sp_results[[ci]]
  arrows(tmp_res[1,1],  ci-0.165, tmp_res[1,3], ci-0.165, length = 0.03, angle = 90,code =3,
         col = "#24d5f7ff", lwd = 2)
  arrows(tmp_res[2,1],  ci-0.5, tmp_res[2,3], ci-0.5, length = 0.03, angle = 90, code = 3,
         col = "#5ee38bff", lwd = 2)
  arrows(tmp_res[3,1],  ci-0.83, tmp_res[3,3], ci-0.83, length = 0.03, angle = 90, code=3,
         col = "#ffb226ff", lwd=2)
  points(tmp_res[1:3,2], c(ci-0.165,ci-0.5,ci-0.83), pch = c(21,22,23), 
         bg = c("#24d5f7ff", "#5ee38bff","#ffb226ff" ), cex = 1)
  }

#abline(v = 0, lty = 2)

axis(1, at = seq(-3,3, 1), labels = F, tck = -0.025)
axis(1, at = seq(-3,3, 0.5), labels = F, tck = -0.025/2)
mtext(text = sprintf("%.0f",seq(-3,3, 1)), 1, 
      line = 0.5, at = seq(-3,3, 1),las = 1)
mtext('Between-city effects on\naverage occupancy',1,at = 0, line = 2.5,cex=0.9)
mtext('(a)', 3, at = u[1]+0.2, line = -0.5 )

par(mar = c(5, 1,0.5,6), xpd = NA)
plot(1~1, type ="n", bty = 'n', xlab = "", ylab = "", xaxt = "n",
     yaxt = "n", ylim = c(0,nspecies), xlim = c(-3,3))
mtext('(b)', 3, at = u[1]+0.2, line = -0.5 )

u <- par("usr")
for(box in seq(2,8,2)){
  rect(xleft = u[1], ybottom = box - 1, xright = u[2], ytop = box,
       col = '#C9C6C3', border = NA)
}
lines(x = c(0,0), y = c(u[3],8), lty = 2)
for(ci in seq_len(8)){
  
  tmp_res <- sp_results[[ci]][4:6,]
  arrows(tmp_res[1,1],  ci-0.165, tmp_res[1,3], ci-0.165, length = 0.03, angle = 90,code =3,
         col = "#24d5f7ff", lwd = 2)
  arrows(tmp_res[2,1],  ci-0.5, tmp_res[2,3], ci-0.5, length = 0.03, angle = 90, code = 3,
         col = "#5ee38bff", lwd = 2)
  arrows(tmp_res[3,1],  ci-0.83, tmp_res[3,3], ci-0.83, length = 0.03, angle = 90, code=3,
         col = "#ffb226ff", lwd = 2)
  points(tmp_res[1:3,2], c(ci-0.165,ci-0.5,ci-0.83), pch = c(21,22,23), 
         bg = c("#24d5f7ff", "#5ee38bff","#ffb226ff" ), cex = 1)
}

axis(1, at = seq(-3,3, 1), labels = F, tck = -0.025)
axis(1, at = seq(-3,3, 0.5), labels = F, tck = -0.025/2)
mtext(text = sprintf("%.0f",seq(-3,3, 1)), 1, 
      line = 0.5, at = seq(-3,3, 1),las = 1)
mtext('Between-city effects on\ninfluence of site-level housing\ndensity on occupancy',
      1,at = 0, line = 3.35, cex = 0.9)

legend(x = 3.3, y = 6.6, legend = c('\nIntercept\n',
                                  '\nProportion\nof habitat\n',
                                  '\nAvg. housing\ndensity\n'),
       pch = 21:23, pt.bg =c("#24d5f7ff", "#5ee38bff","#ffb226ff" ), bty = "n",
       cex = 0.8, pt.cex = 1.4)

dev.off()




