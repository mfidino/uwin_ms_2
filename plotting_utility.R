gen_preds.popdens <- function(mmat = NULL, occ_data = NULL, 
                              city_covs = NULL, species = NULL, intercept = TRUE){
  if(intercept){
    occ_data$city <- factor(as.character(occ_data$city))
    z <- occ_data[,species]
    z[z>1] <- 1
    ob_occ <- table(z, occ_data$city)
  naive <- ob_occ[2,] / colSums(ob_occ)
  nse <- sqrt((ob_occ[2,] * naive * (1 - naive)) / colSums(ob_occ))
  nlow <- naive - nse
  nlow[nlow < 0] <- 0
  nhigh <- naive + nse
  nhigh[nhigh>1] <- 1
  my_se <- cbind(nlow, nhigh)
  
  col_loc <- grep("B\\[\\w,1", colnames(mmat))
  city_mu <- apply(plogis(mmat[,col_loc]), 2, HDIofMCMC)
  
  # calc expected mean from group level regression
  x_range <- seq(100, 2000, by = 5)
  x <- (x_range - unlist(attributes(city_covs))[3]) / 
    unlist(attributes(city_covs))[4]
  
  x <- cbind(1, x)
  col_loc <- grep("G\\[1,1\\]|G\\[1,3\\]", colnames(mmat))
  my_pred <- apply(plogis(mmat[,col_loc] %*% t(x)), 2, HDIofMCMC)
  
  real_vals <- as.numeric(DMwR::unscale(city_covs, city_covs))
  
  return(list(z = naive, mu = my_pred, cmu = city_mu, 
              mu.x = x_range, city.x = real_vals, cse = my_se))
  }else{
    
    col_loc <- grep("B\\[\\w,2", colnames(mmat))
    city_mu <- apply(mmat[,col_loc], 2, HDIofMCMC)
    
    # calc expected mean from group level regression
    x_range <- seq(100, 2000, by = 5)
    x <- (x_range - unlist(attributes(city_covs))[3]) / 
      unlist(attributes(city_covs))[4]
    
    x <- cbind(1, x)
    col_loc <- grep("G\\[2,1\\]|G\\[2,3\\]", colnames(mmat))
    my_pred <- apply(mmat[,col_loc] %*% t(x), 2, HDIofMCMC)
    
    real_vals <- as.numeric(DMwR::unscale(city_covs, city_covs))
    
    return(list(z = NA, mu = my_pred, cmu = city_mu, 
                mu.x = x_range, city.x = real_vals, cse = NA))
    
  }
  
  
}


gen_preds.latitude <- function(mmat = NULL, occ_data = NULL, 
                              city_covs = NULL, species = NULL, intercept = TRUE){
  if(intercept){
    occ_data$city <- factor(as.character(occ_data$city))
    z <- occ_data[,species]
    z[z>1] <- 1
    ob_occ <- table(z, occ_data$city)
    naive <- ob_occ[2,] / colSums(ob_occ)
    nse <- sqrt((ob_occ[2,] * naive * (1 - naive)) / colSums(ob_occ))
    nlow <- naive - nse
    nlow[nlow < 0] <- 0
    nhigh <- naive + nse
    nhigh[nhigh>1] <- 1
    my_se <- cbind(nlow, nhigh)
    
    col_loc <- grep("B\\[\\w,1", colnames(mmat))
    city_mu <- apply(plogis(mmat[,col_loc]), 2, HDIofMCMC)
    
    # calc expected mean from group level regression
    x_range <- seq(30, 44, by = 0.1)
    x <- (x_range - unlist(attributes(city_covs))[3]) / 
      unlist(attributes(city_covs))[4]
    
    x <- cbind(1, x)
    col_loc <- grep("G\\[1,1\\]|G\\[1,4\\]", colnames(mmat))
    my_pred <- apply(plogis(mmat[,col_loc] %*% t(x)), 2, HDIofMCMC)
    
    real_vals <- as.numeric(DMwR::unscale(city_covs, city_covs))
    
    return(list(z = naive, mu = my_pred, cmu = city_mu, 
                mu.x = x_range, city.x = real_vals, cse = my_se))
  }else{
    
    col_loc <- grep("B\\[\\w,2", colnames(mmat))
    city_mu <- apply(mmat[,col_loc], 2, HDIofMCMC)
    
    # calc expected mean from group level regression
    x_range <- seq(30, 44, by = 0.1)
    x <- (x_range - unlist(attributes(city_covs))[3]) / 
      unlist(attributes(city_covs))[4]
    
    x <- cbind(1, x)
    col_loc <- grep("G\\[2,1\\]|G\\[2,4\\]", colnames(mmat))
    my_pred <- apply(mmat[,col_loc] %*% t(x), 2, HDIofMCMC)
    
    real_vals <- as.numeric(DMwR::unscale(city_covs, city_covs))
    
    return(list(z = NA, mu = my_pred, cmu = city_mu, 
                mu.x = x_range, city.x = real_vals, cse = NA))
    
  }
  
  
}


gen_preds.habitat <- function(mmat = NULL, occ_data = NULL, 
                              city_covs = NULL, species = NULL, intercept = TRUE){
  if(intercept){
  occ_data$city <- factor(as.character(occ_data$city))
  z <- occ_data[,species]
  z[z>1] <- 1
  ob_occ <- table(z, occ_data$city)
  naive <- ob_occ[2,] / colSums(ob_occ)
  nse <- sqrt((ob_occ[2,] * naive * (1 - naive)) / colSums(ob_occ))
  nlow <- naive - nse
  nlow[nlow < 0] <- 0
  nhigh <- naive + nse
  nhigh[nhigh>1] <- 1
  my_se <- cbind(nlow, nhigh)
  
  col_loc <- grep("B\\[\\w,1", colnames(mmat))
  city_mu <- apply(plogis(mmat[,col_loc]), 2, HDIofMCMC)
  
  # calc expected mean from group level regression
  x_range <- seq(0.10, 0.70, by = 0.002)
  x <- (x_range - unlist(attributes(city_covs))[3]) / 
    unlist(attributes(city_covs))[4]
  
  x <- cbind(1, x)
  col_loc <- grep("G\\[1,1\\]|G\\[1,2\\]", colnames(mmat))
  my_pred <- apply(plogis(mmat[,col_loc] %*% t(x)), 2, HDIofMCMC)
  
  real_vals <- as.numeric(DMwR::unscale(city_covs, city_covs))
  
  return(list(z = naive, mu = my_pred, cmu = city_mu, 
              mu.x = x_range, city.x = real_vals, cse = my_se))
  }else{
    occ_data$city <- factor(as.character(occ_data$city))
    z <- occ_data[,species]
    z[z>1] <- 1
    ob_occ <- table(z, occ_data$city)
    naive <- ob_occ[2,] / colSums(ob_occ)
    nse <- sqrt((ob_occ[2,] * naive * (1 - naive)) / colSums(ob_occ))
    nlow <- naive - nse
    nlow[nlow < 0] <- 0
    nhigh <- naive + nse
    nhigh[nhigh>1] <- 1
    my_se <- cbind(nlow, nhigh)
    
    col_loc <- grep("B\\[\\w,2", colnames(mmat))
    city_mu <- apply(mmat[,col_loc], 2, HDIofMCMC)
    
    # calc expected mean from group level regression
    x_range <- seq(0.10, 0.70, by = 0.002)
    x <- (x_range - unlist(attributes(city_covs))[3]) / 
      unlist(attributes(city_covs))[4]
    
    x <- cbind(1, x)
    col_loc <- grep("G\\[2,1\\]|G\\[2,2\\]", colnames(mmat))
    my_pred <- apply(mmat[,col_loc] %*% t(x), 2, HDIofMCMC)
    
    real_vals <- as.numeric(DMwR::unscale(city_covs, city_covs))
    
    return(list(z = NA, mu = my_pred, cmu = city_mu, 
                mu.x = x_range, city.x = real_vals, cse = NA))
    
  }
  
}

gen_preds.nocorrelates <- function(mmat = NULL, occ_data = NULL, 
                                   species = NULL, intercept = TRUE){
  if(intercept){
    occ_data$city <- factor(as.character(occ_data$city))
    z <- occ_data[,species]
    z[z>1] <- 1
    ob_occ <- table(z, occ_data$city)
    naive <- ob_occ[2,] / colSums(ob_occ)
    nse <- sqrt((ob_occ[2,] * naive * (1 - naive)) / colSums(ob_occ))
    nlow <- naive - nse
    nlow[nlow < 0] <- 0
    nhigh <- naive + nse
    nhigh[nhigh>1] <- 1
    my_se <- cbind(nlow, nhigh)
    
    col_loc <- grep("B\\[\\w,1", colnames(mmat))
    city_mu <- apply(plogis(mmat[,col_loc]), 2, HDIofMCMC)
    
    # calc expected mean from group level regression
    col_loc <- grep("G\\[1,1\\]", colnames(mmat))
    my_pred <- HDIofMCMC(plogis(mmat[,col_loc]))
    
    
    return(list(z = naive, mu = my_pred, cmu = city_mu, cse = my_se))
  }else{
    occ_data$city <- factor(as.character(occ_data$city))
    z <- occ_data[,species]
    z[z>1] <- 1
    ob_occ <- table(z, occ_data$city)
    naive <- ob_occ[2,] / colSums(ob_occ)
    nse <- sqrt((ob_occ[2,] * naive * (1 - naive)) / colSums(ob_occ))
    nlow <- naive - nse
    nlow[nlow < 0] <- 0
    nhigh <- naive + nse
    nhigh[nhigh>1] <- 1
    my_se <- cbind(nlow, nhigh)
    
    col_loc <- grep("B\\[\\w,2", colnames(mmat))
    city_mu <- apply(mmat[,col_loc], 2, HDIofMCMC)
    
    # calc expected mean from group level regression
    col_loc <- grep("G\\[2,1\\]", colnames(mmat))
    my_pred <- HDIofMCMC(mmat[,col_loc])
    
    
    return(list(z = NA, mu = my_pred, cmu = city_mu, cse = NA))
    
  }
  
}

make_plot.popdens <- function(preds = NULL, species = NULL, add_se = FALSE, 
                              intercept = TRUE){
  
  x_range <- preds$mu.x
  my_pred <- preds$mu
  real_vals <- preds$city.x
  city_mu <- preds$cmu
  naive <- preds$z
  my_se <- preds$cse
  # windows(4,4)
  tiff(paste0("./plots/",species,"/",species, 
              ifelse(intercept, "_intercept_popdens", "_slope_popdens"),".tiff"), 
       height = 4, width = 4, units = "in",
       res = 400, compression = "lzw")
  par(mar = c(3.5,4,0.5,0.5))
  if(intercept){
  plot(1~1, type ="n", bty = 'l', xlab = "", ylab = "", xaxt = "n",
       yaxt = "n", ylim = c(0,1), xlim = c(100,2100))
  }else{
    plot(1~1, type ="n", bty = 'l', xlab = "", ylab = "", xaxt = "n",
         yaxt = "n", ylim = c(-2,2), xlim = c(100,2100))
  }
  x1 <- x_range
  x2 <- rev(x1)
  y1 <- my_pred[1,]
  y2 <- rev(my_pred[3,])
  polygon(c(x1, x2), c(y1, y2), col = scales::alpha("#32DAC3", .20), border = NA)
  if(add_se){
    x1l <- real_vals - 20
    x2r <- real_vals + 20
    y1b <- my_se[,1]
    y2t <- my_se[,2]
    rect(x1l, y1b, x2r, y2t, col = scales::alpha("#FEB600", 0.4), border = NA)
  }
  
  lines(my_pred[2,] ~ x_range, col = "#32DAC3", lwd = 3)
  if(intercept){
  for(i in 1:length(real_vals)){
    lines(x = c(real_vals[i] - 20, real_vals[i ] + 20),
          y = rep(naive[i], 2), col = "#FEB600", lwd = 3 )
    
  }
  }
  # plot predicted occupancy
  for(i in 1:length(real_vals)){
    lines(x = rep(real_vals[i], 2), y = city_mu[-2,i], lwd = 2,
          col = scales::alpha("#424342", 0.8))
  }
  
  points(city_mu[2,] ~ real_vals, pch = 21, 
         bg = scales::alpha("#424342", 0.8), cex = 1)
  axis(1, at = seq(100, 2100, 400), labels = F, tck = -0.025)
  axis(1, at = seq(100, 2100, 200), labels = F, tck = -0.025/2)
  mtext(text = seq(100, 2100, 400), 1, line = 0.35, at = seq(100,2100,400))
  

  if(intercept){
    axis(2, at = seq(0,1, 0.25), labels = F, tck = -0.025)
    axis(2, at = seq(0,1, 0.125), labels = F, tck = -0.025/2)
    mtext(text = sprintf("%.2f",seq(0,1,0.25)), 2, 
          line = 0.5, at = seq(0,1,0.25),las = 1)
  mtext(text = "Average occupancy rate", 2, at = 0.5, cex = 1.2, line = 2.6)
  } else{
    axis(2, at = seq(-2,2, 1), labels = F, tck = -0.025)
    axis(2, at = seq(-2,2, 0.5), labels = F, tck = -0.025/2)
    mtext(text = sprintf("%.2f",seq(-2,2,1)), 2, 
          line = 0.5, at = seq(-2,2,1),las = 1)
    mtext(text = "Response to urbanization", 2, at = 0, cex = 1.2, line = 2.6)
  }
  mtext(text = "Population density in study area", 1, at = 1100, cex = 1.2, line = 1.7)
  #points(naive ~ real_vals, pch = 23, bg = scales::alpha("#FEB600", 0.5), cex = 1.5)
  
  dev.off()
}


make_plot.latitude <- function(preds = NULL, species = NULL, add_se = FALSE, 
                              intercept = TRUE){
  
  x_range <- preds$mu.x
  my_pred <- preds$mu
  real_vals <- preds$city.x
  city_mu <- preds$cmu
  naive <- preds$z
  my_se <- preds$cse
  # windows(4,4)
  tiff(paste0("./plots/",species, "/", species,
              ifelse(intercept, "_intercept_latitude", "_slope_latitude"),".tiff"), 
       height = 4, width = 4, units = "in",
       res = 400, compression = "lzw")
  par(mar = c(3.5,4,0.5,0.5))
  if(intercept){
    plot(1~1, type ="n", bty = 'l', xlab = "", ylab = "", xaxt = "n",
         yaxt = "n", ylim = c(0,1), xlim = c(30,44))
  }else{
    plot(1~1, type ="n", bty = 'l', xlab = "", ylab = "", xaxt = "n",
         yaxt = "n", ylim = c(-2,2), xlim = c(30,44))
  }
  x1 <- x_range
  x2 <- rev(x1)
  y1 <- my_pred[1,]
  y2 <- rev(my_pred[3,])
  polygon(c(x1, x2), c(y1, y2), col = scales::alpha("#32DAC3", .20), border = NA)
  if(add_se){
    x1l <- real_vals - 20
    x2r <- real_vals + 20
    y1b <- my_se[,1]
    y2t <- my_se[,2]
    rect(x1l, y1b, x2r, y2t, col = scales::alpha("#FEB600", 0.4), border = NA)
  }
  
  lines(my_pred[2,] ~ x_range, col = "#32DAC3", lwd = 3)
  if(intercept){
    for(i in 1:length(real_vals)){
      lines(x = c(real_vals[i] - 0.25, real_vals[i ] + 0.25),
            y = rep(naive[i], 2), col = "#FEB600", lwd = 3 )
      
    }
  }
  # plot predicted occupancy
  for(i in 1:length(real_vals)){
    lines(x = rep(real_vals[i], 2), y = city_mu[-2,i], lwd = 2,
          col = scales::alpha("#424342", 0.8))
  }
  
  points(city_mu[2,] ~ real_vals, pch = 21, 
         bg = scales::alpha("#424342", 0.8), cex = 1)
  axis(1, at = seq(30, 44, 1), labels = F, tck = -0.025)
  axis(1, at = seq(30, 44, 0.5), labels = F, tck = -0.025/2)
  mtext(text = seq(30, 44, 2), 1, line = 0.35, at = seq(30, 44, 2))
  
  
  if(intercept){
    axis(2, at = seq(0,1, 0.25), labels = F, tck = -0.025)
    axis(2, at = seq(0,1, 0.125), labels = F, tck = -0.025/2)
    mtext(text = sprintf("%.2f",seq(0,1,0.25)), 2, 
          line = 0.5, at = seq(0,1,0.25),las = 1)
    mtext(text = "Average occupancy rate", 2, at = 0.5, cex = 1.2, line = 2.6)
  } else{
    axis(2, at = seq(-2,2, 1), labels = F, tck = -0.025)
    axis(2, at = seq(-2,2, 0.5), labels = F, tck = -0.025/2)
    mtext(text = sprintf("%.2f",seq(-2,2,1)), 2, 
          line = 0.5, at = seq(-2,2,1),las = 1)
    mtext(text = "Response to urbanization", 2, at = 0, cex = 1.2, line = 2.6)
  }
  mtext(text = "City latitude", 1, at = 37, cex = 1.2, line = 1.7)
  #points(naive ~ real_vals, pch = 23, bg = scales::alpha("#FEB600", 0.5), cex = 1.5)
  
  dev.off()
}

make_plot.habitat <- function(preds = NULL, species = NULL, add_se = FALSE,
                              intercept = TRUE){
  
  x_range <- preds$mu.x
  my_pred <- preds$mu
  real_vals <- preds$city.x
  city_mu <- preds$cmu
  naive <- preds$z
  my_se <- preds$cse
  # windows(4,4)
  tiff(paste0("./plots/",species,"/",species,ifelse(intercept, "_intercept", "_slope"),
              "_habitat.tiff"), height = 4, width = 4, units = "in",
       res = 400, compression = "lzw")
  par(mar = c(3.5,4,0.5,0.5))
  if(intercept){
  plot(1~1, type ="n", bty = 'l', xlab = "", ylab = "", xaxt = "n",
       yaxt = "n", ylim = c(0,1), xlim = c(0.1,0.7))
  }else{
    plot(1~1, type ="n", bty = 'l', xlab = "", ylab = "", xaxt = "n",
         yaxt = "n", ylim = c(-2,2), xlim = c(0.1,0.7))
    
  }
  x1 <- x_range
  x2 <- rev(x1)
  y1 <- my_pred[1,]
  y2 <- rev(my_pred[3,])
  polygon(c(x1, x2), c(y1, y2), col = scales::alpha("#32DAC3", .20), border = NA)
  
  lines(my_pred[2,] ~ x_range, col = "#32DAC3", lwd = 3)
  if(intercept){
  for(i in 1:length(real_vals)){
    lines(x = c(real_vals[i] - 0.01, real_vals[i ] + 0.01),
          y = rep(naive[i], 2), col = "#FEB600", lwd = 3 )
    
  }
  }
  # plot predicted occupancy
  for(i in 1:length(real_vals)){
    lines(x = rep(real_vals[i], 2), y = city_mu[-2,i], lwd = 2,
          col = scales::alpha("#424342", 0.8))
  }
  
  points(city_mu[2,] ~ real_vals, pch = 21, 
         bg = scales::alpha("#424342", 0.8), cex = 1)
  axis(1, at = seq(0.1, 0.7, 0.1), labels = F, tck = -0.025)
  axis(1, at = seq(0.1, 0.7, 0.05), labels = F, tck = -0.025/2)
  mtext(text = seq(0.1, 0.7, 0.1), 1, line = 0.35, at = seq(0.1, 0.7, 0.1))
  

  if(intercept){
    axis(2, at = seq(0,1, 0.25), labels = F, tck = -0.025)
    axis(2, at = seq(0,1, 0.125), labels = F, tck = -0.025/2)
    mtext(text = sprintf("%.2f",seq(0,1,0.25)), 2, 
          line = 0.5, at = seq(0,1,0.25),las = 1)
  mtext(text = "Average occupancy rate", 2, at = 0.5, cex = 1.2, line = 2.6)
  }else{
    axis(2, at = seq(-2,2, 1), labels = F, tck = -0.025)
    axis(2, at = seq(-2,2, 0.5), labels = F, tck = -0.025/2)
    mtext(text = sprintf("%.2f",seq(-2,2,1)), 2, 
          line = 0.5, at = seq(-2,2,1),las = 1)
    mtext(text = "Response to urbanization", 2, at = 0, cex = 1.2, line = 2.6)
  }
  mtext(text = "Proportion habitat in study area", 1, at = 0.4, cex = 1.2, line = 1.7)
  #points(naive ~ real_vals, pch = 23, bg = scales::alpha("#FEB600", 0.5), cex = 1.5)
  
  dev.off()
}

make_plot.nocorrelates <- function(preds = NULL, species = NULL, add_se = FALSE,
                                   intercept = TRUE, city_names = NULL){
  
  my_pred <- preds$mu
  city_mu <- preds$cmu
  naive <- preds$z
  my_se <- preds$cse

  tiff(paste0("./plots/",species,"/",species,ifelse(intercept, "_intercept", "_slope"),
              "_nocorrelates.tiff"), height = 4, width = 4, units = "in",
       res = 400, compression = "lzw")
  par(mar = c(3.5,4,0.5,0.5))
  if(intercept){
    plot(1~1, type ="n", bty = 'l', xlab = "", ylab = "", xaxt = "n",
         yaxt = "n", ylim = c(0,1), xlim = c(0.5,9.5))
  }else{
    plot(1~1, type ="n", bty = 'l', xlab = "", ylab = "", xaxt = "n",
         yaxt = "n", ylim = c(-2,2), xlim = c(0.5,9.5))
    
  }
  x1 <- seq(0.5, 9.5, length.out = length(my_pred))
  x2 <- rev(x1)
  y1 <- rep(my_pred[1],3)
  y2 <- rev(rep(my_pred[3], 3))
  polygon(c(x1, x2), c(y1, y2), col = scales::alpha("#32DAC3", .20), border = NA)
  
  lines(rep(my_pred[2],3) ~ x1, col = "#32DAC3", lwd = 3)
  
  # plot predicted occupancy
  for(i in 1:ncol(city_mu)){
    lines(x = rep(i, 2), y = city_mu[-2,i], lwd = 2,
          col = scales::alpha("#424342", 0.8))
  }
  
  points(city_mu[2,] ~ c(1:9), pch = 21, 
         bg = scales::alpha("#424342", 0.8), cex = 1)
  axis(1, at = seq(1, 9, 1), labels = F, tck = -0.025)
  mtext(text = toupper(city_names), 1, line = 0.35, 
        at = seq(1, 9, 1), cex = 0.65)
  
  
  if(intercept){
    axis(2, at = seq(0,1, 0.25), labels = F, tck = -0.025)
    axis(2, at = seq(0,1, 0.125), labels = F, tck = -0.025/2)
    mtext(text = sprintf("%.2f",seq(0,1,0.25)), 2, 
          line = 0.5, at = seq(0,1,0.25),las = 1)
    mtext(text = "Average occupancy rate", 2, at = 0.5, cex = 1.2, line = 2.6)
  }else{
    axis(2, at = seq(-2,2, 1), labels = F, tck = -0.025)
    axis(2, at = seq(-2,2, 0.5), labels = F, tck = -0.025/2)
    mtext(text = sprintf("%.2f",seq(-2,2,1)), 2, 
          line = 0.5, at = seq(-2,2,1),las = 1)
    mtext(text = "Response to urbanization", 2, at = 0, cex = 1.2, line = 2.6)
  }
  mtext(text = "City", 1, at = 5, cex = 1.2, line = 1.7)
  #points(naive ~ real_vals, pch = 23, bg = scales::alpha("#FEB600", 0.5), cex = 1.5)
  
  dev.off()
}


make_plot.nocorrelates <- function(preds = NULL, species = NULL, add_se = FALSE,
                                   intercept = TRUE){
  
  my_pred <- preds$mu
  city_mu <- preds$cmu
  naive <- preds$z
  my_se <- preds$cse
  tiff(paste0("./plots/",species,"/",species,ifelse(intercept, "_intercept", "_slope"),
              "_nocorrelates.tiff"), height = 4, width = 4, units = "in",
       res = 400, compression = "lzw")
  par(mar = c(3.5,4,0.5,0.5))
  if(intercept){
    plot(1~1, type ="n", bty = 'l', xlab = "", ylab = "", xaxt = "n",
         yaxt = "n", ylim = c(0,1), xlim = c(0.5,9.5))
  }else{
    plot(1~1, type ="n", bty = 'l', xlab = "", ylab = "", xaxt = "n",
         yaxt = "n", ylim = c(-2,2), xlim = c(0.5,9.5))
    
  }
  x1 <- seq(0.5, 9.5, length.out = length(my_pred))
  x2 <- rev(x1)
  y1 <- rep(my_pred[1],3)
  y2 <- rev(rep(my_pred[3], 3))
  polygon(c(x1, x2), c(y1, y2), col = scales::alpha("#32DAC3", .20), border = NA)
  
  lines(rep(my_pred[2],3) ~ x1, col = "#32DAC3", lwd = 3)
  
  # plot predicted occupancy
  for(i in 1:ncol(city_mu)){
    lines(x = rep(i, 2), y = city_mu[-2,i], lwd = 2,
          col = scales::alpha("#424342", 0.8))
  }
  
  points(city_mu[2,] ~ c(1:9), pch = 21, 
         bg = scales::alpha("#424342", 0.8), cex = 1)
  axis(1, at = seq(1, 9, 1), labels = F, tck = -0.025)
  mtext(text = toupper(rownames(preds$cse)), 1, line = 0.35, 
        at = seq(1, 9, 1), cex = 0.65)
  
  
  if(intercept){
    axis(2, at = seq(0,1, 0.25), labels = F, tck = -0.025)
    axis(2, at = seq(0,1, 0.125), labels = F, tck = -0.025/2)
    mtext(text = sprintf("%.2f",seq(0,1,0.25)), 2, 
          line = 0.5, at = seq(0,1,0.25),las = 1)
    mtext(text = "Average occupancy rate", 2, at = 0.5, cex = 1.2, line = 2.6)
  }else{
    axis(2, at = seq(-2,2, 1), labels = F, tck = -0.025)
    axis(2, at = seq(-2,2, 0.5), labels = F, tck = -0.025/2)
    mtext(text = sprintf("%.2f",seq(-2,2,1)), 2, 
          line = 0.5, at = seq(-2,2,1),las = 1)
    mtext(text = "Response to urbanization", 2, at = 0, cex = 1.2, line = 2.6)
  }
  mtext(text = "City", 1, at = 5, cex = 1.2, line = 1.7)
  #points(naive ~ real_vals, pch = 23, bg = scales::alpha("#FEB600", 0.5), cex = 1.5)
  
  dev.off()
}




my_vioplot <- function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, 
          horizontal = FALSE, col = "magenta", border = "black", lty = 1, 
          lwd = 1, rectCol = "black", colMed = "white", pchMed = 19, 
          at, add = FALSE, wex = 1, drawRect = TRUE) 
{
  datas <- list(x, ...)
  n <- length(datas)
  if (missing(at)) 
    at <- 1:n
  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)
  q1 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  med <- vector(mode = "numeric", length = n)
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  baserange <- c(Inf, -Inf)
  args <- list(display = "none")
  if (!(is.null(h))) 
    args <- c(args, h = h)
  for (i in 1:n) {
    data <- datas[[i]]
    data.min <- min(data)
    data.max <- max(data)
    q1[i] <- quantile(data, 0.25)
    q3[i] <- quantile(data, 0.75)
    med[i] <- median(data)
    iqd <- q3[i] - q1[i]
    upper[i] <- min(q3[i] + range * iqd, data.max)
    lower[i] <- max(q1[i] - range * iqd, data.min)
    est.xlim <- c(min(lower[i], data.min), max(upper[i], 
                                               data.max))
    smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
                                     args))
    hscale <- 0.4/max(smout$estimate) * wex
    base[[i]] <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
  }
  if (!add) {
    xlim <- if (n == 1) 
      at + c(-0.5, 0.5)
    else range(at) + min(diff(at))/2 * c(-1, 1)
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  }
  else {
    label <- names
  }
  boxwidth <- 0.05 * wex
  if (!add) 
    plot.new()
  if (!horizontal) {
    if (!add) {
      plot.window(xlim = xlim, ylim = ylim)
      axis(2)
      axis(1, at = at, label = label)
    }
   # box()
    for (i in 1:n) {
      polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
              c(base[[i]], rev(base[[i]])), col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
              lty = lty)
        rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, 
             q3[i], col = rectCol)
        points(at[i], med[i], pch = pchMed, col = colMed)
      }
    }
  }
  else {
    if (!add) {
      plot.window(xlim = ylim, ylim = xlim)
      axis(1)
      axis(2, at = at, label = label)
    }
   # box()
    for (i in 1:n) {
      polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], 
                                              rev(at[i] + height[[i]])), col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
              lty = lty)
        rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + 
               boxwidth/2, col = rectCol)
        points(med[i], at[i], pch = pchMed, col = colMed)
      }
    }
  }
  invisible(list(upper = upper, lower = lower, median = med, 
                 q1 = q1, q3 = q3))
}

