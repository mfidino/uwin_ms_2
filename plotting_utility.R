predict.slope <- function(mmat = NULL, 
                        new_data = NULL, species = NULL,
                        city_data = NULL, 
                        new_city_data = NULL,
                        species_there = NULL){
  
  # check new_data
  if(!class(new_data) == 'data.frame'){
    stop('new_data must be a data.frame')
  }
  if(!all(colnames(new_data) == c('hden', 'prophab_btwn', 'hden_btwn'))){
    err <- paste0('The column names in new_data must be (in this order):\n',
                  '\t(1) hden\n',
                  '\t(2) prophab_btwn\n',
                  '\t(3) hden_btwn\n')
    stop(err)
  }
  
  # for scaling
  s_data <- new_data
  # prepare new data for the posterior
  # step 1. scale the data as needed
  if( !all(s_data$prophab_btwn==0) ){
  s_data$prophab_btwn <- (s_data$prophab_btwn - mean(city_data$habitat)) /sd(city_data$habitat)
  }
  if( !all(s_data$hden_btwn==0) ){
  s_data$hden_btwn <- (s_data$hden_btwn - mean(city_data$hden)) /sd(city_data$hden)
  }
  # construct the model matrix
  mm_between <- model.matrix(~prophab_btwn + hden_btwn, s_data)
  
  slope_preds <- mmat[,grep("b2", colnames(mmat))] %*% t(mm_between)
  
  ci_preds <- t(apply(slope_preds, 2, quantile,
                      probs = c(0.025,0.5,0.975)))
  
  
  return(list( mu = ci_preds))
  
}



predict.intercept <- function(mmat = NULL,
                              new_data = NULL, species = NULL,
                              city_data = NULL, 
                              new_city_data = NULL,
                              species_there = NULL){
  
  # check new_data
  if(!class(new_data) == 'data.frame'){
    stop('new_data must be a data.frame')
  }
  if(!all(colnames(new_data) == c('hden', 'prophab_btwn', 'hden_btwn'))){
    err <- paste0('The column names in new_data must be (in this order):\n',
                  '\t(1) hden\n',
                  '\t(2) prophab_btwn\n',
                  '\t(3) hden_btwn\n')
    stop(err)
  }
  
  # for scaling
  s_data <- new_data
  # prepare new data for the posterior
  # step 1. scale the data as needed
  if( !all(s_data$prophab_btwn==0) ){
    s_data$prophab_btwn <- (s_data$prophab_btwn - mean(city_data$habitat)) /sd(city_data$habitat)
  }
  if( !all(s_data$hden_btwn==0) ){
    s_data$hden_btwn <- (s_data$hden_btwn - mean(city_data$hden)) /sd(city_data$hden)
  }
  # construct the model matrix
  mm_intercept <- model.matrix(~prophab_btwn + hden_btwn, s_data)
  # grab the inter
  
  s_data <- model.matrix(~hden*prophab_btwn + hden*hden_btwn, s_data)
  
  intercept_params <- mmat[,grep("b1", colnames(mmat))]
  
  intercept_preds <- intercept_params %*% t(mm_intercept)
  
  ci_preds <- plogis(t(apply(intercept_preds, 2, quantile, probs = c(0.025,0.5,0.975))))
  
  where_species <- which(species_there == 1)
  
  s_citydata <- data.frame(scale(new_city_data))
  s_citydata <- data.frame(apply(s_citydata, 2, function(x) ifelse(is.nan(x), 0, x)))
  s_citydata <- model.matrix(~prophab_btwn + hden_btwn, s_citydata )
  s_citydata <- s_citydata[where_species,]
  
  city_preds <- matrix(0, ncol = 3, nrow = nrow(s_citydata))
  
  # add random effect
  re <- mmat[,grep('B_diff', colnames(mmat))][,where_species]
  
  for(city in 1:nrow(s_citydata)){
    tmp_params <- intercept_params
    tmp_params[,1] <- tmp_params[,1] + re[,city]
    tmp_preds <- tmp_params %*% s_citydata[city,]
    city_preds[city,] <- plogis(quantile(tmp_preds, c(0.025,0.5,0.975)))
  }
  
  #plot(ci_preds[,2] ~ new_data$hden_btwn, type = 'l', ylim = c(0,1))
  #points(city_preds[,2] ~ new_city_data$hden_btwn[test$cities])  
  
  return(list( mu = ci_preds, cmu = city_preds,
               cities = where_species))
  
}



makeplot.hdens <- function(preds = NULL, species = NULL, x = NULL, 
                              intercept = TRUE, cityx = NULL, 
                              species_there = NULL, window = FALSE, pp = FALSE){
  if(intercept){
  my_pred <- preds$mu
  city_mu <- preds$cmu
  cityx <- cityx[species_there == 1]
  } else {
    my_pred <- preds$mu
    
  }
  if(window){
   windows(4,4, xpos = -350, ypos = 200)
  } else {
    if(pp){
      emf(paste0("./plots/",species,"/",species, 
                 ifelse(intercept, "_intercept_popdens", "_slope_popdens"),".emf"),
          height = 4, width = 4, bg = 'white',
          coordDPI = 500)
      
    } else{
  tiff(paste0("./plots/",species,"/",species, 
                       ifelse(intercept, "_intercept_popdens", "_slope_popdens"),".tiff"),
      height = 4, width = 4, res = 300, units = "in", compression = 'lzw')
    }
  }
  par(mar = c(3.5,4,0.75,0.75))
  if(intercept){
  plot(1~1, type ="n", bty = 'l', xlab = "", ylab = "", xaxt = "n",
       yaxt = "n", ylim = c(0,1), xlim = c(200,1400))
  }else{
    plot(1~1, type ="n", bty = 'l', xlab = "", ylab = "", xaxt = "n",
         yaxt = "n", ylim = c(-4,4), xlim = c(200,1400))
    abline(h = 0, lty = 2)
  }
  # plot predicted occupancy
if(intercept){
  for(i in 1:length(cityx)){
    lines(x = rep(cityx[i], 2), y = city_mu[i,-2], lwd = 2,
          col = scales::alpha("#424342", 0.8))
  }
}
    
    x1 <- x
    x2 <- rev(x1)
    y1 <- my_pred[,1]
    y2 <- rev(my_pred[,3])
    polygon(c(x1, x2), c(y1, y2), col = scales::alpha("#32DAC3", .20), border = NA)
    
    
    
    lines(my_pred[,2] ~ x, col = "#32DAC3", lwd = 3)
    if(intercept){
    points(city_mu[,2] ~ cityx, pch = 21, 
           bg = scales::alpha("#424342", 1), cex = 1)
    }
  

  axis(1, at = seq(200, 1400, 200), labels = F, tck = -0.025)
  axis(1, at = seq(200, 1400, 100), labels = F, tck = -0.025/2)
  mtext(text = seq(200, 1400, 400), 1, line = 0.35, at =seq(200, 1400, 400))
  

  if(intercept){
    axis(2, at = seq(0,1, 0.25), labels = F, tck = -0.025)
    axis(2, at = seq(0,1, 0.125), labels = F, tck = -0.025/2)
    mtext(text = sprintf("%.2f",seq(0,1,0.25)), 2, 
          line = 0.5, at = seq(0,1,0.25),las = 1)
  mtext(text = "Average occupancy in city", 2, at = 0.5, cex = 1.2, line = 2.6)
  } else{
    axis(2, at = seq(-4,4, 1), labels = F, tck = -0.025)
    axis(2, at = seq(-4,4, 0.5), labels = F, tck = -0.025/2)
    mtext(text = sprintf("%.0f",seq(-4,4, 1)), 2, 
          line = 0.5, at = seq(-4,4, 1),las = 1)

    

    mtext('Effect of site-level housing density on occupancy',
          side = 2, at = 0, line = 2, cex = 0.9)
 
  }
  mtext(text = expression("Average housing density of city ( houses" * phantom('h') * km^-2 *")"), 
        1, at = 800, cex = 0.9, line = 1.7)
  
  if(!window){
  dev.off()
  }
}



makeplot.habitat <- function(preds = NULL, species = NULL, x = NULL, 
                           intercept = TRUE, cityx = NULL, window = FALSE,
                           pp = FALSE, species_there = NULL){
  if(intercept){
    my_pred <- preds$mu
    city_mu <- preds$cmu
    cityx <- cityx[species_there == 1]
  } else {
    my_pred <- preds$mu
    
  }
  if(window){
    windows(4,4, xpos = -350, ypos = 200)
  } else {
    if(pp){
      emf(paste0("./plots/",species,"/",species, 
                 ifelse(intercept, "_intercept_habitat", "_slope_habitat"),".emf"),
          height = 4, width = 4, bg = 'white',
          coordDPI = 500)
      
    } else{
    tiff(paste0("./plots/",species,"/",species, 
               ifelse(intercept, "_intercept_habitat", "_slope_habitat"),".tiff"),
        height = 4, width = 4, units = "in" , res = 300, compression = 'lzw')
    }
  }
  par(mar = c(3.5,4,0.75,0.75))
  if(intercept){
    plot(1~1, type ="n", bty = 'l', xlab = "", ylab = "", xaxt = "n",
         yaxt = "n", ylim = c(0,1), xlim = c(0.1,0.7))
  }else{
    plot(1~1, type ="n", bty = 'l', xlab = "", ylab = "", xaxt = "n",
         yaxt = "n", ylim = c(-4,4), xlim = c(0.1,0.7))
    abline(h = 0, lty = 2)
  }
  # plot predicted occupancy
  if(intercept){
    for(i in 1:length(cityx)){
      lines(x = rep(cityx[i], 2), y = city_mu[i,-2], lwd = 2,
            col = scales::alpha("#424342", 0.8))
    }
  }
  
  x1 <- x
  x2 <- rev(x1)
  y1 <- my_pred[,1]
  y2 <- rev(my_pred[,3])
  polygon(c(x1, x2), c(y1, y2), col = scales::alpha("#32DAC3", .20), border = NA)
  
  
  
  lines(my_pred[,2] ~ x, col = "#32DAC3", lwd = 3)
  if(intercept){
    points(city_mu[,2] ~ cityx, pch = 21, 
           bg = scales::alpha("#424342", 1), cex = 1)
  }
  
  
  axis(1, at = seq(0.1, 0.7, 0.1), labels = F, tck = -0.025)
  axis(1, at = seq(0.1, 0.7, 0.05), labels = F, tck = -0.025/2)
  mtext(text = seq(0.1, 0.7, 0.1), 1, line = 0.35, at = seq(0.1, 0.7, 0.1))
  
  
  if(intercept){
    axis(2, at = seq(0,1, 0.25), labels = F, tck = -0.025)
    axis(2, at = seq(0,1, 0.125), labels = F, tck = -0.025/2)
    mtext(text = sprintf("%.2f",seq(0,1,0.25)), 2, 
          line = 0.5, at = seq(0,1,0.25),las = 1)
    mtext(text = "Average occupancy in city", 2, at = 0.5, cex = 1.2, line = 2.6)
  } else{
    axis(2, at = seq(-4,4, 1), labels = F, tck = -0.025)
    axis(2, at = seq(-4,4, 0.5), labels = F, tck = -0.025/2)
    mtext(text = sprintf("%.0f",seq(-4,4, 1)), 2, 
          line = 0.5, at = seq(-4,4, 1),las = 1)
    
    
    
    mtext('Effect of site-level housing density on occupancy',
          side = 2, at = 0, line = 2, cex = 0.9)
    
  }
  mtext(text = "Proportion habitat of city", 
        1, at = 0.4, cex = 0.9, line = 1.7)
  
  if(!window){
    dev.off()
  }
}

# makeplot.habitat <- function(preds = NULL, species = NULL, x = NULL, 
#                              intercept = TRUE, cityx = NULL, 
#                              species_there = NULL){
#   
#   my_pred <- preds$mu
#   city_mu <- preds$cmu
#   cityx <- cityx[species_there == 1]
#   #windows(4,4)
#   #tiff(paste0("./plots/",species,"/",species,ifelse(intercept, "_intercept", "_slope"),
#    #           "_habitat.tiff"), height = 4, width = 4, units = "in",
#     #   res = 400, compression = "lzw")
#   svg(paste0("./plots/",species,"/",species,ifelse(intercept, "_intercept", "_slope"),
#              "_habitat.svg"), height = 4, width = 4)
#   par(mar = c(3.5,4,0.5,0.5))
#   if(intercept){
#     plot(1~1, type ="n", bty = 'l', xlab = "", ylab = "", xaxt = "n",
#          yaxt = "n", ylim = c(0,1), xlim = c(0.1,0.7))
#   }else{
#     plot(1~1, type ="n", bty = 'l', xlab = "", ylab = "", xaxt = "n",
#          yaxt = "n", ylim = c(0,2), xlim = c(0.1,0.7))
#     
#   }
#   # plot predicted occupancy
#   if(intercept){
#   for(i in 1:length(cityx)){
#     lines(x = rep(cityx[i], 2), y = city_mu[i,-2], lwd = 2,
#           col = scales::alpha("#424342", 0.8))
#   }
#   }
#   
#   x1 <- x
#   x2 <- rev(x1)
#   y1 <- my_pred[,1]
#   y2 <- rev(my_pred[,3])
#   polygon(c(x1, x2), c(y1, y2), col = scales::alpha("#32DAC3", .20), border = NA)
#   
#   lines(my_pred[,2] ~ x, col = "#32DAC3", lwd = 3)
#   
#   points(city_mu[,2] ~ cityx, pch = 21, 
#          bg = scales::alpha("#424342", 1), cex = 1)
#   
#   axis(1, at = seq(0.1, 0.7, 0.1), labels = F, tck = -0.025)
#   axis(1, at = seq(0.1, 0.7, 0.05), labels = F, tck = -0.025/2)
#   mtext(text = seq(0.1, 0.7, 0.1), 1, line = 0.35, at = seq(0.1, 0.7, 0.1))
#   
#   
#   if(intercept){
#     axis(2, at = seq(0,1, 0.25), labels = F, tck = -0.025)
#     axis(2, at = seq(0,1, 0.125), labels = F, tck = -0.025/2)
#     mtext(text = sprintf("%.2f",seq(0,1,0.25)), 2, 
#           line = 0.5, at = seq(0,1,0.25),las = 1)
#     mtext(text = "Average occupancy in city", 2, at = 0.5, cex = 1.2, line = 2.6)
#   }else{
#     axis(2, at = seq(0,2, 0.5), labels = F, tck = -0.025)
#     axis(2, at = seq(0,2, 0.25), labels = F, tck = -0.025/2)
#     mtext(text = sprintf("%.2f",seq(0,2, 0.5)), 2, 
#           line = 0.5, at = seq(0,2, 0.5),las = 1)
#     mtext(text = "Response to urbanization", 2, at = 1, cex = 1.2, line = 2.6)
#   }
#   mtext(text = "Proportion habitat in city", 1, at = 0.4, cex = 1.2, line = 1.7)
#   
#   dev.off()
#   
# }

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

