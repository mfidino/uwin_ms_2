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
                       species_there = det_events[,sp_name[species]],
                       model = cpo$best[species])
  if(cpo$best[species] != 'habitat'){
  makeplot.hdens(preds = preds, species = sp_name[species],x = new_data$hden_btwn,
                   cityx = cdat$hden, species_there = det_events[,sp_name[species]],
                 pp = FALSE)
  }
  
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
                                  species_there = det_events[,sp_name[species]],
                             model = cpo$best[species])
  if(cpo$best[species] != 'housing_density'){
  makeplot.habitat(preds = preds, species = sp_name[species],x = new_data$prophab_btwn,
                   cityx = cdat$habitat, species_there = det_events[,sp_name[species]],
                   pp = FALSE)
  }
  
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
                             species_there = det_events[,sp_name[species]],
                         model = cpo$best[species])
  if(cpo$best[species] != 'habitat'){
  makeplot.hdens(preds = preds, species = sp_name[species],x = new_data$hden_btwn,
                 cityx = cdat$hden, species_there = det_events[,sp_name[species]],
                 intercept = FALSE, window = FALSE, pp = FALSE)
  }
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
                             species_there = det_events[,sp_name[species]],
                         model = cpo$best[species])
  if(cpo$best[species] != 'housing_density'){
  makeplot.habitat(preds = preds, species = sp_name[species],x = new_data$prophab_btwn,
                   cityx = cdat$habitat,
                   intercept = FALSE, window = FALSE, pp = FALSE)
  }
  
}




  
  #### plotting for supp mater A
  
  # 
  
  library(vioplot)
  
  supp_plot_data <- patch_covs
  supp_plot_data$hd_1000 <- supp_plot_data$hd_1000 / 1000
  windows(6,8)
  tiff("./plots/supp_mater/hdenrange.tiff", height = 6, width = 8,
       units = "in", res = 600, compression = "lzw")
  par(mar = c(5,7,0.5,0.5), usr =c(0,10,0,10) )
  plot(1~1, type ="n", bty = 'l', xlab = "", ylab = "", xaxt = "n",
       yaxt = "n", ylim = c(0.368,10 * 0.965), xlim = c(0.368,10 * 0.965))
  par("usr")
  
  axis(2, at = seq(0.25, 9.25, 1), labels = F, tck = -0.025)
  
  # get names in the correct order
  
  mtext(text = cplot$pretty[order(cplot$hden)], 
        2, 
        line = 1.25,
        at = seq(0.25,9.5,1),
        las = 1,
        cex = 0.9
  )
  axis(1, at = seq(0,10, 1), labels = F, tck = -0.025)
  axis(1, at = seq(0,10, 1/2), labels = F, tck = -0.025/2)
  mtext(text = sprintf("%.f",seq(0,10,1)),
        1,
        line = 0.75,
        at = 0:10)
  
  mtext(expression("Site-level housing density (1000 houses" * phantom('h') * km^-2 *")"),
        1,
        at = mean(0:10),
        line = 3,
        cex = 1.5
        )
  # plot them my median
  
  to_plot <- cdat$city[order(cdat$hden)]
  u <- par("usr")
  for(i in 0:10){
    lines(x= rep(i, 2),
          y = c(u[1],u[2]),
          col = scales::alpha("#424342", 0.5))
  }
  
  for(i in 1:10){
  my_vioplot(
    supp_plot_data$hd_1000[patch_covs$city== to_plot[i]],
    at = i-1,
    side = "right",
    horizontal = TRUE,
    add = TRUE,
    wex = 2,
    col = "#32DAC3"
  )
  }
  
 dev.off()

