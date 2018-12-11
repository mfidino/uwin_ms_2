# Locate RDS files

source("sourcer.R")

model <- "global"
folder <- paste0("./results/", model,"/")
my_rds <- list.files(folder, pattern = "RDS", full.names = TRUE)
my_rds <- my_rds[-grep("waic", my_rds)]
my_rds_short <- sapply(strsplit(my_rds, "/"), "[", 4)
my_rds_short <- sapply(strsplit(my_rds_short, "\\."), "[", 1)

for(species in 1:length(my_rds)){
  tmp_rds <- readRDS(my_rds[species])
  tmp_mat <- as.matrix(as.mcmc.list(tmp_rds), chains = TRUE)
  # drop the z stuff
  tmp_mat <- tmp_mat[,-c(1,grep("z", colnames(tmp_mat)))]
  tmp_mat <- as.data.frame(tmp_mat)
  # save the file
  data.table::fwrite(tmp_mat, paste0(folder, my_rds_short[species],
                            "_matrix.csv"))
}

my_res <- list.files(folder , pattern = "csv", full.names = TRUE)

# bring in the city data
cdat <- read.csv("data/city_level_data.csv", stringsAsFactors = FALSE)
cdat <- cdat[order(cdat$city),]

# make a plot for each species
for(species in 1:length(my_res)){
  # read in the file
  my_mcmc <- data.table::fread(my_res[species], data.table = FALSE) %>% 
    apply(., 2, HDIofMCMC, .90) %>% t(.)
  # round them 
  rmcmc <- round(my_mcmc, 2)
  # paste them together
  pmcmc <- apply(rmcmc, 1, function(x) {
    paste0(sprintf("%.2f", x[2])," (", 
           sprintf("%.2f", x[1]), " - ", 
           sprintf("%.2f", x[3]), ")")})
  
  g_res <- matrix(pmcmc[grep("^G", names(pmcmc))], ncol = 4, nrow = 2)
  
  # get species name
  g_res <- cbind(my_rds_short[species], c("City intercept", "URB slope"), g_res)
  colnames(g_res) <- c("Species", "City_pars",
                       "Intercept", "Habitat amount", 
                       "Population density", "Latitude")
  
  # standard deviation now
  
  sd_res <- matrix(pmcmc[grep("sigma\\.urb|sigma\\.int|rho_cor", names(pmcmc))], 
                   ncol = 3, 
                   nrow = 1)
  sd_res <- cbind(my_rds_short[species],sd_res)
  colnames(sd_res) <- c("Species","Correlation",  "Intercept", "URB slope")
  
  # city specific now
  
  b_res <- matrix(pmcmc[grep("^B", names(pmcmc))], ncol = 2, nrow = 9)
  
  b_res <- cbind(my_rds_short[species], cdat$city, b_res)
  colnames(b_res) <- c("Species", "City", "Intercept", "URB slope")
  
  
  if(species == 1){
    write.csv(g_res,"./results_summary/city_differences_90.csv", row.names = FALSE)
    write.csv(sd_res, "./results_summary/city_sd_90.csv", row.names = FALSE)
    write.csv(b_res, "./results_summary/within_city_90.csv", row.names = FALSE)
  } else{
    write.table(g_res,"./results_summary/city_differences_90.csv", row.names = FALSE,
                append = TRUE, sep = ",", col.names = FALSE)
    write.table(sd_res,"./results_summary/city_sd_90.csv", row.names = FALSE,
                append = TRUE, sep = ",", col.names = FALSE)
    write.table(b_res, "./results_summary/within_city_90.csv", row.names = FALSE, 
                col.names = FALSE,
                append = TRUE, sep = ",")
  }
  
  
}

g_res <- vector("list", length = length(my_res))
for(species in 1:8 ){
  g_res[[species]] <- fread(my_res[species], data.table = FALSE) 
  

}



# get_mu
my_mu <- matrix(0, ncol = 8, nrow = nrow(g_res[[1]]))

for(i in 1:8){
  my_mu[,i] <- g_res[[i]][,grep("G\\[1,1", colnames(g_res[[i]]))]
}

my_mu <- plogis(my_mu)
mumed <- apply(my_mu, 2, median)

sp_or <- base::order(mumed, decreasing = TRUE)

sp_name_or <- my_rds_short[sp_or]

my_mu <- my_mu[,sp_or]

sp_name_or
fancy_names <- c("Raccoon\n", "Fox\nsquirrel", "Virginia\nopossum" , "Eastern\ncottontail",
                 "Coyote\n", "Gray\nsquirrel", "Striped\nskunk","Red fox\n" )

windows(6,6)
tiff("./plots/avg_occupancy.tiff", height = 6, width = 6, units = "in",
     res = 400, compression = "lzw")
par(mar = c(3,4,0.5,0.5))
plot(1~1, type ="n", bty = 'l', xlab = "", ylab = "", xaxt = "n",
     yaxt = "n", ylim = c(0,1), xlim = c(0.5,8.5))
axis(1, at = seq(1, 9), tck = -0.025/2, labels = FALSE)
mtext(text = fancy_names, side = 1, at = seq(1, 8), line = 1.5, cex = 0.75)
axis(2, at = seq(0,1, 0.25), labels = F, tck = -0.025)
axis(2, at = seq(0,1, 0.125), labels = F, tck = -0.025/2)
mtext(text = sprintf("%.2f",seq(0,1,0.25)), 2, 
      line = 0.8, at = seq(0,1,0.25),las = 1)
mtext(text = "Average occupancy rate", 2, at = 0.5, line = 2.8, cex = 1.3)
my_col <- c('#9c5be8','#ae74ec','#be8ef0','#cea7f4','#ddbff7','#ebd8fb','#f9f2fe',
            '#f9f2fe')
for(i in 1:8){
  to_plot <- my_mu[,i]
  my_hdi <- HDIofMCMC(to_plot)
  to_plot <- to_plot[between(to_plot, my_hdi[1], my_hdi[3])]
my_vioplot(to_plot, at = i, ylim = c(0,1), 
           col = my_col[i], add = TRUE, wex = 0.8, drawRect = FALSE)
points(median(to_plot) ~ i, pch = "-", cex = 3)
}
dev.off()

image <- image_read("./plots/avg_occupancy.tiff")

text(x = seq(1,7)-0.48, par("usr")[3] - 0.07, 
     labels = fancy_names, srt = 35, pos = 1, xpd = TRUE)
  text = sp_name_or, 1, at = seq(1,7))

pinch <- function(wid = NULL, dpi = NULL){
  return( wid * dpi)
  
}

image <- image_scale(image, "1600")

image_write(image, "./plots/avg_occupancy.tiff")
