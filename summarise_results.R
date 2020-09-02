# Locate RDS files

source("sourcer.R")


if(!file.exists("./results_summary")){
  dir.create("./results_summary")
}

model_selection <- TRUE
if(model_selection){
models <- c('global', 'housing_density', 'habitat', 'null')

# get cpo score for each species
cpo_ans <- matrix(NA, ncol = 5, nrow = 8)
colnames(cpo_ans) <- c("species", models)
for(i in 1:length(models)){
  my_cpo <- read.csv(paste0("./results/",models[i],"/cpo.csv"))
  cpo_ans[,1] <- as.character(my_cpo$species)
  cpo_ans[,i+1] <- my_cpo$cpo
}

cpo_ans <- data.frame(cpo_ans)
cpo_ans$best <- apply(cpo_ans[,-1], 1, function(x) which.min(x))
cpo_ans$best <- models[cpo_ans$best]

write.csv(cpo_ans, "cpo_scores.csv", row.names = FALSE, quote = FALSE)


# build up the rds filenames
folder <- "./results/best/" 
my_rds <- paste0("./results/", cpo_ans$best,"/",cpo_ans$species,".RDS")
my_rds_short <- as.character(cpo_ans$species)

for(species in 1:length(my_rds)){
  tmp_rds <- readRDS(my_rds[species])
  tmp_mat <- as.matrix(as.mcmc.list(tmp_rds), chains = TRUE)
  # drop the z stuff
  tmp_matz <- tmp_mat[,grep("z", colnames(tmp_mat))]
  tmp_mat <- tmp_mat[,-c(1,grep("z|lik", colnames(tmp_mat)))]
  tmp_mat <- as.data.frame(tmp_mat)
  tmp_matz <- as.data.frame(tmp_matz)
  # save the file
  data.table::fwrite(tmp_mat, paste0(folder, my_rds_short[species],
                                     "_matrix.csv"))
  data.table::fwrite(tmp_matz, paste0(folder, my_rds_short[species],
                                      "_zed.csv"))
  rm(tmp_mat)
  rm(tmp_matz)
}

my_res <- list.files(folder , pattern = "_matrix.csv$", full.names = TRUE)
my_rds_short <- sapply(strsplit(my_res, "/"), "[", 4)
my_rds_short <- sapply(strsplit(my_rds_short, "\\."), "[", 1)
my_rds_short <- gsub("(\\1)_matrix", "\\1", my_rds_short)

# bring in the city data
cdat <- read.csv("data/city_level_data.csv", stringsAsFactors = FALSE)
cdat <- cdat[order(cdat$city),]




swt <- function(x, cdat){
  ans <- switch(x,
                'global' = c("Bmu", 'between_prophab','between_hden',
                             'within_urb', 'probhab_onurb', 'hden_onurb', 
                             'average_detection'),
                'housing_density' = c("Bmu",'between_hden',
                                      'within_urb', 'hden_onurb', 
                                      'average_detection'),
                'habitat' = c("Bmu", 'between_prophab',
                              'within_urb', 'probhab_onurb',
                              'average_detection'),
                'null' =c("Bmu",
                          'within_urb',
                          'average_detection')
  )
  
  to_ret <- c(paste0(cdat$city,"-intercept"), 
              paste0(cdat$city,"-slope"),
              ans,
              paste0(cdat$city,"-diff")
  )
  return(to_ret)
  
}


# make a plot for each species
for(species in 1:length(my_res)){
  # read in the file
  my_mcmc <- data.table::fread(my_res[species], data.table = FALSE) %>% 
    apply(., 2, HDIofMCMC, .95) %>% t(.)
  # round them 
  rmcmc <- round(my_mcmc, 2)
  # paste them together
  pmcmc <- apply(rmcmc, 1, function(x) {
    paste0(sprintf("%.2f", x[2])," (", 
           sprintf("%.2f", x[1]), " - ", 
           sprintf("%.2f", x[3]), ")")})
  
  
  
  
  # standard deviation now
  
  sd_res <- matrix(pmcmc[grep("sd", names(pmcmc))], 
                   ncol = 2, 
                   nrow = 1)
  sd_res <- cbind(my_rds_short[species],sd_res)
  colnames(sd_res) <- c("Species","psi_sd", "rho_sd")
  
  # city specific now
  npar <- length(grep("^B|^b1|^b2|^Dmu", names(pmcmc)))
  b_res <- matrix(pmcmc[grep("^B|^b1|^b2|^Dmu", names(pmcmc))], ncol = 1, nrow = npar)
  species_is_there <- as.logical(c(rep(det_events[,my_rds_short[species]], 2), 
                                   rep(NA, length(b_res) - 20)))
  
  
 
  
  
  b_res <- cbind(as.character(my_rds_short[species]),
                 swt(cpo_ans$best[order(cpo_ans$species)][species], cdat),
                 species_is_there,
                 b_res
                 )
  colnames(b_res) <- c("Species", "parameter",'species_there', "estimate")
  
  
  if(species == 1){
    write.csv(sd_res, "./results_summary/city_sd_cpo.csv", row.names = FALSE)
    write.csv(b_res, "./results_summary/within_city_cpo.csv", row.names = FALSE)
  } else{
    write.table(sd_res,"./results_summary/city_sd_cpo.csv", row.names = FALSE,
                append = TRUE, sep = ",", col.names = FALSE)
    write.table(b_res, "./results_summary/within_city_cpo.csv", row.names = FALSE, 
                col.names = FALSE,
                append = TRUE, sep = ",")
  }
  
  
}




} else{

model <- "global"
folder <- paste0("./results/", model,"/")
my_rds <- list.files(folder, pattern = "*.RDS", full.names = TRUE)
#my_rds <- my_rds[-grep("waic", my_rds)]
my_rds_short <- sapply(strsplit(my_rds, "/"), "[", 4)
my_rds_short <- sapply(strsplit(my_rds_short, "\\."), "[", 1)
#my_rds_short <- gsub("_ranefshift", "", my_rds_short)
}
for(species in 1:length(my_rds)){
  tmp_rds <- readRDS(my_rds[species])
  tmp_mat <- as.matrix(as.mcmc.list(tmp_rds), chains = TRUE)
  # drop the z stuff
  tmp_matz <- tmp_mat[,grep("z", colnames(tmp_mat))]
  tmp_mat <- tmp_mat[,-c(1,grep("z|lik", colnames(tmp_mat)))]
  tmp_mat <- as.data.frame(tmp_mat)
  tmp_matz <- as.data.frame(tmp_matz)
  # save the file
  data.table::fwrite(tmp_mat, paste0(folder, my_rds_short[species],
                            "_matrix.csv"))
  data.table::fwrite(tmp_matz, paste0(folder, my_rds_short[species],
                                     "_zed.csv"))
  rm(tmp_mat)
  rm(tmp_matz)
}

my_res <- list.files(folder , pattern = "_matrix.csv$", full.names = TRUE)

# bring in the city data
cdat <- read.csv("data/city_level_data.csv", stringsAsFactors = FALSE)
cdat <- cdat[order(cdat$city),]

# make a plot for each species
for(species in 1:length(my_res)){
  # read in the file
  my_mcmc <- data.table::fread(my_res[species], data.table = FALSE) %>% 
    apply(., 2, HDIofMCMC, .95) %>% t(.)
  # round them 
  rmcmc <- round(my_mcmc, 2)
  # paste them together
  pmcmc <- apply(rmcmc, 1, function(x) {
    paste0(sprintf("%.2f", x[2])," (", 
           sprintf("%.2f", x[1]), " - ", 
           sprintf("%.2f", x[3]), ")")})
  

  
  
  # standard deviation now
  
  sd_res <- matrix(pmcmc[grep("sd", names(pmcmc))], 
                   ncol = 2, 
                   nrow = 1)
  sd_res <- cbind(my_rds_short[species],sd_res)
  colnames(sd_res) <- c("Species","psi_sd", "rho_sd")
  
  # city specific now
  
  b_res <- matrix(pmcmc[grep("^B|^b1|^b2|^Dmu", names(pmcmc))], ncol = 1, nrow = 37)
  species_is_there <- as.logical(c(rep(det_events[,my_rds_short[species]], 2), rep(NA, length(b_res) - 20)))
  
  b_res <- cbind(my_rds_short[species], 
                 c(paste0(cdat$city,"-intercept"), paste0(cdat$city,"-slope"), "Bmu", 'between_prophab','between_hden',
                   'within_urb', 'probhab_onurb', 'hden_onurb', 'average_detection', paste0(cdat$city,"-diff")),
                 species_is_there,
                 b_res)
  colnames(b_res) <- c("Species", "parameter",'species_there', "estimate")
  
  
  if(species == 1){
    write.csv(sd_res, "./results_summary/city_sd4.csv", row.names = FALSE)
    write.csv(b_res, "./results_summary/within_city4.csv", row.names = FALSE)
  } else{
    write.table(sd_res,"./results_summary/city_sd4.csv", row.names = FALSE,
                append = TRUE, sep = ",", col.names = FALSE)
    write.table(b_res, "./results_summary/within_city4.csv", row.names = FALSE, 
                col.names = FALSE,
                append = TRUE, sep = ",")
  }
  
  
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

# This is the average occupancy across the cities that have a given 
#  species.

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

