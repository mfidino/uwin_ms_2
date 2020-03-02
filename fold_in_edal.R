# massage in the EDAL data

edal <- read.csv("C:/Users/mfidino/Downloads/OccupancyReport (3).csv",
                 stringsAsFactors = FALSE, skip = 3)


yo <- edal %>% group_by(Species) %>% 
  summarise_if(is.integer, function(x) sum(x, na.rm  = TRUE))
# calculate the number o\

kk <- data.frame(yo$Species, rowSums(yo[,-c(1:2)]))


blank <- data.frame(matrix(0, ncol = ncol(det_data), nrow = length(unique(edal$Site))))

colnames(blank) <- colnames(det_data)
blank$X <- 766:c(766+nrow(blank) -1)
blank$city <- 'edal'
blank$site <- unique(edal$Site)
blank$site_code <- paste0("edal-", unique(edal$Site), "-2")
blank$year <- 2018

blank$coyote <- rowSums(edal[which(edal$Species == "coyote"),grep("Day", colnames(edal))],
                        na.rm = TRUE)
blank$skunk <- rowSums(edal[which(edal$Species == "st. skunk"),grep("Day", colnames(edal))],
                       na.rm = TRUE)
blank$J <- length(grep("Day", colnames(edal))) - 
  rowSums(is.na(edal[which(edal$Species == "coyote"),grep("Day", colnames(edal))]))

# remove J = 0
blank <- blank[-which(blank$J == 0),]

ed_sites <- read.csv("C:/Users/mfidino/Downloads/edal_coords.csv", stringsAsFactors = FALSE)


ed_sites <- ed_sites[which(ed_sites$locationAbbr %in% blank$site),]


det_data %>% group_by(city) %>% 
  summarise_if(is.integer, function(x) sum(x>0) )
