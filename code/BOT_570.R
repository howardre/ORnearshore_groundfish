library(dplyr)
library(fossil)

setwd('/Users/howar/Documents/Oregon State/Thesis/Data visualization')
load("trawl_data")
setwd("C:/Users/howar/Documents/Oregon State/Thesis/Data analysis")
load("triennial_filtered")
load("annual_filtered")


# This creates two matrices, with annual requiring a subset

#--------------------########Annual########-------------------------------------

# Create stratified random sample of annual trawl data to get stations
annual_trawl <- filter(trawl_data,project !="Groundfish Triennial Shelf Survey")
annual_trawl <- annual_trawl[!is.na(annual_trawl$bottom_temp),]
annual_trawl <- annual_trawl[!is.na(annual_trawl$PDO),]
annual_trawl <- annual_trawl[!is.na(annual_trawl$NPGO),]
annual_trawl <- annual_trawl[!is.na(annual_trawl$julian),]
annual_trawl <- annual_trawl[complete.cases(annual_trawl),]
trawl <- annual_trawl
#set.seed(1)
#trawl <- annual_trawl%>% group_by(year) %>% sample_n(65)
trawl_filter <- c(unique(trawl$trawl_id))

#combine datasets
annual_filter <- filter(annual_filter, trawl_id %in% trawl_filter)
match_id<-match(annual_filter$trawl_id,trawl$trawl_id)
annual_filter$PDO<-trawl$PDO[match_id]
annual_filter$NPGO<-trawl$NPGO[match_id]
annual_filter$julian<-trawl$julian[match_id]
annual_filter$bottom_temp<-trawl$bottom_temp[match_id]


#remove unselected samples
subset <- annual_filter[c(1,5,6,12,14,16,19,20,21,22)]
subset <- subset[complete.cases(subset),]

#make matrices
#species
subset <- as.data.frame(subset[order(subset$trawl_id),])
annual_matrix <- create.matrix(subset, tax.name = "common_name",
                               locality = "trawl_id", time.col = NULL, time = NULL,
                               abund = T, abund.col = "lncpue_n")
annual_matrix <- as.data.frame(annual_matrix)
annual_mat <- t(annual_matrix)
species_matrix <- as.data.frame(annual_mat)

#environmental
env_filter <- c(unique(subset$trawl_id))
trawl <- filter(trawl, trawl_id %in% env_filter)
annual_env <- select(trawl, trawl_id, bottom_temp, PDO, NPGO, depth_m, latitude, year, julian)
env_matrix <- annual_env[complete.cases(annual_env),]
env_matrix <- as.data.frame(env_matrix)
rownames(env_matrix) <- env_matrix[,1]
env_matrix <- env_matrix[,-1]

#save as csv
save(species_matrix, file = "species_matrix_a")
write.csv(species_matrix, 'species_matrix_a.csv')
save(env_matrix, file = "env_matrix_a")
write.csv(env_matrix, 'env_matrix_a.csv')

#--------------------#############Triennial############----------------------------

# Create stratified random sample of triennial trawl data to get stations
triennial_trawl <- filter(trawl_data,project !="Groundfish Slope and Shelf Combination Survey")
triennial_trawl <- triennial_trawl[!is.na(triennial_trawl$bottom_temp),]
triennial_trawl <- triennial_trawl[!is.na(triennial_trawl$PDO),]
triennial_trawl <- triennial_trawl[!is.na(triennial_trawl$NPGO),]
triennial_trawl <- triennial_trawl[!is.na(triennial_trawl$julian),]

#combine datasets
match_id<-match(triennial_filter$trawl_id,triennial_trawl$trawl_id)
triennial_filter$PDO<-triennial_trawl$PDO[match_id]
triennial_filter$NPGO<-triennial_trawl$NPGO[match_id]
triennial_filter$julian<-triennial_trawl$julian[match_id]
triennial_filter$bottom_temp<-triennial_trawl$bottom_temp[match_id]

#remove unselected samples
subset <- triennial_filter[c(1,5,6,12,14,16,19,20,21,22)]
subset <- subset[complete.cases(subset),]

#make matrices
#species
subset <- as.data.frame(subset[order(subset$trawl_id),])
triennial_matrix <- create.matrix(subset, tax.name = "common_name",
                               locality = "trawl_id", time.col = NULL, time = NULL,
                               abund = T, abund.col = "lncpue_n")
triennial_matrix <- as.data.frame(triennial_matrix)
triennial_mat <- t(triennial_matrix)
species_matrix <- as.data.frame(triennial_mat)

#environmental
env_filter <- c(unique(subset$trawl_id))
triennial_trawl <- filter(triennial_trawl, trawl_id %in% env_filter)
triennial_env <- select(triennial_trawl, trawl_id, bottom_temp, PDO, NPGO, depth_m, latitude, year, julian)
env_matrix <- na.omit(triennial_env)
env_matrix <- as.data.frame(env_matrix)
rownames(env_matrix) <- env_matrix[,1]
env_matrix <- env_matrix[,-1]

#save as csv
save(species_matrix, file = "species_matrix_t")
write.csv(species_matrix, 'species_matrix_t.csv')
save(env_matrix, file = "env_matrix_t")
write.csv(env_matrix, 'env_matrix_t.csv')
