# Title: NMS Matrices
# This creates two matrices, with annual requiring a subset
# Date created: 10/22/2020

# Load libraries ----
library(dplyr)
library(fossil)

# Load data
setwd('/Users/howar/Documents/Oregon State/ORnearshore_groundfish/code')
options(scipen=999)
trawl_data <- read.delim("../data/NMFS_data/trawl_data.txt", header = T)
triennial_filtered <- read.delim("../data/NMFS_data/triennial_filtered.txt")
annual_filtered <- read.delim("../data/NMFS_data/annual_filtered.txt")

# Annual Survey matrices ----
# Remove the rows with NA's in the environmental variable columns
annual_trawl <- filter(trawl_data, project != "Groundfish Triennial Shelf Survey")
annual_trawl <- annual_trawl[!is.na(annual_trawl$bottom_temp), ]
annual_trawl <- annual_trawl[!is.na(annual_trawl$PDO), ]
annual_trawl <- annual_trawl[!is.na(annual_trawl$NPGO), ]
annual_trawl <- annual_trawl[!is.na(annual_trawl$julian), ]
annual_trawl <- annual_trawl[complete.cases(annual_trawl), ]

# Can create a subsetted matrix if necessary
# set.seed(1)
# trawl <- annual_trawl %>% group_by(year) %>% sample_n(65)
trawl_filter <-  c(unique(annual_trawl$trawl_id)) # Create a way to filter out the trawls with missing data

# Combine datasets
annual_filtered <- filter(annual_filtered, trawl_id %in% trawl_filter) # Filter out the trawls missing environmental data
match_id <- match(annual_filtered$trawl_id, annual_trawl$trawl_id)
annual_filtered$PDO <- annual_trawl$PDO[match_id]
annual_filtered$NPGO <- annual_trawl$NPGO[match_id]
annual_filtered$julian <- annual_trawl$julian[match_id]
annual_filtered$bottom_temp <- annual_trawl$bottom_temp[match_id]

# Remove unselected samples
annual_subset <- annual_filtered[c(1, 5, 6, 12, 14, 16, 19, 20, 21, 22)]
annual_subset <- annual_subset[complete.cases(annual_subset),]

# Make matrices ----
# Species
annual_subset <- as.data.frame(annual_subset[order(annual_subset$trawl_id),])
annual_matrix <- create.matrix(annual_subset,
                               tax.name = "common_name",
                               locality = "trawl_id",
                               time.col = NULL,
                               time = NULL,
                               abund = T,
                               abund.col = "lncpue_n")
species_matrix_annual <- as.data.frame(t(annual_matrix))

# Environmental matrix
env_filter_annual <- c(unique(annual_subset$trawl_id))
annual_trawl <- filter(annual_trawl, trawl_id %in% env_filter_annual)
annual_env <- select(annual_trawl,
                     trawl_id,
                     bottom_temp,
                     PDO,
                     NPGO,
                     depth_m,
                     latitude,
                     year,
                     julian)
env_matrix_annual <- annual_env[complete.cases(annual_env),]
env_matrix_annual <- as.data.frame(env_matrix_annual)
rownames(env_matrix_annual) <- env_matrix_annual[, 1]
env_matrix_annual <- env_matrix_annual[, -1]

# Save as .csv and .R files
save(species_matrix_annual, file = "../data/NMFS_data/species_matrix_a")
write.csv(species_matrix_annual, '../data/NMFS_data/species_matrix_a.csv')
save(env_matrix_annual, file = "../data/NMFS_data/env_matrix_a")
write.csv(env_matrix_annual, '../data/NMFS_data/env_matrix_a.csv')

# Triennial Survey matrices ----
# Remove the rows with NA's in the environmental variable columns
triennial_trawl <- filter(trawl_data, project != "Groundfish Slope and Shelf Combination Survey")
triennial_trawl <- triennial_trawl[!is.na(triennial_trawl$bottom_temp), ]
triennial_trawl <- triennial_trawl[!is.na(triennial_trawl$PDO), ]
triennial_trawl <- triennial_trawl[!is.na(triennial_trawl$NPGO), ]
triennial_trawl <- triennial_trawl[!is.na(triennial_trawl$julian), ]
triennial_trawl <- triennial_trawl[complete.cases(triennial_trawl), ]

# Can create a subsetted matrix if necessary
# set.seed(1)
# trawl <- triennial_trawl %>% group_by(year) %>% sample_n(65)
trawl_filter <-  c(unique(triennial_trawl$trawl_id)) # Create a way to filter out the trawls with missing data

# Combine datasets
triennial_filtered <- filter(triennial_filtered, trawl_id %in% trawl_filter) # Filter out the trawls missing environmental data
match_id <- match(triennial_filtered$trawl_id, triennial_trawl$trawl_id)
triennial_filtered$PDO <- triennial_trawl$PDO[match_id]
triennial_filtered$NPGO <- triennial_trawl$NPGO[match_id]
triennial_filtered$julian <- triennial_trawl$julian[match_id]
triennial_filtered$bottom_temp <- triennial_trawl$bottom_temp[match_id]

# Remove unselected samples
triennial_subset <- triennial_filtered[c(1, 5, 6, 12, 14, 16, 19, 20, 21, 22)]
triennial_subset <- triennial_subset[complete.cases(triennial_subset),]

# Make matrices ----
# Species
triennial_subset <- as.data.frame(triennial_subset[order(triennial_subset$trawl_id),])
triennial_matrix <- create.matrix(triennial_subset,
                               tax.name = "common_name",
                               locality = "trawl_id",
                               time.col = NULL,
                               time = NULL,
                               abund = T,
                               abund.col = "lncpue_n")
species_matrix_triennial <- as.data.frame(t(triennial_matrix))

# Environmental matrix
env_filter_triennial <- c(unique(triennial_subset$trawl_id))
triennial_trawl <- filter(triennial_trawl, trawl_id %in% env_filter_triennial)
triennial_env <- select(triennial_trawl,
                     trawl_id,
                     bottom_temp,
                     PDO,
                     NPGO,
                     depth_m,
                     latitude,
                     year,
                     julian)
env_matrix_triennial <- triennial_env[complete.cases(triennial_env),]
env_matrix_triennial <- as.data.frame(env_matrix_triennial)
rownames(env_matrix_triennial) <- env_matrix_triennial[, 1]
env_matrix_triennial <- env_matrix_triennial[, -1]

# Save as .csv and .R files
save(species_matrix_triennial, file = "../data/NMFS_data/species_matrix_t")
write.csv(species_matrix_triennial, '../data/NMFS_data/species_matrix_t.csv')
save(env_matrix_triennial, file = "../data/NMFS_data/env_matrix_t")
write.csv(env_matrix_triennial, '../data/NMFS_data/env_matrix_t.csv')
