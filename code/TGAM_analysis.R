# Title: Threshold GAM Analysis
# Purpose: Create threshold GAMs for 8 groundfish species to determine regime shifts
# Date Created: 10/18/2020

# Load libraries ----
library(colorRamps)
library(itsadug)
library(fields)
library(sgeostat)
library(dplyr)

# Load data and necessary functions ----
setwd("/Users/howar/Documents/Oregon State/ORnearshore_groundfish/code")
trawl_data <- read.delim("../data/NMFS_data/trawl_data.txt", header = T)
OR_fish <- read.delim("../data/NMFS_data/OR_fish.txt", header = T)
load("../data/bathy.dat")
load("../data/bathy.mat")
source("functions/distance_function.R")
source("functions/vis_gam_COLORS.R")
source("functions/subset_species.R")
source("functions/TGAM_selection.R")
jet_colors<-colorRampPalette(c("#F7FCFD", "#E0ECF4", "#BFD3E6", "#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D", "#6E016B"))

# Subset the data to contain only species of interest ----
# Eight species of interest
arrowtooth_subset <- subset_species("Atheresthes stomias", OR_fish, trawl_data)
english_subset <- subset_species("Parophrys vetulus", OR_fish, trawl_data)
sanddab_subset <- subset_species("Citharichthys sordidus", OR_fish, trawl_data)
dover_subset <- subset_species("Microstomus pacificus", OR_fish, trawl_data)
rex_subset <- subset_species("Glyptocephalus zachirus", OR_fish, trawl_data)
lingcod_subset <- subset_species("Ophiodon elongatus", OR_fish, trawl_data)
petrale_subset <- subset_species("Eopsetta jordani", OR_fish, trawl_data)
sablefish_subset <- subset_species("Anoplopoma fimbria", OR_fish, trawl_data)

# Find thresholds for all 8 species ----
years <- sort(unique(arrowtooth_subset$year))[4:22]
arrowtooth_tgam <- get_tgam(arrowtooth_subset, years)
english_tgam <- get_tgam(english_subset, years)
sanddab_tgam <- get_tgam(sanddab_subset, years)
dover_tgam <- get_tgam_woyear(dover_subset, years)
rex_tgam <- get_tgam_woyear(rex_subset, years)
lingcod_tgam <- get_tgam(lingcod_subset, years)
petrale_tgam <- get_tgam(petrale_subset, years)
sablefish_tgam <- get_tgam(sablefish_subset, years)

# View summary statistics for reference [[1]] and final [[2]] TGAMs ----
summary(arrowtooth_tgam[[1]])
summary(arrowtooth_tgam[[2]])
summary(english_tgam[[1]])
summary(english_tgam[[2]])
summary(sanddab_tgam[[1]])
summary(sanddab_tgam[[2]])
summary(dover_tgam[[1]])
summary(dover_tgam[[2]])
summary(rex_tgam[[1]])
summary(rex_tgam[[2]])
summary(lingcod_tgam[[1]])
summary(lingcod_tgam[[2]])
summary(petrale_tgam[[1]])
summary(petrale_tgam[[2]])
summary(sablefish_tgam[[1]])
summary(sablefish_tgam[[2]])

# Save the results ----
save(arrowtooth_tgam, file = "../results/TGAM/arrowtooth_flounder/arrowtooth_tgam.RData")
save(english_tgam, file = "../results/TGAM/english_sole/english_tgam.RData")
save(sanddab_tgam, file = "../results/TGAM/pacific_sanddab/sanddab_tgam.RData")
save(dover_tgam, file = "../results/TGAM/dover_sole/dover_tgam.RData")
save(rex_tgam, file = "../results/TGAM/rex_sole/rex_tgam.RData")
save(lingcod_tgam, file = "../results/TGAM/lingcod/lingcod_tgam.RData")
save(petrale_tgam, file = "../results/TGAM/petrale_sole/petrale_tgam.RData")
save(sablefish_tgam, file = "../results/TGAM/sablefish/sablefish_tgam.RData")

# Plot AIC for all species ----
windows()
par(family = "serif")
plot_AIC(arrowtooth_tgam, years)
plot_AIC(english_tgam, years)
plot_AIC(sanddab_tgam, years)
plot_AIC(dover_tgam, years)
plot_AIC(rex_tgam, years)
plot_AIC(lingcod_tgam, years)
plot_AIC(petrale_tgam, years)
plot_AIC(sablefish_tgam, years)

