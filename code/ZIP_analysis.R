# Title: ZIP Analysis
# Purpose: Create threshold GAMs for 8 groundfish species to determine regime shifts
# Date Created: 03/18/2021

# Load libraries ----
library(colorRamps)
library(itsadug)
library(fields)
library(sgeostat)
library(dplyr)
library(maps)
library(mapdata)
library(purrr)
library(mgcv)
library(furrr)

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
jet.colors <- colorRampPalette(c("#F7FCFD", "#E0ECF4", "#BFD3E6", "#9EBCDA",
                                 "#8C96C6", "#8C6BB1", "#88419D", "#6E016B"))
trawl_data$program <- recode(trawl_data$project, "Groundfish Slope and Shelf Combination Survey" = "annual",
                             "Groundfish Triennial Shelf Survey" = "triennial")
trawl_data$program[trawl_data$program == "triennial" & trawl_data$year < 1995] <- "triennial1"
trawl_data$program[trawl_data$program == "triennial" & trawl_data$year > 1994] <- "triennial2"

# Subset the data to contain only species of interest ----
# Eight species of interest
arrowtooth_subset <- subset_species_temp("Atheresthes stomias", OR_fish, trawl_data)
english_subset <- subset_species_temp("Parophrys vetulus", OR_fish, trawl_data)
sanddab_subset <- subset_species_temp("Citharichthys sordidus", OR_fish, trawl_data)
dover_subset <- subset_species_temp("Microstomus pacificus", OR_fish, trawl_data)
rex_subset <- subset_species_temp("Glyptocephalus zachirus", OR_fish, trawl_data)
lingcod_subset <- subset_species_temp("Ophiodon elongatus", OR_fish, trawl_data)
petrale_subset <- subset_species_temp("Eopsetta jordani", OR_fish, trawl_data)
sablefish_subset <- subset_species_temp("Anoplopoma fimbria", OR_fish, trawl_data)

# Test zip for Dover
dover_zip <- gam(count ~ factor(year) +
                   s(julian) +
                   s(latitude, longitude) +
                   s(bottom_temp),
                 data = dover_subset,
                 family = ziP())
summary(dover_zip)

# Test ziplss for Dover
# No environmental
dover_ziplss <- gam(
  list(
    count ~ factor(program) +
      s(year) +
      s(julian) +
      s(latitude, longitude),
    ~ factor(program) +
      s(year) +
      s(julian) +
      s(latitude, longitude) +
      s(depth_m)),
  data = dover_subset,
  family = ziplss)
summary(dover_ziplss)
AIC(dover_ziplss)
plot(dover_ziplss)

# Include temperature
dover_ziplss_temp <- gam(
  list(
    count ~ factor(program) +
      s(year) +
      s(julian) +
      s(latitude, longitude),
    ~ factor(program) +
      s(year) +
      s(julian) +
      s(latitude, longitude) +
      s(depth_m) +
      s(bottom_temp)),
  data = dover_subset,
  family = ziplss)
summary(dover_ziplss_temp)
AIC(dover_ziplss_temp)
plot(dover_ziplss_temp)
gam.check(dover_ziplss_temp)


# Remove program variable
dover_ziplss_program <- gam(
  list(
    count ~ s(year) +
      s(julian) +
      s(latitude, longitude),
    ~ s(year) +
      s(julian) +
      s(latitude, longitude) +
      s(depth_m)),
  data = dover_subset,
  family = ziplss)
summary(dover_ziplss_program)
AIC(dover_ziplss_program)
plot(dover_ziplss_program)
