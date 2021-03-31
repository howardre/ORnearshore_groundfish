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
# Load data and necessary functions ----
setwd("/Users/howar/Documents/Oregon State/ORnearshore_groundfish/code")
load('../data/NMFS_data/annual_samples')
load('../data/NMFS_data/annual_tows')
load('../data/NMFS_data/triennial_samples')
load('../data/NMFS_data/triennial_tows')
load("../data/bathy.dat")
load("../data/bathy.mat")
source("functions/subset_species.R")
source("functions/vis_gam_COLORS.R")
jet.colors <- colorRampPalette(rev(c("#b2182b", "#d6604d", "#f4a582", "#fddbc7",
                                     "#f7f7f7", "#d1e5f0", "#92c5de", "#4393c3", "#2166ac" )))
contour_col <- rgb(0, 0, 255, max = 255, alpha = 0, names = "white")

# Subset the data to contain only species of interest for each survey ----
# Eight species of interest
# Annual subsets
arrowtooth_annual <- subset_species_count("Atheresthes stomias", annual, annual_trawl)
english_annual <- subset_species_count("Parophrys vetulus", annual, annual_trawl)
sanddab_annual <- subset_species_count("Citharichthys sordidus", annual, annual_trawl)
dover_annual <- subset_species_count("Microstomus pacificus", annual, annual_trawl)
rex_annual <- subset_species_count("Glyptocephalus zachirus", annual, annual_trawl)
lingcod_annual <- subset_species_count("Ophiodon elongatus", annual, annual_trawl)
petrale_annual <- subset_species_count("Eopsetta jordani", annual, annual_trawl)
sablefish_annual <- subset_species_count("Anoplopoma fimbria", annual, annual_trawl)

# Triennial subsets
arrowtooth_triennial <- subset_species_count("Atheresthes stomias", triennial, triennial_trawl)
english_triennial <- subset_species_count("Parophrys vetulus", triennial, triennial_trawl)
sanddab_triennial <- subset_species_count("Citharichthys sordidus", triennial, triennial_trawl)
dover_triennial <- subset_species_count("Microstomus pacificus", triennial, triennial_trawl)
rex_triennial <- subset_species_count("Glyptocephalus zachirus", triennial, triennial_trawl)
lingcod_triennial <- subset_species_count("Ophiodon elongatus", triennial, triennial_trawl)
petrale_triennial <- subset_species_count("Eopsetta jordani", triennial, triennial_trawl)
sablefish_triennial <- subset_species_count("Anoplopoma fimbria", triennial, triennial_trawl)

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
    count ~ s(year) +
      s(julian) +
      s(latitude, longitude) +
      s(depth_m),
    ~ s(year) +
      s(julian) +
      s(latitude, longitude)),
  data = dover_subset,
  family = ziplss)
summary(dover_ziplss)
AIC(dover_ziplss)
plot(dover_ziplss)

# Include temperature
dover_ziplss_temp <- gam(
  list(
    count ~ s(year) +
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

test <- gam(count ~ s(year) +
  s(julian) +
  s(latitude, longitude) +
  s(depth_m), family = poisson, data = dover_subset[dover_subset$count > 0,])
summary(test)
plot(test)
gam.check(test)

range(dover_subset$count)
range(predict(dover_ziplss))
