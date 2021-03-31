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
source("functions/vis_ziplss_COLORS.R")
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


# Create functions to select the best ZIP for each species ----
ZIP_selection <- function(species_subset){
  year_ziplss <- gam(list(
    count ~ s(year, k = 5) +
      s(julian, k = 5) +
      s(latitude, longitude) +
      s(depth_m),
    ~ s(year, k = 5) +
      s(julian, k = 5) +
      s(latitude, longitude)
  ),
  data = species_subset,
  family = ziplss)
  yeartemp_ziplss <-  gam(list(
    count ~ s(year, k = 5) +
      s(julian, k = 5) +
      s(latitude, longitude) +
      s(depth_m) +
      s(bottom_temp, k = 5),
    ~ s(year, k = 5) +
      s(julian, k = 5) +
      s(latitude, longitude)
  ),
  data = species_subset,
  family = ziplss)
  PDO_ziplss <- gam(list(
    count ~ s(PDO, k = 5) +
      s(julian, k = 5) +
      s(latitude, longitude) +
      s(depth_m),
    ~ s(PDO, k = 5) +
      s(julian, k = 5) +
      s(latitude, longitude)
  ),
  data = species_subset,
  family = ziplss)
  PDOtemp_ziplss <-  gam(list(
    count ~ s(PDO, k = 5) +
      s(julian, k = 5) +
      s(latitude, longitude) +
      s(depth_m) +
      s(bottom_temp, k = 5),
    ~ s(PDO, k = 5) +
      s(julian, k = 5) +
      s(latitude, longitude)
  ),
  data = species_subset,
  family = ziplss)
  NPGO_ziplss <- gam(list(
    count ~ s(NPGO, k = 5) +
      s(julian, k = 5) +
      s(latitude, longitude) +
      s(depth_m),
    ~ s(NPGO, k = 5) +
      s(julian, k = 5) +
      s(latitude, longitude)
  ),
  data = species_subset,
  family = ziplss)
  NPGOtemp_ziplss <-  gam(list(
    count ~ s(NPGO, k = 5) +
      s(julian, k = 5) +
      s(latitude, longitude) +
      s(depth_m) +
      s(bottom_temp, k = 5),
    ~ s(NPGO, k = 5) +
      s(julian, k = 5) +
      s(latitude, longitude)
  ),
  data = species_subset,
  family = ziplss)
  ziplss_list <- list(year_ziplss, yeartemp_ziplss, PDOtemp_ziplss, PDO_ziplss, NPGO_ziplss, NPGOtemp_ziplss)
  best_ziplss <- ziplss_list[[which.min(sapply(1:length(ziplss_list),
                                         function(x) min(ziplss_list[[x]]$aic)))]] # would like to also select by AIC
  return_list <- list(ziplss_list, best_ziplss)
}
ZIP_test <- function(species_subset) {
  test <- gam(count ~ s(year) +
                          s(julian) +
                          s(latitude, longitude) +
                          s(depth_m), family = poisson, data = species_subset[species_subset$count > 0,])
}



# ****Dover sole ----
# Annual
dover_annualzip <- ZIP_selection(dover_annual)
summary(dover_annualzip[[2]])
plot(dover_annualzip[[2]])

windows()
par(mfrow = c(2,2))
gam.check(dover_annualzip[[2]])

# Test to ensure looks like a Poisson distribution
dover_annualcheck <- ZIP_test(dover_annual)
summary(dover_annualcheck)
plot(dover_annualcheck)

windows()
par(mfrow = c(2,2))
gam.check(dover_annualcheck)

range(dover_annual$count)
range(predict(dover_annualzip[[2]]))
range(predict(dover_annualcheck))


# Triennial
dover_triennialzip <- ZIP_selection(dover_triennial)
summary(dover_triennialzip[[2]])
plot(dover_triennialzip[[2]])

windows()
par(mfrow = c(2,2))
gam.check(dover_triennialzip[[2]])

# Test to ensure looks like a Poisson distribution
dover_triennialcheck <- ZIP_test(dover_triennial)
summary(dover_triennialcheck)
plot(dover_triennialcheck)

windows()
par(mfrow = c(2,2))
gam.check(dover_triennialcheck)

range(dover_triennial$count)
range(predict(dover_triennialzip[[2]]))
range(predict(dover_triennialcheck))
