# Title: Stationary GAM Analysis
# Purpose: Create threshold GAMs for 8 groundfish species to determine regime shifts
# Date Created: 11/06/2020

# Load libraries ----
library(maps)
library(mapdata)
library(fields)
library(plotfunctions)
library(mgcv)
library(dplyr)

# Load data and necessary functions ----
setwd("/Users/howar/Documents/Oregon State/ORnearshore_groundfish/code")
load('../data/NMFS_data/annual_samples')
load('../data/NMFS_data/annual_tows')
load('../data/NMFS_data/triennial_samples')
load('../data/NMFS_data/triennial_tows')
load("../data/bathy.dat")
load("../data/bathy.mat")
source("functions/subset_species.R")

# Subset the data to contain only species of interest for each survey ----
# Eight species of interest
# Annual subsets
arrowtooth_annual <- subset_species("Atheresthes stomias", annual, annual_trawl)
english_annual <- subset_species("Parophrys vetulus", annual, annual_trawl)
sanddab_annual <- subset_species("Citharichthys sordidus", annual, annual_trawl)
dover_annual <- subset_species("Microstomus pacificus", annual, annual_trawl)
rex_annual <- subset_species("Glyptocephalus zachirus", annual, annual_trawl)
lingcod_annual <- subset_species("Ophiodon elongatus", annual, annual_trawl)
petrale_annual <- subset_species("Eopsetta jordani", annual, annual_trawl)
sablefish_annual <- subset_species("Anoplopoma fimbria", annual, annual_trawl)

# Triennial subsets
arrowtooth_triennial <- subset_species("Atheresthes stomias", triennial, triennial_trawl)
english_triennial <- subset_species("Parophrys vetulus", triennial, triennial_trawl)
sanddab_triennial <- subset_species("Citharichthys sordidus", triennial, triennial_trawl)
dover_triennial <- subset_species("Microstomus pacificus", triennial, triennial_trawl)
rex_triennial <- subset_species("Glyptocephalus zachirus", triennial, triennial_trawl)
lingcod_triennial <- subset_species("Ophiodon elongatus", triennial, triennial_trawl)
petrale_triennial <- subset_species("Eopsetta jordani", triennial, triennial_trawl)
sablefish_triennial <- subset_species("Anoplopoma fimbria", triennial, triennial_trawl)

# Investigate the data ----
# Properties of CPUE data
CPUE_hist <- function(species_subset){
  windows()
  par(mfrow = c(1,2))
  hist(species_subset$lncpue,
       xlab = "lncpue",
       main = "Raw Data")
  sum(1 * (species_subset$lncpue == 0)) / nrow(species_subset)
  hist(species_subset$lncpue[species_subset$lncpue > 0],
    xlab = "lncpue",
    main = "Normalized Data")
}
CPUE_hist(arrowtooth_annual)
CPUE_hist(english_annual)
CPUE_hist(lingcod_annual)
CPUE_hist(dover_annual)
CPUE_hist(rex_annual)
CPUE_hist(petrale_annual)
CPUE_hist(sablefish_annual)
CPUE_hist(sanddab_annual)
CPUE_hist(arrowtooth_triennial)
CPUE_hist(english_triennial)
CPUE_hist(lingcod_triennial)
CPUE_hist(dover_triennial)
CPUE_hist(rex_triennial)
CPUE_hist(petrale_triennial)
CPUE_hist(sablefish_triennial)
CPUE_hist(sanddab_triennial)

# Variable selection: look at p value, plots, and ultimately AIC and GCV once you remove variable
# GCV: expresses error that model does once you use it for prediction
# AIC: is likelihood of model divided by number of parameters
# GAM on presence/absence ----
gam_check <- function(gam){
  windows()
  par(mfrow = c(2, 2))
  gam.check(gam[[2]])
} # Use to check residuals
gam_plot <- function(gam){
  windows()
  plot(gam[[2]], pages = 1, scale = 0)
  } # View the plots
pa_GAMs <- function(species_subset){
  julian_gam <-  gam(pres ~ factor(year) + s(longitude, latitude) + s(julian),
                     family = binomial,
                     data = species_subset)
  wojulian_gam <-  gam(pres ~ factor(year) + s(longitude, latitude),
                       family = binomial,
                       data = species_subset)
  gam_list <- list(julian_gam, wojulian_gam)
  best_gam <- gam_list[[which.min(sapply(1:length(gam_list),
                                         function(x) min(gam_list[[x]]$aic)))]]
  return_list <- list(gam_list, best_gam)
} # Create the p/a GAMs (2)
paGAM_map <- function(gam, species_subset){
  windows(width = 4, height = 8)
  vis.gam(gam[[2]],
          view = c('longitude', 'latitude'),
          too.far = 0.03, # can change to get better extrapolation
          plot.type = 'contour',
          color = 'topo',
          type = 'response',
          main = paste("Presence - Absence GAM"))
  symbols(species_subset$longitude[species_subset$lncpue > 0],
          species_subset$latitude[species_subset$lncpue > 0],
          circles = species_subset$lncpue[species_subset$lncpue > 0],
          inches = 0.05,
          add = T,
          bg = alpha('purple', 0.2))
  maps::map('worldHires',
            add = T,
            col = 'antiquewhite4',
            fill = T)
} # Use to map the model predictions

# ***Arrowtooth Flounder ----
# Annual
arrowtooth_annualpa <- pa_GAMs(arrowtooth_annual)
summary(arrowtooth_annualpa[[2]]) # view the best model
gam_check(arrowtooth_annualpa)
gam_plot(arrowtooth_annualpa)

# Triennial
arrowtooth_triennialpa <- pa_GAMs(arrowtooth_triennial)
summary(arrowtooth_triennialpa[[2]]) # view the best model
gam_check(arrowtooth_triennialpa)
windows()
plot(arrowtooth_triennialpa[[2]], pages = 1, scale = 0) # View the plots

# Make plot with best models
paGAM_map(arrowtooth_annualpa, arrowtooth_annual)
paGAM_map(arrowtooth_triennialpa, arrowtooth_triennial)

# ***Dover Sole ----
# Annual
dover_annualpa <- pa_GAMs(dover_annual)
summary(dover_annualpa[[2]]) # view the best model
gam_check(dover_annualpa)
gam_plot(dover_annualpa)

# Triennial
dover_triennialpa <- pa_GAMs(dover_triennial)
summary(dover_triennialpa[[2]]) # view the best model
gam_check(dover_triennialpa)
windows()
plot(dover_triennialpa[[2]], pages = 1, scale = 0) # View the plots

# Make plot with best models
paGAM_map(dover_annualpa, dover_annual)
paGAM_map(dover_triennialpa, dover_triennial)
# ***English Sole ----
# Annual
english_annualpa <- pa_GAMs(english_annual)
summary(english_annualpa[[2]]) # view the best model
gam_check(english_annualpa)
gam_plot(english_annualpa)

# Triennial
english_triennialpa <- pa_GAMs(english_triennial)
summary(english_triennialpa[[2]]) # view the best model
gam_check(english_triennialpa)
windows()
plot(english_triennialpa[[2]], pages = 1, scale = 0) # View the plots

# Make plot with best models
paGAM_map(english_annualpa, english_annual)
paGAM_map(english_triennialpa, english_triennial)

# ***Lingcod ----
# Annual
lingcod_annualpa <- pa_GAMs(lingcod_annual)
summary(lingcod_annualpa[[2]]) # view the best model
gam_check(lingcod_annualpa)
gam_plot(lingcod_annualpa)

# Triennial
lingcod_triennialpa <- pa_GAMs(lingcod_triennial)
summary(lingcod_triennialpa[[2]]) # view the best model
gam_check(lingcod_triennialpa)
windows()
plot(lingcod_triennialpa[[2]], pages = 1, scale = 0) # View the plots

# Make plot with best models
paGAM_map(lingcod_annualpa, lingcod_annual)
paGAM_map(lingcod_triennialpa, lingcod_triennial)


# ***Pacific Sanddab ----
# Annual
sanddab_annualpa <- pa_GAMs(sanddab_annual)
summary(sanddab_annualpa[[2]]) # view the best model
gam_check(sanddab_annualpa)
gam_plot(sanddab_annualpa)

# Triennial
sanddab_triennialpa <- pa_GAMs(sanddab_triennial)
summary(sanddab_triennialpa[[2]]) # view the best model
gam_check(sanddab_triennialpa)
windows()
plot(sanddab_triennialpa[[2]], pages = 1, scale = 0) # View the plots

# Make plot with best models
paGAM_map(sanddab_annualpa, sanddab_annual)
paGAM_map(sanddab_triennialpa, sanddab_triennial)

# ***Petrale Sole ----
# Annual
petrale_annualpa <- pa_GAMs(petrale_annual)
summary(petrale_annualpa[[2]]) # view the best model
gam_check(petrale_annualpa)
gam_plot(petrale_annualpa)

# Triennial
petrale_triennialpa <- pa_GAMs(petrale_triennial)
summary(petrale_triennialpa[[2]]) # view the best model
gam_check(petrale_triennialpa)
windows()
plot(petrale_triennialpa[[2]], pages = 1, scale = 0) # View the plots

# Make plot with best models
paGAM_map(petrale_annualpa, petrale_annual)
paGAM_map(petrale_triennialpa, petrale_triennial)
# ***Rex Sole ----
# Annual
rex_annualpa <- pa_GAMs(rex_annual)
summary(rex_annualpa[[2]]) # view the best model
gam_check(rex_annualpa)
gam_plot(rex_annualpa)

# Triennial
rex_triennialpa <- pa_GAMs(rex_triennial)
summary(rex_triennialpa[[2]]) # view the best model
gam_check(rex_triennialpa)
windows()
plot(rex_triennialpa[[2]], pages = 1, scale = 0) # View the plots

# Make plot with best models
paGAM_map(rex_annualpa, rex_annual)
paGAM_map(rex_triennialpa, rex_triennial)
# ***Sablefish ----
# Annual
sablefish_annualpa <- pa_GAMs(sablefish_annual)
summary(sablefish_annualpa[[2]]) # view the best model
gam_check(sablefish_annualpa)
gam_plot(sablefish_annualpa)

# Triennial
sablefish_triennialpa <- pa_GAMs(sablefish_triennial)
summary(sablefish_triennialpa[[2]]) # view the best model
gam_check(sablefish_triennialpa)
windows()
plot(sablefish_triennialpa[[2]], pages = 1, scale = 0) # View the plots

# Make plot with best models
paGAM_map(sablefish_annualpa, sablefish_annual)
paGAM_map(sablefish_triennialpa, sablefish_triennial)

# GAM on CPUE ----
cpue_GAMs <- function(species_subset){
  year_gam <- gam(lncpue ~ s(year, k = 5) + s(longitude, latitude) + s(depth_m) + s(julian, k = 5),
                  data = species_subset[species_subset$lncpue > 0,])
  yeartemp_gam <-  gam(lncpue ~ s(year, k = 5) + s(longitude, latitude) + s(depth_m) +
                         s(julian, k = 5) + s(bottom_temp, k = 5),
                       data = species_subset[species_subset$lncpue > 0,])
  PDO_gam <- gam(lncpue ~ s(PDO, k = 5) + s(longitude, latitude) + s(depth_m) + s(julian, k = 5),
                  data = species_subset[species_subset$lncpue > 0,])
  PDOtemp_gam <-  gam(lncpue ~ s(PDO, k = 5) + s(longitude, latitude) + s(depth_m) +
                         s(julian, k = 5) + s(bottom_temp, k = 5),
                       data = species_subset[species_subset$lncpue > 0,])
  NPGO_gam <- gam(lncpue ~ s(NPGO, k = 5) + s(longitude, latitude) + s(depth_m) + s(julian, k = 5),
                  data = species_subset[species_subset$lncpue > 0,])
  NPGOtemp_gam <-  gam(lncpue ~ s(NPGO, k = 5) + s(longitude, latitude) + s(depth_m) +
                         s(julian, k = 5) + s(bottom_temp, k = 5),
                       data = species_subset[species_subset$lncpue > 0,])
  gam_list <- list(year_gam, yeartemp_gam, PDOtemp_gam, PDO_gam, NPGO_gam, NPGOtemp_gam)
  best_gam <- gam_list[[which.min(sapply(1:length(gam_list),
                                         function(x) min(gam_list[[x]]$aic)))]] # would like to also select by AIC
  return_list <- list(gam_list, best_gam)
} # need to figure out how to add GCV selection criteria
cpueGAM_arrow_a <- cpue_GAMs(arrowtooth_annual)
summary(cpueGAM_arrow_a[[2]]) # view the best model
gam_check(cpueGAM_arrow_a)

windows()
plot(cpueGAM_arrow_a[[2]], )
