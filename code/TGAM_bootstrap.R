# Title: TGAM Bootstrapped Confidence Intervals
# Purpose: Validate the results of the TGAM analyses
# Date created: 10/21/2020

# Load necessary libraries ----
library(mgcv)
library(dplyr)
library(purrr)
library(sp)

# Load data ----
setwd('/Users/howar/Documents/Oregon State/Thesis/Data visualization')
load('trawl_data')
load('OR_fish')
setwd('/Users/howar/Documents/Oregon State/Thesis/Data analysis/Individual Species Analysis/')

# Subset the data to contain only species of interest ----
subset_species <- function(species, catch, tows){
  OR_subset <- catch[catch$scientific_name == species, ]
  match_id <- match(tows$trawl_id, OR_subset$trawl_id)
  tows$lncpue <- OR_subset$lncpue_n[match_id]
  tows$lncpue[is.na(tows$lncpue)] <- 0
  selected_species <- select(tows, julian, year, lncpue, latitude, longitude)
  selected_species <- na.omit(selected_species)
  selected_species$pres <- 1 * (selected_species$lncpue > 0)
  return(selected_species)
}

# Eight species of interest
arrowtooth_subset <- subset_species("Atheresthes stomias", OR_fish, trawl_data)
english_subset <- subset_species("Parophrys vetulus", OR_fish, trawl_data)
sanddab_subset <- subset_species("Citharichthys sordidus", OR_fish, trawl_data)
dover_subset <- subset_species("Microstomus pacificus", OR_fish, trawl_data)
rex_subset <- subset_species("Glyptocephalus zachirus", OR_fish, trawl_data)
lingcod_subset <- subset_species("Ophiodon elongatus", OR_fish, trawl_data)
petrale_subset <- subset_species("Eopsetta jordani", OR_fish, trawl_data)
sablefish_subset <- subset_species("Anoplopoma fimbria", OR_fish, trawl_data)

# Set up number of iterations ----
# Can use any subset if the same size
nsim <- 200 # Can change but running this many GAMs takes a long time
nsamp <- round(0.8 * nrow(arrowtooth_subset)) # Select 80% of the data
years <- sort(unique(arrowtooth_subset$year))[4:22] # Sort out the upper and lower quantiles
# years <- years[1:2] # for testing

# Go through each year to find the threshold year for sample ----
get_tgam_aic <- function(df, yr) {
  df$thr <- ifelse(df$year <= yr, 'before', 'after')
  aic_year <-  gam(
      pres ~ factor(year) + s(longitude, latitude, by = factor(thr)) + s(julian), 
      data = df,
      family = binomial)$aic
  return(aic_year)
}

# Returns AIC difference between TGAM and reference GAM for 1 simulation ----
get_aic_diff <- function(df) {
  samp <- sample_n(df, nsamp)
  ref_gam_aic <- gam(
      pres ~ factor(year) + s(julian) + s(longitude, latitude),
      data = samp,
      family = binomial)$aic
  tgam_all_aic <- rep(0, length(years))
  # tgam_all_aic <- map_dbl(1:length(years), 
  #                         ~ get_tgam_aic(samp, years[.]))
  tgam_all_aic <- map_dbl(years, ~get_tgam_aic(samp, .x))
  best_tgam_aic <- sort(tgam_all_aic)[1]
  diff <- ref_gam_aic - best_tgam_aic
  return(diff)
}

# Rex and Dover version, go through each year to find the threshold year for sample ----
get_tgam_aic2 <- function(df, yr) {
  df$thr <- ifelse(df$year <= yr, 'before', 'after')
  aic_year <-  gam(
    pres ~ factor(thr) + s(julian)+ s(longitude, latitude, by = factor(thr)) - 1, 
    data = df,
    family = binomial)$aic
  return(aic_year)
}

# Rex and Dover version, returns AIC difference between TGAM and reference GAM for 1 simulation ----
get_aic_diff2 <- function(df) {
  samp <- sample_n(df, nsamp)
  ref_gam_aic <- gam(
    pres ~ s(julian) + s(longitude, latitude),
    data = samp,
    family = binomial)$aic
  tgam_all_aic <- rep(0, length(years))
  tgam_all_aic <- map_dbl(years, ~get_tgam_aic2(samp, .x))
  best_tgam_aic <- sort(tgam_all_aic)[1]
  diff <- ref_gam_aic - best_tgam_aic
  return(diff)
}

# Get all the difference in AIC values ----
arrowtooth_differences <- map_dbl(1:nsim, ~ get_aic_diff(arrowtooth_subset))
english_differences <- map_dbl(1:nsim, ~ get_aic_diff(english_subset))
sanddab_differences <- map_dbl(1:nsim, ~ get_aic_diff(sanddab_subset))
lingcod_differences <- map_dbl(1:nsim, ~ get_aic_diff(lingcod_subset))
petrale_differences <- map_dbl(1:nsim, ~ get_aic_diff(petrale_subset))
sablefish_differences <- map_dbl(1:nsim, ~ get_aic_diff(sablefish_subset))
dover_differences <- map_dbl(1:nsim, ~ get_aic_diff2(dover_subset))
rex_differences <- map_dbl(1:nsim, ~ get_aic_diff2(rex_subset))
# Sys.time() # Check how long it takes to run 
plot(arrowtooth_differences) # Plot to see if there are any unusually large values





# Calculate the confidence interval ----
(mean_arrowtooth_differences <- mean(arrowtooth_differences))
CI_half_width_arrowtooth <- 1.96 * sd(arrowtooth_differences) / sqrt(nsim)
mean_arrowtooth_differences - CI_half_width_arrowtooth
mean_arrowtooth_differences + CI_half_width_arrowtooth

(mean_english_differences <- mean(english_differences))
CI_half_width_english <- 1.96 * sd(english_differences) / sqrt(nsim)
mean_english_differences - CI_half_width_english
mean_english_differences + CI_half_width_english

(mean_sanddab_differences <- mean(sanddab_differences))
CI_half_width_sanddab <- 1.96 * sd(sanddab_differences) / sqrt(nsim)
mean_sanddab_differences - CI_half_width_sanddab
mean_sanddab_differences + CI_half_width_sanddab

(mean_dover_differences <- mean(dover_differences))
CI_half_width_dover <- 1.96 * sd(dover_differences) / sqrt(nsim)
mean_dover_differences - CI_half_width_dover
mean_dover_differences + CI_half_width_dover

(mean_rex_differences <- mean(rex_differences))
CI_half_width_rex <- 1.96 * sd(rex_differences) / sqrt(nsim)
mean_rex_differences - CI_half_width_rex
mean_rex_differences + CI_half_width_rex

(mean_lingcod_differences <- mean(lingcod_differences))
CI_half_width_lingcod <- 1.96 * sd(lingcod_differences) / sqrt(nsim)
mean_lingcod_differences - CI_half_width_lingcod
mean_lingcod_differences + CI_half_width_lingcod

(mean_petrale_differences <- mean(petrale_differences))
CI_half_width_petrale <- 1.96 * sd(petrale_differences) / sqrt(nsim)
mean_petrale_differences - CI_half_width_petrale
mean_petrale_differences + CI_half_width_petrale

(mean_sablefish_differences <- mean(sablefish_differences))
CI_half_width_sablefish <- 1.96 * sd(sablefish_differences) / sqrt(nsim)
mean_sablefish_differences - CI_half_width_sablefish
mean_sablefish_differences + CI_half_width_sablefish

