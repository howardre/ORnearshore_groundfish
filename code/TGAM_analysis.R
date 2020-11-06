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
bathy_dat <- load("../data/bathy.dat")
bathy_mat <- load("../data/bathy.mat")
source("functions/distance_function.R")
source("functions/vis_gam_COLORS.R")
source("functions/subset_species.R")
source("functions/TGAM_selection.R")
jet_colors <- colorRampPalette(c("#F7FCFD", "#E0ECF4", "#BFD3E6", "#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D", "#6E016B"))

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
arrowtooth_subset$thr <-  ifelse(arrowtooth_subset$year <= arrowtooth_tgam[[3]], 'before', 'after')
english_tgam <- get_tgam(english_subset, years)
english_subset$thr <-  ifelse(english_subset$year <= english_tgam[[3]], 'before', 'after')
sanddab_tgam <- get_tgam(sanddab_subset, years)
sanddab_subset$thr <-  ifelse(sanddab_subset$year <= sanddab_tgam[[3]], 'before', 'after')
dover_tgam <- get_tgam_woyear(dover_subset, years)
dover_subset$thr <-  ifelse(dover_subset$year <= dover_tgam[[3]], 'before', 'after')
rex_tgam <- get_tgam_woyear(rex_subset, years)
rex_subset$thr <-  ifelse(rex_subset$year <= rex_tgam[[3]], 'before', 'after')
lingcod_tgam <- get_tgam(lingcod_subset, years)
lingcod_subset$thr <-  ifelse(lingcod_subset$year <= lingcod_tgam[[3]], 'before', 'after')
petrale_tgam <- get_tgam(petrale_subset, years)
petrale_subset$thr <-  ifelse(petrale_subset$year <= petrale_tgam[[3]], 'before', 'after')
sablefish_tgam <- get_tgam(sablefish_subset, years)
sablefish_subset$thr <-  ifelse(sablefish_subset$year <= sablefish_tgam[[3]], 'before', 'after')

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
plot_AIC(sablefish_tgam, years) # need to figure out how to select different threshold year

# Validate the results ----
# Calculate distance of each grid point to closest 'positive observation'
# Gives an index, then plot the true and false to get areas of sigificant increase or decrease
subset_distances <- function(tgam, df) {
  nlat = 40 #determine the resolution of the grid
  nlon = 30
  latd = seq(41, 48, length.out = nlat)
  lond = seq(-125,-123.9, length.out = nlon)
  subset_before <- expand.grid(latd, lond)
  names(subset_before) <- c("latitude", "longitude")
  subset_before$year <- 1983 # change depending on species
  subset_before$julian <- 180
  subset_before$thr <- "before"
  subset_before$dist <- NA
  for (k in 1:nrow(subset_before)) {
    dist <- distance_function(
      subset_before$latitude[k],
      subset_before$longitude[k],
      df$latitude,
      df$longitude)
    subset_before$dist[k] <- min(dist)
  }
  return(subset_before)
}

tgam_prediction <- function(tgam, subset_before) {
    predict(tgam[[2]],
            newdata = subset_before,
            se.fit = TRUE,
            type = 'response')
    pred_mean_before <- before_prediction[[1]]
    pred_se_before <- before_prediction[[2]]
    pred_mean_before[subset_before$dist > 10000] <- NA
    pred_se_before[subset_before$dist > 10000] <- NA
    subset_after <- expand.grid(latd, lond)
    names(subset_after) <- c("latitude", "longitude")
    subset_after$year <- 2015
    subset_after$julian <- 180
    subset_after$thr <- "after"
    after_prediction <- predict(tgam[[2]],
                                newdata = subset_after,
                                se.fit = TRUE,
                                type = 'response')
    pred_mean_after <- after_prediction[[1]]
    pred_se_after <- after_prediction[[2]]
    pred_mean_after[subset_before$dist > 10000] <- NA
    pred_se_after[subset_before$dist > 10000] <- NA
    pred_mean_up_before <- pred_mean_before + 1.645 * pred_se_before
    pred_mean_down_before <- pred_mean_before - 1.645 * pred_se_before
    pred_mean_up_after <- pred_mean_after + 1.645 * pred_se_after
    pred_mean_down_after <- pred_mean_after - 1.645 * pred_se_after
    significant_low <- pred_mean_up_after < pred_mean_down_before
    significant_high <- pred_mean_down_after > pred_mean_up_before
    return(list(significant_high, significant_low))
  }

arrowtooth_dist <- subset_distances(arrowtooth_tgam, arrowtooth_subset)
arrowtooth_CI <- tgam_prediction(arrowtooth_tgam, arrowtooth_dist)
english_CI <- tgam_prediction(english_tgam, english_subset)
sanddab_CI <- tgam_prediction(sanddab_tgam, sanddab_subset)
dover_CI <- tgam_prediction(dover_tgam, dover_subset)
rex_CI <- tgam_prediction(rex_tgam, rex_subset)
lingcod_CI <- tgam_prediction(lingcod_tgam, lingcod_subset)
petrale_CI <- tgam_prediction(petrale_tgam, petrale_subset)
sablefish_CI <- tgam_prediction(sablefish_tgam, sablefish_subset)

# Can check the results, will likely get NA
sum(1 * arrowtooth_CI[[2]])
sum(1 * arrowtooth_CI[[1]])

# Make maps of the significant increases and decreases
windows(width = 7, height = 15)
validation_map <- function(df, tgam, conf_int, subset_before, bathy.dat, bathy.mat) {
  significant_high <- conf_int[[1]]
  significant_low <- conf_int[[2]]
  myvis_gam(tgam[[2]],
            view = c('longitude', 'latitude'),
            too.far = 0.025,
            plot.type = 'contour',
            color = 'jet',
            type = 'response',
            cond = list(thr = 'after'),
            contour.col = "gray35",
            xlim = range(df$longitude, na.rm = TRUE) + c(-0.5, 1),
            ylim = range(df$latitude, na.rm = TRUE) + c(-0.5, 0.5),
            cex.lab = 1.5,
            cex.axis = 1.5,
            cex.main = 2,
            xlab = "Longitude °W",
            ylab = "Latitude °N",
            main = deparse(substitute(conf_int)))
  contour(unique(bathy.dat$lon),
          sort(unique(bathy.dat$lat)),
          bathy.mat,
          levels = -seq(200, 200, by = 200),
          labcex = 0.8,
          add = T,
          col = 'black',
          labels = NULL,
          lwd = 1.95)
  points(subset_before$longitude[significant_low],
         subset_before$latitude[significant_low],
         pch = 8,
         cex = 0.5)
  points(subset_before$longitude[significant_high],
         subset_before$latitude[significant_high],
         pch = 16,
         cex = 0.7)
  maps::map('worldHires',
      add = T,
      col = 'peachpuff3',
      fill = T)
}
