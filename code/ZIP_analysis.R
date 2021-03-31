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
    ~ s(year) +
      s(julian) +
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
    ~ s(year) +
      s(julian) +
      s(latitude, longitude)
  ),
  data = species_subset,
  family = ziplss)
  PDO_ziplss <- gam(list(
    count ~ s(PDO, k = 5) +
      s(julian, k = 5) +
      s(latitude, longitude) +
      s(depth_m),
    ~ s(PDO) +
      s(julian) +
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
    ~ s(PDO) +
      s(julian) +
      s(latitude, longitude)
  ),
  data = species_subset,
  family = ziplss)
  NPGO_ziplss <- gam(list(
    count ~ s(NPGO, k = 5) +
      s(julian, k = 5) +
      s(latitude, longitude) +
      s(depth_m),
    ~ s(NPGO) +
      s(julian) +
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
    ~ s(NPGO) +
      s(julian) +
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
  test <- gam(count ~ s(year, k = 5) +
                          s(julian, k = 5) +
                          s(latitude, longitude) +
                          s(depth_m), family = poisson, data = species_subset[species_subset$count > 0,])
}
plot_variable <- function(zip, covariate, bounds, variable, ylabel, yvalues){
  par(
    mar = c(6.4, 7.2, .5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
  plot(zip[[2]],
       pages = 0,
       select = covariate, # 1 = year/PDO/NPGO, 2 = lat/lon, 3 = depth, 4 = julian, 5 = temp
       shade = T,
       shade.col = "lemonchiffon3",
       ylim = bounds,
       xlab = variable,
       ylab = ylabel,
       yaxt = yvalues,
       seWithMean = T,
       scale = 0,
       cex.axis = 3,
       cex.lab = 3,
       family = "serif",
       lwd = 2)
}
location_plot <- function(zip, species_subset) {
  par(
    mar = c(6.4, 7.2, .5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
  myvis_gam(
    zip[[2]],
    view = c('longitude', 'latitude'),
    too.far = 0.03,
    plot.type = 'contour',
    contour.col = contour_col,
    color = "jet" ,
    type = 'response',
    xlim = c(-125.7, -123.6),
    ylim = range(species_subset$latitude, na.rm = TRUE) + c(-.4, .5),
    family = "serif",
    xlab = "Longitude",
    ylab = "Latitude",
    main = " ",
    cex.lab = 2,
    cex.axis = 2)
  maps::map('worldHires',
            add = T,
            col = 'antiquewhite4',
            fill = T)
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

# Create manuscript figures
# Year variable
pdf("../results/ZIP/dover_sole/year_annual.pdf",
    width = 12,
    height = 12)
dover_a_year <- plot_variable(dover_annualzip,
                              covariate = 1,
                              bounds = c(-3.5, 2.1),
                              "Year",
                              "Species Abundance Anomalies",
                              "s")
dev.off()
# Day of year variable
pdf("../results/ZIP/dover_sole/julian_annual.pdf",
    width = 12,
    height = 12)
dover_a_day <- plot_variable(dover_annualzip,
                             covariate = 2,
                             bounds = c(-3.5, 2.1),
                             "Day of Year",
                             " ",
                             "n")
dev.off()
# Depth variable
pdf("../results/ZIP/dover_sole/depth_annual.pdf",
    width = 12,
    height = 12)
dover_a_depth <- plot_variable(dover_annualzip,
                               covariate = 4,
                               bounds = c(-4.2, 2),
                               "Depth (m)",
                               "Species Abundance Anomalies",
                               "s")
dev.off()
# Temperature variable
pdf("../results/ZIP/dover_sole/temp_annual.pdf",
    width = 12,
    height = 12)
dover_a_temp <- plot_variable(dover_annualzip,
                              covariate = 5,
                              bounds = c(-4.2, 2),
                              "Temperature (C)",
                              "",
                              "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/ZIP/dover_sole/location_annual.pdf",
    width = 4,
    height = 11)
location_plot(dover_annualzip, dover_annual)
dev.off()


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

# Create manuscript figures
# Year variable
pdf("../results/ZIP/dover_sole/year_triennial.pdf",
    width = 12,
    height = 12)
dover_a_year <- plot_variable(dover_triennialzip,
                              covariate = 1,
                              bounds = c(-3.5, 2.1),
                              "Year",
                              "Species Abundance Anomalies",
                              "s")
dev.off()
# Day of year variable
pdf("../results/ZIP/dover_sole/julian_triennial.pdf",
    width = 12,
    height = 12)
dover_a_day <- plot_variable(dover_triennialzip,
                             covariate = 2,
                             bounds = c(-3.5, 2.1),
                             "Day of Year",
                             " ",
                             "n")
dev.off()
# Depth variable
pdf("../results/ZIP/dover_sole/depth_triennial.pdf",
    width = 12,
    height = 12)
dover_a_depth <- plot_variable(dover_triennialzip,
                               covariate = 4,
                               bounds = c(-4.2, 2),
                               "Depth (m)",
                               "Species Abundance Anomalies",
                               "s")
dev.off()
# Temperature variable
pdf("../results/ZIP/dover_sole/temp_triennial.pdf",
    width = 12,
    height = 12)
dover_a_temp <- plot_variable(dover_triennialzip,
                              covariate = 5,
                              bounds = c(-4.2, 2),
                              "Temperature (C)",
                              "",
                              "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/ZIP/dover_sole/location_triennial.pdf",
    width = 4,
    height = 11)
location_plot(dover_triennialzip, dover_triennial)
dev.off()
