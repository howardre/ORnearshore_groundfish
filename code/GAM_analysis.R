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
source("functions/vis_gam_COLORS.R")
jet.colors <- colorRampPalette(rev(c("#b2182b", "#d6604d", "#f4a582", "#fddbc7",
                                 "#f7f7f7", "#d1e5f0", "#92c5de", "#4393c3", "#2166ac" )))
contour_col <- rgb(0, 0, 255, max = 255, alpha = 0, names = "white")

# Subset the data to contain only species of interest for each survey ----
# Eight species of interest
# Annual subsets
arrowtooth_annual <- subset_species_temp("Atheresthes stomias", annual, annual_trawl)
english_annual <- subset_species_temp("Parophrys vetulus", annual, annual_trawl)
sanddab_annual <- subset_species_temp("Citharichthys sordidus", annual, annual_trawl)
dover_annual <- subset_species_temp("Microstomus pacificus", annual, annual_trawl)
rex_annual <- subset_species_temp("Glyptocephalus zachirus", annual, annual_trawl)
lingcod_annual <- subset_species_temp("Ophiodon elongatus", annual, annual_trawl)
petrale_annual <- subset_species_temp("Eopsetta jordani", annual, annual_trawl)
sablefish_annual <- subset_species_temp("Anoplopoma fimbria", annual, annual_trawl)

# Triennial subsets
arrowtooth_triennial <- subset_species_temp("Atheresthes stomias", triennial, triennial_trawl)
english_triennial <- subset_species_temp("Parophrys vetulus", triennial, triennial_trawl)
sanddab_triennial <- subset_species_temp("Citharichthys sordidus", triennial, triennial_trawl)
dover_triennial <- subset_species_temp("Microstomus pacificus", triennial, triennial_trawl)
rex_triennial <- subset_species_temp("Glyptocephalus zachirus", triennial, triennial_trawl)
lingcod_triennial <- subset_species_temp("Ophiodon elongatus", triennial, triennial_trawl)
petrale_triennial <- subset_species_temp("Eopsetta jordani", triennial, triennial_trawl)
sablefish_triennial <- subset_species_temp("Anoplopoma fimbria", triennial, triennial_trawl)

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
plot_variable <- function(gam, covariate, bounds, variable, ylabel, yvalues){
  par(
    mar = c(6.4, 7.2, .5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
  plot(gam[[2]],
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
location_plot <- function(gam, species_subset) {
  par(
    mar = c(6.4, 7.2, .5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
  myvis_gam(
    gam[[2]],
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

# ***Arrowtooth Flounder ----
# Annual
cpueGAM_arrow_a <- cpue_GAMs(arrowtooth_annual)
summary(cpueGAM_arrow_a[[2]]) # view the best model
gam_check(cpueGAM_arrow_a)
gam_plot(cpueGAM_arrow_a)
# Year variable
pdf("../results/GAM/arrowtooth_flounder/year_annual.pdf",
    width = 12,
    height = 12)
arrow_a_year <- plot_variable(gam = cpueGAM_arrow_a,
                              covariate = 1,
                              bounds = c(-.6, .7),
                              "Year",
                              "Species Abundance Anomalies",
                              "s")
dev.off()
# Depth variable
pdf("../results/GAM/arrowtooth_flounder/depth_annual.pdf",
     width = 12,
     height = 12)
arrow_a_depth <- plot_variable(cpueGAM_arrow_a,
                               covariate = 3,
                               bounds = c(-3, 2),
                               "Depth (m)",
                               "Species Abundance Anomalies",
                               "s")
dev.off()
# Day of year variable
pdf("../results/GAM/arrowtooth_flounder/julian_annual.pdf",
     width = 12,
     height = 12)
arrow_a_day <- plot_variable(cpueGAM_arrow_a,
                             covariate = 4,
                             bounds = c(-.6, .7),
                             "Day of Year",
                             " ",
                             "n")
dev.off()
# Temperature variable
pdf("../results/GAM/arrowtooth_flounder/temp_annual.pdf",
     width = 12,
     height = 12)
arrow_a_temp <- plot_variable(cpueGAM_arrow_a,
                              covariate = 5,
                              bounds = c(-3, 2),
                              "Temperature (°C)",
                              " ",
                              "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/GAM/arrowtooth_flounder/location_annual.pdf",
    width = 4,
    height = 11)
location_plot(cpueGAM_arrow_a, arrowtooth_annual)
dev.off()

# Triennial
cpueGAM_arrow_t <- cpue_GAMs(arrowtooth_triennial)
summary(cpueGAM_arrow_t[[2]]) # view the best model
gam_check(cpueGAM_arrow_t)
gam_plot(cpueGAM_arrow_t)
pdf("../results/GAM/arrowtooth_flounder/year_triennial.pdf",
    width = 12,
    height = 12)
arrow_t_year <- plot_variable(cpueGAM_arrow_t,
                              covariate = 1,
                              bounds = c(-2, 1),
                              "Year",
                              "Species Abundance Anomalies",
                              "s")
dev.off()
# Depth variable
pdf("../results/GAM/arrowtooth_flounder/depth_triennial.pdf",
    width = 12,
    height = 12)
arrow_t_depth <- plot_variable(cpueGAM_arrow_t,
                               covariate = 3,
                               bounds = c(-2.2, 1.3),
                               "Depth (m)",
                               "Species Abundance Anomalies",
                               "s")
dev.off()
# Day of year variable
pdf("../results/GAM/arrowtooth_flounder/julian_triennial.pdf",
    width = 12,
    height = 12)
arrow_t_day <- plot_variable(cpueGAM_arrow_t,
                             covariate = 4,
                             bounds = c(-2, 1),
                             "Day of Year",
                             " ",
                             "n")
dev.off()
# Temperature variable
pdf("../results/GAM/arrowtooth_flounder/temp_triennial.pdf",
    width = 12,
    height = 12)
arrow_t_temp <- plot_variable(cpueGAM_arrow_t,
                              covariate = 5,
                              bounds = c(-2.2, 1.3),
                              "Temperature (°C)",
                              "",
                              "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/GAM/arrowtooth_flounder/location_triennial.pdf",
    width = 4,
    height = 11)
location_plot(cpueGAM_arrow_t, arrowtooth_triennial)
dev.off()

# ***Dover Sole ----
# Annual
cpueGAM_dover_a <- cpue_GAMs(dover_annual)
summary(cpueGAM_dover_a[[2]]) # view the best model
gam_check(cpueGAM_dover_a)
gam_plot(cpueGAM_dover_a)
# Year variable
pdf("../results/GAM/dover_sole/year_annual.pdf",
    width = 12,
    height = 12)
dover_a_year <- plot_variable(cpueGAM_dover_a,
                              covariate = 1,
                              bounds = c(-.42, .7),
                              "Year",
                              "Species Abundance Anomalies",
                              "s")
dev.off()
# Depth variable
pdf("../results/GAM/dover_sole/depth_annual.pdf",
    width = 12,
    height = 12)
dover_a_depth <- plot_variable(cpueGAM_dover_a,
                               covariate = 3,
                               bounds = c(-3.2, 2.1),
                               "Depth (m)",
                               "Species Abundance Anomalies",
                               "s")
dev.off()
# Day of year variable
pdf("../results/GAM/dover_sole/julian_annual.pdf",
    width = 12,
    height = 12)
dover_a_day <- plot_variable(cpueGAM_dover_a,
                             covariate = 4,
                             bounds = c(-.42, .7),
                             "Day of Year",
                             " ",
                             "n")
dev.off()
# Temperature variable
pdf("../results/GAM/dover_sole/temp_annual.pdf",
    width = 12,
    height = 12)
dover_a_temp <- plot_variable(cpueGAM_dover_a,
                              covariate = 5,
                              bounds = c(-3.2, 2.1),
                              "Temperature (°C)",
                              "",
                              "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/GAM/dover_sole/location_annual.pdf",
    width = 4,
    height = 11)
location_plot(cpueGAM_dover_a, dover_annual)
dev.off()

# Triennial
cpueGAM_dover_t <- cpue_GAMs(dover_triennial)
summary(cpueGAM_dover_t[[2]]) # view the best model
gam_check(cpueGAM_dover_t)
gam_plot(cpueGAM_dover_t)
pdf("../results/GAM/dover_sole/year_triennial.pdf",
    width = 12,
    height = 12)
dover_t_year <- plot_variable(cpueGAM_dover_t,
                              covariate = 1,
                              bounds = c(-1.7, 1.5),
                              "Year",
                              "Species Abundance Anomalies",
                              "s")
dev.off()
# Depth variable
pdf("../results/GAM/dover_sole/depth_triennial.pdf",
    width = 12,
    height = 12)
dover_t_depth <- plot_variable(cpueGAM_dover_t,
                               covariate = 3,
                               bounds = c(-3, 2),
                               "Depth (m)",
                               "Species Abundance Anomalies",
                               "s")
dev.off()
# Day of year variable
pdf("../results/GAM/dover_sole/julian_triennial.pdf",
    width = 12,
    height = 12)
dover_t_day <- plot_variable(cpueGAM_dover_t,
                             covariate = 4,
                             bounds = c(-1.6, 1),
                             "Day of Year",
                             " ",
                             "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/GAM/dover_sole/location_triennial.pdf",
    width = 4,
    height = 11)
location_plot(cpueGAM_dover_t, dover_triennial)
dev.off()

# ***English Sole ----
# Annual
cpueGAM_english_a <- cpue_GAMs(english_annual)
summary(cpueGAM_english_a[[2]]) # view the best model
gam_check(cpueGAM_english_a)
gam_plot(cpueGAM_english_a)
# Year variable
pdf("../results/GAM/english_sole/year_annual.pdf",
    width = 12,
    height = 12)
english_a_year <- plot_variable(cpueGAM_english_a,
                                covariate = 1,
                                bounds = c(-.4, .8),
                                "Year",
                                "Species Abundance Anomalies",
                                "s")
dev.off()
# Depth variable
pdf("../results/GAM/english_sole/depth_annual.pdf",
    width = 12,
    height = 12)
english_a_depth <- plot_variable(cpueGAM_english_a,
                                 covariate = 3,
                                 bounds = c(-1.1, 1.6),
                                 "Depth (m)",
                                 "Species Abundance Anomalies",
                                 "s")
dev.off()
# Day of year variable
pdf("../results/GAM/english_sole/julian_annual.pdf",
    width = 12,
    height = 12)
english_a_day <- plot_variable(cpueGAM_english_a,
                               covariate = 4,
                               bounds = c(-.4, .8),
                               "Day of Year",
                               " ",
                               "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/GAM/english_sole/location_annual.pdf",
    width = 4,
    height = 11)
location_plot(cpueGAM_english_a, english_annual)
dev.off()

# Triennial
cpueGAM_english_t <- cpue_GAMs(english_triennial)
summary(cpueGAM_english_t[[2]]) # view the best model
gam_check(cpueGAM_english_t)
gam_plot(cpueGAM_english_t)
# Year variable
pdf("../results/GAM/english_sole/year_triennial.pdf",
    width = 12,
    height = 12)
english_t_year <- plot_variable(cpueGAM_english_t,
                                covariate = 1,
                                bounds = c(-1.5, 1.5),
                                "Year",
                                "Species Abundance Anomalies",
                                "s")
dev.off()
# Depth variable
pdf("../results/GAM/english_sole/depth_triennial.pdf",
    width = 12,
    height = 12)
english_t_depth <- plot_variable(cpueGAM_english_t,
                                 covariate = 3,
                                 bounds = c(-3, 2),
                                 "Depth (m)",
                                 "Species Abundance Anomalies",
                                 "s")
dev.off()
# Day of year variable
pdf("../results/GAM/english_sole/julian_triennial.pdf",
    width = 12,
    height = 12)
english_t_day <- plot_variable(cpueGAM_english_t,
                               covariate = 4,
                               bounds = c(-1.5, 1.5),
                               "Day of Year",
                               " ",
                               "n")
dev.off()
# Temperature Variable
pdf("../results/GAM/english_sole/temp_triennial.pdf",
    width = 12,
    height = 12)
english_t_temp <- plot_variable(cpueGAM_english_t,
                                covariate = 5,
                                bounds = c(-3, 2),
                                "Temperature (°C)",
                                "",
                                "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/GAM/english_sole/location_triennial.pdf",
    width = 4,
    height = 11)
location_plot(cpueGAM_english_t, english_triennial)
dev.off()

# ***Lingcod ----
# Annual
cpueGAM_lingcod_a <- cpue_GAMs(lingcod_annual)
summary(cpueGAM_lingcod_a[[2]]) # view the best model
gam_check(cpueGAM_lingcod_a)
gam_plot(cpueGAM_lingcod_a)
# Year variable
pdf("../results/GAM/lingcod/year_annual.pdf",
    width = 12,
    height = 12)
lingcod_a_year <- plot_variable(cpueGAM_lingcod_a,
                                covariate = 1,
                                bounds = c(-.5, .4),
                                "Year",
                                "Species Abundance Anomalies",
                                "s")
dev.off()
# Depth variable
pdf("../results/GAM/lingcod/depth_annual.pdf",
    width = 12,
    height = 12)
lingcod_a_depth <- plot_variable(cpueGAM_lingcod_a,
                                 covariate = 3,
                                 bounds = c(-1, 1),
                                 "Depth (m)",
                                 "Species Abundance Anomalies",
                                 "s")
dev.off()
# Day of year variable
pdf("../results/GAM/lingcod/julian_annual.pdf",
    width = 12,
    height = 12)
lingcod_a_day <- plot_variable(cpueGAM_lingcod_a,
                               covariate = 4,
                               bounds = c(-.5, .4),
                               "Day of Year",
                               " ",
                               "n")
dev.off()
# Temperature variable
pdf("../results/GAM/lingcod/temp_annual.pdf",
    width = 12,
    height = 12)
lingcod_a_temp <- plot_variable(cpueGAM_lingcod_a,
                                covariate = 5,
                                bounds = c(-1, 1),
                                "Temperature (°C)",
                                "",
                                "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/GAM/lingcod/location_annual.pdf",
    width = 4,
    height = 11)
location_plot(cpueGAM_lingcod_a, lingcod_annual)
dev.off()

# Triennial
cpueGAM_lingcod_t <- cpue_GAMs(lingcod_triennial)
summary(cpueGAM_lingcod_t[[2]]) # view the best model
gam_check(cpueGAM_lingcod_t)
gam_plot(cpueGAM_lingcod_t)
# NPGO variable
pdf("../results/GAM/lingcod/NPGO_triennial.pdf",
    width = 12,
    height = 12)
lingcod_t_npgo <- plot_variable(cpueGAM_lingcod_t,
                                covariate = 1,
                                bounds = c(-.4, 1),
                                "NPGO",
                                "",
                                "n")
dev.off()
# Depth variable
pdf("../results/GAM/lingcod/depth_triennial.pdf",
    width = 12,
    height = 12)
lingcod_t_depth <- plot_variable(cpueGAM_lingcod_t,
                                 covariate = 3,
                                 bounds = c(-.4, 1),
                                 "Depth (m)",
                                 "Species Abundance Anomalies",
                                 "s")
dev.off()
# Day of year variable
pdf("../results/GAM/lingcod/julian_triennial.pdf",
    width = 12,
    height = 12)
lingcod_t_day <- plot_variable(cpueGAM_lingcod_t,
                               covariate = 4,
                               bounds = c(-.25, .4),
                               "Day of Year",
                               " ",
                               "n")
dev.off()
# Temperature variable
pdf("../results/GAM/lingcod/temp_triennial.pdf",
    width = 12,
    height = 12)
lingcod_t_temp <- plot_variable(cpueGAM_lingcod_t,
                                covariate = 5,
                                bounds = c(-.4, 1),
                                "Temperature (°C)",
                                "",
                                "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/GAM/lingcod/location_triennial.pdf",
    width = 4,
    height = 11)
location_plot(cpueGAM_lingcod_t, lingcod_triennial)
dev.off()

# ***Petrale Sole ----
# Annual
cpueGAM_petrale_a <- cpue_GAMs(petrale_annual)
summary(cpueGAM_petrale_a[[2]]) # view the best model
gam_check(cpueGAM_petrale_a)
gam_plot(cpueGAM_petrale_a)
# Year variable
pdf("../results/GAM/petrale_sole/year_annual.pdf",
    width = 12,
    height = 12)
petrale_a_year <- plot_variable(cpueGAM_petrale_a,
                                covariate = 1,
                                bounds = c(-.7, .5),
                                "Year",
                                "Species Abundance Anomalies",
                                "s")
dev.off()
# Depth variable
pdf("../results/GAM/petrale_sole/depth_annual.pdf",
    width = 12,
    height = 12)
petrale_a_depth <- plot_variable(cpueGAM_petrale_a,
                                 covariate = 3,
                                 bounds = c(-1, .3),
                                 "Depth (m)",
                                 "Species Abundance Anomalies",
                                 "s")
dev.off()
# Day of year variable
pdf("../results/GAM/petrale_sole/julian_annual.pdf",
    width = 12,
    height = 12)
petrale_a_day <- plot_variable(cpueGAM_petrale_a,
                               covariate = 4,
                               bounds = c(-.7, .5),
                               "Day of Year",
                               " ",
                               "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/GAM/petrale_sole/location_annual.pdf",
    width = 4,
    height = 11)
location_plot(cpueGAM_petrale_a, petrale_annual)
dev.off()

# Triennial
cpueGAM_petrale_t <- cpue_GAMs(petrale_triennial)
summary(cpueGAM_petrale_t[[2]]) # view the best model
gam_check(cpueGAM_petrale_t)
gam_plot(cpueGAM_petrale_t)
# Year variable
pdf("../results/GAM/petrale_sole/year_triennial.pdf",
    width = 12,
    height = 12)
petrale_t_year <- plot_variable(cpueGAM_petrale_t,
                                covariate = 1,
                                bounds = c(-.4, .6),
                                "Year",
                                "Species Abundance Anomalies",
                                "s")
dev.off()
# Depth variable
pdf("../results/GAM/petrale_sole/depth_triennial.pdf",
    width = 12,
    height = 12)
petrale_t_depth <- plot_variable(cpueGAM_petrale_t,
                                 covariate = 3,
                                 bounds = c(-.8, .5),
                                 "Depth (m)",
                                 "Species Abundance Anomalies",
                                 "s")
dev.off()
# Day of year variable
pdf("../results/GAM/petrale_sole/julian_triennial.pdf",
    width = 12,
    height = 12)
petrale_t_day <- plot_variable(cpueGAM_petrale_t,
                               covariate = 4,
                               bounds = c(-.4, .6),
                               "Day of Year",
                               " ",
                               "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/GAM/petrale_sole/location_triennial.pdf",
    width = 4,
    height = 11)
location_plot(cpueGAM_petrale_t, petrale_triennial)
dev.off()

# ***Pacific Sanddab ----
# Annual
cpueGAM_sanddab_a <- cpue_GAMs(sanddab_annual)
summary(cpueGAM_sanddab_a[[2]]) # view the best model
gam_check(cpueGAM_sanddab_a)
gam_plot(cpueGAM_sanddab_a)
# Year variable
pdf("../results/GAM/pacific_sanddab/year_annual.pdf",
    width = 12,
    height = 12)
sanddab_a_year <- plot_variable(cpueGAM_sanddab_a,
                                  covariate = 1,
                                  bounds = c(-.5, .6),
                                  "Year",
                                  "Species Abundance Anomalies",
                                  "s")
dev.off()
# Depth variable
pdf("../results/GAM/pacific_sanddab/depth_annual.pdf",
    width = 12,
    height = 12)
sanddab_a_depth <- plot_variable(cpueGAM_sanddab_a,
                                   covariate = 3,
                                   bounds = c(-4.1, 1.1),
                                   "Depth (m)",
                                   "Species Abundance Anomalies",
                                   "s")
dev.off()
# Day of year variable
pdf("../results/GAM/pacific_sanddab/julian_annual.pdf",
    width = 12,
    height = 12)
sanddab_a_day <- plot_variable(cpueGAM_sanddab_a,
                                 covariate = 4,
                                 bounds = c(-.5, .6),
                                 "Day of Year",
                                 " ",
                                 "n")
dev.off()
# Temperature variable
pdf("../results/GAM/pacific_sanddab/temp_annual.pdf",
    width = 12,
    height = 12)
sanddab_a_temp <- plot_variable(cpueGAM_sanddab_a,
                                  covariate = 5,
                                  bounds = c(-4.1, 1.1),
                                  "Temperature (°C)",
                                  "",
                                  "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/GAM/pacific_sanddab/location_annual.pdf",
    width = 4,
    height = 11)
location_plot(cpueGAM_sanddab_a, sanddab_annual)
dev.off()

# Triennial
cpueGAM_sanddab_t <- cpue_GAMs(sanddab_triennial)
summary(cpueGAM_sanddab_t[[2]]) # view the best model
gam_check(cpueGAM_sanddab_t)
gam_plot(cpueGAM_sanddab_t)
# Year variable
pdf("../results/GAM/pacific_sanddab/year_triennial.pdf",
    width = 12,
    height = 12)
sanddab_t_year <- plot_variable(cpueGAM_sanddab_t,
                                  covariate = 1,
                                  bounds = c(-4.2, 1),
                                  "Year",
                                  "Species Abundance Anomalies",
                                  "s")
dev.off()
# Depth variable
pdf("../results/GAM/pacific_sanddab/depth_triennial.pdf",
    width = 12,
    height = 12)
sanddab_t_depth <- plot_variable(cpueGAM_sanddab_t,
                                   covariate = 3,
                                   bounds = c(-4, 2),
                                   "Depth (m)",
                                   "Species Abundance Anomalies",
                                   "s")
dev.off()
# Day of year variable
pdf("../results/GAM/pacific_sanddab/julian_triennial.pdf",
    width = 12,
    height = 12)
sanddab_t_day <- plot_variable(cpueGAM_sanddab_t,
                                 covariate = 4,
                                 bounds = c(-4.2, 1),
                                 "Day of Year",
                                 " ",
                                 "n")
dev.off()
# Temperature variable
pdf("../results/GAM/pacific_sanddab/temp_triennial.pdf",
    width = 12,
    height = 12)
sanddab_t_temp <- plot_variable(cpueGAM_sanddab_t,
                                  covariate = 5,
                                  bounds = c(-4, 2),
                                  "Temperature (°C)",
                                  "",
                                  "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/GAM/pacific_sanddab/location_triennial.pdf",
    width = 4,
    height = 11)
location_plot(cpueGAM_sanddab_t, sanddab_triennial)
dev.off()
# ***Rex Sole ----
# Annual
cpueGAM_rex_a <- cpue_GAMs(rex_annual)
summary(cpueGAM_rex_a[[2]]) # view the best model
gam_check(cpueGAM_rex_a)
gam_plot(cpueGAM_rex_a)
# Year variable
pdf("../results/GAM/rex_sole/year_annual.pdf",
    width = 12,
    height = 12)
rex_a_year <- plot_variable(cpueGAM_rex_a,
                            covariate = 1,
                            bounds = c(-.7, .8),
                            "Year",
                            "Species Abundance Anomalies",
                            "s")
dev.off()
# Depth variable
pdf("../results/GAM/rex_sole/depth_annual.pdf",
    width = 12,
    height = 12)
rex_a_depth <- plot_variable(cpueGAM_rex_a,
                             covariate = 3,
                             bounds = c(-1.8, 1.5),
                             "Depth (m)",
                             "Species Abundance Anomalies",
                             "s")
dev.off()
# Day of year variable
pdf("../results/GAM/rex_sole/julian_annual.pdf",
    width = 12,
    height = 12)
rex_a_day <- plot_variable(cpueGAM_rex_a,
                           covariate = 4,
                           bounds = c(-.7, .8),
                           "Day of Year",
                           " ",
                           "n")
dev.off()
# Temperature Variable
pdf("../results/GAM/rex_sole/temp_annual.pdf",
    width = 12,
    height = 12)
rex_a_temp <- plot_variable(cpueGAM_rex_a,
                            covariate = 5,
                            bounds = c(-1.8, 1.5),
                            "Temperature (°C)",
                            "",
                            "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/GAM/rex_sole/location_annual.pdf",
    width = 4,
    height = 11)
location_plot(cpueGAM_rex_a, rex_annual)
dev.off()

# Triennial
cpueGAM_rex_t <- cpue_GAMs(rex_triennial)
summary(cpueGAM_rex_t[[2]]) # view the best model
gam_check(cpueGAM_rex_t)
gam_plot(cpueGAM_rex_t)
# Year variable
pdf("../results/GAM/rex_sole/year_triennial.pdf",
    width = 12,
    height = 12)
rex_t_year <- plot_variable(cpueGAM_rex_t,
                            covariate = 1,
                            bounds = c(-3, 1),
                            "Year",
                            "Species Abundance Anomalies",
                            "s")
dev.off()
# Depth variable
pdf("../results/GAM/rex_sole/depth_triennial.pdf",
    width = 12,
    height = 12)
rex_t_depth <- plot_variable(cpueGAM_rex_t,
                             covariate = 3,
                             bounds = c(-.7, 1.2),
                             "Depth (m)",
                             "Species Abundance Anomalies",
                             "s")
dev.off()
# Day of year variable
pdf("../results/GAM/rex_sole/julian_triennial.pdf",
    width = 12,
    height = 12)
rex_t_day <- plot_variable(cpueGAM_rex_t,
                           covariate = 4,
                           bounds = c(-3, 1),
                           "Day of Year",
                           " ",
                           "n")
dev.off()
# Temperature Variable
pdf("../results/GAM/rex_sole/temp_triennial.pdf",
    width = 12,
    height = 12)
rex_t_temp <- plot_variable(cpueGAM_rex_t,
                            covariate = 5,
                            bounds = c(-.7, 1.2),
                            "Temperature (°C)",
                            "",
                            "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/GAM/rex_sole/location_triennial.pdf",
    width = 4,
    height = 11)
location_plot(cpueGAM_rex_t, rex_triennial)
dev.off()

# ***Sablefish ----
# Annual
cpueGAM_sablefish_a <- cpue_GAMs(sablefish_annual)
summary(cpueGAM_sablefish_a[[2]]) # view the best model
gam_check(cpueGAM_sablefish_a)
gam_plot(cpueGAM_sablefish_a)
# Year variable
pdf("../results/GAM/sablefish/year_annual.pdf",
    width = 12,
    height = 12)
sablefish_a_year <- plot_variable(cpueGAM_sablefish_a,
                                  covariate = 1,
                                  bounds = c(-.7, .6),
                                  "Year",
                                  "Species Abundance Anomalies",
                                  "s")
dev.off()
# Depth variable
pdf("../results/GAM/sablefish/depth_annual.pdf",
    width = 12,
    height = 12)
sablefish_a_depth <- plot_variable(cpueGAM_sablefish_a,
                                   covariate = 3,
                                   bounds = c(-1.5, 1),
                                   "Depth (m)",
                                   "Species Abundance Anomalies",
                                   "s")
dev.off()
# Day of year variable
pdf("../results/GAM/sablefish/julian_annual.pdf",
    width = 12,
    height = 12)
sablefish_a_day <- plot_variable(cpueGAM_sablefish_a,
                                 covariate = 4,
                                 bounds = c(-.7, .6),
                                 "Day of Year",
                                 " ",
                                 "n")
dev.off()
# Temperature variable
pdf("../results/GAM/sablefish/temp_annual.pdf",
    width = 12,
    height = 12)
sablefish_a_temp <- plot_variable(cpueGAM_sablefish_a,
                                  covariate = 5,
                                  bounds = c(-1.5, 1),
                                  "Temperature (°C)",
                                  "",
                                  "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/GAM/sablefish/location_annual.pdf",
    width = 4,
    height = 11)
location_plot(cpueGAM_sablefish_a, sablefish_annual)
dev.off()

# Triennial
cpueGAM_sablefish_t <- cpue_GAMs(sablefish_triennial)
summary(cpueGAM_sablefish_t[[2]]) # view the best model
gam_check(cpueGAM_sablefish_t)
gam_plot(cpueGAM_sablefish_t)
# Year variable
pdf("../results/GAM/sablefish/year_triennial.pdf",
    width = 12,
    height = 12)
sablefish_t_year <- plot_variable(cpueGAM_sablefish_t,
                                  covariate = 1,
                                  bounds = c(-1.5, 1.5),
                                  "Year",
                                  "Species Abundance Anomalies",
                                  "s")
dev.off()
# Depth variable
pdf("../results/GAM/sablefish/depth_triennial.pdf",
    width = 12,
    height = 12)
sablefish_t_depth <- plot_variable(cpueGAM_sablefish_t,
                                   covariate = 3,
                                   bounds = c(-2.5, 1),
                                   "Depth (m)",
                                   "Species Abundance Anomalies",
                                   "s")
dev.off()
# Day of year variable
pdf("../results/GAM/sablefish/julian_triennial.pdf",
    width = 12,
    height = 12)
sablefish_t_day <- plot_variable(cpueGAM_sablefish_t,
                                 covariate = 4,
                                 bounds = c(-1.5, 1.5),
                                 "Day of Year",
                                 " ",
                                 "n")
dev.off()
# Temperature variable
pdf("../results/GAM/sablefish/temp_triennial.pdf",
    width = 12,
    height = 12)
sablefish_t_temp <- plot_variable(cpueGAM_sablefish_t,
                                  covariate = 5,
                                  bounds = c(-2.5, 1),
                                  "Temperature (°C)",
                                  "",
                                  "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/GAM/sablefish/location_triennial.pdf",
    width = 4,
    height = 11)
location_plot(cpueGAM_sablefish_t, sablefish_triennial)
dev.off()

# ***Proposal Figues----
load("../data/bathy.dat")
load("../data/bathy.mat")

# Edit location function
location_plot2 <- function(gam, species_subset, title) {
  par(
    mar = c(5, 7.2, 1.4, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(4, 1.7, 0))
  myvis_gam(
    gam[[2]],
    view = c('longitude', 'latitude'),
    too.far = 0.03,
    plot.type = 'contour',
    col = contour_col,
    color = "jet" ,
    type = 'response',
    xlim = c(-125.7, -123.6),
    ylim = range(species_subset$latitude, na.rm = TRUE) + c(-.4, .5),
    family = "serif",
    xlab = "Longitude",
    ylab = "Latitude",
    main = title,
    cex.lab = 1.7,
    cex.axis = 1.7,
    cex.main = 1.7)
  maps::map('worldHires',
            add = T,
            col = 'antiquewhite4',
            fill = T)
}

location_plot2 <- function(gam, species_subset, title) {
  par(
    mar = c(5, 7.2, 1.4, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(4, 1.7, 0))
  myvis_gam(
    gam[[2]],
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
    main = title,
    cex.lab = 1.7,
    cex.axis = 1.7,
    cex.main = 1.7)
  maps::map('worldHires',
            add = T,
            col = 'antiquewhite4',
            fill = T)
}


# Sablefish
pdf("../results/GAM/sablefish/proposal.pdf",
    width = 4,
    height = 11,
    family = "serif")
location_plot2(cpueGAM_sablefish_a, sablefish_annual, "Sablefish")
symbols(
  sablefish_annual$longitude[sablefish_annual$lncpue > 0],
  sablefish_annual$latitude[sablefish_annual$lncpue > 0],
  circles = sablefish_annual$lncpue[sablefish_annual$lncpue > 0],
  inches = 0.07,
  add = T,
  fg = 'gray47',
  bg = alpha("gray59", 0.5))
contour(unique(bathy.dat$lon),
        sort(unique(bathy.dat$lat)),
        bathy.mat,
        levels = -seq(200, 200, by = 200),
        labcex = 1.1,
        add = T,
        col = 'black',
        labels = NULL,
        lwd = 1.95)
image.plot(legend.only = T,
           zlim = c(0.6, 2.1),
           col = jet.colors(100),
           legend.shrink = 0.2,
           smallplot = c(.45, .47, .11, .2),
           legend.cex = 1,
           legend.lab = "log(CPUE+1)",
           axis.args = list(cex.axis = 1),
           legend.width = 1,
           legend.line = 2)
dev.off()

# Dover
dover_subset <- dover_annual[sample(nrow(dover_annual), 1000), ]
pdf("../results/GAM/dover_sole/proposal.pdf",
    width = 4,
    height = 11,
    family = "serif")
location_plot2(cpueGAM_dover_a, dover_annual, "Dover Sole")
symbols(
  dover_subset$longitude[dover_subset$lncpue > 0],
  dover_subset$latitude[dover_subset$lncpue > 0],
  circles = dover_subset$lncpue[dover_subset$lncpue > 0],
  inches = 0.07,
  add = T,
  fg = 'gray47',
  bg = alpha("gray59", 0.5))
contour(unique(bathy.dat$lon),
        sort(unique(bathy.dat$lat)),
        bathy.mat,
        levels = -seq(200, 200, by = 200),
        labcex = 1.1,
        add = T,
        col = 'black',
        labels = NULL,
        lwd = 1.95)
image.plot(legend.only = T,
           zlim = c(2, 6.2),
           col = jet.colors(100),
           legend.shrink = 0.2,
           smallplot = c(.45, .47, .11, .2),
           legend.cex = 1,
           legend.lab = "log(CPUE+1)",
           axis.args = list(cex.axis = 1),
           legend.width = 1,
           legend.line = 2)
dev.off()

# Lingcod
pdf("../results/GAM/lingcod/proposal.pdf",
    width = 4,
    height = 11,
    family = "serif")
location_plot2(cpueGAM_lingcod_a, lingcod_annual, "Lingcod")
symbols(
  lingcod_annual$longitude[lingcod_annual$lncpue > 0],
  lingcod_annual$latitude[lingcod_annual$lncpue > 0],
  circles = lingcod_annual$lncpue[lingcod_annual$lncpue > 0],
  inches = 0.07,
  add = T,
  fg = 'gray47',
  bg = alpha("gray59", 0.5))
contour(unique(bathy.dat$lon),
        sort(unique(bathy.dat$lat)),
        bathy.mat,
        levels = -seq(200, 200, by = 200),
        labcex = 1.1,
        add = T,
        col = 'black',
        labels = NULL,
        lwd = 1.95)
image.plot(legend.only = T,
           zlim = c(0.2, 2),
           col = jet.colors(100),
           legend.shrink = 0.2,
           smallplot = c(.45, .47, .11, .2),
           legend.cex = 1,
           legend.lab = "log(CPUE+1)",
           axis.args = list(cex.axis = 1),
           legend.width = 1,
           legend.line = 2)
dev.off()

# Petrale
petrale_subset <- petrale_annual[sample(nrow(petrale_annual), 1000), ]
pdf("../results/GAM/petrale_sole/proposal.pdf",
    width = 4,
    height = 11,
    family = "serif")
location_plot2(cpueGAM_petrale_a, petrale_annual, "Petrale Sole")
symbols(
  petrale_subset$longitude[petrale_subset$lncpue > 0],
  petrale_subset$latitude[petrale_subset$lncpue > 0],
  circles = petrale_subset$lncpue[petrale_subset$lncpue > 0],
  inches = 0.07,
  add = T,
  fg = 'gray47',
  bg = alpha("gray59", 0.5))
contour(unique(bathy.dat$lon),
        sort(unique(bathy.dat$lat)),
        bathy.mat,
        levels = -seq(200, 200, by = 200),
        labcex = 1.1,
        add = T,
        col = 'black',
        labels = NULL,
        lwd = 1.95)
image.plot(legend.only = T,
           zlim = c(2, 3),
           col = jet.colors(100),
           legend.shrink = 0.2,
           smallplot = c(.45, .47, .11, .2),
           legend.cex = 1,
           legend.lab = "log(CPUE+1)",
           axis.args = list(cex.axis = 1),
           legend.width = 1,
           legend.line = 2)
dev.off()
