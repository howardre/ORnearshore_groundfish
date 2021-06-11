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
triennial_trawl <- triennial_trawl[triennial_trawl$bottom_temp < 10, ]
triennial_trawl <- triennial_trawl[triennial_trawl$year < 2004, ]
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


# Create functions to select the best ZIP for each species ----
ZIP_selection <- function(species_subset){
  year_ziplss <- gam(list(
    count ~ offset(log(area_swept)) +
      s(year, k = 5) +
      s(julian, k = 5) +
      s(latitude, longitude) +
      s(depth_m, k = 5) +
      s(grain_size, k = 5),
    ~ offset(log(area_swept)) +
      s(year, k = 5) +
      s(julian, k = 5) +
      s(latitude, longitude)),
  data = species_subset,
  family = ziplss)

  yeartemp_ziplss <-  gam(list(
    count ~ offset(log(area_swept)) +
      s(year, k = 5) +
      s(julian, k = 5) +
      s(latitude, longitude) +
      s(depth_m, k = 5) +
      s(grain_size, k = 5) +
      s(bottom_temp, k = 5),
    ~ offset(log(area_swept)) +
      s(year, k = 5) +
      s(julian, k = 5) +
      s(latitude, longitude)),
  data = species_subset,
  family = ziplss)

  PDO_ziplss <- gam(list(
    count ~ offset(log(area_swept)) +
      s(PDO, k = 5) +
      s(julian, k = 5) +
      s(latitude, longitude) +
      s(depth_m, k = 5) +
      s(grain_size, k = 5),
    ~ offset(log(area_swept)) +
      s(PDO, k = 5) +
      s(julian, k = 5) +
      s(latitude, longitude)),
  data = species_subset,
  family = ziplss)

  PDOtemp_ziplss <-  gam(list(
    count ~ offset(log(area_swept)) +
      s(PDO, k = 5) +
      s(julian, k = 5) +
      s(latitude, longitude) +
      s(depth_m, k = 5) +
      s(grain_size, k = 5) +
      s(bottom_temp, k = 5),
    ~ offset(log(area_swept)) +
      s(PDO, k = 5) +
      s(julian, k = 5) +
      s(latitude, longitude)),
  data = species_subset,
  family = ziplss)

  NPGO_ziplss <- gam(list(
    count ~ offset(log(area_swept)) +
      s(NPGO, k = 5) +
      s(julian, k = 5) +
      s(latitude, longitude) +
      s(depth_m, k = 5) +
      s(grain_size, k = 5),
    ~ offset(log(area_swept)) +
      s(NPGO, k = 5) +
      s(julian, k = 5) +
      s(latitude, longitude)),
  data = species_subset,
  family = ziplss)

  NPGOtemp_ziplss <-  gam(list(
    count ~ offset(log(area_swept)) +
      s(NPGO, k = 5) +
      s(julian, k = 5) +
      s(latitude, longitude) +
      s(depth_m, k = 5) +
      s(grain_size, k = 5) +
      s(bottom_temp, k = 5),
    ~ offset(log(area_swept)) +
      s(NPGO, k = 5) +
      s(julian, k = 5) +
      s(latitude, longitude)),
  data = species_subset,
  family = ziplss)

  ziplss_list <- list(year_ziplss, yeartemp_ziplss, PDO_ziplss, PDOtemp_ziplss, NPGO_ziplss, NPGOtemp_ziplss)
  best_ziplss <- ziplss_list[[which.min(sapply(1:length(ziplss_list),
                                         function(x) min(ziplss_list[[x]]$aic)))]] # would like to also select by AIC
  return_list <- list(ziplss_list, best_ziplss)
}
ZIP_test <- function(species_subset) {
  test <- gam(count ~ offset(log(area_swept)) +
                s(year, k = 5) +
                s(julian, k = 5) +
                s(latitude, longitude) +
                s(depth_m, k = 5) +
                s(grain_size, k = 5),
              family = poisson, data = species_subset[species_subset$count > 0,])
  test <- gam(count ~ offset(log(area_swept)) +
                s(year, k = 5) +
                s(julian, k = 5) +
                s(latitude, longitude) +
                s(depth_m), family = poisson, data = species_subset[species_subset$count > 0,])
}
plot_variable <- function(zip, covariate, bounds, variable, ylabel, yvalues){
  par(
    mar = c(6.4, 7.6, .5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 1.7, 0))
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
       cex.axis = 4,
       cex.lab = 4.5,
       family = "serif",
       lwd = 2)
}
location_plot <- function(zip, species_subset) {
  par(
    mar = c(1, 7.2, .5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 1, 0))
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
    xlab = " ",
    ylab = "Latitude",
    main = " ",
    tck = 0,
    xaxt = "n",
    cex.lab = 3.5,
    cex.axis = 3.5)
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
gam.check(dover_annualcheck) # not a normal distribution
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
                              bounds = c(-1, 1),
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
                             bounds = c(-1, 1),
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
                               bounds = c(-4.2, 1.5),
                               "Depth (m)",
                               "Species Abundance Anomalies",
                               "s")
dev.off()
# Temperature variable
pdf("../results/ZIP/dover_sole/temp_annual.pdf",
    width = 12,
    height = 12)
dover_a_temp <- plot_variable(dover_annualzip,
                              covariate = 6,
                              bounds = c(-4.2, 1.5),
                              "Temperature (C)",
                              "",
                              "n")
dev.off()
# Grain size variable
pdf("../results/ZIP/dover_sole/grain_annual.pdf",
    width = 12,
    height = 12)
dover_a_grain <- plot_variable(dover_annualzip,
                              covariate = 5,
                              bounds = c(-4.2, 1.5),
                              "Grain Size (phi)",
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
dover_t_year <- plot_variable(dover_triennialzip,
                              covariate = 1,
                              bounds = c(-3, 1.2),
                              "Year",
                              "Species Abundance Anomalies",
                              "s")
dev.off()
# Day of year variable
pdf("../results/ZIP/dover_sole/julian_triennial.pdf",
    width = 12,
    height = 12)
dover_t_day <- plot_variable(dover_triennialzip,
                             covariate = 2,
                             bounds = c(-3, 1.2),
                             "Day of Year",
                             " ",
                             "n")
dev.off()
# Depth variable
pdf("../results/ZIP/dover_sole/depth_triennial.pdf",
    width = 12,
    height = 12)
dover_t_depth <- plot_variable(dover_triennialzip,
                               covariate = 4,
                               bounds = c(-1, 1),
                               "Depth (m)",
                               "Species Abundance Anomalies",
                               "s")
dev.off()
# Temperature variable
pdf("../results/ZIP/dover_sole/temp_triennial.pdf",
    width = 12,
    height = 12)
dover_t_temp <- plot_variable(dover_triennialzip,
                              covariate = 6,
                              bounds = c(-1, 1),
                              "Temperature (C)",
                              "",
                              "n")
dev.off()
# Grain size variable
pdf("../results/ZIP/dover_sole/grain_triennial.pdf",
    width = 12,
    height = 12)
dover_t_grain <- plot_variable(dover_triennialzip,
                              covariate = 5,
                              bounds = c(-1, 1),
                              "Grain Size (phi)",
                              "",
                              "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/ZIP/dover_sole/location_triennial.pdf",
    width = 4,
    height = 11)
location_plot(dover_triennialzip, dover_triennial)
dev.off()


# ****Arrowtooth Flounder ----
# Annual
arrowtooth_annualzip <- ZIP_selection(arrowtooth_annual)
summary(arrowtooth_annualzip[[2]])
plot(arrowtooth_annualzip[[2]])

windows()
par(mfrow = c(2,2))
gam.check(arrowtooth_annualzip[[2]])

# Test to ensure looks like a Poisson distribution
arrowtooth_annualcheck <- ZIP_test(arrowtooth_annual)
summary(arrowtooth_annualcheck)
plot(arrowtooth_annualcheck)

windows()
par(mfrow = c(2,2))
gam.check(arrowtooth_annualcheck)

range(arrowtooth_annual$count)
range(predict(arrowtooth_annualzip[[2]]))
range(predict(arrowtooth_annualcheck))

# Create manuscript figures
# Year variable
pdf("../results/ZIP/arrowtooth_flounder/year_annual.pdf",
    width = 12,
    height = 12)
arrowtooth_a_year <- plot_variable(arrowtooth_annualzip,
                              covariate = 1,
                              bounds = c(-1, 1),
                              "Year",
                              "Species Abundance Anomalies",
                              "s")
dev.off()
# Day of year variable
pdf("../results/ZIP/arrowtooth_flounder/julian_annual.pdf",
    width = 12,
    height = 12)
arrowtooth_a_day <- plot_variable(arrowtooth_annualzip,
                             covariate = 2,
                             bounds = c(-1, 1),
                             "Day of Year",
                             " ",
                             "n")
dev.off()
# Depth variable
pdf("../results/ZIP/arrowtooth_flounder/depth_annual.pdf",
    width = 12,
    height = 12)
arrowtooth_a_depth <- plot_variable(arrowtooth_annualzip,
                               covariate = 4,
                               bounds = c(-3.5, 2),
                               "Depth (m)",
                               "Species Abundance Anomalies",
                               "s")
dev.off()
# Temperature variable
pdf("../results/ZIP/arrowtooth_flounder/temp_annual.pdf",
    width = 12,
    height = 12)
arrowtooth_a_temp <- plot_variable(arrowtooth_annualzip,
                              covariate = 6,
                              bounds = c(-3.5, 2),
                              "Temperature (C)",
                              "",
                              "n")
dev.off()
# Grain size variable
pdf("../results/ZIP/arrowtooth_flounder/grain_annual.pdf",
    width = 12,
    height = 12)
arrowtooth_a_grain <- plot_variable(arrowtooth_annualzip,
                              covariate = 5,
                              bounds = c(-3.5, 2),
                              "Grain Size (phi)",
                              "",
                              "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/ZIP/arrowtooth_flounder/location_annual.pdf",
    width = 4,
    height = 11)
location_plot(arrowtooth_annualzip, arrowtooth_annual)
dev.off()


# Triennial
arrowtooth_triennialzip <- ZIP_selection(arrowtooth_triennial)
summary(arrowtooth_triennialzip[[2]])
plot(arrowtooth_triennialzip[[2]])

windows()
par(mfrow = c(2,2))
gam.check(arrowtooth_triennialzip[[2]])

# Test to ensure looks like a Poisson distribution
arrowtooth_triennialcheck <- ZIP_test(arrowtooth_triennial)
summary(arrowtooth_triennialcheck)
plot(arrowtooth_triennialcheck)

windows()
par(mfrow = c(2,2))
gam.check(arrowtooth_triennialcheck)

range(arrowtooth_triennial$count)
range(predict(arrowtooth_triennialzip[[2]]))
range(predict(arrowtooth_triennialcheck))

# Create manuscript figures
# Year variable
pdf("../results/ZIP/arrowtooth_flounder/year_triennial.pdf",
    width = 12,
    height = 12)
arrowtooth_t_year <- plot_variable(arrowtooth_triennialzip,
                              covariate = 1,
                              bounds = c(-3.2, 2),
                              "Year",
                              "Species Abundance Anomalies",
                              "s")
dev.off()
# Day of year variable
pdf("../results/ZIP/arrowtooth_flounder/julian_triennial.pdf",
    width = 12,
    height = 12)
arrowtooth_t_day <- plot_variable(arrowtooth_triennialzip,
                             covariate = 2,
                             bounds = c(-3.2, 2),
                             "Day of Year",
                             " ",
                             "n")
dev.off()
# Depth variable
pdf("../results/ZIP/arrowtooth_flounder/depth_triennial.pdf",
    width = 12,
    height = 12)
arrowtooth_t_depth <- plot_variable(arrowtooth_triennialzip,
                               covariate = 4,
                               bounds = c(-3.5, 2),
                               "Depth (m)",
                               "Species Abundance Anomalies",
                               "s")
dev.off()
# Temperature variable
pdf("../results/ZIP/arrowtooth_flounder/temp_triennial.pdf",
    width = 12,
    height = 12)
arrowtooth_t_temp <- plot_variable(arrowtooth_triennialzip,
                              covariate = 6,
                              bounds = c(-3.5, 2),
                              "Temperature (C)",
                              "",
                              "n")
dev.off()
# Grain size variable
pdf("../results/ZIP/arrowtooth_flounder/grain_triennial.pdf",
    width = 12,
    height = 12)
arrowtooth_t_grain <- plot_variable(arrowtooth_triennialzip,
                              covariate = 5,
                              bounds = c(-3.5, 2),
                              "Grain Size (phi)",
                              "",
                              "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/ZIP/arrowtooth_flounder/location_triennial.pdf",
    width = 4,
    height = 11)
location_plot(arrowtooth_triennialzip, arrowtooth_triennial)
dev.off()

# ****English sole ----
# Annual
english_annualzip <- ZIP_selection(english_annual)
summary(english_annualzip[[2]])
plot(english_annualzip[[2]])

windows()
par(mfrow = c(2,2))
gam.check(english_annualzip[[2]])

# Test to ensure looks like a Poisson distribution
english_annualcheck <- ZIP_test(english_annual)
summary(english_annualcheck)
plot(english_annualcheck)

windows()
par(mfrow = c(2,2))
gam.check(english_annualcheck) # not a normal distribution

range(english_annual$count)
range(predict(english_annualzip[[2]]))
range(predict(english_annualcheck))

# Create manuscript figures
# Year variable
pdf("../results/ZIP/english_sole/year_annual.pdf",
    width = 12,
    height = 12)
english_a_year <- plot_variable(english_annualzip,
                              covariate = 1,
                              bounds = c(-.5, .5),
                              "Year",
                              "Species Abundance Anomalies",
                              "s")
dev.off()
# Day of year variable
pdf("../results/ZIP/english_sole/julian_annual.pdf",
    width = 12,
    height = 12)
english_a_day <- plot_variable(english_annualzip,
                             covariate = 2,
                             bounds = c(-.5, .5),
                             "Day of Year",
                             " ",
                             "n")
dev.off()
# Depth variable
pdf("../results/ZIP/english_sole/depth_annual.pdf",
    width = 12,
    height = 12)
english_a_depth <- plot_variable(english_annualzip,
                               covariate = 4,
                               bounds = c(-1.5, 1),
                               "Depth (m)",
                               "Species Abundance Anomalies",
                               "s")
dev.off()
# Temperature variable
pdf("../results/ZIP/english_sole/temp_annual.pdf",
    width = 12,
    height = 12)
english_a_temp <- plot_variable(english_annualzip,
                              covariate = 6,
                              bounds = c(-1.5, 1),
                              "Temperature (C)",
                              "",
                              "n")
dev.off()
# Grain size variable
pdf("../results/ZIP/english_sole/grain_annual.pdf",
    width = 12,
    height = 12)
english_a_grain <- plot_variable(english_annualzip,
                              covariate = 5,
                              bounds = c(-1.5, 1),
                              "Grain Size (phi)",
                              "",
                              "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/ZIP/english_sole/location_annual.pdf",
    width = 4,
    height = 11)
location_plot(english_annualzip, english_annual)
dev.off()


# Triennial
english_triennialzip <- ZIP_selection(english_triennial)
summary(english_triennialzip[[2]])
plot(english_triennialzip[[2]])

windows()
par(mfrow = c(2,2))
gam.check(english_triennialzip[[2]])

# Test to ensure looks like a Poisson distribution
english_triennialcheck <- ZIP_test(english_triennial)
summary(english_triennialcheck)
plot(english_triennialcheck)

windows()
par(mfrow = c(2,2))
gam.check(english_triennialcheck)

range(english_triennial$count)
range(predict(english_triennialzip[[2]]))
range(predict(english_triennialcheck))

# Create manuscript figures
# Year variable
pdf("../results/ZIP/english_sole/year_triennial.pdf",
    width = 12,
    height = 12)
english_t_year <- plot_variable(english_triennialzip,
                              covariate = 1,
                              bounds = c(-6, 1),
                              "Year",
                              "Species Abundance Anomalies",
                              "s")
dev.off()
# Day of year variable
pdf("../results/ZIP/english_sole/julian_triennial.pdf",
    width = 12,
    height = 12)
english_t_day <- plot_variable(english_triennialzip,
                             covariate = 2,
                             bounds = c(-6, 1),
                             "Day of Year",
                             " ",
                             "n")
dev.off()
# Depth variable
pdf("../results/ZIP/english_sole/depth_triennial.pdf",
    width = 12,
    height = 12)
english_t_depth <- plot_variable(english_triennialzip,
                               covariate = 4,
                               bounds = c(-3, 1),
                               "Depth (m)",
                               "Species Abundance Anomalies",
                               "s")
dev.off()
# Temperature variable
pdf("../results/ZIP/english_sole/temp_triennial.pdf",
    width = 12,
    height = 12)
english_t_temp <- plot_variable(english_triennialzip,
                              covariate = 6,
                              bounds = c(-3, 1),
                              "Temperature (C)",
                              "",
                              "n")
dev.off()
# Grain size variable
pdf("../results/ZIP/english_sole/grain_triennial.pdf",
    width = 12,
    height = 12)
english_t_grain <- plot_variable(english_triennialzip,
                              covariate = 5,
                              bounds = c(-3, 1),
                              "Grain Size (phi)",
                              "",
                              "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/ZIP/english_sole/location_triennial.pdf",
    width = 4,
    height = 11)
location_plot(english_triennialzip, english_triennial)
dev.off()


# ****Lingcod ----
# Annual
lingcod_annualzip <- ZIP_selection(lingcod_annual)
summary(lingcod_annualzip[[2]])
plot(lingcod_annualzip[[2]])

windows()
par(mfrow = c(2,2))
gam.check(lingcod_annualzip[[2]])

# Test to ensure looks like a Poisson distribution
lingcod_annualcheck <- ZIP_test(lingcod_annual)
summary(lingcod_annualcheck)
plot(lingcod_annualcheck)

windows()
par(mfrow = c(2,2))
gam.check(lingcod_annualcheck) # not a normal distribution

range(lingcod_annual$count)
range(predict(lingcod_annualzip[[2]]))
range(predict(lingcod_annualcheck))

# Create manuscript figures
# Year variable
pdf("../results/ZIP/lingcod/year_annual.pdf",
    width = 12,
    height = 12)
lingcod_a_year <- plot_variable(lingcod_annualzip,
                              covariate = 1,
                              bounds = c(-1.2, 1.5),
                              "Year",
                              "Species Abundance Anomalies",
                              "s")
dev.off()
# Day of year variable
pdf("../results/ZIP/lingcod/julian_annual.pdf",
    width = 12,
    height = 12)
lingcod_a_day <- plot_variable(lingcod_annualzip,
                             covariate = 2,
                             bounds = c(-1.2, 1.5),
                             "Day of Year",
                             " ",
                             "n")
dev.off()
# Depth variable
pdf("../results/ZIP/lingcod/depth_annual.pdf",
    width = 12,
    height = 12)
lingcod_a_depth <- plot_variable(lingcod_annualzip,
                               covariate = 4,
                               bounds = c(-1.2, 2.5),
                               "Depth (m)",
                               "Species Abundance Anomalies",
                               "s")
dev.off()
# Temperature variable
pdf("../results/ZIP/lingcod/temp_annual.pdf",
    width = 12,
    height = 12)
lingcod_a_temp <- plot_variable(lingcod_annualzip,
                              covariate = 6,
                              bounds = c(-1.2, 2.5),
                              "Temperature (C)",
                              "",
                              "n")
dev.off()
# Grain size variable
pdf("../results/ZIP/lingcod/grain_annual.pdf",
    width = 12,
    height = 12)
lingcod_a_grain <- plot_variable(lingcod_annualzip,
                              covariate = 5,
                              bounds = c(-1.2, 2.5),
                              "Grain Size (phi)",
                              "",
                              "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/ZIP/lingcod/location_annual.pdf",
    width = 4,
    height = 11)
location_plot(lingcod_annualzip, lingcod_annual)
dev.off()


# Triennial
lingcod_triennialzip <- ZIP_selection(lingcod_triennial)
summary(lingcod_triennialzip[[2]])
plot(lingcod_triennialzip[[2]])

windows()
par(mfrow = c(2,2))
gam.check(lingcod_triennialzip[[2]])

# Test to ensure looks like a Poisson distribution
lingcod_triennialcheck <- ZIP_test(lingcod_triennial)
summary(lingcod_triennialcheck)
plot(lingcod_triennialcheck)

windows()
par(mfrow = c(2,2))
gam.check(lingcod_triennialcheck)

range(lingcod_triennial$count)
range(predict(lingcod_triennialzip[[2]]))
range(predict(lingcod_triennialcheck))

# Create manuscript figures
# NPGO variable
pdf("../results/ZIP/lingcod/NPGO_triennial.pdf",
    width = 12,
    height = 12)
lingcod_t_NPGO <- plot_variable(lingcod_triennialzip,
                              covariate = 1,
                              bounds = c(-1.5, 3.1),
                              "NPGO",
                              " ",
                              "n")
dev.off()
# Day of year variable
pdf("../results/ZIP/lingcod/julian_triennial.pdf",
    width = 12,
    height = 12)
lingcod_t_day <- plot_variable(lingcod_triennialzip,
                             covariate = 2,
                             bounds = c(-1.2, 1.5),
                             "Day of Year",
                             "Species Abundance Anomalies",
                             "s")
dev.off()
# Depth variable
pdf("../results/ZIP/lingcod/depth_triennial.pdf",
    width = 12,
    height = 12)
lingcod_t_depth <- plot_variable(lingcod_triennialzip,
                               covariate = 4,
                               bounds = c(-1.5, 3.1),
                               "Depth (m)",
                               "Species Abundance Anomalies",
                               "s")
dev.off()
# Temperature variable
pdf("../results/ZIP/lingcod/temp_triennial.pdf",
    width = 12,
    height = 12)
lingcod_t_temp <- plot_variable(lingcod_triennialzip,
                              covariate = 6,
                              bounds = c(-1.5, 3.1),
                              "Temperature (C)",
                              "",
                              "n")
dev.off()
# Grain size variable
pdf("../results/ZIP/lingcod/grain_triennial.pdf",
    width = 12,
    height = 12)
lingcod_t_grain <- plot_variable(lingcod_triennialzip,
                              covariate = 5,
                              bounds = c(-1.5, 3.1),
                              "Grain Size (phi)",
                              "",
                              "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/ZIP/lingcod/location_triennial.pdf",
    width = 4,
    height = 11)
location_plot(lingcod_triennialzip, lingcod_triennial)
dev.off()


# ****Pacific Sanddab ----
# Annual
sanddab_annualzip <- ZIP_selection(sanddab_annual)
summary(sanddab_annualzip[[2]])
plot(sanddab_annualzip[[2]])

windows()
par(mfrow = c(2,2))
gam.check(sanddab_annualzip[[2]])

# Test to ensure looks like a Poisson distribution
sanddab_annualcheck <- ZIP_test(sanddab_annual)
summary(sanddab_annualcheck)
plot(sanddab_annualcheck)

windows()
par(mfrow = c(2,2))
gam.check(sanddab_annualcheck) # not a normal distribution

range(sanddab_annual$count)
range(predict(sanddab_annualzip[[2]]))
range(predict(sanddab_annualcheck))

# Create manuscript figures
# NPGO variable
pdf("../results/ZIP/pacific_sanddab/NPGO_annual.pdf",
    width = 12,
    height = 12)
sanddab_a_NPGO <- plot_variable(sanddab_annualzip,
                              covariate = 1,
                              bounds = c(-5.6, 1.2),
                              "NPGO",
                              " ",
                              "n")
dev.off()
# Day of year variable
pdf("../results/ZIP/pacific_sanddab/julian_annual.pdf",
    width = 12,
    height = 12)
sanddab_a_day <- plot_variable(sanddab_annualzip,
                             covariate = 2,
                             bounds = c(-0.5, 1),
                             "Day of Year",
                             "Species Abundance Anomalies",
                             "s")
dev.off()
# Depth variable
pdf("../results/ZIP/pacific_sanddab/depth_annual.pdf",
    width = 12,
    height = 12)
sanddab_a_depth <- plot_variable(sanddab_annualzip,
                               covariate = 4,
                               bounds = c(-5.6, 1.2),
                               "Depth (m)",
                               "Species Abundance Anomalies",
                               "s")
dev.off()
# Temperature variable
pdf("../results/ZIP/pacific_sanddab/temp_annual.pdf",
    width = 12,
    height = 12)
sanddab_a_temp <- plot_variable(sanddab_annualzip,
                              covariate = 6,
                              bounds = c(-5.6, 1.2),
                              "Temperature (C)",
                              "",
                              "n")
dev.off()
# Grain size variable
pdf("../results/ZIP/pacific_sanddab/grain_annual.pdf",
    width = 12,
    height = 12)
sanddab_a_grain <- plot_variable(sanddab_annualzip,
                              covariate = 5,
                              bounds = c(-5.6, 1.2),
                              "Grain Size (phi)",
                              "",
                              "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/ZIP/pacific_sanddab/location_annual.pdf",
    width = 4,
    height = 11)
location_plot(sanddab_annualzip, sanddab_annual)
dev.off()


# Triennial
sanddab_triennialzip <- ZIP_selection(sanddab_triennial)
summary(sanddab_triennialzip[[2]])
plot(sanddab_triennialzip[[2]])

windows()
par(mfrow = c(2,2))
gam.check(sanddab_triennialzip[[2]])

# Test to ensure looks like a Poisson distribution
sanddab_triennialcheck <- ZIP_test(sanddab_triennial)
summary(sanddab_triennialcheck)
plot(sanddab_triennialcheck)

windows()
par(mfrow = c(2,2))
gam.check(sanddab_triennialcheck)

range(sanddab_triennial$count)
range(predict(sanddab_triennialzip[[2]]))
range(predict(sanddab_triennialcheck))

# Create manuscript figures
# Year variable
pdf("../results/ZIP/pacific_sanddab/year_triennial.pdf",
    width = 12,
    height = 12)
sanddab_t_year <- plot_variable(sanddab_triennialzip,
                              covariate = 1,
                              bounds = c(-4, 1),
                              "Year",
                              "Species Abundance Anomalies",
                              "s")
dev.off()
# Day of year variable
pdf("../results/ZIP/pacific_sanddab/julian_triennial.pdf",
    width = 12,
    height = 12)
sanddab_t_day <- plot_variable(sanddab_triennialzip,
                             covariate = 2,
                             bounds = c(-4, 1),
                             "Day of Year",
                             " ",
                             "n")
dev.off()
# Depth variable
pdf("../results/ZIP/pacific_sanddab/depth_triennial.pdf",
    width = 12,
    height = 12)
sanddab_t_depth <- plot_variable(sanddab_triennialzip,
                               covariate = 4,
                               bounds = c(-6.7, 2),
                               "Depth (m)",
                               "Species Abundance Anomalies",
                               "s")
dev.off()
# Temperature variable
pdf("../results/ZIP/pacific_sanddab/temp_triennial.pdf",
    width = 12,
    height = 12)
sanddab_t_temp <- plot_variable(sanddab_triennialzip,
                              covariate = 6,
                              bounds = c(-6.7, 2),
                              "Temperature (C)",
                              "",
                              "n")
dev.off()
# Grain size variable
pdf("../results/ZIP/pacific_sanddab/grain_triennial.pdf",
    width = 12,
    height = 12)
sanddab_t_grain <- plot_variable(sanddab_triennialzip,
                              covariate = 5,
                              bounds = c(-6.7, 2),
                              "Grain Size (phi)",
                              "",
                              "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/ZIP/pacific_sanddab/location_triennial.pdf",
    width = 4,
    height = 11)
location_plot(sanddab_triennialzip, sanddab_triennial)
dev.off()


# ****Petrale sole ----
# Annual
petrale_annualzip <- ZIP_selection(petrale_annual)
summary(petrale_annualzip[[2]])
plot(petrale_annualzip[[2]])

windows()
par(mfrow = c(2,2))
gam.check(petrale_annualzip[[2]])

# Test to ensure looks like a Poisson distribution
petrale_annualcheck <- ZIP_test(petrale_annual)
summary(petrale_annualcheck)
plot(petrale_annualcheck)

windows()
par(mfrow = c(2,2))
gam.check(petrale_annualcheck) # not a normal distribution

range(petrale_annual$count)
range(predict(petrale_annualzip[[2]]))
range(predict(petrale_annualcheck))

# Create manuscript figures
# Year variable
pdf("../results/ZIP/petrale_sole/year_annual.pdf",
    width = 12,
    height = 12)
petrale_a_year <- plot_variable(petrale_annualzip,
                              covariate = 1,
                              bounds = c(-1, 1),
                              "Year",
                              "Species Abundance Anomalies",
                              "s")
dev.off()
# Day of year variable
pdf("../results/ZIP/petrale_sole/julian_annual.pdf",
    width = 12,
    height = 12)
petrale_a_day <- plot_variable(petrale_annualzip,
                             covariate = 2,
                             bounds = c(-1, 1),
                             "Day of Year",
                             " ",
                             "n")
dev.off()
# Depth variable
pdf("../results/ZIP/petrale_sole/depth_annual.pdf",
    width = 12,
    height = 12)
petrale_a_depth <- plot_variable(petrale_annualzip,
                               covariate = 4,
                               bounds = c(-1.5, 0.5),
                               "Depth (m)",
                               "Species Abundance Anomalies",
                               "s")
dev.off()
# Temperature variable
pdf("../results/ZIP/petrale_sole/temp_annual.pdf",
    width = 12,
    height = 12)
petrale_a_temp <- plot_variable(petrale_annualzip,
                              covariate = 6,
                              bounds = c(-1.5, 0.5),
                              "Temperature (C)",
                              "",
                              "n")
dev.off()
# Grain size variable
pdf("../results/ZIP/petrale_sole/grain_annual.pdf",
    width = 12,
    height = 12)
petrale_a_grain <- plot_variable(petrale_annualzip,
                              covariate = 5,
                              bounds = c(-1.5, 0.5),
                              "Grain Size (phi)",
                              "",
                              "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/ZIP/petrale_sole/location_annual.pdf",
    width = 4,
    height = 11)
location_plot(petrale_annualzip, petrale_annual)
dev.off()


# Triennial
petrale_triennialzip <- ZIP_selection(petrale_triennial)
summary(petrale_triennialzip[[2]])
plot(petrale_triennialzip[[2]])

windows()
par(mfrow = c(2,2))
gam.check(petrale_triennialzip[[2]])

# Test to ensure looks like a Poisson distribution
petrale_triennialcheck <- ZIP_test(petrale_triennial)
summary(petrale_triennialcheck)
plot(petrale_triennialcheck)

windows()
par(mfrow = c(2,2))
gam.check(petrale_triennialcheck)

range(petrale_triennial$count)
range(predict(petrale_triennialzip[[2]]))
range(predict(petrale_triennialcheck))

# Create manuscript figures
# Year variable
pdf("../results/ZIP/petrale_sole/year_triennial.pdf",
    width = 12,
    height = 12)
petrale_t_year <- plot_variable(petrale_triennialzip,
                              covariate = 1,
                              bounds = c(-1.5, 1),
                              "Year",
                              "Species Abundance Anomalies",
                              "s")
dev.off()
# Day of year variable
pdf("../results/ZIP/petrale_sole/julian_triennial.pdf",
    width = 12,
    height = 12)
petrale_t_day <- plot_variable(petrale_triennialzip,
                             covariate = 2,
                             bounds = c(-1.5, 1),
                             "Day of Year",
                             " ",
                             "n")
dev.off()
# Depth variable
pdf("../results/ZIP/petrale_sole/depth_triennial.pdf",
    width = 12,
    height = 12)
petrale_t_depth <- plot_variable(petrale_triennialzip,
                               covariate = 4,
                               bounds = c(-1, 1.5),
                               "Depth (m)",
                               "Species Abundance Anomalies",
                               "s")
dev.off()
# Temperature variable
pdf("../results/ZIP/petrale_sole/temp_triennial.pdf",
    width = 12,
    height = 12)
petrale_t_temp <- plot_variable(petrale_triennialzip,
                              covariate = 6,
                              bounds = c(-1, 1.5),
                              "Temperature (C)",
                              "",
                              "n")
dev.off()
# Grain size variable
pdf("../results/ZIP/petrale_sole/grain_triennial.pdf",
    width = 12,
    height = 12)
petrale_t_grain <- plot_variable(petrale_triennialzip,
                              covariate = 5,
                              bounds = c(-1, 1.5),
                              "Grain Size (phi)",
                              "",
                              "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/ZIP/petrale_sole/location_triennial.pdf",
    width = 4,
    height = 11)
location_plot(petrale_triennialzip, petrale_triennial)
dev.off()


# ****Rex sole ----
# Annual
rex_annualzip <- ZIP_selection(rex_annual)
summary(rex_annualzip[[2]])
plot(rex_annualzip[[2]])

windows()
par(mfrow = c(2,2))
gam.check(rex_annualzip[[2]])

# Test to ensure looks like a Poisson distribution
rex_annualcheck <- ZIP_test(rex_annual)
summary(rex_annualcheck)
plot(rex_annualcheck)

windows()
par(mfrow = c(2,2))
gam.check(rex_annualcheck) # not a normal distribution

range(rex_annual$count)
range(predict(rex_annualzip[[2]]))
range(predict(rex_annualcheck))

# Create manuscript figures
# Year variable
pdf("../results/ZIP/rex_sole/year_annual.pdf",
    width = 12,
    height = 12)
rex_a_year <- plot_variable(rex_annualzip,
                              covariate = 1,
                              bounds = c(-1.2, 0.7),
                              "Year",
                              "Species Abundance Anomalies",
                              "s")
dev.off()
# Day of year variable
pdf("../results/ZIP/rex_sole/julian_annual.pdf",
    width = 12,
    height = 12)
rex_a_day <- plot_variable(rex_annualzip,
                             covariate = 2,
                             bounds = c(-1.2, 0.7),
                             "Day of Year",
                             " ",
                             "n")
dev.off()
# Depth variable
pdf("../results/ZIP/rex_sole/depth_annual.pdf",
    width = 12,
    height = 12)
rex_a_depth <- plot_variable(rex_annualzip,
                               covariate = 4,
                               bounds = c(-1, 1),
                               "Depth (m)",
                               "Species Abundance Anomalies",
                               "s")
dev.off()
# Temperature variable
pdf("../results/ZIP/rex_sole/temp_annual.pdf",
    width = 12,
    height = 12)
rex_a_temp <- plot_variable(rex_annualzip,
                              covariate = 6,
                              bounds = c(-1, 1),
                              "Temperature (C)",
                              "",
                              "n")
dev.off()
# Grain size variable
pdf("../results/ZIP/rex_sole/grain_annual.pdf",
    width = 12,
    height = 12)
rex_a_grain <- plot_variable(rex_annualzip,
                              covariate = 5,
                              bounds = c(-1, 1),
                              "Grain Size (phi)",
                              "",
                              "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/ZIP/rex_sole/location_annual.pdf",
    width = 4,
    height = 11)
location_plot(rex_annualzip, rex_annual)
dev.off()


# Triennial
rex_triennialzip <- ZIP_selection(rex_triennial)
summary(rex_triennialzip[[2]])
plot(rex_triennialzip[[2]])

windows()
par(mfrow = c(2,2))
gam.check(rex_triennialzip[[2]])

# Test to ensure looks like a Poisson distribution
rex_triennialcheck <- ZIP_test(rex_triennial)
summary(rex_triennialcheck)
plot(rex_triennialcheck)

windows()
par(mfrow = c(2,2))
gam.check(rex_triennialcheck)

range(rex_triennial$count)
range(predict(rex_triennialzip[[2]]))
range(predict(rex_triennialcheck))

# Create manuscript figures
# Year variable
pdf("../results/ZIP/rex_sole/year_triennial.pdf",
    width = 12,
    height = 12)
rex_t_year <- plot_variable(rex_triennialzip,
                              covariate = 1,
                              bounds = c(-4, 1.2),
                              "Year",
                              "Species Abundance Anomalies",
                              "s")
dev.off()
# Day of year variable
pdf("../results/ZIP/rex_sole/julian_triennial.pdf",
    width = 12,
    height = 12)
rex_t_day <- plot_variable(rex_triennialzip,
                             covariate = 2,
                             bounds = c(-4, 1.2),
                             "Day of Year",
                             " ",
                             "n")
dev.off()
# Depth variable
pdf("../results/ZIP/rex_sole/depth_triennial.pdf",
    width = 12,
    height = 12)
rex_t_depth <- plot_variable(rex_triennialzip,
                               covariate = 4,
                               bounds = c(-1.2, 1.5),
                               "Depth (m)",
                               "Species Abundance Anomalies",
                               "s")
dev.off()
# Temperature variable
pdf("../results/ZIP/rex_sole/temp_triennial.pdf",
    width = 12,
    height = 12)
rex_t_temp <- plot_variable(rex_triennialzip,
                              covariate = 6,
                              bounds = c(-1.2, 1.5),
                              "Temperature (C)",
                              "",
                              "n")
dev.off()
# Grain size variable
pdf("../results/ZIP/rex_sole/grain_triennial.pdf",
    width = 12,
    height = 12)
rex_t_grain <- plot_variable(rex_triennialzip,
                              covariate = 5,
                              bounds = c(-1.2, 1.5),
                              "Grain Size (phi)",
                              "",
                              "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/ZIP/rex_sole/location_triennial.pdf",
    width = 4,
    height = 11)
location_plot(rex_triennialzip, rex_triennial)
dev.off()

# ****Sablefish ----
# Annual
sablefish_annualzip <- ZIP_selection(sablefish_annual)
summary(sablefish_annualzip[[2]])
plot(sablefish_annualzip[[2]])

windows()
par(mfrow = c(2,2))
gam.check(sablefish_annualzip[[2]])

# Test to ensure looks like a Poisson distribution
sablefish_annualcheck <- ZIP_test(sablefish_annual)
summary(sablefish_annualcheck)
plot(sablefish_annualcheck)

windows()
par(mfrow = c(2,2))
gam.check(sablefish_annualcheck) # not a normal distribution

range(sablefish_annual$count)
range(predict(sablefish_annualzip[[2]]))
range(predict(sablefish_annualcheck))

# Create manuscript figures
# Year variable
pdf("../results/ZIP/sablefish/year_annual.pdf",
    width = 12,
    height = 12)
sablefish_a_year <- plot_variable(sablefish_annualzip,
                              covariate = 1,
                              bounds = c(-1, 1.5),
                              "Year",
                              "Species Abundance Anomalies",
                              "s")
dev.off()
# Day of year variable
pdf("../results/ZIP/sablefish/julian_annual.pdf",
    width = 12,
    height = 12)
sablefish_a_day <- plot_variable(sablefish_annualzip,
                             covariate = 2,
                             bounds = c(-1, 1.5),
                             "Day of Year",
                             " ",
                             "n")
dev.off()
# Depth variable
pdf("../results/ZIP/sablefish/depth_annual.pdf",
    width = 12,
    height = 12)
sablefish_a_depth <- plot_variable(sablefish_annualzip,
                               covariate = 4,
                               bounds = c(-6, 1),
                               "Depth (m)",
                               "Species Abundance Anomalies",
                               "s")
dev.off()
# Temperature variable
pdf("../results/ZIP/sablefish/temp_annual.pdf",
    width = 12,
    height = 12)
sablefish_a_temp <- plot_variable(sablefish_annualzip,
                              covariate = 6,
                              bounds = c(-6, 1),
                              "Temperature (C)",
                              "",
                              "n")
dev.off()
# Grain size variable
pdf("../results/ZIP/sablefish/grain_annual.pdf",
    width = 12,
    height = 12)
sablefish_a_grain <- plot_variable(sablefish_annualzip,
                              covariate = 5,
                              bounds = c(-6, 1),
                              "Grain Size (phi)",
                              "",
                              "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/ZIP/sablefish/location_annual.pdf",
    width = 4,
    height = 11)
location_plot(sablefish_annualzip, sablefish_annual)
dev.off()


# Triennial
sablefish_triennialzip <- ZIP_selection(sablefish_triennial)
summary(sablefish_triennialzip[[2]])
plot(sablefish_triennialzip[[2]])

windows()
par(mfrow = c(2,2))
gam.check(sablefish_triennialzip[[2]])

# Test to ensure looks like a Poisson distribution
sablefish_triennialcheck <- ZIP_test(sablefish_triennial)
summary(sablefish_triennialcheck)
plot(sablefish_triennialcheck)

windows()
par(mfrow = c(2,2))
gam.check(sablefish_triennialcheck)

range(sablefish_triennial$count)
range(predict(sablefish_triennialzip[[2]]))
range(predict(sablefish_triennialcheck))

# Create manuscript figures
# NPGO variable
pdf("../results/ZIP/sablefish/NPGO_triennial.pdf",
    width = 12,
    height = 12)
sablefish_t_NPGO <- plot_variable(sablefish_triennialzip,
                              covariate = 1,
                              bounds = c(-25, 2.3),
                              "NPGO",
                              " ",
                              "n")
dev.off()
# Day of year variable
pdf("../results/ZIP/sablefish/julian_triennial.pdf",
    width = 12,
    height = 12)
sablefish_t_day <- plot_variable(sablefish_triennialzip,
                             covariate = 2,
                             bounds = c(-1, 2),
                             "Day of Year",
                             "Species Abundance Anomalies",
                             "s")
dev.off()
# Depth variable
pdf("../results/ZIP/sablefish/depth_triennial.pdf",
    width = 12,
    height = 12)
sablefish_t_depth <- plot_variable(sablefish_triennialzip,
                               covariate = 4,
                               bounds = c(-25, 2.3),
                               "Depth (m)",
                               "Species Abundance Anomalies",
                               "s")
dev.off()
# Temperature variable
pdf("../results/ZIP/sablefish/temp_triennial.pdf",
    width = 12,
    height = 12)
sablefish_t_temp <- plot_variable(sablefish_triennialzip,
                              covariate = 6,
                              bounds = c(-25, 2.3),
                              "Temperature (C)",
                              "",
                              "n")
dev.off()
# Grain size variable
pdf("../results/ZIP/sablefish/grain_triennial.pdf",
    width = 12,
    height = 12)
sablefish_t_grain <- plot_variable(sablefish_triennialzip,
                              covariate = 5,
                              bounds = c(-25, 2.3),
                              "Grain Size (phi)",
                              "",
                              "n")
dev.off()
# Latitude/Longitude Map
pdf("../results/ZIP/sablefish/location_triennial.pdf",
    width = 4,
    height = 11)
location_plot(sablefish_triennialzip, sablefish_triennial)
dev.off()
