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

# Find thresholds for all 8 species
years <- sort(unique(arrowtooth_subset$year))[4:22]
arrowtooth_tgam <- get_tgam(arrowtooth_subset, years)
english_tgam <- get_tgam(english_subset, years)
sanddab_tgam <- get_tgam(sanddab_subset, years)
dover_tgam <- get_tgam_woyear(dover_subset, years)
rex_tgam <- get_tgam_woyear(rex_subset, years)
lingcod_tgam <- get_tgam(lingcod_subset, years)
petrale_tgam <- get_tgam(petrale_subset, years)
sablefish_tgam <- get_tgam(sablefish_subset, years)

# Plot AIC for all species
windows()
par(family = "serif")
plot(
  years,
  arrowtooth_tgam[[5]], # need y axis to be the AIC range for all GAMs
  type = 'b',
  xlab = 'Year',
  ylab = 'AIC',
  main = "Arrowtooth Flounder",
  cex.main = 1.4,
  cex.lab = 1.2,
  cex.axis = 1.2
)
abline(v = arrowtooth_tgam[[3]], lty = 2)
text(2003, 3465, 'Before', cex = 1.2)
text(2011, 3465, 'After', cex = 1.2)
abline(h = AIC(arrowtooth_tgam[[1]]), lty = 2)

plot_AIC <- function(tgam, years) {
  plot(years,
       tgam[[5]], # need y axis to be the AIC range for all GAMs
       type = 'b',
       xlab = 'Year',
       ylab = 'AIC',
       main = deparse(substitute(tgam)),
       cex.main = 1.4,
       cex.lab = 1.2,
       cex.axis = 1.2)
  abline(v = tgam[[3]], lty = 2)
  abline(h = AIC(tgam[[1]]), lty = 2)
}
