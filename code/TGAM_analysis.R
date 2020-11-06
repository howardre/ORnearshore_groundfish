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

# Validate the results and map the final product ----
# Calculate distance of each grid point to closest 'positive observation'
# Gives an index, then plot the true and false to get areas of sigificant increase or decrease
subset_distances <- function(tgam, df) {
  nlat = 40 # determine the resolution of the grid
  nlon = 30
  latd = seq(41, 48, length.out = nlat)
  lond = seq(-125,-123.9, length.out = nlon)
  subset_before <- expand.grid(latd, lond)
  names(subset_before) <- c("latitude", "longitude")
  subset_before$year <- 2003 # change depending on species
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
    before_prediction <- predict(tgam[[2]],
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
    return(list(significant_high, significant_low, before_prediction, after_prediction))
  }

arrowtooth_dist <- subset_distances(arrowtooth_tgam, arrowtooth_subset)
arrowtooth_CI <- tgam_prediction(arrowtooth_tgam, arrowtooth_dist)
english_dist <- subset_distances(english_tgam, english_subset)
english_CI <- tgam_prediction(english_tgam, english_subset)
sanddab_dist <- subset_distances(sanddab_tgam, sanddab_subset)
sanddab_CI <- tgam_prediction(sanddab_tgam, sanddab_subset)
dover_dist <- subset_distances(dover_tgam, dover_subset)
dover_CI <- tgam_prediction(dover_tgam, dover_subset)
rex_dist <- subset_distances(rex_tgam, rex_subset)
rex_CI <- tgam_prediction(rex_tgam, rex_subset)
lingcod_dist <- subset_distances(lingcod_tgam, lingcod_subset)
lingcod_CI <- tgam_prediction(lingcod_tgam, lingcod_subset)
petrale_dist <- subset_distances(petrale_tgam, petrale_subset)
petrale_CI <- tgam_prediction(petrale_tgam, petrale_subset)
sablefish_dist <- subset_distances(sablefish_tgam, sablefish_subset)
sablefish_CI <- tgam_prediction(sablefish_tgam, sablefish_subset)

# Can check the results, will likely get NA
sum(1 * arrowtooth_CI[[2]])
sum(1 * arrowtooth_CI[[1]])

# Make maps of the significant increases and decreases
# Polygons drawn manually to calculate differences in the real data based on the predictions
validation_map <- function(species_subset, tgam, conf_int, species_dist, bathy.dat, bathy.mat){
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
            xlim = range(species_subset$longitude, na.rm = TRUE) + c(-0.5, 1),
            ylim = range(species_subset$latitude, na.rm = TRUE) + c(-0.5, 0.5),
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
  points(species_dist$longitude[significant_low],
         species_dist$latitude[significant_low],
         pch = 8,
         cex = 0.5)
  points(species_dist$longitude[significant_high],
         species_dist$latitude[significant_high],
         pch = 16,
         cex = 0.7)
  maps::map('worldHires',
            add = T,
            col = 'peachpuff3',
            fill = T)
}

# Arrowtooth Flounder ----
# ***Validate the results ----
windows(width = 7, height = 15)
validation_map(arrowtooth_subset, arrowtooth_tgam, arrowtooth_CI, arrowtooth_dist, bathy.dat, bathy.mat)
# savePol = locator(40, type = "o")
# arrowtooth_poly_n1 = data.frame(x = savePol$x, y = savePol$y) # northern portion of decrease part 1
# arrowtooth_poly_n2 = data.frame(x = savePol$x, y = savePol$y) # northern portion of decrease part 2
# arrowtooth_poly_c = data.frame(x = savePol$x, y = savePol$y) # central
# arrowtooth_poly_s = data.frame(x = savePol$x, y = savePol$y) # southern
# arrowtooth_poly_n <- rbind(arrowtooth_poly_n1, arrowtooth_poly_n2)

# Can use this index to check if there are an appropriate number of data points in a polygon
data_check <- function(species_subset, polygon){
  in.chull(species_subset$longitude,
           species_subset$latitude,
           polygon$x,
           polygon$y)
}

# Subset data to only pick up those that are inside the polygon for real data
polygon_subset <- function(species_subset, polygon){
  polygon_differences <- species_subset[in.chull(species_subset$longitude,
                                                 species_subset$latitude,
                                                 polygon$x,
                                                 polygon$y), ]
  return(polygon_differences)
}

arrowtooth_n_decrease1 <- polygon_subset(arrowtooth_subset, arrowtooth_poly_n1)
arrowtooth_n_decrease2 <- polygon_subset(arrowtooth_subset, arrowtooth_poly_n2)
arrowtooth_c_decrease <- polygon_subset(arrowtooth_subset, arrowtooth_poly_c)
arrowtooth_s_decrease <- polygon_subset(arrowtooth_subset, arrowtooth_poly_s)
arrowtooth_n_decrease <- rbind(arrowtooth_n_decrease1, arrowtooth_n_decrease2)

# Add the subset of data to the map to check if polygon is in right spot
polygon_map_check <- function(species_subset, difference, polygon){
  plot(species_subset$longitude,
       species_subset$latitude,
       pch = 16,
       ylim = range(species_subset$latitude),
       xlim = range(species_subset$longitude))
  points(difference$longitude,
         difference$latitude,
         col = 'red')
  polygon(polygon$x, polygon$y)
  maps::map('worldHires',
            add = T,
            col = 'peachpuff3',
            fill = T)
}

windows(width = 7, height = 15)
polygon_map_check(arrowtooth_subset, arrowtooth_n_decrease1, arrowtooth_poly_n1)

windows(width = 7, height = 15)
polygon_map_check(arrowtooth_subset, arrowtooth_n_decrease2, arrowtooth_poly_n2)

windows(width = 7, height = 15)
polygon_map_check(arrowtooth_subset, arrowtooth_c_decrease, arrowtooth_poly_c)

windows(width = 7, height = 15)
polygon_map_check(arrowtooth_subset, arrowtooth_s_decrease, arrowtooth_poly_s)

# Plot prediction grid
prediction_map <- function(species_dist, species_CI){
  nlat = 40
  nlon = 30
  latd = seq(41, 48, length.out = nlat)
  lond = seq(-125,-123.9, length.out = nlon)
  image.plot(
    lond, latd,
    t(matrix(species_dist$diff,
             nrow = length(latd),
             ncol = length(lond),
             byrow = T)),
    col = tim.colors(100),
    ylab = "",
    xlab = "",
    xlim = range(species_dist$longitude),
    ylim = range(species_dist$latitude),
    main = 'Predictions',
    cex.main = 1.5,
    cex.lab = 1.4,
    cex.axis = 1.4)
  points(species_CI$longitude[significant_low],
         species_CI$latitude[significant_low],
         pch = 8,
         cex = 0.5)
  points(species_CI$longitude[significant_high],
         species_CI$latitude[significant_high],
         pch = 16,
         cex = 0.7)
  maps::map('worldHires',
            add = T,
            col = 'grey',
            fill = T)
}
windows(width = 7, height = 15)
prediction_map(arrowtooth_dist, arrowtooth_CI)

# Subset data to only pick up those that are inside the polygon for real data
polygon_pred <- function(species_dist, polygon){
  species_distr[in.chull(
    species_dist$longitude,
    species_dist$latitude,
    polygon$x,
    polygon$y), ]
}

arrowtooth_north <- polygon_subset(arrowtooth_dist, arrowtooth_poly_n)
arrowtooth_central <- polygon_subset(arrowtooth_dist, arrowtooth_poly_c)
arrowtooth_south <- polygon_subset(arrowtooth_dist, arrowtooth_poly_s)

# Calculate difference before and after the threshold year in the real data
avg_pres_change <- function(difference){
  before <- difference[difference$thr == 'before', ]
  after <- difference[difference$thr == 'after',]
  average <- sum(after$pres == 1) / nrow(after) -
    (sum(before$pres == 1) / nrow(before))
  return(average)
}

arrowtooth_n_avg <- avg_pres_change(arrowtooth_n_decrease)
arrowtooth_c_avg <- avg_pres_change(arrowtooth_c_decrease)
arrowtooth_s_avg <- avg_pres_change(arrowtooth_s_decrease)

# Calculate difference before and after the threshold year for the predictions
arrowtooth_n_pred <- sum(arrowtooth_north$diff) / nrow(arrowtooth_north)
arrowtooth_c_pred <- sum(arrowtooth_central$diff) / nrow(arrowtooth_central)
arrowtooth_s_pred <- sum(arrowtooth_south$diff) / nrow(arrowtooth_south)

# ***Create the final maps ----
# Calculate difference before and after
arrowtooth_dist$mean_after <- arrowtooth_CI[[3]]$fit
arrowtooth_dist$mean_before <- arrowtooth_CI[[4]]$fit
arrowtooth_dist$diff <-  arrowtooth_dist$mean_after - arrowtooth_dist$mean_before
arrowtooth_dist$diff[is.na(arrowtooth_dist$diff)] <- 0
arrowtooth_dist$diff[arrowtooth_dist$dist > 10000] <- NA

tgam_map <- function(species_subset, tgam, threshold, title, bathy.dat, bathy.mat){
  myvis_gam(tgam[[2]],
            view = c('longitude', 'latitude'),
            too.far = 0.025,
            plot.type = 'contour',
            color = 'jet',
            type = 'response',
            cond = list(thr = threshold),
            main = title,
            contour.col = "gray35",
            xlim = range(species_subset$longitude, na.rm = TRUE) + c(-0.5, 1),
            ylim = range(species_subset$latitude, na.rm = TRUE) + c(-0.5, 0.5),
            cex.lab = 2,
            cex.axis = 2,
            cex.main = 2.4,
            xlab = "Longitude °W",
            ylab = "Latitude °N")
  contour(unique(bathy.dat$lon),
          sort(unique(bathy.dat$lat)),
          bathy.mat,
          levels = -seq(200, 200, by = 200),
          labcex = 0.8,
          add = T,
          col = 'black',
          labels = NULL,
          lwd = 1.95)
  symbols(species_subset$longitude[species_subset$pres > 0 &
                                     species_subset$thr == threshold],
          species_subset$latitude[species_subset$pres > 0 &
                                    species_subset$thr == threshold],
  circles = species_subset$pres[species_subset$pres > 0 &
                                  species_subset$thr == threshold],
  inches = 0.035,
  add = T,
  bg = alpha('gray', 0.4),
  fg = alpha('gray', 0.4))
  maps::map('worldHires',
            add = T,
            col = 'peachpuff3',
            fill = T)
  points(-124.0535,44.6368, pch=20)
  text(-124.0535,44.6368,"Newport", adj=c(0,1.2), cex = 2)
  points(-123.8313,46.1879, pch=20)
  text(-123.88028,46.13361,"Astoria", adj=c(0,1.2), cex = 2)
  points(-124.3,43.3, pch=20)
  text(-124.3,43.3,"Charleston", adj=c(0,1.2), cex = 2)
}
pred_map <- function(species_subset, species_dist, bathy.dat, bathy.mat){
  nlat = 40
  nlon = 30
  latd = seq(41, 48, length.out = nlat)
  lond = seq(-125,-123.9, length.out = nlon)
  image(lond, latd,
        t(matrix(species_dist$diff,
                 nrow = length(latd),
                 ncol = length(lond),
        byrow = F)),
        col = jet.colors(100),
        ylab = expression(paste("Latitude (" ^ 0, 'N)')),
        xlab = expression(paste("Longitude (" ^ 0, 'W)')),
        xlim = range(species_subset$longitude, na.rm = TRUE) + c(-0.5, 1),
        ylim = range(species_subset$latitude, na.rm = TRUE) + c(-0.5, 0.5),
        main = 'Predicted Change',
        cex.main = 2.4,
        cex.lab = 2,
        cex.axis = 2)
  contour(unique(bathy.dat$lon),
          sort(unique(bathy.dat$lat)),
          bathy.mat,
          levels = -seq(200, 200, by = 200),
          labcex = 0.8,
          add = T,
          col = 'black',
          labels = NULL,
          lwd = 1.95)
  points(species_dist$longitude[significant_low],
         species_dist$latitude[significant_low],
         pch = 2,
         cex = 0.9)
  points(species_dist$longitude[significant_high],
         species_dist$latitude[significant_high],
         pch = 16,
         cex = 0.9)
  maps::map('worldHires',
            add = T,
            col = 'peachpuff3',
            fill = T)
  points(-124.0535,44.6368, pch=20)
  text(-124.0535,44.6368,"Newport", adj=c(0,1.2), cex = 2)
  points(-123.8313,46.1879, pch=20)
  text(-123.88028,46.13361,"Astoria", adj=c(0,1.2), cex = 2)
  points(-124.3,43.3, pch=20)
  text(-124.3,43.3,"Charleston", adj=c(0,1.2), cex = 2)
  image.plot(legend.only = T,
             zlim = c(0, 1),
             col = jet.colors(100),
             legend.shrink = 0.2,
             smallplot = c(.16, .18, .06, .21),
             legend.cex = 1.5,
             legend.lab = "change",
             axis.args = list(cex.axis = 1.5),
             legend.width = 1,
             legend.line = 3)
  legend("bottomleft",
         legend = c("Decrease", "Increase"),
         pch = c(2, 16),
         bty = "n",
         pt.cex = 2,
         cex = 1.9,
         inset = c(0.03, 0.19))
  }

pdf("../results/TGAM/arrowtooth_flounder/arrowtooth_threshold.pdf",
    width = 17,
    height = 12)
par(
  mfrow = c(1, 3),
  family = 'serif',
  mar = c(4, 5, 3, 0.2) + 0.1,
  oma = c(1.5, 2, 4, 2))
jet.colors<-colorRampPalette(c("#F7FCFD", "#E0ECF4", "#BFD3E6", "#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D", "#6E016B"))
tgam_map(arrowtooth_subset, arrowtooth_tgam, threshold = "before", title = "Before 2007", bathy.dat, bathy.mat)
tgam_map(arrowtooth_subset, arrowtooth_tgam, threshold = "after", title = "After 2007", bathy.dat, bathy.mat)
jet.colors <- colorRampPalette(c("#F7FCF090", "#E0F3DB90", "#CCEBC590", "#A8DDB590", "#7BCCC490",
                                 "#4EB3D390", "#2B8CBE90", "#0868AC90",  "#08408190"), alpha = T) # Create new palette for predictions
pred_map(arrowtooth_subset, arrowtooth_dist, bathy.dat, bathy.mat)
mtext("Arrowtooth Flounder",
      line = 1,
      outer = T,
      cex = 2.5)
dev.off()

# English Sole ----
# ***Validate the results ----
windows(width = 7, height = 15)
validation_map(english_subset, english_tgam, english_CI, english_dist, bathy.dat, bathy.mat)
# savePol = locator(40, type = "o")
# english_poly_c = data.frame(x = savePol$x, y = savePol$y) # central decrease
# english_poly_s = data.frame(x = savePol$x, y = savePol$y) # southern increase

# Can use the data_check() function to see if there are appropriate number of data points in a polygon
# Subset data to only pick up those that are inside the polygon for real data
english_c_decrease <- polygon_subset(english_subset, english_poly_c)
english_s_increase <- polygon_subset(english_subset, english_poly_s)

# Add the subset of data to the map to check if polygon is in right spot
windows(width = 7, height = 15)
polygon_map_check(english_subset, english_c_decrease, english_poly_c)

windows(width = 7, height = 15)
polygon_map_check(english_subset, english_s_increase, english_poly_s)

# Plot prediction grid
windows(width = 7, height = 15)
prediction_map(english_dist, english_CI)

# Subset data to only pick up those that are inside the polygon for real data
english_central <- polygon_subset(english_dist, english_poly_c)
english_south <- polygon_subset(english_dist, english_poly_s)

# Calculate difference before and after the threshold year in the real data
english_c_avg <- avg_pres_change(english_c_decrease)
english_s_avg <- avg_pres_change(english_s_increase)

# Calculate difference before and after the threshold year for the predictions
english_c_pred <- sum(english_central$diff) / nrow(english_central)
english_s_pred <- sum(english_south$diff) / nrow(english_south)

# ***Create the final maps ----
# Calculate difference before and after
english_dist$mean_after <- english_CI[[3]]$fit
english_dist$mean_before <- english_CI[[4]]$fit
english_dist$diff <-  english_dist$mean_after - english_dist$mean_before
english_dist$diff[is.na(english_dist$diff)] <- 0
english_dist$diff[english_dist$dist > 10000] <- NA

pdf("../results/TGAM/english_sole/english_threshold.pdf",
    width = 17,
    height = 12)
par(
  mfrow = c(1, 3),
  family = 'serif',
  mar = c(4, 5, 3, 0.2) + 0.1,
  oma = c(1.5, 2, 4, 2))
jet.colors<-colorRampPalette(c("#F7FCFD", "#E0ECF4", "#BFD3E6", "#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D", "#6E016B"))
tgam_map(english_subset, english_tgam, threshold = "before", title = "Before 2007", bathy.dat, bathy.mat)
tgam_map(english_subset, english_tgam, threshold = "after", title = "After 2007", bathy.dat, bathy.mat)
jet.colors <- colorRampPalette(c("#F7FCF090", "#E0F3DB90", "#CCEBC590", "#A8DDB590", "#7BCCC490",
                                 "#4EB3D390", "#2B8CBE90", "#0868AC90",  "#08408190"), alpha = T) # Create new palette for predictions
pred_map(english_subset, english_dist, bathy.dat, bathy.mat)
mtext("English Sole",
      line = 1,
      outer = T,
      cex = 2.5)
dev.off()

# Pacific Sanddab ----
# ***Validate the results ----
windows(width = 7, height = 15)
validation_map(sanddab_subset, sanddab_tgam, sanddab_CI, sanddab_dist, bathy.dat, bathy.mat)
# savePol = locator(40, type = "o")
# sanddab_poly_n = data.frame(x = savePol$x, y = savePol$y) # northern decrease
# sanddab_poly_c = data.frame(x = savePol$x, y = savePol$y) # central decrease
# sanddab_poly_s = data.frame(x = savePol$x, y = savePol$y) # southern decrease

# Can use the data_check() function to see if there are appropriate number of data points in a polygon
# Subset data to only pick up those that are inside the polygon for real data
sanddab_n_decrease <- polygon_subset(sanddab_subset, sanddab_poly_n)
sanddab_c_decrease <- polygon_subset(sanddab_subset, sanddab_poly_c)
sanddab_s_decrease <- polygon_subset(sanddab_subset, sanddab_poly_s)

# Add the subset of data to the map to check if polygon is in right spot
windows(width = 7, height = 15)
polygon_map_check(sanddab_subset, sanddab_n_decrease, sanddab_poly_n)

windows(width = 7, height = 15)
polygon_map_check(sanddab_subset, sanddab_c_decrease, sanddab_poly_c)

windows(width = 7, height = 15)
polygon_map_check(sanddab_subset, sanddab_s_decrease, sanddab_poly_s)

# Plot prediction grid
windows(width = 7, height = 15)
prediction_map(sanddab_dist, sanddab_CI)

# Subset data to only pick up those that are inside the polygon for real data
sanddab_north <- polygon_subset(sanddab_dist, sanddab_poly_n)
sanddab_central <- polygon_subset(sanddab_dist, sanddab_poly_c)
sanddab_south <- polygon_subset(sanddab_dist, sanddab_poly_s)

# Calculate difference before and after the threshold year in the real data
sanddab_n_avg <- avg_pres_change(sanddab_n_decrease)
sanddab_c_avg <- avg_pres_change(sanddab_c_decrease)
sanddab_s_avg <- avg_pres_change(sanddab_s_decrease)

# Calculate difference before and after the threshold year for the predictions
sanddab_n_pred <- sum(sanddab_north$diff) / nrow(sanddab_north)
sanddab_c_pred <- sum(sanddab_central$diff) / nrow(sanddab_central)
sanddab_s_pred <- sum(sanddab_south$diff) / nrow(sanddab_south)

# ***Create the final maps ----
# Calculate difference before and after
sanddab_dist$mean_after <- sanddab_CI[[3]]$fit
sanddab_dist$mean_before <- sanddab_CI[[4]]$fit
sanddab_dist$diff <-  sanddab_dist$mean_after - sanddab_dist$mean_before
sanddab_dist$diff[is.na(sanddab_dist$diff)] <- 0
sanddab_dist$diff[sanddab_dist$dist > 10000] <- NA

pdf("../results/TGAM/pacific_sanddab/sanddab_threshold.pdf",
    width = 17,
    height = 12)
par(
  mfrow = c(1, 3),
  family = 'serif',
  mar = c(4, 5, 3, 0.2) + 0.1,
  oma = c(1.5, 2, 4, 2))
jet.colors<-colorRampPalette(c("#F7FCFD", "#E0ECF4", "#BFD3E6", "#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D", "#6E016B"))
tgam_map(sanddab_subset, sanddab_tgam, threshold = "before", title = "Before 2007", bathy.dat, bathy.mat)
tgam_map(sanddab_subset, sanddab_tgam, threshold = "after", title = "After 2007", bathy.dat, bathy.mat)
jet.colors <- colorRampPalette(c("#F7FCF090", "#E0F3DB90", "#CCEBC590", "#A8DDB590", "#7BCCC490",
                                 "#4EB3D390", "#2B8CBE90", "#0868AC90",  "#08408190"), alpha = T) # Create new palette for predictions
pred_map(sanddab_subset, sanddab_dist, bathy.dat, bathy.mat)
mtext("Pacific Sanddab",
      line = 1,
      outer = T,
      cex = 2.5)
dev.off()

# Dover Sole ----
# ***Validate the results ----
windows(width = 7, height = 15)
validation_map(dover_subset, dover_tgam, dover_CI, dover_dist, bathy.dat, bathy.mat)
# savePol = locator(40, type = "o")
# dover_poly_n1 = data.frame(x = savePol$x, y = savePol$y) # northern increase part 1
# dover_poly_n2 = data.frame(x = savePol$x, y = savePol$y) # northern increase part 2
# dover_poly_c1 = data.frame(x = savePol$x, y = savePol$y) # central increase
# dover_poly_c2 = data.frame(x = savePol$x, y = savePol$y) # central increase
# dover_poly_s = data.frame(x = savePol$x, y = savePol$y) # southern increase
# dover_poly_n <- rbind(dover_poly_n1, dover_poly_n2) # merge two northern polygons
# dover_poly_c <- rbind(dover_poly_c1, dover_poly_c2) # merge two northern polygons

# Can use the data_check() function to see if there are appropriate number of data points in a polygon
# Subset data to only pick up those that are inside the polygon for real data
dover_n_increase1 <- polygon_subset(dover_subset, dover_poly_n1)
dover_n_increase2 <- polygon_subset(dover_subset, dover_poly_n2)
dover_c_increase1 <- polygon_subset(dover_subset, dover_poly_c1)
dover_c_increase2 <- polygon_subset(dover_subset, dover_poly_c2)
dover_s_increase <- polygon_subset(dover_subset, dover_poly_s)
dover_n_increase <- rbind(dover_n_increase1, dover_n_increase2)
dover_c_increase <- rbind(dover_c_increase1, dover_c_increase2)

# Add the subset of data to the map to check if polygon is in right spot
windows(width = 7, height = 15)
polygon_map_check(dover_subset, dover_n_increase, dover_poly_n)

windows(width = 7, height = 15)
polygon_map_check(dover_subset, dover_c_increase, dover_poly_c)

windows(width = 7, height = 15)
polygon_map_check(dover_subset, dover_s_increase, dover_poly_s)

# Plot prediction grid
windows(width = 7, height = 15)
prediction_map(dover_dist, dover_CI)

# Subset data to only pick up those that are inside the polygon for real data
dover_north <- polygon_subset(dover_dist, dover_poly_n)
dover_central <- polygon_subset(dover_dist, dover_poly_c)
dover_south <- polygon_subset(dover_dist, dover_poly_s)

# Calculate difference before and after the threshold year in the real data
dover_n_avg <- avg_pres_change(dover_n_increase)
dover_c_avg <- avg_pres_change(dover_c_increase)
dover_s_avg <- avg_pres_change(dover_s_increase)

# Calculate difference before and after the threshold year for the predictions
dover_n_pred <- sum(dover_north$diff) / nrow(dover_north)
dover_c_pred <- sum(dover_central$diff) / nrow(dover_central)
dover_s_pred <- sum(dover_south$diff) / nrow(dover_south)

# ***Create the final maps ----
# Calculate difference before and after
dover_dist$mean_after <- dover_CI[[3]]$fit
dover_dist$mean_before <- dover_CI[[4]]$fit
dover_dist$diff <-  dover_dist$mean_after - dover_dist$mean_before
dover_dist$diff[is.na(dover_dist$diff)] <- 0
dover_dist$diff[dover_dist$dist > 10000] <- NA

pdf("../results/TGAM/dover_sole/dover_threshold.pdf",
    width = 17,
    height = 12)
par(
  mfrow = c(1, 3),
  family = 'serif',
  mar = c(4, 5, 3, 0.2) + 0.1,
  oma = c(1.5, 2, 4, 2))
jet.colors<-colorRampPalette(c("#F7FCFD", "#E0ECF4", "#BFD3E6", "#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D", "#6E016B"))
tgam_map(dover_subset, dover_tgam, threshold = "before", title = "Before 2007", bathy.dat, bathy.mat)
tgam_map(dover_subset, dover_tgam, threshold = "after", title = "After 2007", bathy.dat, bathy.mat)
jet.colors <- colorRampPalette(c("#F7FCF090", "#E0F3DB90", "#CCEBC590", "#A8DDB590", "#7BCCC490",
                                 "#4EB3D390", "#2B8CBE90", "#0868AC90",  "#08408190"), alpha = T) # Create new palette for predictions
pred_map(dover_subset, dover_dist, bathy.dat, bathy.mat)
mtext("Dover Sole",
      line = 1,
      outer = T,
      cex = 2.5)
dev.off()
