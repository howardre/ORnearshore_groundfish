# Libraries
library(maps)
library(mapdata)
library(fields)
library(marmap)
library(colorRamps)
library(itsadug)
library(RColorBrewer)
library(purrr)
library(mgcv)
library(furrr)

# Load data and functions
setwd('/Users/howar/Documents/Oregon State/ORnearshore_groundfish/code')
source("functions/vis_gam_COLORS.R")
source("functions/distance_function.R")
trawl_data <- read.delim("../data/NMFS_data/trawl_data.txt", header = T)
load("../data/bathy.dat")
load("../data/bathy.mat")
jet.colors <- colorRampPalette(c("#F7FCFD", "#E5F5F9", "#CCECE6", "#99D8C9",
                                 "#66C2A4", "#41AE76", "#238B45", "#006D2C", "#00441B"))
blues <- c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
greys <- c(grey(0.6), grey(0.93), grey(0.99))

# Create bathymetry dataset
OR_bathy <- getNOAA.bathy(lon1= -127, lon2= -121, lat1= 49, lat2= 39, resolution=1)

# Survey two panel figure ----
survey_map <- function(title, bathymetry, survey, years, letter){
plot.bathy(
  bathymetry,
  image = T,
  axes = T,
  lwd = 0.03,
  land = T,
  n = 0,
  bpal = list(c(0, max(bathymetry), greys), c(min(bathymetry), 0, blues)),
  ylim = c(40.3, 48.2),
  xlim = c(-125.5,-123.5),
  ylab = "Latitude °N",
  xlab = "Longitude °W",
  main = title,
  cex.lab = 2,
  cex.main = 2.5,
  cex.axis = 1.8)
  points(-124.0535, 44.6368, pch = 20)
  text(-124.0535,
       44.6368,
       "Newport",
       adj = c(0, 1.2),
       cex = 1.8)
  points(-123.8313, 46.1879, pch = 20)
  text(-123.88028,
       46.13361,
       "Astoria",
       adj = c(0, 1.2),
       cex = 1.8)
  points(-124.3, 43.3, pch = 20)
  text(-124.3,
       43.3,
       "Charleston",
       adj = c(0, 1.2),
       cex = 1.8)
  plot(bathymetry,
       deep = 0,
       shallow = 0,
       lwd = 1,
       add =  T)
  plot(bathymetry,
       deep = -50,
       shallow = -50,
       lwd = 0.4,
       drawlabels = T,
       add = T,
       col = "slategrey")
  plot(bathymetry,
       deep = -200,
       shallow = -200,
       lwd = 0.4,
       drawlabels = T,
       add = T,
       col =  "slategrey")
map.scale(-126., 40.3, cex = 1)
points(survey$longitude[survey$year == years],
       survey$latitude[survey$year == years],
       pch = 18,
       col = 'black',
       cex = 1.3)
title(outer = outer,
      adj = .025,
      main = letter,
      cex.main = 2,
      col = "black",
      font = 2,
      line = -2)
}

tiff("../final_figs/Figure_1.tiff", width = 11, height = 12, units = "in", res = 300)
par(mfrow = c(1, 2),
    family = 'serif',
    mar = c(4, 5, 3, .2) + .15)
survey_map("AFSC Triennial Survey", OR_bathy, trawl_data, 1995, "(A)")
survey_map("NWFSC Annual Survey", OR_bathy, trawl_data, 2018, "(B)")
dev.off()


intro_map <- function(bathymetry){

}

OR_bathy1 <- getNOAA.bathy(lon1= -127, lon2= -121, lat1= 48, lat2= 39, resolution=1)
tiff("../final_figs/manuscript2_fig_tables/Figure_1.tiff",
     width = 5.5,
     height = 12,
     units = "in",
     res = 300)
par(mfrow = c(1, 1),
    family = 'serif',
    mar = c(4, 5, 3, .2) + .15)
plot.bathy(
  OR_bathy1,
  image = T,
  axes = T,
  lwd = 0.03,
  land = T,
  n = 0,
  bpal = list(c(0,
                max(OR_bathy1),
                greys),
              c(min(OR_bathy1),
                0,
                blues)),
  ylim = c(40.3, 47.2),
  xlim = c(-125.5,-123.5),
  ylab = "Latitude °N",
  xlab = "Longitude °W",
  main = "",
  cex.lab = 2,
  cex.main = 2.5,
  cex.axis = 1.8)
points(-124.0535, 44.6368, pch = 20)
text(-124.0535,
     44.6368,
     "Newport",
     adj = c(0, 1.2),
     cex = 1.8)
points(-123.8313, 46.1879, pch = 20)
text(-123.88028,
     46.13361,
     "Astoria",
     adj = c(0, 1.2),
     cex = 1.8)
points(-124.3, 43.3, pch = 20)
text(-124.3,
     43.3,
     "Charleston",
     adj = c(0, 1.2),
     cex = 1.8)
plot(OR_bathy1,
     deep = 0,
     shallow = 0,
     lwd = 1,
     add =  T)
plot(OR_bathy1,
     deep = -50,
     shallow = -50,
     lwd = 0.7,
     drawlabels = T,
     add = T,
     col = "gray19")
plot(OR_bathy1,
     deep = -200,
     shallow = -200,
     lwd = 0.7,
     drawlabels = T,
     add = T,
     col =  "gray19")
map.scale(-126., 40.3, cex = 1)
dev.off()

###########################################################################################################################################################
# Make temperature TGAM ----
trawl_data$bottom_temp[is.na(trawl_data$bottom_temp)] <- 0
trawl_data <- trawl_data[trawl_data$bottom_temp < 10, ]
trawl_data <- trawl_data[trawl_data$bottom_temp > 4, ]

get_yr_aic <- function(df, yr) {
  df$thr <- ifelse(df$year <= yr, 'before', 'after')
  ctrl <- list(nthreads = 6)
  aic_year <-  gam(bottom_temp ~ factor(year) +
                     s(longitude, latitude, by = factor(thr)) +
                     s(julian),
                   control = ctrl,
                   data = df)$aic
  return_list <- data.frame(aic_year, yr)
}
get_tgam <- function(df, years) {
  ctrl <- list(nthreads = 6)
  ref_gam <- gam(bottom_temp ~ factor(year) +
                   s(julian) +
                   s(longitude, latitude),
                 control = ctrl,
                 data = df)
  all_tgams <- data.frame(aic_year = rep(0, length(years)),
                          year = rep(0, length(years)))
  all_tgams <- future_map(years, ~ get_yr_aic(df, .x))
  best_tgam <- all_tgams[[which.min(sapply(1:length(all_tgams),
                                           function(x) min(all_tgams[[x]]$aic_year)))]]
  diff <- ref_gam$aic - best_tgam$aic_year
  df$thr <- ifelse(df$year <= best_tgam$yr, 'before', 'after')
  final_tgam <-  gam(bottom_temp ~ factor(year) +
                       s(longitude, latitude, by = factor(thr)) +
                       s(julian),
                     control = ctrl,
                     data = df)
  all_aic <- sapply(all_tgams, function(x) {as.numeric(x[1])})
  return_list <- list(ref_gam, final_tgam, best_tgam$yr, diff, all_aic)
}
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

years <- sort(unique(trawl_data$year))[4:22]
temp_tgam <- get_tgam(trawl_data, years)
plot_AIC(temp_tgam, years)

trawl_data$thr <-  ifelse(trawl_data$year <= trawl_data[[3]], 'before', 'after')
subset_distances <- function(tgam, df, year) {
  nlat = 40 # determine the resolution of the grid
  nlon = 30
  latd = seq(41, 48, length.out = nlat)
  lond = seq(-125, -123.9, length.out = nlon)
  subset_before <- expand.grid(latd, lond)
  names(subset_before) <- c("latitude", "longitude")
  subset_before$year <- year # change depending on species
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

tgam_prediction <- function(tgam, subset_before, year) {
  before_prediction <- predict(tgam[[2]],
                               newdata = subset_before,
                               se.fit = TRUE,
                               type = 'response')
  nlat = 40
  nlon = 30
  latd = seq(41, 48, length.out = nlat)
  lond = seq(-125,-123.9, length.out = nlon)
  pred_mean_before <- before_prediction[[1]]
  pred_se_before <- before_prediction[[2]]
  pred_mean_before[subset_before$dist > 10000] <- NA
  pred_se_before[subset_before$dist > 10000] <- NA
  subset_after <- expand.grid(latd, lond)
  names(subset_after) <- c("latitude", "longitude")
  subset_after$year <- year
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
  pred_mean_up_before <- pred_mean_before + 1.96 * pred_se_before
  pred_mean_down_before <- pred_mean_before - 1.96 * pred_se_before
  pred_mean_up_after <- pred_mean_after + 1.96 * pred_se_after
  pred_mean_down_after <- pred_mean_after - 1.96 * pred_se_after
  significant_low <- pred_mean_up_after < pred_mean_down_before
  significant_high <- pred_mean_down_after > pred_mean_up_before
  return(list(significant_high, significant_low, before_prediction, after_prediction))
}
temp_dist <- subset_distances(temp_tgam, trawl_data, 1986)
temp_CI <- tgam_prediction(temp_tgam, temp_dist, 2014)
significant_high <- temp_CI[[1]]
significant_low <- temp_CI[[2]]

# TGAM Maps
tgam_map <- function(species_subset, tgam, latitude, threshold, title, bathy.dat, bathy.mat){
  myvis_gam(tgam[[2]],
            view = c('longitude', 'latitude'),
            too.far = 0.025,
            plot.type = 'contour',
            color = 'jet',
            type = 'response',
            cond = list(thr = threshold),
            main = title,
            contour.col = "gray35",
            xlim = range(species_subset$longitude, na.rm = TRUE) + c(-0.3, .8),
            ylim = range(species_subset$latitude, na.rm = TRUE) + c(-0.5, 0.5),
            cex.lab = 2,
            cex.main = 2.5,
            cex.axis = 1.8,
            ylab = latitude,
            xlab = "Longitude °W")
  contour(unique(bathy.dat$lon),
          sort(unique(bathy.dat$lat)),
          bathy.mat,
          levels = -seq(200, 200, by = 200),
          labcex = 1.1,
          add = T,
          col = 'black',
          labels = NULL,
          lwd = 1.95)
  maps::map('worldHires',
            add = T,
            col = 'peachpuff3',
            fill = T)
  points(-124.0535, 44.6368, pch = 20)
  text(-124.0535,
       44.6368,
       "Newport",
       adj = c(0, 1.2),
       cex = 1.8)
  points(-123.8313, 46.1879, pch = 20)
  text(-123.88028,
       46.13361,
       "Astoria",
       adj = c(0, 1.2),
       cex = 1.8)
  points(-124.3, 43.3, pch = 20)
  text(-124.3,
       43.3,
       "Charleston",
       adj = c(0, 1.2),
       cex = 1.8)
}

tiff("../final_figs/Figure_3.tiff", width = 11, height = 12, units = "in", res = 300)
par(mfrow=c(1, 2),
    family = 'serif',
    mar=c(4, 5, 3, .2) + .1)
tgam_map(trawl_data, temp_tgam, "Latitude °N", "before", "Before 1992", bathy.dat, bathy.mat)
image.plot(legend.only = T,
           zlim = c(4, 10),
           col = jet.colors(100),
           legend.shrink = 0.2,
           smallplot = c(.25, .27, .1, .25),
           legend.cex = 2,
           legend.lab = "temperature (°C)",
           axis.args = list(cex.axis = 1.8),
           legend.width = 1,
           legend.line = 3.3)
tgam_map(trawl_data, temp_tgam, " ", "after", "After 1992", bathy.dat, bathy.mat)
points(temp_dist$longitude[significant_low],
       temp_dist$latitude[significant_low],
       pch = 2,
       cex = .9)
points(temp_dist$longitude[significant_high],
       temp_dist$latitude[significant_high],
       pch = 16,
       cex = 1.1)
legend("bottomleft",
       legend = c("Decrease", "Increase"),
       pch = c(2, 16),
       bty = "n",
       pt.cex = 2.5,
       cex = 1.7,
       inset = c(0.01, 0.001))
dev.off()

# Presentation Location map ----
browns <-
  c(
    "navajowhite4",
    "navajowhite3",
    "navajowhite2",
    "navajowhite1",
    "navajowhite",
    "moccasin"
  )
blues <-
  c("lightsteelblue4",
    "lightsteelblue3",
    "lightsteelblue2",
    "lightsteelblue1")

westcoast_bathy <-
  getNOAA.bathy(
    lon1 = -129,
    lon2 = -117,
    lat1 = 49,
    lat2 = 30,
    resolution = 1
  )

windows(width = 8, height = 12.8)
par(family = "serif")
plot.bathy(
  westcoast_bathy,
  image = T,
  axes = T,
  lwd = 0.03,
  land = T,
  n = 0,
  bpal = list(c(0, max(westcoast_bathy), browns), c(min(westcoast_bathy), 0, blues)),
  ylim = c(30, 48.5),
  xlim = c(-128,-117),
  xlab = "Longitude °W",
  ylab = "Latitude °N",
  main = "Northern California Current",
  cex.lab = 1,
  cex.main = 1,
  cex.axis = 1.1
)
points(-118.2437, 34.0522, pch = 20)
text(-120.2, 35.3, "Los Angeles", adj = c(0, 1.2))
points(-122.6750, 45.5051, pch = 20, cex = 0.8)
text(-122.6750, 45.5051, "Portland", adj = c(0, 1.2))
points(-122.3321, 47.6062, pch = 20)
text(-122.32, 47.6062, "Seattle", adj = c(0, 1.2))
points(-122.4194, 37.7749, pch = 20)
text(-121.4, 37.7749, "San Francisco", adj = c(0, 1.2))
points(-121.4944, 38.5816, pch = 20)
text(-121.4944, 39.4, "Sacramento", adj = c(0, 1.2))
points(-123.2620, 44.5646, pch = 20)
text(-123.2620, 44.5646, "Corvallis", adj = c(0, 1.2))
