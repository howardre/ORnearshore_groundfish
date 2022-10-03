### Title: Visualization of effort
### Purpose: Compare sampling and fishing effort
### Date Created: 06/02/2020

## Load libraries and data ----
setwd("C:/Users/howar/Documents/Oregon State/ORnearshore_groundfish/code/")
load("../data/ODFW_data/logbooks_corrected")
load('../data/NMFS_data/trawl_data')
survey_data <- trawl_data
source("functions/grid_fill.R")
source("functions/grid_data.R")
source("functions/grid_map.R")
source("functions/grid_pdf.R")
library(sgeostat)
library(maps)
library(mapdata)
library(fields)
library(marmap)
library(dplyr)
library(viridis)

# For depth, import data and show contour on a map
# .xyz option no longer available for download
bathy_dat <- read.table("../data/etopo1.xyz", sep = '')
names(bathy_dat) <- c('lon', 'lat', 'depth')
bathy_dat$depth[bathy_dat$depth > 0] <- NA # Avoid points above water
head(bathy_dat)
bathy_mat <- matrix(bathy_dat$depth,
                    nrow = length(unique(bathy_dat$lon)),
                    ncol = length(unique(bathy_dat$lat)))[, order(unique(bathy_dat$lat))]

## Reduce logbook data set ----
logbooks_final$month_day <- as.numeric(format(logbooks_final$TOWDATE, '%m%d'))
filtered <- logbooks_final[logbooks_final$month_day >= 517 &
                             logbooks_final$month_day <= 929 &
                             logbooks_final$depth <= -5, ]

# Generate a column to count each trawl
filtered$hauls <- filtered$GEAR[filtered$GEAR > 0] <- 1
trawl_counts <- filtered[filtered$species == 'PTRL_ADJ', ] # use one species since just counting effort

## Logbooks
### Make a regular grid and count stations within each grid cell ----
# These demonstrate the filling of grid cells
# grid_fill just gives the plot, grid_data gives the matrix and plot
# All years
nlat = 20 # determine resolution of grid
nlon = 15
latd = seq(42, 47, length.out = nlat)
lond = seq(-125, -123.9, length.out = nlon)
years_data <- length(unique(trawl_counts$year))
zlat <- (latd[1:(length(latd) - 1)] + latd[2:length(latd)]) / 2
zlon <- (lond[1:(length(lond) - 1)] + lond[2:length(lond)]) / 2
z_matrix <- matrix(nstations,
                   ncol = length(zlat),
                   nrow = length(zlon),
    byrow = F) / years_data # n of tows/year

# 80s
dev.new(width = 4, height = 10)
plot(trawl_counts$lon,
     trawl_counts$lat,
     pch = '.',
     ylim = c(42,47),
     xlim = c(-125, -123.9))
grid_fill(30, 25, 1980, 1990, trawl_counts)


# 90s
dev.new(width = 4, height = 10)
plot(trawl_counts$lon,
     trawl_counts$lat,
     pch = '.',
     ylim = c(42,47),
     xlim = c(-125, -123.9))
grid_fill(30, 25, 1989, 2002, trawl_counts)

# 00s
dev.new(width = 4, height = 10)
plot(trawl_counts$lon,
     trawl_counts$lat,
     pch = '.',
     ylim = c(42,47),
     xlim = c(-125, -123.9))
grid_fill(30, 25, 2001, 2010, trawl_counts)

# 10s
dev.new(width = 4, height = 10)
plot(trawl_counts$lon,
     trawl_counts$lat,
     pch = '.',
     ylim = c(42,47),
     xlim = c(-125, -123.9))
grid_fill(30, 25, 2009, 2019, trawl_counts)

## Decades
# This is the average effort per year for each decade
eighties_logbook <- grid_data(20, 15, 1980, 1990, trawl_counts)
nineties_logbook <- grid_data(20, 15, 1989, 2002, trawl_counts)
thousands_logbook <- grid_data(20, 15, 2001, 2010, trawl_counts)
teens_logbook <- grid_data(20, 15, 2009, 2018, trawl_counts)

save(eighties_logbook, file = "../data/ODFW_data/eighties_logbook")
save(nineties_logbook, file = "../data/ODFW_data/nineties_logbook")
save(thousands_logbook, file = "../data/ODFW_data/thousands_logbook")
save(teens_logbook, file = "../data/ODFW_data/teens_logbook")


### Map the grids ----
windows(width = 15, height = 9)
par(mfrow = c(1, 4),
    family = 'serif',
    mar = c(4, 5, 3, .3) + .1)
grid_map(20, 15, 1980, 1990, 230, trawl_counts, eighties_logbook, "Fishing Effort 1980s", bathy_dat, bathy_mat)
grid_map(20, 15, 1989, 2002, 190, trawl_counts, nineties_logbook, "Fishing Effort 1990s", bathy_dat, bathy_mat)
grid_map(20, 15, 2001, 2010, 180, trawl_counts, thousands_logbook, "Fishing Effort 2000s", bathy_dat, bathy_mat)
grid_map(20, 15, 2009, 2018, 180, trawl_counts, teens_logbook, "Fishing Effort 2010s", bathy_dat, bathy_mat)

## Survey ----
### Make a regular grid and count stations within each grid cell ----
survey_data$hauls <- trawl_data$month_day[survey_data$year_month > 0] <- 1

# 80s
dev.new(width = 4, height = 10)
plot(survey_data$lon,
     survey_data$lat,
     pch = '.',
     ylim = c(42,47),
     xlim = c(-125, -123.9))
grid_fill(30, 25, 1980, 1990, survey_data)


# 90s
dev.new(width = 4, height = 10)
plot(survey_data$lon,
     survey_data$lat,
     pch = '.',
     ylim = c(42,47),
     xlim = c(-125, -123.9))
grid_fill(30, 25, 1989, 2002, survey_data)

# 00s
dev.new(width = 4, height = 10)
plot(survey_data$lon,
     survey_data$lat,
     pch = '.',
     ylim = c(42,47),
     xlim = c(-125, -123.9))
grid_fill(30, 25, 2001, 2010, survey_data)

# 10s
dev.new(width = 4, height = 10)
plot(survey_data$lon,
     survey_data$lat,
     pch = '.',
     ylim = c(42,47),
     xlim = c(-125, -123.9))
grid_fill(30, 25, 2009, 2019, survey_data)

# All years
years_data <- length(unique(survey_data$year))
zlat <- (latd[1:(length(latd) - 1)] + latd[2:length(latd)]) / 2
zlon <- (lond[1:(length(lond) - 1)] + lond[2:length(lond)]) / 2

# Decades
# This is the average effort per year for each decade
eighties_survey <- grid_data(20, 15, 1980, 1990, survey_data)
nineties_survey <- grid_data(20, 15, 1989, 2002, survey_data)
thousands_survey <- grid_data(20, 15, 2001, 2010, survey_data)
teens_survey <- grid_data(20, 15, 2009, 2018, survey_data)

### Map the grids ----
windows(width = 15, height = 9)
par(mfrow = c(1, 4),
    family = 'serif',
    mar = c(4, 5, 3, .3) + .1)
grid_map(20, 15, 1980, 1990, 3, survey_data, eighties_survey, "Survey Effort 1980s", bathy_dat, bathy_mat)
grid_map(20, 15, 1989, 2002, 3, survey_data, nineties_survey, "Survey Effort 1990s", bathy_dat, bathy_mat)
grid_map(20, 15, 2001, 2010, 3, survey_data, thousands_survey, "Survey Effort 2000s", bathy_dat, bathy_mat)
grid_map(20, 15, 2009, 2018, 3, survey_data, teens_survey, "Survey Effort 2010s", bathy_dat, bathy_mat)

### Final figure ----
pdf("../results/Gear/four_panel.pdf",
    width = 7.5,
    height = 18)
par(mfrow = c(4, 4),
    family = "serif",
    mar = c(4, 5, 3, .3) + .1)
grid_pdf(20, 15, 1980, 1990, 230, trawl_counts, eighties_logbook,
         "Fishing Effort 1980s", bathy_dat, bathy_mat, "Latitude")
image.plot(legend.only = T,
           col = viridis(100, option = "A", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.76, .82, .09, .24),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.6),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, 230),
           legend.args = list("avg. annual number \n of tows",
                              side = 2, cex = 1.4))
grid_pdf(20, 15, 2009, 2018, 230, trawl_counts, teens_logbook, "Fishing Effort 2010s", bathy_dat, bathy_mat)
grid_pdf(20, 15, 1980, 1990, 4.5, survey_data, eighties_survey, "Survey Effort 1980s", bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "A", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.76, .82, .09, .24),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.6),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, 4.5),
           legend.args = list("avg. annual number \n of tows",
                              side = 2, cex = 1.4))
grid_pdf(20, 15, 2009, 2018, 3, survey_data, teens_survey, "Survey Effort 2010s", bathy_dat, bathy_mat)
dev.off()


tiff("../results/Gear/eight_panel.tiff",
    width = 15,
    height = 17.5,
    units = "in",
    compression = "lzw",
    res = 300)
par(mfrow = c(2, 4),
    family = "serif",
    mar = c(4, 5, 3, .3) + .1)
grid_pdf(20, 15, 1980, 1990, 230, trawl_counts, eighties_logbook, "Fishing Effort 1980s", bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "A", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.76, .82, .09, .24),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.6),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, 230),
           legend.args = list("avg. annual \n number of tows",
                              side = 2, cex = 1.4))
grid_pdf(20, 15, 1989, 2002, 230, trawl_counts, nineties_logbook, "Fishing Effort 1990s", bathy_dat, bathy_mat)
grid_pdf(20, 15, 2001, 2010, 230, trawl_counts, thousands_logbook, "Fishing Effort 2000s", bathy_dat, bathy_mat)
grid_pdf(20, 15, 2009, 2018, 230, trawl_counts, teens_logbook, "Fishing Effort 2010s", bathy_dat, bathy_mat)
grid_pdf(20, 15, 1980, 1990, 4.5, survey_data, eighties_survey, "Survey Effort 1980s", bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "A", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.76, .82, .09, .24),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.6),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, 4.5),
           legend.args = list("avg. annual \n number of tows",
                              side = 2, cex = 1.4))
grid_pdf(20, 15, 1989, 2002, 3, survey_data, nineties_survey, "Survey Effort 1990s", bathy_dat, bathy_mat)
grid_pdf(20, 15, 2001, 2010, 3, survey_data, thousands_survey, "Survey Effort 2000s", bathy_dat, bathy_mat)
grid_pdf(20, 15, 2009, 2018, 3, survey_data, teens_survey, "Survey Effort 2010s", bathy_dat, bathy_mat)
mtext(expression(paste("Latitude ("^0,'N)')),
      side = 2,
      line = -2.3,
      cex = 1.5,
      outer = TRUE)
mtext(expression(paste("Longitude ("^0, 'W)')),
      side = 1,
      line = -1.2,
      cex = 1.5,
      outer = TRUE)
dev.off()
