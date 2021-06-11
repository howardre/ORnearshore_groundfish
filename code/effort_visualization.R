### Title: Visualization of effort
### Purpose: Compare sampling and fishing effort
### Date Created: 06/02/2020

## Load libraries and data ----
setwd("C:/Users/howar/Documents/Oregon State/ORnearshore_groundfish/code/")
load("../data/ODFW_data/logbooks_corrected")
source("functions/grid_fill.R")
library(sgeostat)
library(maps)
library(mapdata)
library(fields)
library(marmap)
library(dplyr)

# For depth, import data and show contour on a map
# .xyz option no longer available for download
bathy_dat <- read.table("../data/etopo1.xyz", sep = '')
names(bathy_dat) <- c('lon', 'lat', 'depth')
bathy_dat$depth[bathy_dat$depth > 0] <- NA # Avoid points above water
head(bathy_dat)
bathy_mat <- matrix(bathy_dat$depth,
                    nrow = length(unique(bathy_dat$lon)),
                    ncol = length(unique(bathy_dat$lat)))[, order(unique(bathy_dat$lat))]

## Reduce data set ----
logbooks_final$month_day <- as.numeric(format(logbooks_final$TOWDATE, '%m%d'))
filtered <- logbooks_final[logbooks_final$month_day >= 517 &
                             logbooks_final$month_day <= 929 &
                             logbooks_final$depth <= -5, ]

# generate a column to count each trawl
filtered$hauls <- filtered$GEAR[filtered$GEAR > 0] <- 1
trawl_counts <- filtered[filtered$species == 'PTRL_ADJ', ] # use one species since just counting effort

## Make a regular grid and count stations within each grid cell ----
# All years
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
grid_fill(30, 25, 1989, 2000, trawl_counts)

# 00s
dev.new(width = 4, height = 10)
plot(trawl_counts$lon,
     trawl_counts$lat,
     pch = '.',
     ylim = c(42,47),
     xlim = c(-125, -123.9))
grid_fill(30, 25, 1999, 2010, trawl_counts)

# 10s
dev.new(width = 4, height = 10)
plot(trawl_counts$lon,
     trawl_counts$lat,
     pch = '.',
     ylim = c(42,47),
     xlim = c(-125, -123.9))
grid_fill(30, 25, 2009, 2019, trawl_counts)

## decades
zlat <- (latd[1:(length(latd) - 1)] + latd[2:length(latd)]) / 2
zlon <- (lond[1:(length(lond) - 1)] + lond[2:length(lond)]) / 2
eighties.data <- length(unique(trawl_counts$year[trawl_counts$year > 1980 &
                                trawl_counts$year < 1990]))
z_matrix1 <- matrix(nstations1,
                    ncol = length(zlat),
                    nrow = length(zlon),
                    byrow = F) / eighties.data
z_matrix1 <- ifelse(z_matrix1 == 0, NA, z_matrix1)
nineties.data <- length(unique(trawl_counts$year[trawl_counts$year > 1989 &
                                trawl_counts$year < 2000]))
z_matrix2 <- matrix(nstations2,
                    ncol = length(zlat),
                    nrow = length(zlon),
                    byrow = F) / nineties.data
z_matrix2 <- ifelse(z_matrix2 == 0, NA, z_matrix2)
thousands.data <- length(unique(trawl_counts$year[trawl_counts$year > 1999 &
                                trawl_counts$year < 2010]))
z_matrix3 <- matrix(nstations3,
                    ncol = length(zlat),
                    nrow = length(zlon),
                    byrow = F) / thousands.data
z_matrix3 <- ifelse(z_matrix3 == 0, NA, z_matrix3)
tens.data <- length(unique(trawl_counts$year[trawl_counts$year > 2009 &
                                trawl_counts$year < 2018]))
z_matrix4 <- matrix(nstations4,
                    ncol = length(zlat),
                    nrow = length(zlon),
                    byrow = F) / tens.data
z_matrix4 <- ifelse(z_matrix4 == 0, NA, z_matrix4)

## Map the grids ----
# windows(width=11,height=5)
# par(mfrow=c(1,2))
windows(width = 15, height = 9)
par(mfrow = c(1, 4),
    family = 'serif',
    mar = c(4, 5, 3, .3) + .1)

image(
  zlon,
  zlat,
  z_matrix1,
  col = hcl.colors(40, "RdYlBu", rev = T),
  xlim = c(-125, -123.6),
  ylim = c(42, 47),
  main = "Fishing Effort 1980s",
  ylab = expression(paste("Latitude (" ^ 0, 'N)')),
  xlab = expression(paste("Longitude (" ^ 0, 'W)')),
  zlim = c(0, 230)
)
map("worldHires",
    fill = T,
    col = "grey",
    add = T)
image.plot(
  legend.only = T,
  col = hcl.colors(40, "RdYlBu", rev = T),
  legend.shrink = 0.2,
  smallplot = c(.76, .81, .09, .25),
  legend.cex = 0.7,
  axis.args = list(cex.axis = 0.9),
  legend.width = 1,
  zlim = c(0, 230),
  legend.lab = "avg. number of tows"
)
contour(
  unique(bathy_dat$lon),
  sort(unique(bathy_dat$lat)),
  bathy_mat,
  lwd = 1,
  levels = -100,
  col = "gray19",
  add = T
)
contour(
  unique(bathy_dat$lon),
  sort(unique(bathy_dat$lat)),
  bathy_mat,
  lwd = 1,
  levels = -200,
  col = "gray19",
  add = T
)
points(-124.0535, 44.6368, pch = 20)
text(-124.0535, 44.6368, "Newport", adj = c(0, 1.2))
points(-123.8313, 46.1879, pch = 20)
text(-123.88028, 46.13361, "Astoria", adj = c(0, 1.2))
points(-124.3, 43.3, pch = 20)
text(-124.3, 43.3, "Charleston", adj = c(0, 1.2))

image(
  zlon,
  zlat,
  z_matrix2,
  col = hcl.colors(40, "RdYlBu", rev = T),
  xlim = c(-125, -123.6),
  ylim = c(42, 47),
  main = "Fishing Effort 1990s",
  ylab = expression(paste("Latitude (" ^ 0, 'N)')),
  xlab = expression(paste("Longitude (" ^ 0, 'W)')),
  zlim = c(0, 230)
)
map("worldHires",
    fill = T,
    col = "grey",
    add = T)
image.plot(
  legend.only = T,
  col = hcl.colors(40, "RdYlBu", rev = T),
  legend.shrink = 0.2,
  smallplot = c(.76, .81, .09, .25),
  legend.cex = 0.7,
  axis.args = list(cex.axis = 0.9),
  legend.width = 1,
  zlim = c(0, 230),
  legend.lab = "avg. number of tows"
)
contour(
  unique(bathy_dat$lon),
  sort(unique(bathy_dat$lat)),
  bathy_mat,
  lwd = 1,
  levels = -100,
  col = "gray19",
  add = T
)
contour(
  unique(bathy_dat$lon),
  sort(unique(bathy_dat$lat)),
  bathy_mat,
  lwd = 1,
  levels = -200,
  col = "gray19",
  add = T
)
points(-124.0535, 44.6368, pch = 20)
text(-124.0535, 44.6368, "Newport", adj = c(0, 1.2))
points(-123.8313, 46.1879, pch = 20)
text(-123.88028, 46.13361, "Astoria", adj = c(0, 1.2))
points(-124.3, 43.3, pch = 20)
text(-124.3, 43.3, "Charleston", adj = c(0, 1.2))

image(
  zlon,
  zlat,
  z_matrix3,
  col = hcl.colors(40, "RdYlBu", rev = T),
  xlim = c(-125, -123.6),
  ylim = c(42, 47),
  main = "Fishing Effort 2000s",
  ylab = expression(paste("Latitude (" ^ 0, 'N)')),
  xlab = expression(paste("Longitude (" ^ 0, 'W)')),
  zlim = c(0, 230)
)
map("worldHires",
    fill = T,
    col = "grey",
    add = T)
image.plot(
  legend.only = T,
  col = hcl.colors(40, "RdYlBu", rev = T),
  legend.shrink = 0.2,
  smallplot = c(.76, .81, .09, .25),
  legend.cex = 0.7,
  axis.args = list(cex.axis = 0.9),
  legend.width = 1,
  zlim = c(0, 230),
  legend.lab = "avg. number of tows"
)
contour(
  unique(bathy_dat$lon),
  sort(unique(bathy_dat$lat)),
  bathy_mat,
  lwd = 1,
  levels = -100,
  col = "gray19",
  add = T
)
contour(
  unique(bathy_dat$lon),
  sort(unique(bathy_dat$lat)),
  bathy_mat,
  lwd = 1,
  levels = -200,
  col = "gray19",
  add = T
)
points(-124.0535, 44.6368, pch = 20)
text(-124.0535, 44.6368, "Newport", adj = c(0, 1.2))
points(-123.8313, 46.1879, pch = 20)
text(-123.88028, 46.13361, "Astoria", adj = c(0, 1.2))
points(-124.3, 43.3, pch = 20)
text(-124.3, 43.3, "Charleston", adj = c(0, 1.2))

image(
  zlon,
  zlat,
  z_matrix4,
  col = hcl.colors(40, "RdYlBu", rev = T),
  xlim = c(-125, -123.6),
  ylim = c(42, 47),
  main = "Fishing Effort 2010s",
  ylab = expression(paste("Latitude (" ^ 0, 'N)')),
  xlab = expression(paste("Longitude (" ^ 0, 'W)')),
  zlim = c(0, 230)
)
map("worldHires",
    fill = T,
    col = "grey",
    add = T)
image.plot(
  legend.only = T,
  col = hcl.colors(40, "RdYlBu", rev = T),
  legend.shrink = 0.2,
  smallplot = c(.76, .81, .09, .25),
  legend.cex = 0.7,
  axis.args = list(cex.axis = 0.9),
  legend.width = 1,
  zlim = c(0, 230),
  legend.lab = "avg. number of tows"
)
contour(
  unique(bathy_dat$lon),
  sort(unique(bathy_dat$lat)),
  bathy_mat,
  lwd = 1,
  levels = -100,
  col = "gray19",
  add = T
)
contour(
  unique(bathy_dat$lon),
  sort(unique(bathy_dat$lat)),
  bathy_mat,
  lwd = 1,
  levels = -200,
  col = "gray19",
  add = T
)
points(-124.0535, 44.6368, pch = 20)
text(-124.0535, 44.6368, "Newport", adj = c(0, 1.2))
points(-123.8313, 46.1879, pch = 20)
text(-123.88028, 46.13361, "Astoria", adj = c(0, 1.2))
points(-124.3, 43.3, pch = 20)
text(-124.3, 43.3, "Charleston", adj = c(0, 1.2))

## Survey
load('/Users/howar/Documents/Oregon State/Thesis/Data Visualization/trawl_data')
trawl_data$hauls <-
  trawl_data$month_day[trawl_data$year_month > 0] <- 1

nlat = 20#determine resolution of grid
nlon = 15
latd = seq(42, 47, length.out = nlat)
lond = seq(-125, -123.9, length.out = nlon)


grid_lon = data.frame(
  lon1 = rep(lond[-length(lond)], (nlat - 1)),
  lon2 = rep(lond[-1], (nlat - 1)),
  lon3 = rep(lond[-1], (nlat - 1)),
  lon4 = rep(lond[-length(lond)], (nlat - 1))
)#make dataframe of just longitude

grid_lat = data.frame(
  lat1 = sort(rep(latd[-length(latd)], (nlon - 1))),
  lat2 = sort(rep(latd[-length(latd)], (nlon - 1))),
  lat3 = sort(rep(latd[-1], (nlon - 1))),
  lat4 = sort(rep(latd[-1], (nlon - 1)))
)#lat dataframe

#80s
dev.new(width = 4, height = 10)
plot(
  trawl_data$longitude[trawl_data$year > 1979 &
                         trawl_data$year < 1990],
  trawl_data$latitude[trawl_data$year > 1979 &
                        trawl_data$year < 1990],
  pch = '.',
  ylim = c(42, 47),
  xlim = c(-125, -123.9)
)
nstations2 = NA * (1:nrow(grid_lon))

for (i in 1:length(nstations1)) {
  tmp = in.chull(trawl_data$longitude[trawl_data$year > 1979 &
                                        trawl_data$year < 1990],
                 trawl_data$latitude[trawl_data$year > 1979 &
                                       trawl_data$year < 1990],
                 grid_lon[i, ],
                 grid_lat[i, ])
  nstations1[i] = sum(trawl_data$hauls[trawl_data$year > 1979 &
                                          trawl_data$year < 1990] * tmp)#This decides what goes into each grid pixel
  points(trawl_data$longitude[tmp],
         trawl_data$latitude[tmp],
         col = i,
         pch = 16)
  polygon(grid_lon[i, ], grid_lat[i, ])
}

##90s
dev.new(width = 4, height = 10)
plot(
  trawl_data$longitude[trawl_data$year > 1989 &
                         trawl_data$year < 2000],
  trawl_data$latitude[trawl_data$year > 1989 &
                        trawl_data$year < 2000],
  pch = '.',
  ylim = c(42, 47),
  xlim = c(-125, -123.9)
)
nstations2 = NA * (1:nrow(grid_lon))

for (i in 1:length(nstations2)) {
  tmp = in.chull(trawl_data$longitude[trawl_data$year > 1989 &
                                        trawl_data$year < 2000],
                 trawl_data$latitude[trawl_data$year > 1989 &
                                       trawl_data$year < 2000],
                 grid_lon[i, ],
                 grid_lat[i, ])
  nstations2[i] = sum(trawl_data$hauls[trawl_data$year > 1989 &
                                          trawl_data$year < 2000] * tmp)#This decides what goes into each grid pixel
  points(trawl_data$longitude[tmp],
         trawl_data$latitude[tmp],
         col = i,
         pch = 16)
  polygon(grid_lon[i, ], grid_lat[i, ])
}


##00s
dev.new(width = 4, height = 10)
plot(
  trawl_data$longitude[trawl_data$year >= 2000 &
                         trawl_data$year <= 2009],
  trawl_data$latitude[trawl_data$year >= 2000 &
                        trawl_data$year <= 2009],
  pch = '.',
  ylim = c(42, 47),
  xlim = c(-125, -123.9)
)
nstations3 = NA * (1:nrow(grid_lon))

for (i in 1:length(nstations3)) {
  tmp = in.chull(trawl_data$longitude[trawl_data$year >= 2000 &
                                        trawl_data$year <= 2009],
                 trawl_data$latitude[trawl_data$year >= 2000 &
                                       trawl_data$year <= 2009],
                 grid_lon[i, ],
                 grid_lat[i, ])
  nstations3[i] = sum(trawl_data$hauls[trawl_data$year >= 2000 &
                                          trawl_data$year <= 2009] * tmp)#This decides what goes into each grid pixel
  points(trawl_data$longitude[tmp],
         trawl_data$latitude[tmp],
         col = i,
         pch = 16)
  polygon(grid_lon[i, ], grid_lat[i, ])
}

##10s
dev.new(width = 4, height = 10)
plot(
  trawl_data$longitude[trawl_data$year >= 2010 &
                         trawl_data$year <= 2017],
  trawl_data$latitude[trawl_data$year >= 2010 &
                        trawl_data$year <= 2017],
  pch = '.',
  ylim = c(42, 47),
  xlim = c(-125, -123.9)
)
nstations4 = NA * (1:nrow(grid_lon))

for (i in 1:length(nstations4)) {
  tmp = in.chull(trawl_data$longitude[trawl_data$year >= 2010 &
                                        trawl_data$year <= 2017],
                 trawl_data$latitude[trawl_data$year >= 2010 &
                                       trawl_data$year <= 2017],
                 grid_lon[i, ],
                 grid_lat[i, ])
  nstations4[i] = sum(trawl_data$hauls[trawl_data$year >= 2010 &
                                          trawl_data$year <= 2017] * tmp)#This decides what goes into each grid pixel
  points(trawl_data$longitude[tmp],
         trawl_data$latitude[tmp],
         col = i,
         pch = 16)
  polygon(grid_lon[i, ], grid_lat[i, ])
}

##ALL YEARS
years_data <- length(unique(trawl_data$year))
zlat <- (latd[1:(length(latd) - 1)] + latd[2:length(latd)]) / 2
zlon <- (lond[1:(length(lond) - 1)] + lond[2:length(lond)]) / 2

## decades
eighties.data <-
  length(unique(trawl_data$year[trawl_data$year >= 1980 &
                                  trawl_data$year <= 1989]))
z_matrix1 <-
  matrix(
    nstations1,
    ncol = length(zlat),
    nrow = length(zlon),
    byrow = F
  ) / eighties.data
z_matrix1 <- ifelse(z_matrix1 == 0, NA, z_matrix1)
nineties.data <-
  length(unique(trawl_data$year[trawl_data$year > 1989 &
                                  trawl_data$year < 2000]))
z_matrix2 <-
  matrix(
    nstations2,
    ncol = length(zlat),
    nrow = length(zlon),
    byrow = F
  ) / nineties.data
z_matrix2 <- ifelse(z_matrix2 == 0, NA, z_matrix2)
thousands.data <-
  length(unique(trawl_data$year[trawl_data$year > 1999 &
                                  trawl_data$year < 2010]))
z_matrix3 <-
  matrix(
    nstations3,
    ncol = length(zlat),
    nrow = length(zlon),
    byrow = F
  ) / thousands.data
z_matrix3 <- ifelse(z_matrix3 == 0, NA, z_matrix3)
tens.data <-
  length(unique(trawl_data$year[trawl_data$year > 2009 &
                                  trawl_data$year < 2018]))
z_matrix4 <-
  matrix(
    nstations4,
    ncol = length(zlat),
    nrow = length(zlon),
    byrow = F
  ) / tens.data
z_matrix4 <- ifelse(z_matrix4 == 0, NA, z_matrix4)

###MAKE MAPS

#windows(width=11,height=5)
#par(mfrow=c(1,2))
windows(width = 15, height = 9)
par(mfrow = c(1, 4),
    family = 'serif',
    mar = c(4, 5, 3, .3) + .1)

image(
  zlon,
  zlat,
  z_matrix1,
  col = hcl.colors(40, "RdYlBu", rev = T),
  xlim = c(-125, -123.6),
  ylim = c(42, 47),
  main = "Survey Effort 1980s",
  ylab = expression(paste("Latitude (" ^ 0, 'N)')),
  xlab = expression(paste("Longitude (" ^ 0, 'W)')),
  zlim = c(0, 4.5)
)
map("worldHires",
    fill = T,
    col = "grey",
    add = T)
image.plot(
  legend.only = T,
  col = hcl.colors(40, "RdYlBu", rev = T),
  legend.shrink = 0.2,
  smallplot = c(.76, .81, .09, .25),
  legend.cex = 0.7,
  axis.args = list(cex.axis = 0.9),
  legend.width = 1,
  zlim = c(0, 4.5),
  legend.lab = "avg. number of tows"
)
contour(
  unique(bathy_dat$lon),
  sort(unique(bathy_dat$lat)),
  bathy_mat,
  lwd = 1,
  levels = -100,
  col = "gray19",
  add = T
)
contour(
  unique(bathy_dat$lon),
  sort(unique(bathy_dat$lat)),
  bathy_mat,
  lwd = 1,
  levels = -200,
  col = "gray19",
  add = T
)
points(-124.0535, 44.6368, pch = 20)
text(-124.0535, 44.6368, "Newport", adj = c(0, 1.2))
points(-123.8313, 46.1879, pch = 20)
text(-123.88028, 46.13361, "Astoria", adj = c(0, 1.2))
points(-124.3, 43.3, pch = 20)
text(-124.3, 43.3, "Charleston", adj = c(0, 1.2))

image(
  zlon,
  zlat,
  z_matrix2,
  col = hcl.colors(40, "RdYlBu", rev = T),
  xlim = c(-125, -123.6),
  zlim = c(0, 4.5),
  main = "Survey Effort 1990s",
  ylab = expression(paste("Latitude (" ^ 0, 'N)')),
  xlab = expression(paste("Longitude (" ^ 0, 'W)'))
)
map("worldHires",
    fill = T,
    col = "grey",
    add = T)
image.plot(
  legend.only = T,
  col = hcl.colors(40, "RdYlBu", rev = T),
  legend.shrink = 0.2,
  smallplot = c(.76, .81, .09, .25),
  legend.cex = 0.7,
  axis.args = list(cex.axis = 0.9),
  legend.width = 1,
  zlim = c(0, 4.5),
  legend.lab = "avg. number of tows"
)
contour(
  unique(bathy_dat$lon),
  sort(unique(bathy_dat$lat)),
  bathy_mat,
  lwd = 1,
  levels = -100,
  col = "gray19",
  add = T
)
contour(
  unique(bathy_dat$lon),
  sort(unique(bathy_dat$lat)),
  bathy_mat,
  lwd = 1,
  levels = -200,
  col = "gray19",
  add = T
)
points(-124.0535, 44.6368, pch = 20)
text(-124.0535, 44.6368, "Newport", adj = c(0, 1.2))
points(-123.8313, 46.1879, pch = 20)
text(-123.88028, 46.13361, "Astoria", adj = c(0, 1.2))
points(-124.3, 43.3, pch = 20)
text(-124.3, 43.3, "Charleston", adj = c(0, 1.2))

image(
  zlon,
  zlat,
  z_matrix3,
  col = hcl.colors(40, "RdYlBu", rev = T),
  xlim = c(-125, -123.6),
  ylim = c(42, 47),
  main = "Survey Effort 2000s",
  ylab = expression(paste("Latitude (" ^ 0, 'N)')),
  xlab = expression(paste("Longitude (" ^ 0, 'W)')),
  zlim = c(0, 4.5)
)
map("worldHires",
    fill = T,
    col = "grey",
    add = T)
image.plot(
  legend.only = T,
  col = hcl.colors(40, "RdYlBu", rev = T),
  legend.shrink = 0.2,
  smallplot = c(.76, .81, .09, .25),
  legend.cex = 0.7,
  axis.args = list(cex.axis = 0.9),
  legend.width = 1,
  zlim = c(0, 4.5),
  legend.lab = "avg. number of tows"
)
contour(
  unique(bathy_dat$lon),
  sort(unique(bathy_dat$lat)),
  bathy_mat,
  lwd = 1,
  levels = -100,
  col = "gray19",
  add = T
)
contour(
  unique(bathy_dat$lon),
  sort(unique(bathy_dat$lat)),
  bathy_mat,
  lwd = 1,
  levels = -200,
  col = "gray19",
  add = T
)
points(-124.0535, 44.6368, pch = 20)
text(-124.0535, 44.6368, "Newport", adj = c(0, 1.2))
points(-123.8313, 46.1879, pch = 20)
text(-123.88028, 46.13361, "Astoria", adj = c(0, 1.2))
points(-124.3, 43.3, pch = 20)
text(-124.3, 43.3, "Charleston", adj = c(0, 1.2))

image(
  zlon,
  zlat,
  z_matrix4,
  col = hcl.colors(40, "RdYlBu", rev = T),
  xlim = c(-125, -123.6),
  ylim = c(42, 47),
  main = "Survey Effort 2010s",
  ylab = expression(paste("Latitude (" ^ 0, 'N)')),
  xlab = expression(paste("Longitude (" ^ 0, 'W)')),
  zlim = c(0, 4.5)
)
map("worldHires",
    fill = T,
    col = "grey",
    add = T)
image.plot(
  legend.only = T,
  col = hcl.colors(40, "RdYlBu", rev = T),
  legend.shrink = 0.2,
  smallplot = c(.76, .81, .09, .25),
  legend.cex = 0.7,
  axis.args = list(cex.axis = 0.9),
  legend.width = 1,
  zlim = c(0, 4.5),
  legend.lab = "avg. number of tows"
)
contour(
  unique(bathy_dat$lon),
  sort(unique(bathy_dat$lat)),
  bathy_mat,
  lwd = 1,
  levels = -100,
  col = "gray19",
  add = T
)
contour(
  unique(bathy_dat$lon),
  sort(unique(bathy_dat$lat)),
  bathy_mat,
  lwd = 1,
  levels = -200,
  col = "gray19",
  add = T
)
points(-124.0535, 44.6368, pch = 20)
text(-124.0535, 44.6368, "Newport", adj = c(0, 1.2))
points(-123.8313, 46.1879, pch = 20)
text(-123.88028, 46.13361, "Astoria", adj = c(0, 1.2))
points(-124.3, 43.3, pch = 20)
text(-124.3, 43.3, "Charleston", adj = c(0, 1.2))
