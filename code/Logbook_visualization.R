### Title: Visualization of species catches
### Purpose: Compare catches of 6 species
### Date Created: 06/02/2020

## Load libraries, functions, and data ----
library(dplyr)
library(tidyr)
library(reshape)
library(reshape2)
library(maps)
library(mapdata)
library(sgeostat)
library(fields)
library(ggplot2)
library(ggthemes)
library(ineq)
library(reldist)
library(gridExtra)

setwd("C:/Users/howar/Documents/Oregon State/ORnearshore_groundfish/code/")
load("../data/ODFW_data/logbooks_corrected")
load("../data/ODFW_data/fish_tickets_final")
load('../data/NMFS_data/trawl_data')
load('../data/NMFS_data/OR_fish')
survey_data <- trawl_data
remove(trawl_data)
source("functions/exploration_panels.R")
source("functions/cpue_fillpts.R")
source("functions/cpue_grid.R")
source("functions/cpue_map.R")
source("functions/depth_contours.R")
source("functions/species_grid_pdf.R")

# For depth, import data and show contour on a map
# .xyz option no longer available for download
bathy_dat <- read.table("../data/etopo1.xyz", sep = '')
names(bathy_dat) <- c('lon', 'lat', 'depth')
bathy_dat$depth[bathy_dat$depth > 0] <- NA # Avoid points above water
bathy_mat <- matrix(bathy_dat$depth,
                    nrow = length(unique(bathy_dat$lon)),
                    ncol = length(unique(bathy_dat$lat)))[, order(unique(bathy_dat$lat))]

# Logbooks: convert CPUE to kg/hr and ln(x+1), add other necessary columns
logbooks_final$kg_caught <- logbooks_final$species_weight * 0.4535924
logbooks_final <- logbooks_final[!logbooks_final$DURATION == 0, ] # remove tows with no trawl duration
logbooks_final$CPUE <- logbooks_final$kg_caught / logbooks_final$DURATION
logbooks_final <- logbooks_final[!is.na(logbooks_final$CPUE), ]
logbooks_final$lncpue <- log(logbooks_final$CPUE + 1)
logbooks_final <- logbooks_final[logbooks_final$depth <= -5, ] # remove unreasonably shallow tows
logbooks_final$pres <- 1 * (logbooks_final$species_weight > 0)
logbooks_final$month_day <- as.numeric(format(logbooks_final$TOWDATE, '%m%d'))

# Survey: fill any NAs with 0s, rename columns to match logbooks, filter out columns
OR_fish$lncpue[is.na(OR_fish$lncpue)] <- 0
OR_fish$cpue_kg[is.na(OR_fish$cpue_kg)] <- 0
OR_fish$total_catch_wt_kg[is.na(OR_fish$total_catch_wt_kg)] <- 0
names(OR_fish)[names(OR_fish) == "total_catch_wt_kg"] <- "catch"
names(OR_fish)[names(OR_fish) == "cpue_kg"] <- "CPUE"
names(OR_fish)[names(OR_fish) == "longitude"] <- "lon"
names(OR_fish)[names(OR_fish) == "latitude"] <- "lat"
names(OR_fish)[names(OR_fish) == "depth_m"] <- "depth"
OR_fish <- select(OR_fish, scientific_name, depth, lat, lon, catch, CPUE, lncpue, year)
OR_fish$pres <- 1 * (OR_fish$catch > 0)
OR_fish$depth <- -abs(OR_fish$depth)

# Make grid for the map panels
nlat = 20 # determine resolution of grid
nlon = 15
latd = seq(42, 47, length.out = nlat)
lond = seq(-125, -123.9, length.out = nlon)

grid_lon = data.frame(
  lon1 = rep(lond[-length(lond)], (nlat - 1)),
  lon2 = rep(lond[-1], (nlat - 1)),
  lon3 = rep(lond[-1], (nlat - 1)),
  lon4 = rep(lond[-length(lond)], (nlat - 1))) # make dataframe of just longitude

grid_lat = data.frame(
  lat1 = sort(rep(latd[-length(latd)], (nlon - 1))),
  lat2 = sort(rep(latd[-length(latd)], (nlon - 1))),
  lat3 = sort(rep(latd[-1], (nlon - 1))),
  lat4 = sort(rep(latd[-1], (nlon - 1)))) # lat dataframe

zlat <- (latd[1:(length(latd) - 1)] + latd[2:length(latd)]) / 2
zlon <- (lond[1:(length(lond) - 1)] + lond[2:length(lond)]) / 2

## Petrale sole ----
### Logbooks ----
# Filter to just petrale and survey months
subset_petrale_logbook <- logbooks_final[logbooks_final$species == 'PTRL_ADJ', ]
subset_petrale_logbook <- subset_petrale_logbook[subset_petrale_logbook$month_day >= 517 &
                                   subset_petrale_logbook$month_day <= 929, ]

# Create data exploration map
windows(width = 28, height = 18)
par(mfrow = c(1, 4))
exploration_panels(subset_petrale_logbook, bathy_dat, bathy_mat, 1981, 1989, "Petrale Sole 1980s")
exploration_panels(subset_petrale_logbook, bathy_dat, bathy_mat, 1990, 1999, "Petrale Sole 1990s")
exploration_panels(subset_petrale_logbook, bathy_dat, bathy_mat, 2000, 2009, "Petrale Sole 2000s")
exploration_panels(subset_petrale_logbook, bathy_dat, bathy_mat, 2010, 2017, "Petrale Sole 2010s")

## Decade maps
# Make sure grid has been created (above)
dev.new(width = 4, height = 10)
cpue_fillpts(subset_petrale_logbook, 1981, 1989)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_petrale_logbook, 1990, 1999)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_petrale_logbook, 2000, 2009)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_petrale_logbook, 2010, 2017)

## Decade data
# Make sure grid chunk has been run
eighties_logbooks_ptrl <- cpue_grid(subset_petrale_logbook, 1981, 1989)
nineties_logbooks_ptrl <- cpue_grid(subset_petrale_logbook, 1990, 1999)
thousands_logbooks_ptrl <- cpue_grid(subset_petrale_logbook, 2000, 2009)
tens_logbooks_ptrl <- cpue_grid(subset_petrale_logbook, 2010, 2017)

## Make maps
# Change legend depending on maximum lncpue - use matrix with highest CPUE
max(eighties_logbooks_ptrl, na.rm = T)
max(nineties_logbooks_ptrl, na.rm = T)
max(thousands_logbooks_ptrl, na.rm = T)
max(tens_logbooks_ptrl, na.rm = T)

windows(width = 15, height = 9)
par(mfrow = c(1, 4),
    family = 'serif',
    mar = c(4, 5, 3, .3) + .1)
cpue_map(eighties_logbooks_ptrl, tens_logbooks_ptrl, "PurpOr", "Logbook Petrale Sole 1980s", bathy_dat, bathy_mat)
cpue_map(nineties_logbooks_ptrl, tens_logbooks_ptrl, "PurpOr", "Logbook Petrale Sole 1990s", bathy_dat, bathy_mat)
cpue_map(thousands_logbooks_ptrl, tens_logbooks_ptrl, "PurpOr", "Logbook Petrale Sole 2000s", bathy_dat, bathy_mat)
cpue_map(tens_logbooks_ptrl, tens_logbooks_ptrl, "PurpOr", "Logbook Petrale Sole 2010s", bathy_dat, bathy_mat)

## Depth distribution by year
logs_eighties_depth_ptrl <- depth_contours(subset_petrale_logbook, 1989, "BuPu", "1980s")
logs_nineties_depth_ptrl <- depth_contours(subset_petrale_logbook, 1999, "BuPu", "1990s")
logs_thousands_depth_ptrl <- depth_contours(subset_petrale_logbook, 2009, "BuPu", "2000s")
logs_tens_depth_ptrl <- depth_contours(subset_petrale_logbook, 2017, "BuPu", "2010s")

windows()
grid.arrange(logs_eighties_depth_ptrl,
             logs_nineties_depth_ptrl,
             logs_thousands_depth_ptrl,
             logs_tens_depth_ptrl,
             ncol = 2,
             top = textGrob("Petrale Sole Logbook Depth Distribution",
                            gp = gpar(fontfamily = "serif",
                                      cex = 1,
                                      fontface = "bold")))

## Catch over time
catch_petrale <- tickets_final %>% filter(common_name == 'Petrale Sole') %>%
  group_by(YEAR) %>%
  summarise(species_weight_sum = sum(TIK_LBS))

# Plot
windows(width = 20, height = 10)
par(mfrow = c(1, 1))
ggplot(data = catch_petrale,
       aes(x = YEAR, y = species_weight_sum / 1000)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
  geom_col(fill = "blue4") +
  theme_tufte() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "Year",
       y = "Total Catch (1000s of lbs)",
       title = "Petrale Sole Total Nearshore Catch")

## Average CPUE over time
logbook_petrale_cpue <- subset_petrale_logbook %>% group_by(year) %>%
  summarise(cpue_mean = mean(CPUE))

# Plot
windows(width = 20, height = 10)
par(mfrow = c(1, 1))
ggplot(data = logbook_petrale_cpue, aes(x = year, y = cpue_mean)) +
  geom_path() +
  geom_point() +
  theme_tufte() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "Year",
       y = "Cpue (kg/hr)",
       title = "Mean CPUE of Nearshore Petrale Sole Caught in Groundfish Fishery")

## Seasonality
# Create winter subset
winter_petrale <- filter(logbooks_final, !between(month_day, 516, 930),
                         species == 'PTRL_ADJ')

## Winter Decades
dev.new(width = 4, height = 10)
cpue_fillpts(winter_petrale, 1981, 1989)

dev.new(width = 4, height = 10)
cpue_fillpts(winter_petrale, 1990, 1999)

dev.new(width = 4, height = 10)
cpue_fillpts(winter_petrale, 2000, 2009)

dev.new(width = 4, height = 10)
cpue_fillpts(winter_petrale, 2010, 2017)

## Decade data
eighties_logbooks_ptrl_winter <- cpue_grid(winter_petrale, 1981, 1989)
nineties_logbooks_ptrl_winter <- cpue_grid(winter_petrale, 1990, 1999)
thousands_logbooks_ptrl_winter <- cpue_grid(winter_petrale, 2000, 2009)
tens_logbooks_ptrl_winter <- cpue_grid(winter_petrale, 2010, 2017)

## Make maps
# Change legend depending on maximum lncpue - use matrix with highest CPUE
max(eighties_logbooks_ptrl_winter, na.rm = T)
max(nineties_logbooks_ptrl_winter, na.rm = T)
max(thousands_logbooks_ptrl_winter, na.rm = T)
max(tens_logbooks_ptrl_winter, na.rm = T)

windows(width = 15, height = 9)
par(mfrow = c(1, 4),
    family = 'serif',
    mar = c(4, 5, 3, .3) + .1)
cpue_map(eighties_logbooks_ptrl_winter, tens_logbooks_ptrl_winter, "PurpOr", "Logbook Petrale Sole 1980s", bathy_dat, bathy_mat)
cpue_map(nineties_logbooks_ptrl_winter, tens_logbooks_ptrl_winter, "PurpOr", "Logbook Petrale Sole 1990s", bathy_dat, bathy_mat)
cpue_map(thousands_logbooks_ptrl_winter, tens_logbooks_ptrl_winter, "PurpOr", "Logbook Petrale Sole 2000s", bathy_dat, bathy_mat)
cpue_map(tens_logbooks_ptrl_winter, tens_logbooks_ptrl_winter, "PurpOr", "Logbook Petrale Sole 2010s", bathy_dat, bathy_mat)

## Depth distribution by year (winter)
logs_eighties_depth_ptrlw <- depth_contours(subset_petrale_logbook, 1989, "BuPu", "1980s")
logs_nineties_depth_ptrlw <- depth_contours(subset_petrale_logbook, 1999, "BuPu", "1990s")
logs_thousands_depth_ptrlw <- depth_contours(subset_petrale_logbook, 2009, "BuPu", "2000s")
logs_tens_depth_ptrlw <- depth_contours(subset_petrale_logbook, 2017, "BuPu", "2010s")

windows()
grid.arrange(logs_eighties_depth_ptrlw,
             logs_nineties_depth_ptrlw,
             logs_thousands_depth_ptrlw,
             logs_tens_depth_ptrlw,
             ncol = 2,
             top = textGrob("Petrale Sole Logbook Depth Distribution",
                            gp = gpar(fontfamily = "serif",
                                      cex = 1,
                                      fontface = "bold")))

### Survey Maps ----
subset_petrale_survey <- OR_fish[OR_fish$scientific_name == 'Eopsetta jordani', ]
subset_petrale_survey <- subset_petrale_survey[subset_petrale_survey$catch <= 300, ] # filter out outlier large hauls

# Create data exploration map
windows(width = 28, height = 18)
par(mfrow = c(1, 4))
exploration_panels(subset_petrale_survey, bathy_dat, bathy_mat, 1980, 1989, "Petrale Sole 1980s")
exploration_panels(subset_petrale_survey, bathy_dat, bathy_mat, 1990, 1999, "Petrale Sole 1990s")
exploration_panels(subset_petrale_survey, bathy_dat, bathy_mat, 2000, 2009, "Petrale Sole 2000s")
exploration_panels(subset_petrale_survey, bathy_dat, bathy_mat, 2010, 2018, "Petrale Sole 2010s")

## Decade maps
# Make sure grid has been created (above)
dev.new(width = 4, height = 10)
cpue_fillpts(subset_petrale_survey, 1980, 1989)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_petrale_survey, 1990, 1999)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_petrale_survey, 2000, 2009)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_petrale_survey, 2010, 2018)

## Decade data
# Make sure grid chunk has been run
eighties_surveys_ptrl <- cpue_grid(subset_petrale_survey, 1980, 1989)
nineties_surveys_ptrl <- cpue_grid(subset_petrale_survey, 1990, 1999)
thousands_surveys_ptrl <- cpue_grid(subset_petrale_survey, 2000, 2009)
tens_surveys_ptrl <- cpue_grid(subset_petrale_survey, 2010, 2018)

## Make maps
# Change legend depending on maximum lncpue - use matrix with highest CPUE
max(eighties_surveys_ptrl, na.rm = T)
max(nineties_surveys_ptrl, na.rm = T)
max(thousands_surveys_ptrl, na.rm = T)
max(tens_surveys_ptrl, na.rm = T)

windows(width = 15, height = 9)
par(mfrow = c(1, 4),
    family = 'serif',
    mar = c(4, 5, 3, .3) + .1)
cpue_map(eighties_surveys_ptrl, eighties_surveys_ptrl, "PurpOr", "Survey Petrale Sole 1980s", bathy_dat, bathy_mat)
cpue_map(nineties_surveys_ptrl, eighties_surveys_ptrl, "PurpOr", "Survey Petrale Sole 1990s", bathy_dat, bathy_mat)
cpue_map(thousands_surveys_ptrl, eighties_surveys_ptrl, "PurpOr", "Survey Petrale Sole 2000s", bathy_dat, bathy_mat)
cpue_map(tens_surveys_ptrl, eighties_surveys_ptrl, "PurpOr", "Survey Petrale Sole 2010s", bathy_dat, bathy_mat)

## Depth distribution by year
survey_eighties_depth_ptrl <- depth_contours(subset_petrale_survey, 1989, "BuPu", "1980s")
survey_nineties_depth_ptrl <- depth_contours(subset_petrale_survey, 1999, "BuPu", "1990s")
survey_thousands_depth_ptrl <- depth_contours(subset_petrale_survey, 2009, "BuPu", "2000s")
survey_tens_depth_ptrl <- depth_contours(subset_petrale_survey, 2018, "BuPu", "2010s")

windows()
grid.arrange(survey_eighties_depth_ptrl,
             survey_nineties_depth_ptrl,
             survey_thousands_depth_ptrl,
             survey_tens_depth_ptrl,
             ncol = 2,
             top = textGrob("Petrale Sole Survey Depth Distribution",
                            gp = gpar(fontfamily = "serif",
                                      cex = 1,
                                      fontface = "bold")))

## Average CPUE over time
survey_petrale_cpue <- subset_petrale_survey %>% group_by(year) %>%
  summarise(cpue_mean = mean(CPUE))

# Plot
windows(width = 20, height = 10)
par(mfrow = c(1, 1))
ggplot(data = survey_petrale_cpue, aes(x = year, y = cpue_mean)) +
  geom_path() +
  geom_point() +
  theme_tufte() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "Year",
       y = "Cpue (kg/hr)",
       title = "Mean CPUE of Nearshore Petrale Sole Caught in Groundfish Fishery")

## Manuscript Maps
# Manuscript Figure
pdf("../results/logbook_survey_maps/petrale_sole_maps.pdf",
    width = 7.5,
    height = 18)
par(mfrow = c(2, 2),
    family = "serif",
    mar = c(4, 5, 3, .3) + .1)
species_grid_pdf(1980, 1990, eighties_logbooks_ptrl,
                 tens_logbooks_ptrl, "Logbook Petrale Sole 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "E", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.76, .82, .09, .24),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.6),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(tens_logbooks_ptrl, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.4))
species_grid_pdf(2009, 2018, tens_logbooks_ptrl,
         tens_logbooks_ptrl, "Logbook Petrale Sole 2010s",
         bathy_dat, bathy_mat)
species_grid_pdf(1980, 1990, eighties_surveys_ptrl,
         eighties_surveys_ptrl, "Survey Petrale Sole 1980s",
         bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "E", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.76, .82, .09, .24),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.6),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, 4.5),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.4))
species_grid_pdf(2009, 2018, tens_surveys_ptrl,
         eighties_surveys_ptrl, "Survey Petrale Sole 2010s", bathy_dat, bathy_mat)
dev.off()
