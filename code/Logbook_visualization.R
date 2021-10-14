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
source("functions/species_grid_sup.R")

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
OR_fish <- select(OR_fish, scientific_name, catch, CPUE, lncpue, year, trawl_id)
OR_fish$pres <- 1 * (OR_fish$catch > 0)

names(survey_data)[names(survey_data) == "longitude"] <- "lon"
names(survey_data)[names(survey_data) == "latitude"] <- "lat"
names(survey_data)[names(survey_data) == "depth_m"] <- "depth"
survey_data$depth <- -abs(survey_data$depth)

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
max(tens_logbooks_ptrl, na.rm = T) # max

windows(width = 15, height = 9)
par(mfrow = c(1, 4),
    family = 'serif',
    mar = c(4, 5, 3, .3) + .1)
cpue_map(eighties_logbooks_ptrl, tens_logbooks_ptrl, "PurpOr", "Logbook Petrale Sole 1980s", bathy_dat, bathy_mat)
cpue_map(nineties_logbooks_ptrl, tens_logbooks_ptrl, "PurpOr", "Logbook Petrale Sole 1990s", bathy_dat, bathy_mat)
cpue_map(thousands_logbooks_ptrl, tens_logbooks_ptrl, "PurpOr", "Logbook Petrale Sole 2000s", bathy_dat, bathy_mat)
cpue_map(tens_logbooks_ptrl, tens_logbooks_ptrl, "PurpOr", "Logbook Petrale Sole 2010s", bathy_dat, bathy_mat)
dev.copy(tiff, "../results/visualization/petrale_sole_fourpanel_logs.tiff",
         width = 15, height = 9, units = "in", res = 200)
dev.off()

## Depth distribution by year
logs_eighties_depth_ptrl <- depth_contours(subset_petrale_logbook, 1989, "1980s")
logs_nineties_depth_ptrl <- depth_contours(subset_petrale_logbook, 1999, "1990s")
logs_thousands_depth_ptrl <- depth_contours(subset_petrale_logbook, 2009, "2000s")
logs_tens_depth_ptrl <- depth_contours(subset_petrale_logbook, 2017, "2010s")

windows()
grid.arrange(logs_eighties_depth_ptrl,
             logs_nineties_depth_ptrl,
             logs_thousands_depth_ptrl,
             logs_tens_depth_ptrl,
             ncol = 2,
             top = textGrob("Petrale Sole Logbook Depth Distribution",
                            gp = gpar(fontfamily = "serif",
                                      cex = 0.8,
                                      fontface = "bold")),
             bottom = textGrob("Depth (m)",
                               gp = gpar(fontfamily = "serif",
                                         cex = 0.8,
                                         fontface = "bold")),
             left = textGrob("Latitude",
                             rot = 90,
                             gp = gpar(fontfamily = "serif",
                                       cex = 0.8,
                                       fontface = "bold",
                                       vjust = 1)))
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/petrale_sole_depth.tiff",
         width = 4, height = 4, units = "in", res = 200)
dev.off()

## Catch over time
catch_petrale <- tickets_final %>% filter(common_name == 'Petrale Sole') %>%
  group_by(YEAR) %>%
  summarise(species_weight_sum = sum(TIK_LBS))

# Plot
windows(width = 200, height = 100)
ggplot(data = catch_petrale,
       aes(x = YEAR, y = species_weight_sum / 1000)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
  scale_x_continuous(breaks = c(1980, 1990, 2000, 2010)) +
  geom_col(fill = "orangered4") +
  theme_tufte() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10)) +
  labs(x = "",
       y = "Total Catch (1000s of lbs)",
       title = "")
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/petrale_sole_catch.tiff",
         width = 4, height = 2, units = "in", res = 200)
dev.off()

## Average CPUE over time
logbook_petrale_cpue <- subset_petrale_logbook %>% group_by(year) %>%
  summarise(cpue_mean = mean(CPUE))

# Plot
windows(width = 20, height = 10)
par(mfrow = c(1, 1))
ggplot(data = logbook_petrale_cpue, aes(x = year, y = cpue_mean)) +
  geom_path(color = "orangered4", size = 0.5) +
  geom_point(size = 0.9) +
  scale_x_continuous(breaks = c(1980, 1990, 2000, 2010)) +
  theme_tufte() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 8),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10)) +
  labs(x = "",
       y = "CPUE (kg/hr)",
       title = "")
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/petrale_sole_change.tiff",
         width = 4, height = 2, units = "in", res = 200)
dev.off()

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
max(tens_logbooks_ptrl_winter, na.rm = T) # max

windows(width = 15, height = 9)
par(mfrow = c(1, 4),
    family = 'serif',
    mar = c(4, 5, 3, .3) + .1)
cpue_map(eighties_logbooks_ptrl_winter, tens_logbooks_ptrl_winter, "PurpOr", "Logbook Petrale Sole 1980s", bathy_dat, bathy_mat)
cpue_map(nineties_logbooks_ptrl_winter, tens_logbooks_ptrl_winter, "PurpOr", "Logbook Petrale Sole 1990s", bathy_dat, bathy_mat)
cpue_map(thousands_logbooks_ptrl_winter, tens_logbooks_ptrl_winter, "PurpOr", "Logbook Petrale Sole 2000s", bathy_dat, bathy_mat)
cpue_map(tens_logbooks_ptrl_winter, tens_logbooks_ptrl_winter, "PurpOr", "Logbook Petrale Sole 2010s", bathy_dat, bathy_mat)
dev.copy(tiff, "../results/visualization/petrale_sole_fourpanel_winter.tiff",
         width = 15, height = 9, units = "in", res = 200)
dev.off()

## Depth distribution by year (winter)
logs_eighties_depth_ptrlw <- depth_contours(winter_petrale, 1989, "1980s")
logs_nineties_depth_ptrlw <- depth_contours(winter_petrale, 1999, "1990s")
logs_thousands_depth_ptrlw <- depth_contours(winter_petrale, 2009, "2000s")
logs_tens_depth_ptrlw <- depth_contours(winter_petrale, 2017, "2010s")

windows()
grid.arrange(logs_eighties_depth_ptrlw,
             logs_nineties_depth_ptrlw,
             logs_thousands_depth_ptrlw,
             logs_tens_depth_ptrlw,
             ncol = 2,
             top = textGrob("Petrale Sole Logbook Winter Depth Distribution",
                            gp = gpar(fontfamily = "serif",
                                      cex = 0.8,
                                      fontface = "bold")),
             bottom = textGrob("Depth (m)",
                               gp = gpar(fontfamily = "serif",
                                         cex = 0.8,
                                         fontface = "bold")),
             left = textGrob("Latitude",
                             rot = 90,
                             gp = gpar(fontfamily = "serif",
                                       cex = 0.8,
                                       fontface = "bold",
                                       vjust = 1)))
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/petrale_sole_depth_winter.tiff",
         width = 4, height = 4, units = "in", res = 200)
dev.off()

### Survey ----
subset_petrale <- OR_fish[OR_fish$scientific_name == 'Eopsetta jordani', ]
subset_petrale <- subset_petrale[subset_petrale$catch <= 300, ] # filter out outlier large hauls
match_id <- match(survey_data$trawl_id, subset_petrale$trawl_id)
survey_data$lncpue <- subset_petrale$lncpue[match_id]
survey_data$CPUE <- subset_petrale$CPUE[match_id]
survey_data$catch <- subset_petrale$catch[match_id]
survey_data$pres <- subset_petrale$pres[match_id]
subset_petrale_survey <- survey_data
subset_petrale_survey$lncpue[is.na(subset_petrale_survey$lncpue)] <- 0
subset_petrale_survey$CPUE[is.na(subset_petrale_survey$CPUE)] <- 0
subset_petrale_survey$catch[is.na(subset_petrale_survey$catch)] <- 0
subset_petrale_survey$pres[is.na(subset_petrale_survey$pres)] <- 0

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
max(eighties_surveys_ptrl, na.rm = T) # max
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
dev.copy(tiff, "../results/visualization/petrale_sole_fourpanel_survey.tiff",
         width = 15, height = 9, units = "in", res = 200)
dev.off()

## Depth distribution by year
survey_eighties_depth_ptrl <- depth_contours(subset_petrale_survey, 1989, "1980s")
survey_nineties_depth_ptrl <- depth_contours(subset_petrale_survey, 1999, "1990s")
survey_thousands_depth_ptrl <- depth_contours(subset_petrale_survey, 2009, "2000s")
survey_tens_depth_ptrl <- depth_contours(subset_petrale_survey, 2018, "2010s")

windows()
grid.arrange(survey_eighties_depth_ptrl,
             survey_nineties_depth_ptrl,
             survey_thousands_depth_ptrl,
             survey_tens_depth_ptrl,
             ncol = 2,
             top = textGrob("Petrale Sole Survey Depth Distribution",
                            gp = gpar(fontfamily = "serif",
                                      cex = 0.8,
                                      fontface = "bold")),
             bottom = textGrob("Depth (m)",
                               gp = gpar(fontfamily = "serif",
                                         cex = 0.8,
                                         fontface = "bold")),
             left = textGrob("Latitude",
                             rot = 90,
                             gp = gpar(fontfamily = "serif",
                                       cex = 0.8,
                                       fontface = "bold",
                                       vjust = 1)))
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/petrale_sole_depth_survey.tiff",
         width = 4, height = 4, units = "in", res = 200)
dev.off()

## Average CPUE over time
survey_petrale_cpue <- subset_petrale_survey %>% group_by(year) %>%
  summarise(cpue_mean = mean(CPUE))

# Plot
windows(width = 20, height = 10)
par(mfrow = c(1, 1))
ggplot(data = survey_petrale_cpue, aes(x = year, y = cpue_mean)) +
  geom_path(color = "orangered4", size = 0.5) +
  geom_point(size = 0.9) +
  scale_x_continuous(breaks = c(1980, 1990, 2000, 2010)) +
  theme_tufte() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 8),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10)) +
  labs(x = "",
       y = "CPUE (kg/ha)",
       title = "")
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/petrale_change_survey.tiff",
         width = 4, height = 2, units = "in", res = 200)
dev.off()

### Manuscript Maps ----
# Four panel maps
pdf("../final_figs/manuscript2_fig_tables/petrale_sole_maps.pdf",
    width = 7.5,
    height = 18)
par(mfrow = c(2, 2),
    family = "serif",
    mar = c(4, 5, 3, .3) + .1)
species_grid_pdf(eighties_logbooks_ptrl,
                 tens_logbooks_ptrl, "Logbook Petrale Sole 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.76, .82, .09, .24),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.6),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(tens_logbooks_ptrl, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.4))
species_grid_pdf(tens_logbooks_ptrl,
         tens_logbooks_ptrl, "Logbook Petrale Sole 2010s",
         bathy_dat, bathy_mat)
species_grid_pdf(eighties_surveys_ptrl,
         eighties_surveys_ptrl, "Survey Petrale Sole 1980s",
         bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.76, .82, .09, .24),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.6),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(eighties_surveys_ptrl, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.4))
species_grid_pdf(tens_surveys_ptrl,
         eighties_surveys_ptrl, "Survey Petrale Sole 2010s", bathy_dat, bathy_mat)
dev.off()

# Winter maps (logbooks only)
pdf("../final_figs/manuscript2_fig_tables/petrale_sole_winter.pdf",
    width = 7.5,
    height = 18)
par(mfrow = c(2, 2),
    family = "serif",
    mar = c(4, 5, 3, .3) + .1)
species_grid_pdf(eighties_logbooks_ptrl_winter,
                 tens_logbooks_ptrl_winter, "Winter Petrale Sole 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.76, .82, .09, .24),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.6),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(tens_logbooks_ptrl_winter, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.4))
species_grid_pdf(tens_logbooks_ptrl_winter,
                 tens_logbooks_ptrl_winter, "Winter Petrale Sole 2010s",
                 bathy_dat, bathy_mat)
species_grid_pdf(eighties_logbooks_ptrl,
                 tens_logbooks_ptrl_winter, "Summer Petrale Sole 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.76, .82, .09, .24),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.6),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(tens_logbooks_ptrl_winter, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.4))
species_grid_pdf(tens_logbooks_ptrl,
                 tens_logbooks_ptrl_winter, "Summer Petrale Sole 2010s", bathy_dat, bathy_mat)
dev.off()

### Supplement Maps ----
# Eight panel logbooks and surveys
pdf("../final_figs/manuscript2_fig_tables/petrale_sole_supplement.pdf",
    width = 15,
    height = 17)
par(mfrow = c(2, 4),
    family = "serif",
    mar = c(5.5, 6, 3, .3) + .1)
species_grid_sup(eighties_logbooks_ptrl,
                 tens_logbooks_ptrl,
                 "Logbook Petrale Sole 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.79, .85, .11, .25),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.8),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(tens_logbooks_ptrl, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.2))
species_grid_sup(nineties_logbooks_ptrl,
                 tens_logbooks_ptrl,
                 "Logbook Petrale Sole 1990s",
                 bathy_dat, bathy_mat)
species_grid_sup(thousands_logbooks_ptrl,
                 tens_logbooks_ptrl,
                 "Logbook Petrale Sole 2000s",
                 bathy_dat, bathy_mat)
species_grid_sup(tens_logbooks_ptrl,
                 tens_logbooks_ptrl,
                 "Logbook Petrale Sole 2010s",
                 bathy_dat, bathy_mat)
species_grid_sup(eighties_surveys_ptrl,
                 eighties_surveys_ptrl,
                 "Survey Petrale Sole 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.79, .85, .11, .25),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.8),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(eighties_surveys_ptrl, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.2))
species_grid_sup(nineties_surveys_ptrl,
                 eighties_surveys_ptrl,
                 "Survey Petrale Sole 1990s",
                 bathy_dat, bathy_mat)
species_grid_sup(thousands_surveys_ptrl,
                 eighties_surveys_ptrl,
                 "Survey Petrale Sole 2000s",
                 bathy_dat, bathy_mat)
species_grid_sup(tens_surveys_ptrl,
                 eighties_surveys_ptrl,
                 "Survey Petrale Sole 2010s",
                 bathy_dat, bathy_mat)
dev.off()

# Eight panel seasonality
pdf("../final_figs/manuscript2_fig_tables/petrale_seasonal_supplement.pdf",
    width = 15,
    height = 17)
par(mfrow = c(2, 4),
    family = "serif",
    mar = c(5.5, 6, 3, .3) + .1)
species_grid_sup(eighties_logbooks_ptrl_winter,
                 tens_logbooks_ptrl_winter,
                 "Winter Petrale Sole 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.79, .85, .11, .25),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.8),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(tens_logbooks_ptrl_winter, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.2))
species_grid_sup(nineties_logbooks_ptrl_winter,
                 tens_logbooks_ptrl_winter,
                 "Winter Petrale Sole 1990s",
                 bathy_dat, bathy_mat)
species_grid_sup(thousands_logbooks_ptrl_winter,
                 tens_logbooks_ptrl_winter,
                 "Winter Petrale Sole 2000s",
                 bathy_dat, bathy_mat)
species_grid_sup(tens_logbooks_ptrl_winter,
                 tens_logbooks_ptrl_winter,
                 "Winter Petrale Sole 2010s",
                 bathy_dat, bathy_mat)
species_grid_sup(eighties_logbooks_ptrl,
                 tens_logbooks_ptrl_winter,
                 "Summer Petrale Sole 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.79, .85, .11, .25),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.8),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(tens_logbooks_ptrl_winter, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.2))
species_grid_sup(nineties_logbooks_ptrl,
                 tens_logbooks_ptrl_winter,
                 "Summer Petrale Sole 1990s",
                 bathy_dat, bathy_mat)
species_grid_sup(thousands_logbooks_ptrl,
                 tens_logbooks_ptrl_winter,
                 "Summer Petrale Sole 2000s",
                 bathy_dat, bathy_mat)
species_grid_sup(tens_logbooks_ptrl,
                 tens_logbooks_ptrl_winter,
                 "Summer Petrale Sole 2010s",
                 bathy_dat, bathy_mat)
dev.off()


## Dover Sole ----
### Logbooks ----
# Filter to just dover and survey months
subset_dover_logbook <- logbooks_final[logbooks_final$species == 'DOVR_ADJ', ]
subset_dover_logbook <- subset_dover_logbook[subset_dover_logbook$month_day >= 517 &
                                                   subset_dover_logbook$month_day <= 929, ]

# Create data exploration map
windows(width = 28, height = 18)
par(mfrow = c(1, 4))
exploration_panels(subset_dover_logbook, bathy_dat, bathy_mat, 1981, 1989, "Dover Sole 1980s")
exploration_panels(subset_dover_logbook, bathy_dat, bathy_mat, 1990, 1999, "Dover Sole 1990s")
exploration_panels(subset_dover_logbook, bathy_dat, bathy_mat, 2000, 2009, "Dover Sole 2000s")
exploration_panels(subset_dover_logbook, bathy_dat, bathy_mat, 2010, 2017, "Dover Sole 2010s")

## Decade maps
# Make sure grid has been created (above)
dev.new(width = 4, height = 10)
cpue_fillpts(subset_dover_logbook, 1981, 1989)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_dover_logbook, 1990, 1999)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_dover_logbook, 2000, 2009)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_dover_logbook, 2010, 2017)

## Decade data
# Make sure grid chunk has been run
eighties_logbooks_dovr <- cpue_grid(subset_dover_logbook, 1981, 1989)
nineties_logbooks_dovr <- cpue_grid(subset_dover_logbook, 1990, 1999)
thousands_logbooks_dovr <- cpue_grid(subset_dover_logbook, 2000, 2009)
tens_logbooks_dovr <- cpue_grid(subset_dover_logbook, 2010, 2017)

## Make maps
# Change legend depending on maximum lncpue - use matrix with highest CPUE
max(eighties_logbooks_dovr, na.rm = T)
max(nineties_logbooks_dovr, na.rm = T)
max(thousands_logbooks_dovr, na.rm = T)
max(tens_logbooks_dovr, na.rm = T) # max

windows(width = 15, height = 9)
par(mfrow = c(1, 4),
    family = 'serif',
    mar = c(4, 5, 3, .3) + .1)
cpue_map(eighties_logbooks_dovr, tens_logbooks_dovr, "PurpOr", "Logbook Dover Sole 1980s", bathy_dat, bathy_mat)
cpue_map(nineties_logbooks_dovr, tens_logbooks_dovr, "PurpOr", "Logbook Dover Sole 1990s", bathy_dat, bathy_mat)
cpue_map(thousands_logbooks_dovr, tens_logbooks_dovr, "PurpOr", "Logbook Dover Sole 2000s", bathy_dat, bathy_mat)
cpue_map(tens_logbooks_dovr, tens_logbooks_dovr, "PurpOr", "Logbook Dover Sole 2010s", bathy_dat, bathy_mat)
dev.copy(tiff, "../results/visualization/dover_sole_fourpanel_logs.tiff",
         width = 15, height = 9, units = "in", res = 200)
dev.off()

## Depth distribution by year
logs_eighties_depth_dovr <- depth_contours(subset_dover_logbook, 1989, "1980s")
logs_nineties_depth_dovr <- depth_contours(subset_dover_logbook, 1999, "1990s")
logs_thousands_depth_dovr <- depth_contours(subset_dover_logbook, 2009, "2000s")
logs_tens_depth_dovr <- depth_contours(subset_dover_logbook, 2017, "2010s")

windows()
grid.arrange(logs_eighties_depth_dovr,
             logs_nineties_depth_dovr,
             logs_thousands_depth_dovr,
             logs_tens_depth_dovr,
             ncol = 2,
             top = textGrob("Dover Sole Logbook Depth Distribution",
                            gp = gpar(fontfamily = "serif",
                                      cex = 0.8,
                                      fontface = "bold")),
             bottom = textGrob("Depth (m)",
                               gp = gpar(fontfamily = "serif",
                                         cex = 0.8,
                                         fontface = "bold")),
             left = textGrob("Latitude",
                             rot = 90,
                             gp = gpar(fontfamily = "serif",
                                       cex = 0.8,
                                       fontface = "bold",
                                       vjust = 1)))
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/dover_sole_depth.tiff",
         width = 4, height = 4, units = "in", res = 200)
dev.off()

## Catch over time
catch_dover <- tickets_final %>% filter(common_name == 'Dover Sole') %>%
  group_by(YEAR) %>%
  summarise(species_weight_sum = sum(TIK_LBS))

# Plot
windows(width = 200, height = 100)
ggplot(data = catch_dover,
       aes(x = YEAR, y = species_weight_sum / 1000)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
  scale_x_continuous(breaks = c(1980, 1990, 2000, 2010)) +
  geom_col(fill = "orangered4") +
  theme_tufte() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10)) +
  labs(x = "",
       y = "Total Catch (1000s of lbs)",
       title = "")
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/dover_sole_catch.tiff",
         width = 4, height = 2, units = "in", res = 200)
dev.off()

## Average CPUE over time
logbook_dover_cpue <- subset_dover_logbook %>% group_by(year) %>%
  summarise(cpue_mean = mean(CPUE))

# Plot
windows(width = 20, height = 10)
par(mfrow = c(1, 1))
ggplot(data = logbook_dover_cpue, aes(x = year, y = cpue_mean)) +
  geom_path(color = "orangered4", size = 0.5) +
  geom_point(size = 0.9) +
  scale_x_continuous(breaks = c(1980, 1990, 2000, 2010)) +
  theme_tufte() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 8),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10)) +
  labs(x = "",
       y = "CPUE (kg/hr)",
       title = "")
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/dover_sole_change.tiff",
         width = 4, height = 2, units = "in", res = 200)
dev.off()

## Seasonality
# Create winter subset
winter_dover <- filter(logbooks_final, !between(month_day, 516, 930),
                         species == 'DOVR_ADJ')

## Winter Decades
dev.new(width = 4, height = 10)
cpue_fillpts(winter_dover, 1981, 1989)

dev.new(width = 4, height = 10)
cpue_fillpts(winter_dover, 1990, 1999)

dev.new(width = 4, height = 10)
cpue_fillpts(winter_dover, 2000, 2009)

dev.new(width = 4, height = 10)
cpue_fillpts(winter_dover, 2010, 2017)

## Decade data
eighties_logbooks_dovr_winter <- cpue_grid(winter_dover, 1981, 1989)
nineties_logbooks_dovr_winter <- cpue_grid(winter_dover, 1990, 1999)
thousands_logbooks_dovr_winter <- cpue_grid(winter_dover, 2000, 2009)
tens_logbooks_dovr_winter <- cpue_grid(winter_dover, 2010, 2017)

## Make maps
# Change legend depending on maximum lncpue - use matrix with highest CPUE
max(eighties_logbooks_dovr_winter, na.rm = T)
max(nineties_logbooks_dovr_winter, na.rm = T)
max(thousands_logbooks_dovr_winter, na.rm = T)
max(tens_logbooks_dovr_winter, na.rm = T) # max

windows(width = 15, height = 9)
par(mfrow = c(1, 4),
    family = 'serif',
    mar = c(4, 5, 3, .3) + .1)
cpue_map(eighties_logbooks_dovr_winter, tens_logbooks_dovr_winter, "PurpOr", "Logbook Dover Sole 1980s", bathy_dat, bathy_mat)
cpue_map(nineties_logbooks_dovr_winter, tens_logbooks_dovr_winter, "PurpOr", "Logbook Dover Sole 1990s", bathy_dat, bathy_mat)
cpue_map(thousands_logbooks_dovr_winter, tens_logbooks_dovr_winter, "PurpOr", "Logbook Dover Sole 2000s", bathy_dat, bathy_mat)
cpue_map(tens_logbooks_dovr_winter, tens_logbooks_dovr_winter, "PurpOr", "Logbook Dover Sole 2010s", bathy_dat, bathy_mat)
dev.copy(tiff, "../results/visualization/dover_sole_fourpanel_winter.tiff",
         width = 15, height = 9, units = "in", res = 200)
dev.off()

## Depth distribution by year (winter)
logs_eighties_depth_dovrw <- depth_contours(winter_dover, 1989, "1980s")
logs_nineties_depth_dovrw <- depth_contours(winter_dover, 1999, "1990s")
logs_thousands_depth_dovrw <- depth_contours(winter_dover, 2009, "2000s")
logs_tens_depth_dovrw <- depth_contours(winter_dover, 2017, "2010s")

windows()
grid.arrange(logs_eighties_depth_dovrw,
             logs_nineties_depth_dovrw,
             logs_thousands_depth_dovrw,
             logs_tens_depth_dovrw,
             ncol = 2,
             top = textGrob("Dover Sole Logbook Winter Depth Distribution",
                            gp = gpar(fontfamily = "serif",
                                      cex = 0.8,
                                      fontface = "bold")),
             bottom = textGrob("Depth (m)",
                               gp = gpar(fontfamily = "serif",
                                         cex = 0.8,
                                         fontface = "bold")),
             left = textGrob("Latitude",
                             rot = 90,
                             gp = gpar(fontfamily = "serif",
                                       cex = 0.8,
                                       fontface = "bold",
                                       vjust = 1)))
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/dover_sole_depth_winter.tiff",
         width = 4, height = 4, units = "in", res = 200)
dev.off()

### Survey ----
subset_dover <- OR_fish[OR_fish$scientific_name == 'Microstomus pacificus', ]
match_id <- match(survey_data$trawl_id, subset_dover$trawl_id)
survey_data$lncpue <- subset_dover$lncpue[match_id]
survey_data$CPUE <- subset_dover$CPUE[match_id]
survey_data$catch <- subset_dover$catch[match_id]
survey_data$pres <- subset_dover$pres[match_id]
subset_dover_survey <- survey_data
subset_dover_survey$lncpue[is.na(subset_dover_survey$lncpue)] <- 0
subset_dover_survey$CPUE[is.na(subset_dover_survey$CPUE)] <- 0
subset_dover_survey$catch[is.na(subset_dover_survey$catch)] <- 0
subset_dover_survey$pres[is.na(subset_dover_survey$pres)] <- 0

# Create data exploration map
windows(width = 28, height = 18)
par(mfrow = c(1, 4))
exploration_panels(subset_dover_survey, bathy_dat, bathy_mat, 1980, 1989, "Dover Sole 1980s")
exploration_panels(subset_dover_survey, bathy_dat, bathy_mat, 1990, 1999, "Dover Sole 1990s")
exploration_panels(subset_dover_survey, bathy_dat, bathy_mat, 2000, 2009, "Dover Sole 2000s")
exploration_panels(subset_dover_survey, bathy_dat, bathy_mat, 2010, 2018, "Dover Sole 2010s")

## Decade maps
# Make sure grid has been created (above)
dev.new(width = 4, height = 10)
cpue_fillpts(subset_dover_survey, 1980, 1989)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_dover_survey, 1990, 1999)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_dover_survey, 2000, 2009)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_dover_survey, 2010, 2018)

## Decade data
# Make sure grid chunk has been run
eighties_surveys_dovr <- cpue_grid(subset_dover_survey, 1980, 1989)
nineties_surveys_dovr <- cpue_grid(subset_dover_survey, 1990, 1999)
thousands_surveys_dovr <- cpue_grid(subset_dover_survey, 2000, 2009)
tens_surveys_dovr <- cpue_grid(subset_dover_survey, 2010, 2018)

## Make maps
# Change legend depending on maximum lncpue - use matrix with highest CPUE
max(eighties_surveys_dovr, na.rm = T)
max(nineties_surveys_dovr, na.rm = T)
max(thousands_surveys_dovr, na.rm = T) # max
max(tens_surveys_dovr, na.rm = T)

windows(width = 15, height = 9)
par(mfrow = c(1, 4),
    family = 'serif',
    mar = c(4, 5, 3, .3) + .1)
cpue_map(eighties_surveys_dovr, thousands_surveys_dovr, "PurpOr", "Survey Dover Sole 1980s", bathy_dat, bathy_mat)
cpue_map(nineties_surveys_dovr, thousands_surveys_dovr, "PurpOr", "Survey Dover Sole 1990s", bathy_dat, bathy_mat)
cpue_map(thousands_surveys_dovr, thousands_surveys_dovr, "PurpOr", "Survey Dover Sole 2000s", bathy_dat, bathy_mat)
cpue_map(tens_surveys_dovr, thousands_surveys_dovr, "PurpOr", "Survey Dover Sole 2010s", bathy_dat, bathy_mat)
dev.copy(tiff, "../results/visualization/dover_sole_fourpanel_survey.tiff",
         width = 15, height = 9, units = "in", res = 200)
dev.off()

## Depth distribution by year
survey_eighties_depth_dovr <- depth_contours(subset_dover_survey, 1989, "1980s")
survey_nineties_depth_dovr <- depth_contours(subset_dover_survey, 1999, "1990s")
survey_thousands_depth_dovr <- depth_contours(subset_dover_survey, 2009, "2000s")
survey_tens_depth_dovr <- depth_contours(subset_dover_survey, 2018, "2010s")

windows()
grid.arrange(survey_eighties_depth_dovr,
             survey_nineties_depth_dovr,
             survey_thousands_depth_dovr,
             survey_tens_depth_dovr,
             ncol = 2,
             top = textGrob("Dover Sole Survey Depth Distribution",
                            gp = gpar(fontfamily = "serif",
                                      cex = 0.8,
                                      fontface = "bold")),
             bottom = textGrob("Depth (m)",
                               gp = gpar(fontfamily = "serif",
                                         cex = 0.8,
                                         fontface = "bold")),
             left = textGrob("Latitude",
                             rot = 90,
                             gp = gpar(fontfamily = "serif",
                                       cex = 0.8,
                                       fontface = "bold",
                                       vjust = 1)))
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/dover_sole_depth_survey.tiff",
         width = 4, height = 4, units = "in", res = 200)
dev.off()

## Average CPUE over time
survey_dover_cpue <- subset_dover_survey %>% group_by(year) %>%
  summarise(cpue_mean = mean(CPUE))

# Plot
windows(width = 20, height = 10)
par(mfrow = c(1, 1))
ggplot(data = survey_dover_cpue, aes(x = year, y = cpue_mean)) +
  geom_path(color = "orangered4", size = 0.5) +
  geom_point(size = 0.9) +
  scale_x_continuous(breaks = c(1980, 1990, 2000, 2010)) +
  theme_tufte() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 8),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10)) +
  labs(x = "",
       y = "CPUE (kg/ha)",
       title = "")
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/dover_change_survey.tiff",
         width = 4, height = 2, units = "in", res = 200)
dev.off()

### Manuscript Maps ----
# Four panel maps
pdf("../final_figs/manuscript2_fig_tables/dover_sole_maps.pdf",
    width = 7.5,
    height = 18)
par(mfrow = c(2, 2),
    family = "serif",
    mar = c(4, 5, 3, .3) + .1)
species_grid_pdf(eighties_logbooks_dovr,
                 tens_logbooks_dovr, "Logbook Dover Sole 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.76, .82, .09, .24),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.6),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(tens_logbooks_dovr, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.4))
species_grid_pdf(tens_logbooks_dovr,
                 tens_logbooks_dovr, "Logbook Dover Sole 2010s",
                 bathy_dat, bathy_mat)
species_grid_pdf(eighties_surveys_dovr,
                 thousands_surveys_dovr, "Survey Dover Sole 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.76, .82, .09, .24),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.6),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(thousands_surveys_dovr, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.4))
species_grid_pdf(tens_surveys_dovr,
                 thousands_surveys_dovr, "Survey Dover Sole 2010s", bathy_dat, bathy_mat)
dev.off()

# Winter maps (logbooks only)
pdf("../final_figs/manuscript2_fig_tables/dover_sole_winter.pdf",
    width = 7.5,
    height = 18)
par(mfrow = c(2, 2),
    family = "serif",
    mar = c(4, 5, 3, .3) + .1)
species_grid_pdf(eighties_logbooks_dovr_winter,
                 tens_logbooks_dovr_winter, "Winter Dover Sole 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.76, .82, .09, .24),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.6),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(tens_logbooks_dovr_winter, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.4))
species_grid_pdf(tens_logbooks_dovr_winter,
                 tens_logbooks_dovr_winter, "Winter Dover Sole 2010s",
                 bathy_dat, bathy_mat)
species_grid_pdf(eighties_logbooks_dovr,
                 tens_logbooks_dovr_winter, "Summer Dover Sole 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.76, .82, .09, .24),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.6),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(tens_logbooks_dovr_winter, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.4))
species_grid_pdf(tens_logbooks_dovr,
                 tens_logbooks_dovr_winter, "Summer Dover Sole 2010s", bathy_dat, bathy_mat)
dev.off()

### Supplement Maps ----
# Eight panel logbooks and surveys
pdf("../final_figs/manuscript2_fig_tables/dover_sole_supplement.pdf",
    width = 15,
    height = 17)
par(mfrow = c(2, 4),
    family = "serif",
    mar = c(5.5, 6, 3, .3) + .1)
species_grid_sup(eighties_logbooks_dovr,
                 tens_logbooks_dovr,
                 "Logbook Dover Sole 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.79, .85, .11, .25),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.8),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(tens_logbooks_dovr, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.2))
species_grid_sup(nineties_logbooks_dovr,
                 tens_logbooks_dovr,
                 "Logbook Dover Sole 1990s",
                 bathy_dat, bathy_mat)
species_grid_sup(thousands_logbooks_dovr,
                 tens_logbooks_dovr,
                 "Logbook Dover Sole 2000s",
                 bathy_dat, bathy_mat)
species_grid_sup(tens_logbooks_dovr,
                 tens_logbooks_dovr,
                 "Logbook Dover Sole 2010s",
                 bathy_dat, bathy_mat)
species_grid_sup(eighties_surveys_dovr,
                 thousands_surveys_dovr,
                 "Survey Dover Sole 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.79, .85, .11, .25),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.8),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(thousands_surveys_dovr, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.2))
species_grid_sup(nineties_surveys_dovr,
                 thousands_surveys_dovr,
                 "Survey Dover Sole 1990s",
                 bathy_dat, bathy_mat)
species_grid_sup(thousands_surveys_dovr,
                 thousands_surveys_dovr,
                 "Survey Dover Sole 2000s",
                 bathy_dat, bathy_mat)
species_grid_sup(tens_surveys_dovr,
                 thousands_surveys_dovr,
                 "Survey Dover Sole 2010s",
                 bathy_dat, bathy_mat)
dev.off()

# Eight panel seasonality
pdf("../final_figs/manuscript2_fig_tables/dover_seasonal_supplement.pdf",
    width = 15,
    height = 17)
par(mfrow = c(2, 4),
    family = "serif",
    mar = c(5.5, 6, 3, .3) + .1)
species_grid_sup(eighties_logbooks_dovr_winter,
                 tens_logbooks_dovr_winter,
                 "Winter Dover Sole 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.79, .85, .11, .25),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.8),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(tens_logbooks_dovr_winter, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.2))
species_grid_sup(nineties_logbooks_dovr_winter,
                 tens_logbooks_dovr_winter,
                 "Winter Dover Sole 1990s",
                 bathy_dat, bathy_mat)
species_grid_sup(thousands_logbooks_dovr_winter,
                 tens_logbooks_dovr_winter,
                 "Winter Dover Sole 2000s",
                 bathy_dat, bathy_mat)
species_grid_sup(tens_logbooks_dovr_winter,
                 tens_logbooks_dovr_winter,
                 "Winter Dover Sole 2010s",
                 bathy_dat, bathy_mat)
species_grid_sup(eighties_logbooks_dovr,
                 tens_logbooks_dovr_winter,
                 "Summer Dover Sole 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.79, .85, .11, .25),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.8),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(tens_logbooks_dovr_winter, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.2))
species_grid_sup(nineties_logbooks_dovr,
                 tens_logbooks_dovr_winter,
                 "Summer Dover Sole 1990s",
                 bathy_dat, bathy_mat)
species_grid_sup(thousands_logbooks_dovr,
                 tens_logbooks_dovr_winter,
                 "Summer Dover Sole 2000s",
                 bathy_dat, bathy_mat)
species_grid_sup(tens_logbooks_dovr,
                 tens_logbooks_dovr_winter,
                 "Summer Dover Sole 2010s",
                 bathy_dat, bathy_mat)
dev.off()

## Pacific Sanddab ----
### Logbooks ----
# Filter to just sanddab and survey months
subset_sanddab_logbook <- logbooks_final[logbooks_final$species == 'SDAB_ADJ', ]
subset_sanddab_logbook <- subset_sanddab_logbook[subset_sanddab_logbook$month_day >= 517 &
                                               subset_sanddab_logbook$month_day <= 929, ]

# Create data exploration map
windows(width = 28, height = 18)
par(mfrow = c(1, 4))
exploration_panels(subset_sanddab_logbook, bathy_dat, bathy_mat, 1981, 1989, "Pacific Sanddab 1980s")
exploration_panels(subset_sanddab_logbook, bathy_dat, bathy_mat, 1990, 1999, "Pacific Sanddab 1990s")
exploration_panels(subset_sanddab_logbook, bathy_dat, bathy_mat, 2000, 2009, "Pacific Sanddab 2000s")
exploration_panels(subset_sanddab_logbook, bathy_dat, bathy_mat, 2010, 2017, "Pacific Sanddab 2010s")

## Decade maps
# Make sure grid has been created (above)
dev.new(width = 4, height = 10)
cpue_fillpts(subset_sanddab_logbook, 1981, 1989)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_sanddab_logbook, 1990, 1999)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_sanddab_logbook, 2000, 2009)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_sanddab_logbook, 2010, 2017)

## Decade data
# Make sure grid chunk has been run
eighties_logbooks_sdab <- cpue_grid(subset_sanddab_logbook, 1981, 1989)
nineties_logbooks_sdab <- cpue_grid(subset_sanddab_logbook, 1990, 1999)
thousands_logbooks_sdab <- cpue_grid(subset_sanddab_logbook, 2000, 2009)
tens_logbooks_sdab <- cpue_grid(subset_sanddab_logbook, 2010, 2017)

## Make maps
# Change legend depending on maximum lncpue - use matrix with highest CPUE
max(eighties_logbooks_sdab, na.rm = T)
max(nineties_logbooks_sdab, na.rm = T) # max
max(thousands_logbooks_sdab, na.rm = T)
max(tens_logbooks_sdab, na.rm = T)

windows(width = 15, height = 9)
par(mfrow = c(1, 4),
    family = 'serif',
    mar = c(4, 5, 3, .3) + .1)
cpue_map(eighties_logbooks_sdab, nineties_logbooks_sdab, "PurpOr", "Logbook Pacific Sanddab 1980s", bathy_dat, bathy_mat)
cpue_map(nineties_logbooks_sdab, nineties_logbooks_sdab, "PurpOr", "Logbook Pacific Sanddab 1990s", bathy_dat, bathy_mat)
cpue_map(thousands_logbooks_sdab, nineties_logbooks_sdab, "PurpOr", "Logbook Pacific Sanddab 2000s", bathy_dat, bathy_mat)
cpue_map(tens_logbooks_sdab, nineties_logbooks_sdab, "PurpOr", "Logbook Pacific Sanddab 2010s", bathy_dat, bathy_mat)
dev.copy(tiff, "../results/visualization/sanddab_fourpanel_logs.tiff",
         width = 15, height = 9, units = "in", res = 200)
dev.off()

## Depth distribution by year
logs_eighties_depth_sdab <- depth_contours(subset_sanddab_logbook, 1989, "1980s")
logs_nineties_depth_sdab <- depth_contours(subset_sanddab_logbook, 1999, "1990s")
logs_thousands_depth_sdab <- depth_contours(subset_sanddab_logbook, 2009, "2000s")
logs_tens_depth_sdab <- depth_contours(subset_sanddab_logbook, 2017, "2010s")

windows()
grid.arrange(logs_eighties_depth_sdab,
             logs_nineties_depth_sdab,
             logs_thousands_depth_sdab,
             logs_tens_depth_sdab,
             ncol = 2,
             top = textGrob("Pacific Sanddab Logbook Depth Distribution",
                            gp = gpar(fontfamily = "serif",
                                      cex = 0.8,
                                      fontface = "bold")),
             bottom = textGrob("Depth (m)",
                               gp = gpar(fontfamily = "serif",
                                         cex = 0.8,
                                         fontface = "bold")),
             left = textGrob("Latitude",
                             rot = 90,
                             gp = gpar(fontfamily = "serif",
                                       cex = 0.8,
                                       fontface = "bold",
                                       vjust = 1)))
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/sanddab_depth.tiff",
         width = 4, height = 4, units = "in", res = 200)
dev.off()

## Catch over time
catch_sanddab <- tickets_final %>% filter(common_name == 'Pacific Sanddab') %>%
  group_by(YEAR) %>%
  summarise(species_weight_sum = sum(TIK_LBS))

# Plot
windows(width = 200, height = 100)
ggplot(data = catch_sanddab,
       aes(x = YEAR, y = species_weight_sum / 1000)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
  geom_col(fill = "orangered4") +
  scale_x_continuous(breaks = c(1980, 1990, 2000, 2010)) +
  theme_tufte() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10)) +
  labs(x = "",
       y = "Total Catch (1000s of lbs)",
       title = "")
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/sanddab_catch.tiff",
         width = 4, height = 2, units = "in", res = 200)
dev.off()

## Average CPUE over time
logbook_sanddab_cpue <- subset_sanddab_logbook %>% group_by(year) %>%
  summarise(cpue_mean = mean(CPUE))

# Plot
windows(width = 20, height = 10)
par(mfrow = c(1, 1))
ggplot(data = logbook_sanddab_cpue, aes(x = year, y = cpue_mean)) +
  geom_path(color = "orangered4", size = 0.5) +
  geom_point(size = 0.9) +
  scale_x_continuous(breaks = c(1980, 1990, 2000, 2010)) +
  theme_tufte() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 8),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10)) +
  labs(x = "",
       y = "CPUE (kg/hr)",
       title = "")
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/sanddab_change.tiff",
         width = 4, height = 2, units = "in", res = 200)
dev.off()

### Survey ----
subset_sanddab <- OR_fish[OR_fish$scientific_name == 'Citharichthys sordidus', ]
match_id <- match(survey_data$trawl_id, subset_sanddab$trawl_id)
survey_data$lncpue <- subset_sanddab$lncpue[match_id]
survey_data$CPUE <- subset_sanddab$CPUE[match_id]
survey_data$catch <- subset_sanddab$catch[match_id]
survey_data$pres <- subset_sanddab$pres[match_id]
subset_sanddab_survey <- survey_data
subset_sanddab_survey$lncpue[is.na(subset_sanddab_survey$lncpue)] <- 0
subset_sanddab_survey$CPUE[is.na(subset_sanddab_survey$CPUE)] <- 0
subset_sanddab_survey$catch[is.na(subset_sanddab_survey$catch)] <- 0
subset_sanddab_survey$pres[is.na(subset_sanddab_survey$pres)] <- 0

# Create data exploration map
windows(width = 28, height = 18)
par(mfrow = c(1, 4))
exploration_panels(subset_sanddab_survey, bathy_dat, bathy_mat, 1980, 1989, "Pacific Sanddab 1980s")
exploration_panels(subset_sanddab_survey, bathy_dat, bathy_mat, 1990, 1999, "Pacific Sanddab 1990s")
exploration_panels(subset_sanddab_survey, bathy_dat, bathy_mat, 2000, 2009, "Pacific Sanddab 2000s")
exploration_panels(subset_sanddab_survey, bathy_dat, bathy_mat, 2010, 2018, "Pacific Sanddab 2010s")

## Decade maps
# Make sure grid has been created (above)
dev.new(width = 4, height = 10)
cpue_fillpts(subset_sanddab_survey, 1980, 1989)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_sanddab_survey, 1990, 1999)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_sanddab_survey, 2000, 2009)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_sanddab_survey, 2010, 2018)

## Decade data
# Make sure grid chunk has been run
eighties_surveys_sdab <- cpue_grid(subset_sanddab_survey, 1980, 1989)
nineties_surveys_sdab <- cpue_grid(subset_sanddab_survey, 1990, 1999)
thousands_surveys_sdab <- cpue_grid(subset_sanddab_survey, 2000, 2009)
tens_surveys_sdab <- cpue_grid(subset_sanddab_survey, 2010, 2018)

## Make maps
# Change legend depending on maximum lncpue - use matrix with highest CPUE
max(eighties_surveys_sdab, na.rm = T)
max(nineties_surveys_sdab, na.rm = T)# max
max(thousands_surveys_sdab, na.rm = T)
max(tens_surveys_sdab, na.rm = T)

windows(width = 15, height = 9)
par(mfrow = c(1, 4),
    family = 'serif',
    mar = c(4, 5, 3, .3) + .1)
cpue_map(eighties_surveys_sdab, nineties_surveys_sdab, "PurpOr", "Survey Pacific Sanddab 1980s", bathy_dat, bathy_mat)
cpue_map(nineties_surveys_sdab, nineties_surveys_sdab, "PurpOr", "Survey Pacific Sanddab 1990s", bathy_dat, bathy_mat)
cpue_map(thousands_surveys_sdab, nineties_surveys_sdab, "PurpOr", "Survey Pacific Sanddab 2000s", bathy_dat, bathy_mat)
cpue_map(tens_surveys_sdab, nineties_surveys_sdab, "PurpOr", "Survey Pacific Sanddab 2010s", bathy_dat, bathy_mat)
dev.copy(tiff, "../results/visualization/sanddab_fourpanel_survey.tiff",
         width = 15, height = 9, units = "in", res = 200)
dev.off()

## Depth distribution by year
survey_eighties_depth_sdab <- depth_contours(subset_sanddab_survey, 1989, "1980s")
survey_nineties_depth_sdab <- depth_contours(subset_sanddab_survey, 1999, "1990s")
survey_thousands_depth_sdab <- depth_contours(subset_sanddab_survey, 2009, "2000s")
survey_tens_depth_sdab <- depth_contours(subset_sanddab_survey, 2018, "2010s")

windows()
grid.arrange(survey_eighties_depth_sdab,
             survey_nineties_depth_sdab,
             survey_thousands_depth_sdab,
             survey_tens_depth_sdab,
             ncol = 2,
             top = textGrob("Pacific Sanddab Survey Depth Distribution",
                            gp = gpar(fontfamily = "serif",
                                      cex = 0.8,
                                      fontface = "bold")),
             bottom = textGrob("Depth (m)",
                               gp = gpar(fontfamily = "serif",
                                         cex = 0.8,
                                         fontface = "bold")),
             left = textGrob("Latitude",
                             rot = 90,
                             gp = gpar(fontfamily = "serif",
                                       cex = 0.8,
                                       fontface = "bold",
                                       vjust = 1)))
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/sanddab_depth_survey.tiff",
         width = 4, height = 4, units = "in", res = 200)
dev.off()

## Average CPUE over time
survey_sanddab_cpue <- subset_sanddab_survey %>% group_by(year) %>%
  summarise(cpue_mean = mean(CPUE))

# Plot
windows(width = 20, height = 10)
par(mfrow = c(1, 1))
ggplot(data = survey_sanddab_cpue, aes(x = year, y = cpue_mean)) +
  geom_path(color = "orangered4", size = 0.5) +
  geom_point(size = 0.9) +
  scale_x_continuous(breaks = c(1980, 1990, 2000, 2010)) +
  theme_tufte() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 8),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10)) +
  labs(x = "",
       y = "CPUE (kg/ha)",
       title = "")
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/sanddab_change_survey.tiff",
         width = 4, height = 2, units = "in", res = 200)
dev.off()

### Manuscript Maps ----
# Four panel maps
pdf("../final_figs/manuscript2_fig_tables/sanddab_maps.pdf",
    width = 7.5,
    height = 18)
par(mfrow = c(2, 2),
    family = "serif",
    mar = c(4, 5, 3, .3) + .1)
species_grid_pdf(eighties_logbooks_sdab,
                 nineties_logbooks_sdab, "Logbook Pacific Sanddab \n 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.76, .82, .09, .24),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.6),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(nineties_logbooks_sdab, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.4))
species_grid_pdf(tens_logbooks_sdab,
                 nineties_logbooks_sdab, "Logbook Pacific Sanddab \n 2010s",
                 bathy_dat, bathy_mat)
species_grid_pdf(eighties_surveys_sdab,
                 nineties_surveys_sdab, "Survey Pacific Sanddab \n 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.76, .82, .09, .24),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.6),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(nineties_surveys_sdab, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.4))
species_grid_pdf(tens_surveys_sdab,
                 nineties_surveys_sdab, "Survey Pacific Sanddab \n 2010s", bathy_dat, bathy_mat)
dev.off()

### Supplement Maps ----
# Eight panel logbooks and surveys
pdf("../final_figs/manuscript2_fig_tables/sanddab_supplement.pdf",
    width = 15,
    height = 17)
par(mfrow = c(2, 4),
    family = "serif",
    mar = c(5.5, 6, 3, .3) + .1)
species_grid_sup(eighties_logbooks_sdab,
                 nineties_logbooks_sdab,
                 "Logbook Pacific Sanddab 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.79, .85, .11, .25),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.8),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(nineties_logbooks_sdab, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.2))
species_grid_sup(nineties_logbooks_sdab,
                 nineties_logbooks_sdab,
                 "Logbook Pacific Sanddab 1990s",
                 bathy_dat, bathy_mat)
species_grid_sup(thousands_logbooks_sdab,
                 nineties_logbooks_sdab,
                 "Logbook Pacific Sanddab 2000s",
                 bathy_dat, bathy_mat)
species_grid_sup(tens_logbooks_sdab,
                 nineties_logbooks_sdab,
                 "Logbook Pacific Sanddab 2010s",
                 bathy_dat, bathy_mat)
species_grid_sup(eighties_surveys_sdab,
                 nineties_surveys_sdab,
                 "Survey Pacific Sanddab 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.79, .85, .11, .25),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.8),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(nineties_surveys_sdab, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.2))
species_grid_sup(nineties_surveys_sdab,
                 nineties_surveys_sdab,
                 "Survey Pacific Sanddab 1990s",
                 bathy_dat, bathy_mat)
species_grid_sup(thousands_surveys_sdab,
                 nineties_surveys_sdab,
                 "Survey Pacific Sanddab 2000s",
                 bathy_dat, bathy_mat)
species_grid_sup(tens_surveys_sdab,
                 nineties_surveys_sdab,
                 "Survey Pacific Sanddab 2010s",
                 bathy_dat, bathy_mat)
dev.off()

## English Sole ----
### Logbooks ----
# Filter to just English sole and survey months
subset_english_logbook <- logbooks_final[logbooks_final$species == 'EGLS_ADJ', ]
subset_english_logbook <- subset_english_logbook[subset_english_logbook$month_day >= 517 &
                                                       subset_english_logbook$month_day <= 929, ]

# Create data exploration map
windows(width = 28, height = 18)
par(mfrow = c(1, 4))
exploration_panels(subset_english_logbook, bathy_dat, bathy_mat, 1981, 1989, "English Sole 1980s")
exploration_panels(subset_english_logbook, bathy_dat, bathy_mat, 1990, 1999, "English Sole 1990s")
exploration_panels(subset_english_logbook, bathy_dat, bathy_mat, 2000, 2009, "English Sole 2000s")
exploration_panels(subset_english_logbook, bathy_dat, bathy_mat, 2010, 2017, "English Sole 2010s")

## Decade maps
# Make sure grid has been created (above)
dev.new(width = 4, height = 10)
cpue_fillpts(subset_english_logbook, 1981, 1989)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_english_logbook, 1990, 1999)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_english_logbook, 2000, 2009)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_english_logbook, 2010, 2017)

## Decade data
# Make sure grid chunk has been run
eighties_logbooks_engl <- cpue_grid(subset_english_logbook, 1981, 1989)
nineties_logbooks_engl <- cpue_grid(subset_english_logbook, 1990, 1999)
thousands_logbooks_engl <- cpue_grid(subset_english_logbook, 2000, 2009)
tens_logbooks_engl <- cpue_grid(subset_english_logbook, 2010, 2017)

## Make maps
# Change legend depending on maximum lncpue - use matrix with highest CPUE
max(eighties_logbooks_engl, na.rm = T)
max(nineties_logbooks_engl, na.rm = T)
max(thousands_logbooks_engl, na.rm = T)
max(tens_logbooks_engl, na.rm = T)# max

windows(width = 15, height = 9)
par(mfrow = c(1, 4),
    family = 'serif',
    mar = c(4, 5, 3, .3) + .1)
cpue_map(eighties_logbooks_engl, tens_logbooks_engl, "PurpOr", "Logbook English Sole 1980s", bathy_dat, bathy_mat)
cpue_map(nineties_logbooks_engl, tens_logbooks_engl, "PurpOr", "Logbook English Sole 1990s", bathy_dat, bathy_mat)
cpue_map(thousands_logbooks_engl, tens_logbooks_engl, "PurpOr", "Logbook English Sole 2000s", bathy_dat, bathy_mat)
cpue_map(tens_logbooks_engl, tens_logbooks_engl, "PurpOr", "Logbook English Sole 2010s", bathy_dat, bathy_mat)
dev.copy(tiff, "../results/visualization/english_fourpanel_logs.tiff",
         width = 15, height = 9, units = "in", res = 200)
dev.off()

## Depth distribution by year
logs_eighties_depth_engl <- depth_contours(subset_english_logbook, 1989, "1980s")
logs_nineties_depth_engl <- depth_contours(subset_english_logbook, 1999, "1990s")
logs_thousands_depth_engl <- depth_contours(subset_english_logbook, 2009, "2000s")
logs_tens_depth_engl <- depth_contours(subset_english_logbook, 2017, "2010s")

windows()
grid.arrange(logs_eighties_depth_engl,
             logs_nineties_depth_engl,
             logs_thousands_depth_engl,
             logs_tens_depth_engl,
             ncol = 2,
             top = textGrob("English Sole Logbook Depth Distribution",
                            gp = gpar(fontfamily = "serif",
                                      cex = 0.8,
                                      fontface = "bold")),
             bottom = textGrob("Depth (m)",
                               gp = gpar(fontfamily = "serif",
                                         cex = 0.8,
                                         fontface = "bold")),
             left = textGrob("Latitude",
                             rot = 90,
                             gp = gpar(fontfamily = "serif",
                                       cex = 0.8,
                                       fontface = "bold",
                                       vjust = 1)))
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/english_depth.tiff",
         width = 4, height = 4, units = "in", res = 200)
dev.off()

## Catch over time
catch_sanddab <- tickets_final %>% filter(common_name == 'English Sole') %>%
  group_by(YEAR) %>%
  summarise(species_weight_sum = sum(TIK_LBS))

# Plot
windows(width = 200, height = 100)
ggplot(data = catch_sanddab,
       aes(x = YEAR, y = species_weight_sum / 1000)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
  geom_col(fill = "orangered4") +
  scale_x_continuous(breaks = c(1980, 1990, 2000, 2010)) +
  theme_tufte() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10)) +
  labs(x = "",
       y = "Total Catch (1000s of lbs)",
       title = "")
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/english_catch.tiff",
         width = 4, height = 2, units = "in", res = 200)
dev.off()

## Average CPUE over time
logbook_english_cpue <- subset_english_logbook %>% group_by(year) %>%
  summarise(cpue_mean = mean(CPUE))

# Plot
windows(width = 20, height = 10)
par(mfrow = c(1, 1))
ggplot(data = logbook_english_cpue, aes(x = year, y = cpue_mean)) +
  geom_path(color = "orangered4", size = 0.5) +
  geom_point(size = 0.9) +
  scale_x_continuous(breaks = c(1980, 1990, 2000, 2010)) +
  theme_tufte() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 8),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10)) +
  labs(x = "",
       y = "CPUE (kg/hr)",
       title = "")
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/english_change.tiff",
         width = 4, height = 2, units = "in", res = 200)
dev.off()

### Survey ----
subset_english <- OR_fish[OR_fish$scientific_name == 'Parophrys vetulus', ]
match_id <- match(survey_data$trawl_id, subset_english$trawl_id)
survey_data$lncpue <- subset_english$lncpue[match_id]
survey_data$CPUE <- subset_english$CPUE[match_id]
survey_data$catch <- subset_english$catch[match_id]
survey_data$pres <- subset_english$pres[match_id]
subset_english_survey <- survey_data
subset_english_survey$lncpue[is.na(subset_english_survey$lncpue)] <- 0
subset_english_survey$CPUE[is.na(subset_english_survey$CPUE)] <- 0
subset_english_survey$catch[is.na(subset_english_survey$catch)] <- 0
subset_english_survey$pres[is.na(subset_english_survey$pres)] <- 0

# Create data exploration map
windows(width = 28, height = 18)
par(mfrow = c(1, 4))
exploration_panels(subset_english_survey, bathy_dat, bathy_mat, 1980, 1989, "English Sole 1980s")
exploration_panels(subset_english_survey, bathy_dat, bathy_mat, 1990, 1999, "English Sole 1990s")
exploration_panels(subset_english_survey, bathy_dat, bathy_mat, 2000, 2009, "English Sole 2000s")
exploration_panels(subset_english_survey, bathy_dat, bathy_mat, 2010, 2018, "English Sole 2010s")

## Decade maps
# Make sure grid has been created (above)
dev.new(width = 4, height = 10)
cpue_fillpts(subset_english_survey, 1980, 1989)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_english_survey, 1990, 1999)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_english_survey, 2000, 2009)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_english_survey, 2010, 2018)

## Decade data
# Make sure grid chunk has been run
eighties_surveys_engl <- cpue_grid(subset_english_survey, 1980, 1989)
nineties_surveys_engl <- cpue_grid(subset_english_survey, 1990, 1999)
thousands_surveys_engl <- cpue_grid(subset_english_survey, 2000, 2009)
tens_surveys_engl <- cpue_grid(subset_english_survey, 2010, 2018)

## Make maps
# Change legend depending on maximum lncpue - use matrix with highest CPUE
max(eighties_surveys_engl, na.rm = T) # max
max(nineties_surveys_engl, na.rm = T)
max(thousands_surveys_engl, na.rm = T)
max(tens_surveys_engl, na.rm = T)

windows(width = 15, height = 9)
par(mfrow = c(1, 4),
    family = 'serif',
    mar = c(4, 5, 3, .3) + .1)
cpue_map(eighties_surveys_engl, eighties_surveys_engl, "PurpOr", "Survey English Sole 1980s", bathy_dat, bathy_mat)
cpue_map(nineties_surveys_engl, eighties_surveys_engl, "PurpOr", "Survey English Sole 1990s", bathy_dat, bathy_mat)
cpue_map(thousands_surveys_engl, eighties_surveys_engl, "PurpOr", "Survey English Sole 2000s", bathy_dat, bathy_mat)
cpue_map(tens_surveys_engl, eighties_surveys_engl, "PurpOr", "Survey English Sole 2010s", bathy_dat, bathy_mat)
dev.copy(tiff, "../results/visualization/english_fourpanel_survey.tiff",
         width = 15, height = 9, units = "in", res = 200)
dev.off()

## Depth distribution by year
survey_eighties_depth_engl <- depth_contours(subset_english_survey, 1989, "1980s")
survey_nineties_depth_engl <- depth_contours(subset_english_survey, 1999, "1990s")
survey_thousands_depth_engl <- depth_contours(subset_english_survey, 2009, "2000s")
survey_tens_depth_engl <- depth_contours(subset_english_survey, 2018, "2010s")

windows()
grid.arrange(survey_eighties_depth_engl,
             survey_nineties_depth_engl,
             survey_thousands_depth_engl,
             survey_tens_depth_engl,
             ncol = 2,
             top = textGrob("English Sole Survey Depth Distribution",
                            gp = gpar(fontfamily = "serif",
                                      cex = 0.8,
                                      fontface = "bold")),
             bottom = textGrob("Depth (m)",
                               gp = gpar(fontfamily = "serif",
                                         cex = 0.8,
                                         fontface = "bold")),
             left = textGrob("Latitude",
                             rot = 90,
                             gp = gpar(fontfamily = "serif",
                                       cex = 0.8,
                                       fontface = "bold",
                                       vjust = 1)))
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/english_depth_survey.tiff",
         width = 4, height = 4, units = "in", res = 200)
dev.off()

## Average CPUE over time
survey_english_cpue <- subset_english_survey %>% group_by(year) %>%
  summarise(cpue_mean = mean(CPUE))

# Plot
windows(width = 20, height = 10)
par(mfrow = c(1, 1))
ggplot(data = survey_english_cpue, aes(x = year, y = cpue_mean)) +
  geom_path(color = "orangered4", size = 0.5) +
  geom_point(size = 0.9) +
  scale_x_continuous(breaks = c(1980, 1990, 2000, 2010)) +
  theme_tufte() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 8),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10)) +
  labs(x = "",
       y = "CPUE (kg/ha)",
       title = "")
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/english_change_survey.tiff",
         width = 4, height = 2, units = "in", res = 200)
dev.off()

### Manuscript Maps ----
# Four panel maps
pdf("../final_figs/manuscript2_fig_tables/english_maps.pdf",
    width = 7.5,
    height = 18)
par(mfrow = c(2, 2),
    family = "serif",
    mar = c(4, 5, 3, .3) + .1)
species_grid_pdf(eighties_logbooks_engl,
                 tens_logbooks_engl, "Logbook English Sole 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.76, .82, .09, .24),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.6),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(tens_logbooks_engl, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.4))
species_grid_pdf(tens_logbooks_engl,
                 tens_logbooks_engl, "Logbook English Sole 2010s",
                 bathy_dat, bathy_mat)
species_grid_pdf(eighties_surveys_engl,
                 eighties_surveys_engl, "Survey English Sole 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.76, .82, .09, .24),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.6),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(eighties_surveys_engl, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.4))
species_grid_pdf(tens_surveys_engl,
                 eighties_surveys_engl, "Survey English Sole 2010s", bathy_dat, bathy_mat)
dev.off()

### Supplement Maps ----
# Eight panel logbooks and surveys
pdf("../final_figs/manuscript2_fig_tables/english_supplement.pdf",
    width = 15,
    height = 17)
par(mfrow = c(2, 4),
    family = "serif",
    mar = c(5.5, 6, 3, .3) + .1)
species_grid_sup(eighties_logbooks_engl,
                 tens_logbooks_engl,
                 "Logbook English Sole 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.79, .85, .11, .25),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.8),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(tens_logbooks_engl, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.2))
species_grid_sup(nineties_logbooks_engl,
                 tens_logbooks_engl,
                 "Logbook English Sole 1990s",
                 bathy_dat, bathy_mat)
species_grid_sup(thousands_logbooks_engl,
                 tens_logbooks_engl,
                 "Logbook English Sole 2000s",
                 bathy_dat, bathy_mat)
species_grid_sup(tens_logbooks_engl,
                 tens_logbooks_engl,
                 "Logbook English Sole 2010s",
                 bathy_dat, bathy_mat)
species_grid_sup(eighties_surveys_engl,
                 eighties_surveys_engl,
                 "Survey English Sole 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.79, .85, .11, .25),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.8),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(eighties_surveys_engl, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.2))
species_grid_sup(nineties_surveys_engl,
                 eighties_surveys_engl,
                 "Survey English Sole 1990s",
                 bathy_dat, bathy_mat)
species_grid_sup(thousands_surveys_engl,
                 eighties_surveys_engl,
                 "Survey English Sole 2000s",
                 bathy_dat, bathy_mat)
species_grid_sup(tens_surveys_engl,
                 eighties_surveys_engl,
                 "Survey English Sole 2010s",
                 bathy_dat, bathy_mat)
dev.off()

## Sand Sole ----
### Logbooks ----
# Filter to just sand sole and survey months
subset_sand_sole_logbook <- logbooks_final[logbooks_final$species == 'SSOL_ADJ', ]
subset_sand_sole_logbook <- subset_sand_sole_logbook[subset_sand_sole_logbook$month_day >= 517 &
                                                   subset_sand_sole_logbook$month_day <= 929, ]

# Create data exploration map
windows(width = 28, height = 18)
par(mfrow = c(1, 4))
exploration_panels(subset_sand_sole_logbook, bathy_dat, bathy_mat, 1981, 1989, "Sand Sole 1980s")
exploration_panels(subset_sand_sole_logbook, bathy_dat, bathy_mat, 1990, 1999, "Sand Sole 1990s")
exploration_panels(subset_sand_sole_logbook, bathy_dat, bathy_mat, 2000, 2009, "Sand Sole 2000s")
exploration_panels(subset_sand_sole_logbook, bathy_dat, bathy_mat, 2010, 2017, "Sand Sole 2010s")

## Decade maps
# Make sure grid has been created (above)
dev.new(width = 4, height = 10)
cpue_fillpts(subset_sand_sole_logbook, 1981, 1989)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_sand_sole_logbook, 1990, 1999)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_sand_sole_logbook, 2000, 2009)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_sand_sole_logbook, 2010, 2017)

## Decade data
# Make sure grid chunk has been run
eighties_logbooks_ssol <- cpue_grid(subset_sand_sole_logbook, 1981, 1989)
nineties_logbooks_ssol <- cpue_grid(subset_sand_sole_logbook, 1990, 1999)
thousands_logbooks_ssol <- cpue_grid(subset_sand_sole_logbook, 2000, 2009)
tens_logbooks_ssol <- cpue_grid(subset_sand_sole_logbook, 2010, 2017)

## Make maps
# Change legend depending on maximum lncpue - use matrix with highest CPUE
max(eighties_logbooks_ssol, na.rm = T)
max(nineties_logbooks_ssol, na.rm = T)
max(thousands_logbooks_ssol, na.rm = T)
max(tens_logbooks_ssol, na.rm = T)# max

windows(width = 15, height = 9)
par(mfrow = c(1, 4),
    family = 'serif',
    mar = c(4, 5, 3, .3) + .1)
cpue_map(eighties_logbooks_ssol, tens_logbooks_ssol, "PurpOr", "Logbook Sand Sole 1980s", bathy_dat, bathy_mat)
cpue_map(nineties_logbooks_ssol, tens_logbooks_ssol, "PurpOr", "Logbook Sand Sole 1990s", bathy_dat, bathy_mat)
cpue_map(thousands_logbooks_ssol, tens_logbooks_ssol, "PurpOr", "Logbook Sand Sole 2000s", bathy_dat, bathy_mat)
cpue_map(tens_logbooks_ssol, tens_logbooks_ssol, "PurpOr", "Logbook Sand Sole 2010s", bathy_dat, bathy_mat)
dev.copy(tiff, "../results/visualization/sand_sole_fourpanel_logs.tiff",
         width = 15, height = 9, units = "in", res = 200)
dev.off()

## Depth distribution by year
logs_eighties_depth_ssol <- depth_contours(subset_sand_sole_logbook, 1989, "1980s")
logs_nineties_depth_ssol <- depth_contours(subset_sand_sole_logbook, 1999, "1990s")
logs_thousands_depth_ssol <- depth_contours(subset_sand_sole_logbook, 2009, "2000s")
logs_tens_depth_ssol <- depth_contours(subset_sand_sole_logbook, 2017, "2010s")

windows()
grid.arrange(logs_eighties_depth_ssol,
             logs_nineties_depth_ssol,
             logs_thousands_depth_ssol,
             logs_tens_depth_ssol,
             ncol = 2,
             top = textGrob("Sand Sole Logbook Depth Distribution",
                            gp = gpar(fontfamily = "serif",
                                      cex = 0.8,
                                      fontface = "bold")),
             bottom = textGrob("Depth (m)",
                               gp = gpar(fontfamily = "serif",
                                         cex = 0.8,
                                         fontface = "bold")),
             left = textGrob("Latitude",
                             rot = 90,
                             gp = gpar(fontfamily = "serif",
                                       cex = 0.8,
                                       fontface = "bold",
                                       vjust = 1)))
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/sand_sole_depth.tiff",
         width = 4, height = 4, units = "in", res = 200)
dev.off()

## Catch over time
catch_sanddab <- tickets_final %>% filter(common_name == 'Sand Sole') %>%
  group_by(YEAR) %>%
  summarise(species_weight_sum = sum(TIK_LBS))

# Plot
windows(width = 200, height = 100)
ggplot(data = catch_sanddab,
       aes(x = YEAR, y = species_weight_sum / 1000)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
  geom_col(fill = "orangered4") +
  theme_tufte() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7)) +
  labs(x = "Year",
       y = "Total Catch (1000s of lbs)",
       title = "Sand Sole Total Nearshore Catch")
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/sand_sole_catch.tiff",
         width = 4, height = 2, units = "in", res = 200)
dev.off()

## Average CPUE over time
logbook_sand_sole_cpue <- subset_sand_sole_logbook %>% group_by(year) %>%
  summarise(cpue_mean = mean(CPUE))

# Plot
windows(width = 20, height = 10)
par(mfrow = c(1, 1))
ggplot(data = logbook_sand_sole_cpue, aes(x = year, y = cpue_mean)) +
  geom_path(color = "orangered4", size = 0.5) +
  geom_point(size = 0.9) +
  theme_tufte() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 8),
        axis.title = element_text(size = 7.5),
        axis.text = element_text(size = 7)) +
  labs(x = "Year",
       y = "CPUE (kg/hr)",
       title = "Mean CPUE of Nearshore Sand Sole Caught in Groundfish Fishery")
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/sand_sole_change.tiff",
         width = 4, height = 2, units = "in", res = 200)
dev.off()

### Survey ----
subset_sandsole <- OR_fish[OR_fish$scientific_name == 'Psettichthys melanostictus', ]
match_id <- match(survey_data$trawl_id, subset_sandsole$trawl_id)
survey_data$lncpue <- subset_sandsole$lncpue[match_id]
survey_data$CPUE <- subset_sandsole$CPUE[match_id]
survey_data$catch <- subset_sandsole$catch[match_id]
survey_data$pres <- subset_sandsole$pres[match_id]
subset_sand_sole_survey <- survey_data
subset_sand_sole_survey$lncpue[is.na(subset_sand_sole_survey$lncpue)] <- 0
subset_sand_sole_survey$CPUE[is.na(subset_sand_sole_survey$CPUE)] <- 0
subset_sand_sole_survey$catch[is.na(subset_sand_sole_survey$catch)] <- 0
subset_sand_sole_survey$pres[is.na(subset_sand_sole_survey$pres)] <- 0

# Create data exploration map
windows(width = 28, height = 18)
par(mfrow = c(1, 4))
exploration_panels(subset_sand_sole_survey, bathy_dat, bathy_mat, 1980, 1989, "Sand Sole 1980s")
exploration_panels(subset_sand_sole_survey, bathy_dat, bathy_mat, 1990, 1999, "Sand Sole 1990s")
exploration_panels(subset_sand_sole_survey, bathy_dat, bathy_mat, 2000, 2009, "Sand Sole 2000s")
exploration_panels(subset_sand_sole_survey, bathy_dat, bathy_mat, 2010, 2018, "Sand Sole 2010s")

## Decade maps
# Make sure grid has been created (above)
dev.new(width = 4, height = 10)
cpue_fillpts(subset_sand_sole_survey, 1980, 1989)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_sand_sole_survey, 1990, 1999)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_sand_sole_survey, 2000, 2009)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_sand_sole_survey, 2010, 2018)

## Decade data
# Make sure grid chunk has been run
eighties_surveys_ssol <- cpue_grid(subset_sand_sole_survey, 1980, 1989)
nineties_surveys_ssol <- cpue_grid(subset_sand_sole_survey, 1990, 1999)
thousands_surveys_ssol <- cpue_grid(subset_sand_sole_survey, 2000, 2009)
tens_surveys_ssol <- cpue_grid(subset_sand_sole_survey, 2010, 2018)

## Make maps
# Change legend depending on maximum lncpue - use matrix with highest CPUE
max(eighties_surveys_ssol, na.rm = T)
max(nineties_surveys_ssol, na.rm = T)
max(thousands_surveys_ssol, na.rm = T) # max
max(tens_surveys_ssol, na.rm = T)

windows(width = 15, height = 9)
par(mfrow = c(1, 4),
    family = 'serif',
    mar = c(4, 5, 3, .3) + .1)
cpue_map(eighties_surveys_ssol, thousands_surveys_ssol, "PurpOr", "Survey Sand Sole 1980s", bathy_dat, bathy_mat)
cpue_map(nineties_surveys_ssol, thousands_surveys_ssol, "PurpOr", "Survey Sand Sole 1990s", bathy_dat, bathy_mat)
cpue_map(thousands_surveys_ssol, thousands_surveys_ssol, "PurpOr", "Survey Sand Sole 2000s", bathy_dat, bathy_mat)
cpue_map(tens_surveys_ssol, thousands_surveys_ssol, "PurpOr", "Survey Sand Sole 2010s", bathy_dat, bathy_mat)
dev.copy(tiff, "../results/visualization/sand_sole_fourpanel_survey.tiff",
         width = 15, height = 9, units = "in", res = 200)
dev.off()

## Depth distribution by year
survey_eighties_depth_ssol <- depth_contours(subset_sand_sole_survey, 1989, "1980s")
survey_nineties_depth_ssol <- depth_contours(subset_sand_sole_survey, 1999, "1990s")
survey_thousands_depth_ssol <- depth_contours(subset_sand_sole_survey, 2009, "2000s")
survey_tens_depth_ssol <- depth_contours(subset_sand_sole_survey, 2018, "2010s")

windows()
grid.arrange(survey_eighties_depth_ssol,
             survey_nineties_depth_ssol,
             survey_thousands_depth_ssol,
             survey_tens_depth_ssol,
             ncol = 2,
             top = textGrob("Sand Sole Survey Depth Distribution",
                            gp = gpar(fontfamily = "serif",
                                      cex = 0.8,
                                      fontface = "bold")),
             bottom = textGrob("Depth (m)",
                               gp = gpar(fontfamily = "serif",
                                         cex = 0.8,
                                         fontface = "bold")),
             left = textGrob("Latitude",
                             rot = 90,
                             gp = gpar(fontfamily = "serif",
                                       cex = 0.8,
                                       fontface = "bold",
                                       vjust = 1)))
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/sand_sole_depth_survey.tiff",
         width = 4, height = 4, units = "in", res = 200)
dev.off()

## Average CPUE over time
survey_sand_sole_cpue <- subset_sand_sole_survey %>% group_by(year) %>%
  summarise(cpue_mean = mean(CPUE))

# Plot
windows(width = 20, height = 10)
par(mfrow = c(1, 1))
ggplot(data = survey_sand_sole_cpue, aes(x = year, y = cpue_mean)) +
  geom_path(color = "orangered4", size = 0.5) +
  geom_point(size = 0.9) +
  theme_tufte() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 8),
        axis.title = element_text(size = 7.5),
        axis.text = element_text(size = 7)) +
  labs(x = "Year",
       y = "CPUE (kg/ha)",
       title = "Mean CPUE of Nearshore Sand Sole Caught in Groundfish Fishery")
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/sandsole_change_survey.tiff",
         width = 4, height = 2, units = "in", res = 200)
dev.off()

### Manuscript Maps ----
# Four panel maps
pdf("../final_figs/manuscript2_fig_tables/sandsole_maps.pdf",
    width = 7.5,
    height = 18)
par(mfrow = c(2, 2),
    family = "serif",
    mar = c(4, 5, 3, .3) + .1)
species_grid_pdf(eighties_logbooks_ssol,
                 tens_logbooks_ssol, "Logbook Sand Sole 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.76, .82, .09, .24),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.6),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(tens_logbooks_ssol, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.4))
species_grid_pdf(tens_logbooks_ssol,
                 tens_logbooks_ssol, "Logbook Sand Sole 2010s",
                 bathy_dat, bathy_mat)
species_grid_pdf(eighties_surveys_ssol,
                 thousands_surveys_ssol, "Survey Sand Sole 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.76, .82, .09, .24),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.6),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(thousands_surveys_ssol, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.4))
species_grid_pdf(tens_surveys_ssol,
                 thousands_surveys_ssol, "Survey Sand Sole 2010s", bathy_dat, bathy_mat)
dev.off()

### Supplement Maps ----
# Eight panel logbooks and surveys
pdf("../final_figs/manuscript2_fig_tables/sandsole_supplement.pdf",
    width = 15,
    height = 17)
par(mfrow = c(2, 4),
    family = "serif",
    mar = c(5.5, 6, 3, .3) + .1)
species_grid_sup(eighties_logbooks_ssol,
                 tens_logbooks_ssol,
                 "Logbook Sand Sole 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.79, .85, .11, .25),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.8),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(tens_logbooks_ssol, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.2))
species_grid_sup(nineties_logbooks_ssol,
                 tens_logbooks_ssol,
                 "Logbook Sand Sole 1990s",
                 bathy_dat, bathy_mat)
species_grid_sup(thousands_logbooks_ssol,
                 tens_logbooks_ssol,
                 "Logbook Sand Sole 2000s",
                 bathy_dat, bathy_mat)
species_grid_sup(tens_logbooks_ssol,
                 tens_logbooks_ssol,
                 "Logbook Sand Sole 2010s",
                 bathy_dat, bathy_mat)
species_grid_sup(eighties_surveys_ssol,
                 thousands_surveys_ssol,
                 "Survey Sand Sole 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.79, .85, .11, .25),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.8),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(thousands_surveys_ssol, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.2))
species_grid_sup(nineties_surveys_ssol,
                 thousands_surveys_ssol,
                 "Survey Sand Sole 1990s",
                 bathy_dat, bathy_mat)
species_grid_sup(thousands_surveys_ssol,
                 thousands_surveys_ssol,
                 "Survey Sand Sole 2000s",
                 bathy_dat, bathy_mat)
species_grid_sup(tens_surveys_ssol,
                 thousands_surveys_ssol,
                 "Survey Sand Sole 2010s",
                 bathy_dat, bathy_mat)
dev.off()


## Starry Flounder ----
### Logbooks ----
# Filter to just Starry Flounder and survey months
subset_starry_logbook <- logbooks_final[logbooks_final$species == 'STRY_ADJ', ]
subset_starry_logbook <- subset_starry_logbook[subset_starry_logbook$month_day >= 517 &
                                                       subset_starry_logbook$month_day <= 929, ]

# Create data exploration map
windows(width = 28, height = 18)
par(mfrow = c(1, 4))
exploration_panels(subset_starry_logbook, bathy_dat, bathy_mat, 1981, 1989, "Starry Flounder 1980s")
exploration_panels(subset_starry_logbook, bathy_dat, bathy_mat, 1990, 1999, "Starry Flounder 1990s")
exploration_panels(subset_starry_logbook, bathy_dat, bathy_mat, 2000, 2009, "Starry Flounder 2000s")
exploration_panels(subset_starry_logbook, bathy_dat, bathy_mat, 2010, 2017, "Starry Flounder 2010s")

## Decade maps
# Make sure grid has been created (above)
dev.new(width = 4, height = 10)
cpue_fillpts(subset_starry_logbook, 1981, 1989)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_starry_logbook, 1990, 1999)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_starry_logbook, 2000, 2009)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_starry_logbook, 2010, 2017)

## Decade data
# Make sure grid chunk has been run
eighties_logbooks_stry <- cpue_grid(subset_starry_logbook, 1981, 1989)
nineties_logbooks_stry <- cpue_grid(subset_starry_logbook, 1990, 1999)
thousands_logbooks_stry <- cpue_grid(subset_starry_logbook, 2000, 2009)
tens_logbooks_stry <- cpue_grid(subset_starry_logbook, 2010, 2017)

## Make maps
# Change legend depending on maximum lncpue - use matrix with highest CPUE
max(eighties_logbooks_stry, na.rm = T)
max(nineties_logbooks_stry, na.rm = T)
max(thousands_logbooks_stry, na.rm = T) # max
max(tens_logbooks_stry, na.rm = T)

windows(width = 15, height = 9)
par(mfrow = c(1, 4),
    family = 'serif',
    mar = c(4, 5, 3, .3) + .1)
cpue_map(eighties_logbooks_stry, thousands_logbooks_stry, "PurpOr", "Logbook Starry Flounder 1980s", bathy_dat, bathy_mat)
cpue_map(nineties_logbooks_stry, thousands_logbooks_stry, "PurpOr", "Logbook Starry Flounder 1990s", bathy_dat, bathy_mat)
cpue_map(thousands_logbooks_stry, thousands_logbooks_stry, "PurpOr", "Logbook Starry Flounder 2000s", bathy_dat, bathy_mat)
cpue_map(tens_logbooks_stry, thousands_logbooks_stry, "PurpOr", "Logbook Starry Flounder 2010s", bathy_dat, bathy_mat)
dev.copy(tiff, "../results/visualization/starry_fourpanel_logs.tiff",
         width = 15, height = 9, units = "in", res = 200)
dev.off()

## Depth distribution by year
logs_eighties_depth_stry <- depth_contours(subset_starry_logbook, 1989, "1980s")
logs_nineties_depth_stry <- depth_contours(subset_starry_logbook, 1999, "1990s")
logs_thousands_depth_stry <- depth_contours(subset_starry_logbook, 2009, "2000s")
logs_tens_depth_stry <- depth_contours(subset_starry_logbook, 2017, "2010s")

windows()
grid.arrange(logs_eighties_depth_stry,
             logs_nineties_depth_stry,
             logs_thousands_depth_stry,
             logs_tens_depth_stry,
             ncol = 2,
             top = textGrob("Starry Flounder Logbook Depth Distribution",
                            gp = gpar(fontfamily = "serif",
                                      cex = 0.8,
                                      fontface = "bold")),
             bottom = textGrob("Depth (m)",
                               gp = gpar(fontfamily = "serif",
                                         cex = 0.8,
                                         fontface = "bold")),
             left = textGrob("Latitude",
                             rot = 90,
                             gp = gpar(fontfamily = "serif",
                                       cex = 0.8,
                                       fontface = "bold",
                                       vjust = 1)))
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/starry_depth.tiff",
         width = 4, height = 4, units = "in", res = 200)
dev.off()

## Catch over time
catch_sanddab <- tickets_final %>% filter(common_name == 'Starry Flounder') %>%
  group_by(YEAR) %>%
  summarise(species_weight_sum = sum(TIK_LBS))

# Plot
windows(width = 200, height = 100)
ggplot(data = catch_sanddab,
       aes(x = YEAR, y = species_weight_sum / 1000)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
  geom_col(fill = "orangered4") +
  theme_tufte() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7)) +
  labs(x = "Year",
       y = "Total Catch (1000s of lbs)",
       title = "Starry Flounder Total Nearshore Catch")
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/starry_catch.tiff",
         width = 4, height = 2, units = "in", res = 200)
dev.off()

## Average CPUE over time
logbook_starry_cpue <- subset_starry_logbook %>% group_by(year) %>%
  summarise(cpue_mean = mean(CPUE))

# Plot
windows(width = 20, height = 10)
par(mfrow = c(1, 1))
ggplot(data = logbook_starry_cpue, aes(x = year, y = cpue_mean)) +
  geom_path(color = "orangered4", size = 0.5) +
  geom_point(size = 0.9) +
  theme_tufte() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 8),
        axis.title = element_text(size = 7.5),
        axis.text = element_text(size = 7)) +
  labs(x = "Year",
       y = "CPUE (kg/hr)",
       title = "Mean CPUE of Nearshore Starry Flounder Caught in Groundfish Fishery")
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/starry_change.tiff",
         width = 4, height = 2, units = "in", res = 200)
dev.off()

### Survey ----
subset_starry <- OR_fish[OR_fish$scientific_name == 'Platichthys stellatus', ]
match_id <- match(survey_data$trawl_id, subset_starry$trawl_id)
survey_data$lncpue <- subset_starry$lncpue[match_id]
survey_data$CPUE <- subset_starry$CPUE[match_id]
survey_data$catch <- subset_starry$catch[match_id]
survey_data$pres <- subset_starry$pres[match_id]
subset_starry_survey <- survey_data
subset_starry_survey$lncpue[is.na(subset_starry_survey$lncpue)] <- 0
subset_starry_survey$CPUE[is.na(subset_starry_survey$CPUE)] <- 0
subset_starry_survey$catch[is.na(subset_starry_survey$catch)] <- 0
subset_starry_survey$pres[is.na(subset_starry_survey$pres)] <- 0

# Create data exploration map
windows(width = 28, height = 18)
par(mfrow = c(1, 4))
exploration_panels(subset_starry_survey, bathy_dat, bathy_mat, 1980, 1989, "Starry Flounder 1980s")
exploration_panels(subset_starry_survey, bathy_dat, bathy_mat, 1990, 1999, "Starry Flounder 1990s")
exploration_panels(subset_starry_survey, bathy_dat, bathy_mat, 2000, 2009, "Starry Flounder 2000s")
exploration_panels(subset_starry_survey, bathy_dat, bathy_mat, 2010, 2018, "Starry Flounder 2010s")

## Decade maps
# Make sure grid has been created (above)
dev.new(width = 4, height = 10)
cpue_fillpts(subset_starry_survey, 1980, 1989)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_starry_survey, 1990, 1999)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_starry_survey, 2000, 2009)

dev.new(width = 4, height = 10)
cpue_fillpts(subset_starry_survey, 2010, 2018)

## Decade data
# Make sure grid chunk has been run
eighties_surveys_stry <- cpue_grid(subset_starry_survey, 1980, 1989)
nineties_surveys_stry <- cpue_grid(subset_starry_survey, 1990, 1999)
thousands_surveys_stry <- cpue_grid(subset_starry_survey, 2000, 2009)
tens_surveys_stry <- cpue_grid(subset_starry_survey, 2010, 2018)

## Make maps
# Change legend depending on maximum lncpue - use matrix with highest CPUE
max(eighties_surveys_stry, na.rm = T)
max(nineties_surveys_stry, na.rm = T)
max(thousands_surveys_stry, na.rm = T) # max
max(tens_surveys_stry, na.rm = T)

windows(width = 15, height = 9)
par(mfrow = c(1, 4),
    family = 'serif',
    mar = c(4, 5, 3, .3) + .1)
cpue_map(eighties_surveys_stry, thousands_surveys_stry, "PurpOr", "Survey Starry Flounder 1980s", bathy_dat, bathy_mat)
cpue_map(nineties_surveys_stry, thousands_surveys_stry, "PurpOr", "Survey Starry Flounder 1990s", bathy_dat, bathy_mat)
cpue_map(thousands_surveys_stry, thousands_surveys_stry, "PurpOr", "Survey Starry Flounder 2000s", bathy_dat, bathy_mat)
cpue_map(tens_surveys_stry, thousands_surveys_stry, "PurpOr", "Survey Starry Flounder 2010s", bathy_dat, bathy_mat)
dev.copy(tiff, "../results/visualization/starry_fourpanel_survey.tiff",
         width = 15, height = 9, units = "in", res = 200)
dev.off()

## Depth distribution by year
survey_eighties_depth_stry <- depth_contours(subset_starry_survey, 1989, "1980s")
survey_nineties_depth_stry <- depth_contours(subset_starry_survey, 1999, "1990s")
survey_thousands_depth_stry <- depth_contours(subset_starry_survey, 2009, "2000s")
survey_tens_depth_stry <- depth_contours(subset_starry_survey, 2018, "2010s")

windows()
grid.arrange(survey_eighties_depth_stry,
             survey_nineties_depth_stry,
             survey_thousands_depth_stry,
             survey_tens_depth_stry,
             ncol = 2,
             top = textGrob("Starry Flounder Survey Depth Distribution",
                            gp = gpar(fontfamily = "serif",
                                      cex = 0.8,
                                      fontface = "bold")),
             bottom = textGrob("Depth (m)",
                               gp = gpar(fontfamily = "serif",
                                         cex = 0.8,
                                         fontface = "bold")),
             left = textGrob("Latitude",
                             rot = 90,
                             gp = gpar(fontfamily = "serif",
                                       cex = 0.8,
                                       fontface = "bold",
                                       vjust = 1)))
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/starry_depth_survey.tiff",
         width = 4, height = 4, units = "in", res = 200)
dev.off()

## Average CPUE over time
survey_starry_cpue <- subset_starry_survey %>% group_by(year) %>%
  summarise(cpue_mean = mean(CPUE))

# Plot
windows(width = 20, height = 10)
par(mfrow = c(1, 1))
ggplot(data = survey_starry_cpue, aes(x = year, y = cpue_mean)) +
  geom_path(color = "orangered4", size = 0.5) +
  geom_point(size = 0.9) +
  theme_tufte() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 8),
        axis.title = element_text(size = 7.5),
        axis.text = element_text(size = 7)) +
  labs(x = "Year",
       y = "CPUE (kg/ha)",
       title = "Mean CPUE of Nearshore Starry Flounder Caught in Groundfish Fishery")
dev.copy(tiff, "../final_figs/manuscript2_fig_tables/starry_change_survey.tiff",
         width = 4, height = 2, units = "in", res = 200)
dev.off()

### Manuscript Maps ----
# Four panel maps
pdf("../final_figs/manuscript2_fig_tables/starry_maps.pdf",
    width = 7.5,
    height = 18)
par(mfrow = c(2, 2),
    family = "serif",
    mar = c(4, 5, 3, .3) + .1)
species_grid_pdf(eighties_logbooks_stry,
                 thousands_logbooks_stry, "Logbook Starry Flounder \n 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.76, .82, .09, .24),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.6),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(thousands_logbooks_stry, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.4))
species_grid_pdf(tens_logbooks_stry,
                 thousands_logbooks_stry, "Logbook Starry Flounder \n 2010s",
                 bathy_dat, bathy_mat)
species_grid_pdf(eighties_surveys_stry,
                 thousands_surveys_stry, "Survey Starry Flounder \n 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.76, .82, .09, .24),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.6),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(thousands_surveys_stry, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.4))
species_grid_pdf(tens_surveys_stry,
                 thousands_surveys_stry, "Survey Starry Flounder \n 2010s", bathy_dat, bathy_mat)
dev.off()

### Supplement Maps ----
# Eight panel logbooks and surveys
pdf("../final_figs/manuscript2_fig_tables/starry_supplement.pdf",
    width = 15,
    height = 17)
par(mfrow = c(2, 4),
    family = "serif",
    mar = c(5.5, 6, 3, .3) + .1)
species_grid_sup(eighties_logbooks_stry,
                 thousands_logbooks_stry,
                 "Logbook Starry Flounder 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.79, .85, .11, .25),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.8),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(thousands_logbooks_stry, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.2))
species_grid_sup(nineties_logbooks_stry,
                 thousands_logbooks_stry,
                 "Logbook Starry Flounder 1990s",
                 bathy_dat, bathy_mat)
species_grid_sup(thousands_logbooks_stry,
                 thousands_logbooks_stry,
                 "Logbook Starry Flounder 2000s",
                 bathy_dat, bathy_mat)
species_grid_sup(tens_logbooks_stry,
                 thousands_logbooks_stry,
                 "Logbook Starry Flounder 2010s",
                 bathy_dat, bathy_mat)
species_grid_sup(eighties_surveys_stry,
                 thousands_surveys_stry,
                 "Survey Starry Flounder 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = viridis(100, option = "F", direction = -1),
           legend.shrink = 0.2,
           smallplot = c(.79, .85, .11, .25),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.8),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(thousands_surveys_stry, na.rm = T)),
           legend.args = list("ln(CPUE+1)",
                              side = 2, cex = 1.2))
species_grid_sup(nineties_surveys_stry,
                 thousands_surveys_stry,
                 "Survey Starry Flounder 1990s",
                 bathy_dat, bathy_mat)
species_grid_sup(thousands_surveys_stry,
                 thousands_surveys_stry,
                 "Survey Starry Flounder 2000s",
                 bathy_dat, bathy_mat)
species_grid_sup(tens_surveys_stry,
                 thousands_surveys_stry,
                 "Survey Starry Flounder 2010s",
                 bathy_dat, bathy_mat)
dev.off()
