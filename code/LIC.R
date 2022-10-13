### Title: Local Index of Collocation
### Purpose: Determine overlap of the data sets
### Date Created: 07/09/2020

## Load libraries, functions, and data ----
library(dplyr)
library(tidyr)
library(reshape)
library(reshape2)
library(maps)
library(mapdata)
library(sgeostat)
library(fields)
library(here)

# Data
load(here('data/ODFW_data', 'logbooks_corrected'))
load(here ('data/ODFW_data', 'fish_tickets_final'))
load(here('data/NMFS_data', 'trawl_data'))
load(here('data/NMFS_data', 'OR_fish'))
load(here('data/ODFW_data', 'eighties_logbook'))
load(here('data/ODFW_data', 'nineties_logbook'))
load(here('data/ODFW_data', 'thousands_logbook'))
load(here('data/ODFW_data', 'teens_logbook'))
survey_data <- trawl_data
rm(trawl_data)

# Functions
# local index of collocation function (just gives singular value)
loc_collocfn <- function(prey, pred) {
  p_prey <- prey / sum(prey, na.rm = T)
  p_pred <- pred / sum(pred, na.rm = T)
  sum(p_prey * p_pred, na.rm = T) / (sqrt(sum(p_prey ^ 2, na.rm = T) * sum(p_pred ^
                                                                             2, na.rm = T)))
}
source(here("code/functions", "biomass_fillpts.R"))
source(here("code/functions", "biomass_grid.R"))
source(here("code/functions", "lic_map.R"))
source(here("code/functions", "lic_sup.R"))
source(here("code/functions", "spatial_lic.R"))

# For depth, import data and show contour on a map
# .xyz option no longer available for download
bathy_dat <- read.table(here("data", "etopo1.xyz"), sep = '')
names(bathy_dat) <- c('lon', 'lat', 'depth')
bathy_dat$depth[bathy_dat$depth > 0] <- NA # Avoid points above water
bathy_mat <- matrix(bathy_dat$depth,
                    nrow = length(unique(bathy_dat$lon)),
                    ncol = length(unique(bathy_dat$lat)))[, order(unique(bathy_dat$lat))]

## Prepare data and grids ----
# Logbooks: convert CPUE to kg/hr and ln(x+1), add other necessary columns
logbooks_final$kg_caught <- logbooks_final$species_weight * 0.4535924
logbooks_final <- logbooks_final[!logbooks_final$DURATION == 0, ] # remove tows with no trawl duration
logbooks_final$CPUE <- logbooks_final$kg_caught / logbooks_final$DURATION
logbooks_final <- logbooks_final[!is.na(logbooks_final$CPUE), ]
logbooks_final$lncpue <- log(logbooks_final$CPUE + 1)
logbooks_final <- logbooks_final[logbooks_final$depth <= -5, ] # remove unreasonably shallow tows
logbooks_final$pres <- 1 * (logbooks_final$species_weight > 0)
logbooks_final$month_day <- as.numeric(format(logbooks_final$TOWDATE, '%m%d'))
logbooks_final <- logbooks_final[logbooks_final$month_day >= 517 &
                                 logbooks_final$month_day <= 929, ]

# Create index matrices from average tow counts
# Use if reducing to only high effort cells
# eighties_index <- replace(eighties_logbook, eighties_logbook < 0.1 * mean(eighties_logbook, na.rm = T), NA) != TRUE
# eighties_index[is.na(eighties_index)] <- FALSE
# nineties_index <- replace(nineties_logbook, nineties_logbook < 0.1 * mean(nineties_logbook, na.rm = T), NA) != TRUE
# nineties_index[is.na(nineties_index)] <- FALSE
# thousands_index <- replace(thousands_logbook, thousands_logbook < 0.1 * mean(thousands_logbook, na.rm = T), NA) != TRUE
# thousands_index[is.na(thousands_index)] <- FALSE
# teens_index <- replace(teens_logbook, teens_logbook < 0.1 * mean(teens_logbook, na.rm = T), NA) != TRUE
# teens_index[is.na(teens_index)] <- FALSE

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

## Petrale Sole ----
# Filter to just petrale
logbook_petrale <- logbooks_final[logbooks_final$species == 'PTRL_ADJ', ]

subset_petrale <- OR_fish[OR_fish$scientific_name == 'Eopsetta jordani', ]
match_id <- match(survey_data$trawl_id, subset_petrale$trawl_id)
survey_data$CPUE <- subset_petrale$CPUE[match_id]
survey_petrale <- survey_data
survey_petrale$CPUE[is.na(survey_petrale$CPUE)] <- 0

### Logbooks
# Decade grids fill with proportion of biomass in each cell for corresponding points
dev.new(width = 4, height = 10)
biomass_fillpts(logbook_petrale, 1981, 1989)
dev.new(width = 4, height = 10)
biomass_fillpts(logbook_petrale, 1990, 2001)
dev.new(width = 4, height = 10)
biomass_fillpts(logbook_petrale, 2002, 2009)
dev.new(width = 4, height = 10)
biomass_fillpts(logbook_petrale, 2010, 2017)

# Fill decade grids with data
eighties_logbooks_ptrl <- biomass_grid(logbook_petrale, 1981, 1989)
nineties_logbooks_ptrl <- biomass_grid(logbook_petrale, 1990, 2001)
thousands_logbooks_ptrl <- biomass_grid(logbook_petrale, 2002, 2009)
tens_logbooks_ptrl <- biomass_grid(logbook_petrale, 2010, 2017)

### Survey
dev.new(width = 4, height = 10)
biomass_fillpts(survey_petrale, 1981, 1989)
dev.new(width = 4, height = 10)
biomass_fillpts(survey_petrale, 1990, 2001)
dev.new(width = 4, height = 10)
biomass_fillpts(survey_petrale, 2002, 2009)
dev.new(width = 4, height = 10)
biomass_fillpts(survey_petrale, 2010, 2017)

# Fill decade grids with data
eighties_survey_ptrl <- biomass_grid(survey_petrale, 1980, 1989)
nineties_survey_ptrl <- biomass_grid(survey_petrale, 1990, 2001)
thousands_survey_ptrl <- biomass_grid(survey_petrale, 2002, 2009)
tens_survey_ptrl <- biomass_grid(survey_petrale, 2010, 2018)

### Spatial indicators
# local index of collocation per decade for petrale sole overall
loc_collocfn(eighties_logbooks_ptrl, eighties_survey_ptrl)
loc_collocfn(nineties_logbooks_ptrl, nineties_survey_ptrl)
loc_collocfn(thousands_logbooks_ptrl, thousands_survey_ptrl)
loc_collocfn(tens_logbooks_ptrl, tens_survey_ptrl)

# Fill cell values for LIC
eighties_petrale <- spatial_lic(eighties_logbooks_ptrl, logbook_petrale,
                               eighties_survey_ptrl, survey_petrale, 1980, 1989)
nineties_petrale <- spatial_lic(nineties_logbooks_ptrl, logbook_petrale,
                               nineties_survey_ptrl, survey_petrale, 1990, 2001)
thousands_petrale <- spatial_lic(thousands_logbooks_ptrl, logbook_petrale,
                                thousands_survey_ptrl, survey_petrale, 2002, 2009)
tens_petrale <- spatial_lic(tens_logbooks_ptrl, logbook_petrale,
                           tens_survey_ptrl, survey_petrale, 2010, 2018)

# Determine maximum value to scale appropriately
max(eighties_petrale, na.rm = T)
max(nineties_petrale, na.rm = T)
max(thousands_petrale, na.rm = T) # max
max(tens_petrale, na.rm = T)

# Mean and SD
mean(eighties_petrale, na.rm = T)
mean(nineties_petrale, na.rm = T)
mean(thousands_petrale, na.rm = T)
mean(tens_petrale, na.rm = T)

sd(eighties_petrale, na.rm = T)
sd(nineties_petrale, na.rm = T)
sd(thousands_petrale, na.rm = T)
sd(tens_petrale, na.rm = T)

# Create map
pdf("../final_figs/fish_res_fig_tables/petrale_overlap.pdf",
    width = 7.5,
    height = 9)
par(mfrow = c(1, 2),
    family = "serif",
    mar = c(4, 5, 3, .3) + .1)
lic_map(eighties_petrale, thousands_petrale, "Petrale Sole 1980s",
        bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = hcl.colors(100, "YlOrRd", rev = F),
           legend.shrink = 0.2,
           smallplot = c(.73, .79, .11, .25),
           legend.cex = 1.4,
           axis.args = list(cex.axis = 1.2),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(thousands_petrale, na.rm = T)),
           legend.args = list("LIC",
                              side = 2, cex = 1.2))
lic_map(tens_petrale, thousands_petrale, "Petrale Sole 2010s",
        bathy_dat, bathy_mat)
dev.off()

# Supplement figure
pdf("../final_figs/fish_res_fig_tables/petrale_overlap_supplement.pdf",
    width = 15,
    height = 8.5)
par(mfrow = c(1, 4),
    family = "serif",
    mar = c(5.5, 6, 3, .3) + .1)
lic_sup(eighties_petrale, eighties_petrale, "Petrale Sole 1980s",
        bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = hcl.colors(100, "YlOrRd", rev = F),
           legend.shrink = 0.2,
           smallplot = c(.76, .82, .11, .25),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.8),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(eighties_petrale, na.rm = T)),
           legend.args = list("LIC",
                              side = 2, cex = 1.2))
lic_sup(nineties_petrale, eighties_petrale, "Petrale Sole 1990s",
        bathy_dat, bathy_mat)
lic_sup(thousands_petrale, eighties_petrale, "Petrale Sole 2000s",
        bathy_dat, bathy_mat)
lic_sup(tens_petrale, eighties_petrale, "Petrale Sole 2010s",
        bathy_dat, bathy_mat)
dev.off()

## Dover Sole ----
# Filter to just dover
logbook_dover <- logbooks_final[logbooks_final$species == 'DOVR_ADJ', ]

subset_dover <- OR_fish[OR_fish$scientific_name == 'Microstomus pacificus', ]
match_id <- match(survey_data$trawl_id, subset_dover$trawl_id)
survey_data$CPUE <- subset_dover$CPUE[match_id]
survey_dover <- survey_data
survey_dover$CPUE[is.na(survey_dover$CPUE)] <- 0

### Logbooks
# Decade grids fill with proportion of biomass in each cell for corresponding points
dev.new(width = 4, height = 10)
biomass_fillpts(logbook_dover, 1981, 1989)
dev.new(width = 4, height = 10)
biomass_fillpts(logbook_dover, 1990, 2001)
dev.new(width = 4, height = 10)
biomass_fillpts(logbook_dover, 2002, 2009)
dev.new(width = 4, height = 10)
biomass_fillpts(logbook_dover, 2010, 2017)

# Fill decade grids with data
eighties_logbooks_dovr <- biomass_grid(logbook_dover, 1981, 1989)
nineties_logbooks_dovr <- biomass_grid(logbook_dover, 1990, 2001)
thousands_logbooks_dovr <- biomass_grid(logbook_dover, 2002, 2009)
tens_logbooks_dovr <- biomass_grid(logbook_dover, 2010, 2017)

### Survey
dev.new(width = 4, height = 10)
biomass_fillpts(survey_dover, 1981, 1989)
dev.new(width = 4, height = 10)
biomass_fillpts(survey_dover, 1990, 2001)
dev.new(width = 4, height = 10)
biomass_fillpts(survey_dover, 2002, 2009)
dev.new(width = 4, height = 10)
biomass_fillpts(survey_dover, 2010, 2017)

# Fill decade grids with data
eighties_survey_dovr <- biomass_grid(survey_dover, 1980, 1989)
nineties_survey_dovr <- biomass_grid(survey_dover, 1990, 2001)
thousands_survey_dovr <- biomass_grid(survey_dover, 2002, 2009)
tens_survey_dovr <- biomass_grid(survey_dover, 2010, 2018)

### Spatial indicators
# local index of collocation per decade for dover sole overall
loc_collocfn(eighties_logbooks_dovr, eighties_survey_dovr)
loc_collocfn(nineties_logbooks_dovr, nineties_survey_dovr)
loc_collocfn(thousands_logbooks_dovr, thousands_survey_dovr)
loc_collocfn(tens_logbooks_dovr, tens_survey_dovr)

# Fill cell values for LIC
eighties_dover <- spatial_lic(eighties_logbooks_dovr, logbook_dover,
                                eighties_survey_dovr, survey_dover, 1980, 1989)
nineties_dover <- spatial_lic(nineties_logbooks_dovr, logbook_dover,
                                nineties_survey_dovr, survey_dover, 1990, 2001)
thousands_dover <- spatial_lic(thousands_logbooks_dovr, logbook_dover,
                                 thousands_survey_dovr, survey_dover, 2002, 2009)
tens_dover <- spatial_lic(tens_logbooks_dovr, logbook_dover,
                            tens_survey_dovr, survey_dover, 2010, 2018)

# Determine maximum value to scale appropriately
max(eighties_dover, na.rm = T) # max
max(nineties_dover, na.rm = T)
max(thousands_dover, na.rm = T)
max(tens_dover, na.rm = T)

# Mean and SD
mean(eighties_dover, na.rm = T)
mean(nineties_dover, na.rm = T)
mean(thousands_dover, na.rm = T)
mean(tens_dover, na.rm = T)

sd(eighties_dover, na.rm = T)
sd(nineties_dover, na.rm = T)
sd(thousands_dover, na.rm = T)
sd(tens_dover, na.rm = T)

# Create map
# Manuscript figure
pdf("../final_figs/fish_res_fig_tables/dover_overlap.pdf",
    width = 7.5,
    height = 9)
par(mfrow = c(1, 2),
    family = "serif",
    mar = c(4, 5, 3, .3) + .1)
lic_map(eighties_dover, eighties_dover, "Dover Sole 1980s",
        bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = hcl.colors(100, "YlOrRd", rev = F),
           legend.shrink = 0.2,
           smallplot = c(.73, .79, .11, .25),
           legend.cex = 1.4,
           axis.args = list(cex.axis = 1.2),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(eighties_dover, na.rm = T)),
           legend.args = list("LIC",
                              side = 2, cex = 1.2))
lic_map(tens_dover, eighties_dover, "Dover Sole 2010s",
        bathy_dat, bathy_mat)
dev.off()

# Supplement figure
pdf("../final_figs/fish_res_fig_tables/dover_overlap_supplement.pdf",
    width = 15,
    height = 8.5)
par(mfrow = c(1, 4),
    family = "serif",
    mar = c(5.5, 6, 3, .3) + .1)
lic_sup(eighties_dover, eighties_dover, "Dover Sole 1980s",
        bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = hcl.colors(100, "YlOrRd", rev = F),
           legend.shrink = 0.2,
           smallplot = c(.76, .82, .11, .25),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.8),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(eighties_dover, na.rm = T)),
           legend.args = list("LIC",
                              side = 2, cex = 1.2))
lic_sup(nineties_dover, eighties_dover, "Dover Sole 1990s",
        bathy_dat, bathy_mat)
lic_sup(thousands_dover, eighties_dover, "Dover Sole 2000s",
        bathy_dat, bathy_mat)
lic_sup(tens_dover, eighties_dover, "Dover Sole 2010s",
        bathy_dat, bathy_mat)
dev.off()

## English Sole ----
# Filter to just English
logbook_english <- logbooks_final[logbooks_final$species == 'EGLS_ADJ', ]

subset_english <- OR_fish[OR_fish$scientific_name == 'Parophrys vetulus', ]
match_id <- match(survey_data$trawl_id, subset_english$trawl_id)
survey_data$CPUE <- subset_english$CPUE[match_id]
survey_english <- survey_data
survey_english$CPUE[is.na(survey_english$CPUE)] <- 0

### Logbooks
# Decade grids fill with proportion of biomass in each cell for corresponding points
dev.new(width = 4, height = 10)
biomass_fillpts(logbook_english, 1981, 1989)
dev.new(width = 4, height = 10)
biomass_fillpts(logbook_english, 1990, 2001)
dev.new(width = 4, height = 10)
biomass_fillpts(logbook_english, 2002, 2009)
dev.new(width = 4, height = 10)
biomass_fillpts(logbook_english, 2010, 2017)

# Fill decade grids with data
eighties_logbooks_egls <- biomass_grid(logbook_english, 1981, 1989)
nineties_logbooks_egls <- biomass_grid(logbook_english, 1990, 2001)
thousands_logbooks_egls <- biomass_grid(logbook_english, 2002, 2009)
tens_logbooks_egls <- biomass_grid(logbook_english, 2010, 2017)

### Survey
dev.new(width = 4, height = 10)
biomass_fillpts(survey_english, 1981, 1989)
dev.new(width = 4, height = 10)
biomass_fillpts(survey_english, 1990, 2001)
dev.new(width = 4, height = 10)
biomass_fillpts(survey_english, 2002, 2009)
dev.new(width = 4, height = 10)
biomass_fillpts(survey_english, 2010, 2017)

# Fill decade grids with data
eighties_survey_egls <- biomass_grid(survey_english, 1980, 1989)
nineties_survey_egls <- biomass_grid(survey_english, 1990, 2001)
thousands_survey_egls <- biomass_grid(survey_english, 2002, 2009)
tens_survey_egls <- biomass_grid(survey_english, 2010, 2018)

### Spatial indicators
# local index of collocation per decade for English Sole overall
loc_collocfn(eighties_logbooks_egls, eighties_survey_egls)
loc_collocfn(nineties_logbooks_egls, nineties_survey_egls)
loc_collocfn(thousands_logbooks_egls, thousands_survey_egls)
loc_collocfn(tens_logbooks_egls, tens_survey_egls)

# Fill cell values for LIC
eighties_english <- spatial_lic(eighties_logbooks_egls, logbook_english,
                              eighties_survey_egls, survey_english, 1980, 1989)
nineties_english <- spatial_lic(nineties_logbooks_egls, logbook_english,
                              nineties_survey_egls, survey_english, 1990, 2001)
thousands_english <- spatial_lic(thousands_logbooks_egls, logbook_english,
                               thousands_survey_egls, survey_english, 2002, 2009)
tens_english <- spatial_lic(tens_logbooks_egls, logbook_english,
                          tens_survey_egls, survey_english, 2010, 2018)

# Determine maximum value to scale appropriately
max(eighties_english, na.rm = T) # max
max(nineties_english, na.rm = T)
max(thousands_english, na.rm = T)
max(tens_english, na.rm = T)

# Mean and SD
mean(eighties_english, na.rm = T)
mean(nineties_english, na.rm = T)
mean(thousands_english, na.rm = T)
mean(tens_english, na.rm = T)

sd(eighties_english, na.rm = T)
sd(nineties_english, na.rm = T)
sd(thousands_english, na.rm = T)
sd(tens_english, na.rm = T)

# Create map
# Manuscript figure
pdf("../final_figs/fish_res_fig_tables/english_overlap.pdf",
    width = 7.5,
    height = 9)
par(mfrow = c(1, 2),
    family = "serif",
    mar = c(4, 5, 3, .3) + .1)
lic_map(eighties_english, eighties_english, "English Sole 1980s",
        bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = hcl.colors(100, "YlOrRd", rev = F),
           legend.shrink = 0.2,
           smallplot = c(.73, .79, .11, .25),
           legend.cex = 1.4,
           axis.args = list(cex.axis = 1.2),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(eighties_english, na.rm = T)),
           legend.args = list("LIC",
                              side = 2, cex = 1.2))
lic_map(tens_english, eighties_english, "English Sole 2010s",
        bathy_dat, bathy_mat)
dev.off()

# Supplement figure
pdf("../final_figs/fish_res_fig_tables/english_overlap_supplement.pdf",
    width = 15,
    height = 8.5)
par(mfrow = c(1, 4),
    family = "serif",
    mar = c(5.5, 6, 3, .3) + .1)
lic_sup(eighties_english, eighties_english, "English Sole 1980s",
        bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = hcl.colors(100, "YlOrRd", rev = F),
           legend.shrink = 0.2,
           smallplot = c(.76, .82, .11, .25),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.8),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(eighties_english, na.rm = T)),
           legend.args = list("LIC",
                              side = 2, cex = 1.2))
lic_sup(nineties_english, eighties_english, "English Sole 1990s",
        bathy_dat, bathy_mat)
lic_sup(thousands_english, eighties_english, "English Sole 2000s",
        bathy_dat, bathy_mat)
lic_sup(tens_english, eighties_english, "English Sole 2010s",
        bathy_dat, bathy_mat)
dev.off()

## Pacific Sanddab ----
# Filter to just sanddab
logbook_sanddab <- logbooks_final[logbooks_final$species == 'SDAB_ADJ', ]

subset_sanddab <- OR_fish[OR_fish$scientific_name == 'Citharichthys sordidus', ]
match_id <- match(survey_data$trawl_id, subset_sanddab$trawl_id)
survey_data$CPUE <- subset_sanddab$CPUE[match_id]
survey_sanddab <- survey_data
survey_sanddab$CPUE[is.na(survey_sanddab$CPUE)] <- 0

### Logbooks
# Decade grids fill with proportion of biomass in each cell for corresponding points
dev.new(width = 4, height = 10)
biomass_fillpts(logbook_sanddab, 1981, 1989)
dev.new(width = 4, height = 10)
biomass_fillpts(logbook_sanddab, 1990, 2001)
dev.new(width = 4, height = 10)
biomass_fillpts(logbook_sanddab, 2002, 2009)
dev.new(width = 4, height = 10)
biomass_fillpts(logbook_sanddab, 2010, 2017)

# Fill decade grids with data
eighties_logbooks_sdab <- biomass_grid(logbook_sanddab, 1981, 1989)
nineties_logbooks_sdab <- biomass_grid(logbook_sanddab, 1990, 2001)
thousands_logbooks_sdab <- biomass_grid(logbook_sanddab, 2002, 2009)
tens_logbooks_sdab <- biomass_grid(logbook_sanddab, 2010, 2017)

### Survey
dev.new(width = 4, height = 10)
biomass_fillpts(survey_sanddab, 1981, 1989)
dev.new(width = 4, height = 10)
biomass_fillpts(survey_sanddab, 1990, 2001)
dev.new(width = 4, height = 10)
biomass_fillpts(survey_sanddab, 2002, 2009)
dev.new(width = 4, height = 10)
biomass_fillpts(survey_sanddab, 2010, 2017)

# Fill decade grids with data
eighties_survey_sdab <- biomass_grid(survey_sanddab, 1980, 1989)
nineties_survey_sdab <- biomass_grid(survey_sanddab, 1990, 2001)
thousands_survey_sdab <- biomass_grid(survey_sanddab, 2002, 2009)
tens_survey_sdab <- biomass_grid(survey_sanddab, 2010, 2018)

### Spatial indicators
# local index of collocation per decade for Pacific Sanddab overall
loc_collocfn(eighties_logbooks_sdab, eighties_survey_sdab)
loc_collocfn(nineties_logbooks_sdab, nineties_survey_sdab)
loc_collocfn(thousands_logbooks_sdab, thousands_survey_sdab)
loc_collocfn(tens_logbooks_sdab, tens_survey_sdab)

# Fill cell values for LIC
eighties_sanddab <- spatial_lic(eighties_logbooks_sdab, logbook_sanddab,
                                eighties_survey_sdab, survey_sanddab, 1980, 1989)
nineties_sanddab <- spatial_lic(nineties_logbooks_sdab, logbook_sanddab,
                                nineties_survey_sdab, survey_sanddab, 1990, 2001)
thousands_sanddab <- spatial_lic(thousands_logbooks_sdab, logbook_sanddab,
                                 thousands_survey_sdab, survey_sanddab, 2002, 2009)
tens_sanddab <- spatial_lic(tens_logbooks_sdab, logbook_sanddab,
                            tens_survey_sdab, survey_sanddab, 2010, 2018)

# Determine maximum value to scale appropriately
max(eighties_sanddab, na.rm = T) # max
max(nineties_sanddab, na.rm = T)
max(thousands_sanddab, na.rm = T)
max(tens_sanddab, na.rm = T)

# Mean and SD
mean(eighties_sanddab, na.rm = T)
mean(nineties_sanddab, na.rm = T)
mean(thousands_sanddab, na.rm = T)
mean(tens_sanddab, na.rm = T)

sd(eighties_sanddab, na.rm = T)
sd(nineties_sanddab, na.rm = T)
sd(thousands_sanddab, na.rm = T)
sd(tens_sanddab, na.rm = T)

# Create map
# Manuscript figure
pdf("../final_figs/fish_res_fig_tables/sanddab_overlap.pdf",
    width = 7.5,
    height = 9)
par(mfrow = c(1, 2),
    family = "serif",
    mar = c(4, 5, 3, .3) + .1)
lic_map(eighties_sanddab, eighties_sanddab, "Pacific Sanddab 1980s",
        bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = hcl.colors(100, "YlOrRd", rev = F),
           legend.shrink = 0.2,
           smallplot = c(.73, .79, .11, .25),
           legend.cex = 1.4,
           axis.args = list(cex.axis = 1.2),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(eighties_sanddab, na.rm = T)),
           legend.args = list("LIC",
                              side = 2, cex = 1.2))
lic_map(tens_sanddab, eighties_sanddab, "Pacific Sanddab 2010s",
        bathy_dat, bathy_mat)
dev.off()

# Supplement figure
pdf("../final_figs/fish_res_fig_tables/sanddab_overlap_supplement.pdf",
    width = 15,
    height = 8.5)
par(mfrow = c(1, 4),
    family = "serif",
    mar = c(5.5, 6, 3, .3) + .1)
lic_sup(eighties_sanddab, eighties_sanddab, "Pacific Sanddab 1980s",
        bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = hcl.colors(100, "YlOrRd", rev = F),
           legend.shrink = 0.2,
           smallplot = c(.76, .82, .11, .25),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.8),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(eighties_sanddab, na.rm = T)),
           legend.args = list("LIC",
                              side = 2, cex = 1.2))
lic_sup(nineties_sanddab, eighties_sanddab, "Pacific Sanddab 1990s",
        bathy_dat, bathy_mat)
lic_sup(thousands_sanddab, eighties_sanddab, "Pacific Sanddab 2000s",
        bathy_dat, bathy_mat)
lic_sup(tens_sanddab, eighties_sanddab, "Pacific Sanddab 2010s",
        bathy_dat, bathy_mat)
dev.off()

## Sand Sole ----
# Filter to just sand sole
logbook_sandsole <- logbooks_final[logbooks_final$species == 'SSOL_ADJ', ]

subset_sandsole <- OR_fish[OR_fish$scientific_name == 'Psettichthys melanostictus', ]
match_id <- match(survey_data$trawl_id, subset_sandsole$trawl_id)
survey_data$CPUE <- subset_sandsole$CPUE[match_id]
survey_sandsole <- survey_data
survey_sandsole$CPUE[is.na(survey_sandsole$CPUE)] <- 0

### Logbooks
# Decade grids fill with proportion of biomass in each cell for corresponding points
dev.new(width = 4, height = 10)
biomass_fillpts(logbook_sandsole, 1981, 1989)
dev.new(width = 4, height = 10)
biomass_fillpts(logbook_sandsole, 1990, 2001)
dev.new(width = 4, height = 10)
biomass_fillpts(logbook_sandsole, 2002, 2009)
dev.new(width = 4, height = 10)
biomass_fillpts(logbook_sandsole, 2010, 2017)

# Fill decade grids with data
eighties_logbooks_ssol <- biomass_grid(logbook_sandsole, 1981, 1989)
nineties_logbooks_ssol <- biomass_grid(logbook_sandsole, 1990, 2001)
thousands_logbooks_ssol <- biomass_grid(logbook_sandsole, 2002, 2009)
tens_logbooks_ssol <- biomass_grid(logbook_sandsole, 2010, 2017)

### Survey
dev.new(width = 4, height = 10)
biomass_fillpts(survey_sandsole, 1981, 1989)
dev.new(width = 4, height = 10)
biomass_fillpts(survey_sandsole, 1990, 2001)
dev.new(width = 4, height = 10)
biomass_fillpts(survey_sandsole, 2002, 2009)
dev.new(width = 4, height = 10)
biomass_fillpts(survey_sandsole, 2010, 2017)

# Fill decade grids with data
eighties_survey_ssol <- biomass_grid(survey_sandsole, 1980, 1989)
nineties_survey_ssol <- biomass_grid(survey_sandsole, 1990, 2001)
thousands_survey_ssol <- biomass_grid(survey_sandsole, 2002, 2009)
tens_survey_ssol <- biomass_grid(survey_sandsole, 2010, 2018)

### Spatial indicators
# local index of collocation per decade for Sand Sole overall
loc_collocfn(eighties_logbooks_ssol, eighties_survey_ssol)
loc_collocfn(nineties_logbooks_ssol, nineties_survey_ssol)
loc_collocfn(thousands_logbooks_ssol, thousands_survey_ssol)
loc_collocfn(tens_logbooks_ssol, tens_survey_ssol)

# Fill cell values for LIC
eighties_sandsole <- spatial_lic(eighties_logbooks_ssol, logbook_sandsole,
                                eighties_survey_ssol, survey_sandsole, 1980, 1989)
nineties_sandsole <- spatial_lic(nineties_logbooks_ssol, logbook_sandsole,
                                nineties_survey_ssol, survey_sandsole, 1990, 2001)
thousands_sandsole <- spatial_lic(thousands_logbooks_ssol, logbook_sandsole,
                                 thousands_survey_ssol, survey_sandsole, 2002, 2009)
tens_sandsole <- spatial_lic(tens_logbooks_ssol, logbook_sandsole,
                            tens_survey_ssol, survey_sandsole, 2010, 2018)

# Determine maximum value to scale appropriately
max(eighties_sandsole, na.rm = T) # max
max(nineties_sandsole, na.rm = T)
max(thousands_sandsole, na.rm = T)
max(tens_sandsole, na.rm = T)

# Mean and SD
mean(eighties_sandsole, na.rm = T)
mean(nineties_sandsole, na.rm = T)
mean(thousands_sandsole, na.rm = T)
mean(tens_sandsole, na.rm = T)

sd(eighties_sandsole, na.rm = T)
sd(nineties_sandsole, na.rm = T)
sd(thousands_sandsole, na.rm = T)
sd(tens_sandsole, na.rm = T)

# Create map
# Manuscript figure
pdf("../final_figs/fish_res_fig_tables/sandsole_overlap.pdf",
    width = 7.5,
    height = 9)
par(mfrow = c(1, 2),
    family = "serif",
    mar = c(4, 5, 3, .3) + .1)
lic_map(eighties_sandsole, eighties_sandsole, "Sand Sole 1980s",
        bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = hcl.colors(100, "YlOrRd", rev = F),
           legend.shrink = 0.2,
           smallplot = c(.73, .79, .11, .25),
           legend.cex = 1.4,
           axis.args = list(cex.axis = 1.2),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(eighties_sandsole, na.rm = T)),
           legend.args = list("LIC",
                              side = 2, cex = 1.2))
lic_map(tens_sandsole, eighties_sandsole, "Sand Sole 2010s",
        bathy_dat, bathy_mat)
dev.off()

# Supplement figure
pdf("../final_figs/fish_res_fig_tables/sandsole_overlap_supplement.pdf",
    width = 15,
    height = 8.5)
par(mfrow = c(1, 4),
    family = "serif",
    mar = c(5.5, 6, 3, .3) + .1)
lic_sup(eighties_sandsole, eighties_sandsole, "Sand Sole 1980s",
        bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = hcl.colors(100, "YlOrRd", rev = F),
           legend.shrink = 0.2,
           smallplot = c(.76, .82, .11, .25),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.8),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(eighties_sandsole, na.rm = T)),
           legend.args = list("LIC",
                              side = 2, cex = 1.2))
lic_sup(nineties_sandsole, eighties_sandsole, "Sand Sole 1990s",
        bathy_dat, bathy_mat)
lic_sup(thousands_sandsole, eighties_sandsole, "Sand Sole 2000s",
        bathy_dat, bathy_mat)
lic_sup(tens_sandsole, eighties_sandsole, "Sand Sole 2010s",
        bathy_dat, bathy_mat)
dev.off()

## Starry Flounder ----
# Filter to just starry
logbook_starry <- logbooks_final[logbooks_final$species == 'STRY_ADJ', ]

subset_starry <- OR_fish[OR_fish$scientific_name == 'Platichthys stellatus', ]
match_id <- match(survey_data$trawl_id, subset_starry$trawl_id)
survey_data$CPUE <- subset_starry$CPUE[match_id]
survey_starry <- survey_data
survey_starry$CPUE[is.na(survey_starry$CPUE)] <- 0

### Logbooks
# Decade grids fill with proportion of biomass in each cell for corresponding points
dev.new(width = 4, height = 10)
biomass_fillpts(logbook_starry, 1981, 1989)
dev.new(width = 4, height = 10)
biomass_fillpts(logbook_starry, 1990, 2001)
dev.new(width = 4, height = 10)
biomass_fillpts(logbook_starry, 2002, 2009)
dev.new(width = 4, height = 10)
biomass_fillpts(logbook_starry, 2010, 2017)

# Fill decade grids with data
eighties_logbooks_stry <- biomass_grid(logbook_starry, 1981, 1989)
nineties_logbooks_stry <- biomass_grid(logbook_starry, 1990, 2001)
thousands_logbooks_stry <- biomass_grid(logbook_starry, 2002, 2009)
tens_logbooks_stry <- biomass_grid(logbook_starry, 2010, 2017)

### Survey
dev.new(width = 4, height = 10)
biomass_fillpts(survey_starry, 1981, 1989)
dev.new(width = 4, height = 10)
biomass_fillpts(survey_starry, 1990, 2001)
dev.new(width = 4, height = 10)
biomass_fillpts(survey_starry, 2002, 2009)
dev.new(width = 4, height = 10)
biomass_fillpts(survey_starry, 2010, 2017)

# Fill decade grids with data
eighties_survey_stry <- biomass_grid(survey_starry, 1980, 1989)
nineties_survey_stry <- biomass_grid(survey_starry, 1990, 2001)
thousands_survey_stry <- biomass_grid(survey_starry, 2002, 2009)
tens_survey_stry <- biomass_grid(survey_starry, 2010, 2018)

### Spatial indicators
# local index of collocation per decade for Starry Flounder overall
loc_collocfn(eighties_logbooks_stry, eighties_survey_stry)
loc_collocfn(nineties_logbooks_stry, nineties_survey_stry)
loc_collocfn(thousands_logbooks_stry, thousands_survey_stry)
loc_collocfn(tens_logbooks_stry, tens_survey_stry)

# Fill cell values for LIC
eighties_starry <- spatial_lic(eighties_logbooks_stry, logbook_starry,
                                 eighties_survey_stry, survey_starry, 1980, 1989)
nineties_starry <- spatial_lic(nineties_logbooks_stry, logbook_starry,
                                 nineties_survey_stry, survey_starry, 1990, 2001)
thousands_starry <- spatial_lic(thousands_logbooks_stry, logbook_starry,
                                  thousands_survey_stry, survey_starry, 2002, 2009)
tens_starry <- spatial_lic(tens_logbooks_stry, logbook_starry,
                             tens_survey_stry, survey_starry, 2010, 2018)

# Determine maximum value to scale appropriately
max(eighties_starry, na.rm = T) # max
max(nineties_starry, na.rm = T)
max(thousands_starry, na.rm = T)
max(tens_starry, na.rm = T)

# Mean and SD
mean(eighties_starry, na.rm = T)
mean(nineties_starry, na.rm = T)
mean(thousands_starry, na.rm = T)
mean(tens_starry, na.rm = T)

sd(eighties_starry, na.rm = T)
sd(nineties_starry, na.rm = T)
sd(thousands_starry, na.rm = T)
sd(tens_starry, na.rm = T)

# Create map
pdf("../final_figs/fish_res_fig_tables/starry_overlap.pdf",
    width = 7.5,
    height = 9)
par(mfrow = c(1, 2),
    family = "serif",
    mar = c(4, 5, 3, .3) + .1)
lic_map(eighties_starry, eighties_starry, "Starry Flounder 1980s",
                 bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = hcl.colors(100, "YlOrRd", rev = F),
           legend.shrink = 0.2,
           smallplot = c(.73, .79, .11, .25),
           legend.cex = 1.4,
           axis.args = list(cex.axis = 1.2),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(eighties_starry, na.rm = T)),
           legend.args = list("LIC",
                              side = 2, cex = 1.2))
lic_map(tens_starry, eighties_starry, "Starry Flounder 2010s",
                 bathy_dat, bathy_mat)
dev.off()

# Supplement figure
pdf("../final_figs/fish_res_fig_tables/starry_overlap_supplement.pdf",
    width = 15,
    height = 8.5)
par(mfrow = c(1, 4),
    family = "serif",
    mar = c(5.5, 6, 3, .3) + .1)
lic_sup(eighties_starry, eighties_starry, "Starry Flounder 1980s",
        bathy_dat, bathy_mat)
image.plot(legend.only = T,
           col = hcl.colors(100, "YlOrRd", rev = F),
           legend.shrink = 0.2,
           smallplot = c(.76, .82, .11, .25),
           legend.cex = 1.6,
           axis.args = list(cex.axis = 1.8),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(0, max(eighties_starry, na.rm = T)),
           legend.args = list("LIC",
                              side = 2, cex = 1.2))
lic_sup(nineties_starry, eighties_starry, "Starry Flounder 1990s",
        bathy_dat, bathy_mat)
lic_sup(thousands_starry, eighties_starry, "Starry Flounder 2000s",
        bathy_dat, bathy_mat)
lic_sup(tens_starry, eighties_starry, "Starry Flounder 2010s",
        bathy_dat, bathy_mat)
dev.off()
