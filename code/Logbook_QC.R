library(dplyr)
library(tidyr)
library(reshape)
library(reshape2)
library(MASS)
library(marmap)

##############################################################################################################
# Load data for filtering, cleaning, etc.
setwd("C:/Users/howar/Documents/Oregon State/ORnearshore_groundfish/code/")
load("../data/ODFW_data/logbooks")
load("../data/ODFW_data/fish_tickets")
load("../data/ODFW_data/vessel_data")

##############################################################################################################
# Reduce the size of the logbook dataset by removing unnecessary columns
logbooks_reduced <- logbooks[-c(4, 5, 10, 17:20, 22, 23, 25:28, 30, 31, 33, 34, 38)] # Removes DOCNUM, FUEL, RETURNTIME, SETTIME through LORAN locations, other set and up types, MINDEPTH

##############################################################################################################
# Change name of latitude and longitude to lat and lon, add negatives for longitude
logbooks_reduced$lat <- logbooks_reduced$SET_LAT
logbooks_reduced$lon <- -abs((logbooks_reduced$SET_LONG)) # Add negs

##############################################################################################################
# Filter logbooks to locations within approximately Oregon shelf fishing grounds (can look at minimum latitudes to get latitude values)
logbooks_shelf <- filter(logbooks_reduced,
               lat >= 42.0000,
               lat <= 47.0000,
               lon <= -123.9,
               lon >= -125)

##############################################################################################################
# Reformat dates and add a year column
logbooks_shelf$year <- format(as.Date(logbooks_shelf$TOWDATE, format = "%Y-%m-%d"), '%Y')
logbooks_shelf$year <- as.numeric(as.character((logbooks_shelf$year)))
logbooks_shelf$year[is.na(logbooks_shelf$year)] <- 0
logbooks_shelf <- logbooks_shelf[!(logbooks_shelf$year == 0), ]

##############################################################################################################
# Build LOESS model: play with SPAN and DEGREE of the LOESS function to find the best fit
bathy.dat <- read.table("../data/etopo1.xyz", sep = '') # Load NOAA ETOPO1
names(bathy.dat) <- c('lon', 'lat', 'depth')
bathy.dat$depth[bathy.dat$depth > 0] <- NA # Make points on land NA
depth.loess <- loess(depth ~ lon * lat,
                     span = 0.01,
                     degree = 2,
                     data = bathy.dat) # Build LOESS
# summary(lm(depth.loess$fitted ~ bathy.dat$depth)) # Will produce an error unless NA's are removed I think

# Predict depth on the new grid
logbooks_shelf$depth.pred <- predict(depth.loess, newdata = logbooks_shelf)

##############################################################################################################
# Filter dataset by depth, gear
# Remove points on land and in extremely unreasonably shallow water
logbooks_depth <- logbooks_shelf[logbooks_shelf$depth.pred <= -10, ]

# Remove depths above 200 meters
logbooks_depth <- logbooks_depth[logbooks_depth$depth.pred >= -200, ]

# Remove the large footrope gear
logbooks_gear <- logbooks_depth[!(logbooks_depth$GEAR == 391), ]
logbooks_gear <- logbooks_gear[!is.na(logbooks_gear$GEAR), ]

##############################################################################################################
# Create unique ID for each haul
logbooks_gear$Trawl_ID <-paste(logbooks_gear$TICKET,
                              logbooks_gear$NTOW,
                              sep = ".")
logbooks_trawls <- logbooks_gear[!is.na(logbooks_gear$Trawl_ID), ] # Remove NA's just in case

##############################################################################################################
# Reshape to get a species column for each haul
logbooks_species <- logbooks_trawls %>% dplyr::select(Trawl_ID,
                                                      c(23:131)) # Dataframe with just the species data
logbooks_characteristics <- logbooks_trawls %>% dplyr::select(Trawl_ID,
                                                              c(1:22),
                                                              c(132:141)) # Dataframe with just the haul characteristics
logbooks_characteristics <- logbook_characteristics %>% dplyr::select(Trawl_ID,
                                        TICKET,
                                        lon,
                                        lat,
                                        depth.pred,
                                        TOWDATE,
                                        year,
                                        DURATION,
                                        GEAR,
                                        DOCNUM)
logbooks_species <- melt(logbooks_species, id = "Trawl_ID")

###############################################################################################################
# Put the two datasets back together
match_id <- match(logbooks_species$trawl_ID,
                  logbooks_characteristics$trawl_ID)
logbooks_species$TripID <- logbooks_characteristics$TripID[match_id]
logbooks_expanded <- merge(logbooks_characteristics,
                           logbooks_species,
                           by = "Trawl_ID",
                           all = TRUE)
colnames(logbooks_expanded)[5] <- "depth"
colnames(logbooks_expanded)[11] <- "species"
colnames(logbooks_expanded)[12] <- "species_weight"

###############################################################################################################
#### Raw data map to filter out any extra points on land
windows(width = 28, height = 18)
par(mfrow = c(1, 4))
plot(1, 1,
     xlim = range(logbooks_expanded$lon, na.rm = TRUE) + c(-.5, .2),
     ylim = range(logbooks_expanded$lat, na.rm = TRUE) + c(-0.2, .2),
     ylab = "latitude °N",
     xlab = "longitude °W",
     main = '1980s')
map("worldHires",
    fill = T,
    col = "grey",
    add = T)
points(logbooks_expanded$lon[logbooks_expanded$year >= 1980 & logbooks_expanded$year <= 1989],
       logbooks_expanded$lat[logbooks_expanded$year >= 1981 & logbooks_expanded$year <= 1989],
       pch = ".",
       col = 'purple')
#contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=c(-10,-200),
#        labcex=0.7,add=T,col='black', labels = NULL, lwd = 2)

plot(1, 1, xlim = range(logbooks_expanded$lon, na.rm = TRUE) + c(-.5, .2),
     ylim = range(logbooks_expanded$lat, na.rm = TRUE) + c(-.2, .2),
     ylab = expression(paste("latitude (" ^ 0, 'N)')),
     xlab = expression(paste("longitude (" ^ 0, 'E)')),
     main = paste('Petrale Sole 1990s'))
map("worldHires",
    fill = T,
    col = "grey",
    add = T)
points(logbooks_expanded$lon[logbooks_expanded$year >= 1990 & logbooks_expanded$year <= 1999],
       logbooks_expanded$lat[logbooks_expanded$year >= 1990 & logbooks_expanded$year <= 1999],
       pch = ".",
       col = 'purple')
contour( unique(bathy.dat$lon),
         sort(unique(bathy.dat$lat)),
         bathy.mat,
         levels = c(-10, -200),
         labcex = 0.7,
         add = T,
         col = 'black',
         labels = NULL,
         lwd = 2)

plot(1, 1, xlim = range(logbooks_expanded$lon, na.rm = TRUE) + c(-.5, .2),
     ylim = range(logbooks_expanded$lat, na.rm = TRUE) + c(-.2, .2),
     ylab = expression(paste("latitude (" ^ 0, 'N)')),
     xlab = expression(paste("longitude (" ^ 0, 'E)')),
     main = paste('Petrale Sole 2000s'))
map("worldHires",
    fill = T,
    col = "grey",
    add = T)
points(logbooks_expanded$lon[logbooks_expanded$year >= 2000 & logbooks_expanded$year <= 2009],
       logbooks_expanded$lat[logbooks_expanded$year >= 2000 & logbooks_expanded$year <= 2009],
       pch = ".",
       col = 'purple')

contour(unique(bathy.dat$lon),
        sort(unique(bathy.dat$lat)),
        bathy.mat,
        levels = c(-10, -200),
        labcex = 0.7,
        add = T,
        col = 'black',
        labels = NULL,
        lwd = 2)

plot(1, 1, xlim = range(logbooks_expanded$lon, na.rm = TRUE) + c(-.5, .2),
     ylim = range(logbooks_expanded$lat, na.rm = TRUE) + c(-.2, .2),
     ylab = expression(paste("latitude (" ^ 0, 'N)')),
     xlab = expression(paste("longitude (" ^ 0, 'E)')),
     main = paste('Petrale Sole 2010s'))
map("worldHires",
    fill = T,
    col = "grey",
    add = T)
points(logbooks_expanded$lon[logbooks_expanded$year >= 2010 & logbooks_expanded$year <= 2017],
       logbooks_expanded$lat[logbooks_expanded$year >= 2010 & logbooks_expanded$year <= 2017],
       pch = ".",
       col = 'purple')
contour( unique(bathy.dat$lon),
         sort(unique(bathy.dat$lat)),
         bathy.mat,
         levels = c(-10, -200),
         labcex = 0.7,
         add = T,
         col = 'black',
         labels = NULL,
         lwd = 2)
###############################################################################################################
# Filter out trawls still on land (identified through maps of each decade)
logbooks_final <- filter(logbooks_expanded,
                Trawl_ID != 5000406.1 &
                Trawl_ID != 5000406.3 &
                Trawl_ID != 5000422.1 &
                Trawl_ID != 5000454.1 &
                Trawl_ID != 5000480.1 &
                Trawl_ID != 5000500.1 &
                Trawl_ID != 5000522.1 &
                Trawl_ID != 5000563.1 &
                Trawl_ID != 5000563.3 &
                Trawl_ID != 5000589.1 &
                Trawl_ID != 5000589.3 &
                Trawl_ID != 5000589.5 &
                Trawl_ID != 5000605.1 &
                Trawl_ID != 5000605.3 &
                Trawl_ID != 5000679.5 &
                Trawl_ID != 5000679.6 &
                Trawl_ID != 5000679.7 &
                Trawl_ID != 5000691.1 &
                Trawl_ID != 5000708.11 &
                Trawl_ID != 5000708.5 &
                Trawl_ID != 5000708.6 &
                Trawl_ID != 5000708.7 &
                Trawl_ID != 5000708.9 &
                Trawl_ID != 5000736.1 &
                Trawl_ID != 5000798.1 &
                Trawl_ID != 5000819.1 &
                Trawl_ID != 5000819.2 &
                Trawl_ID != 5000819.3 &
                Trawl_ID != 5000819.4 &
                Trawl_ID != 5000819.5 &
                Trawl_ID != 5000845.1 &
                Trawl_ID != 5000845.2 &
                Trawl_ID != 5000845.3 &
                Trawl_ID != 5000845.4 &
                Trawl_ID != 5000845.5 &
                Trawl_ID != 5000845.6 &
                Trawl_ID != 5000845.7 &
                Trawl_ID != 5000871.1 &
                Trawl_ID != 5000874.1 &
                Trawl_ID != 5000874.3 &
                Trawl_ID != 5000893.10 &
                Trawl_ID != 5000893.11 &
                Trawl_ID != 5000893.12 &
                Trawl_ID != 5000893.8 &
                Trawl_ID != 5000893.9 &
                Trawl_ID != 5000904.1 &
                Trawl_ID != 5000904.2 &
                Trawl_ID != 5000904.3 &
                Trawl_ID != 5000970.1 &
                Trawl_ID != 5000970.3 &
                Trawl_ID != 5000976.1 &
                Trawl_ID != 5000984.1 &
                Trawl_ID != 5000984.3 &
                Trawl_ID != 5001000.1 &
                Trawl_ID != 5001000.3 &
                Trawl_ID != 5001005.1 &
                Trawl_ID != 5001005.3 &
                Trawl_ID != 5001035.1 &
                Trawl_ID != 5001035.3 &
                Trawl_ID != 5001056.1 &
                Trawl_ID != 5001059.1 &
                Trawl_ID != 5001070.1 &
                Trawl_ID != 5001114.3 &
                Trawl_ID != 5001114.4 &
                Trawl_ID != 5001114.5 &
                Trawl_ID != 5001139.11 &
                Trawl_ID != 5001139.6 &
                Trawl_ID != 5001139.7 &
                Trawl_ID != 5001201.1 &
                Trawl_ID != 5001213.1 &
                Trawl_ID != 5001213.9 &
                Trawl_ID != 5001224.1 &
                Trawl_ID != 5001254.1 &
                Trawl_ID != 5001254.3 &
                Trawl_ID != 5001266.10 &
                Trawl_ID != 5001266.4 &
                Trawl_ID != 5001266.5 &
                Trawl_ID != 5001266.6 &
                Trawl_ID != 5001301.1 &
                Trawl_ID != 5001301.2 &
                Trawl_ID != 5001301.3 &
                Trawl_ID != 5001323.1 &
                Trawl_ID != 5001323.2 &
                Trawl_ID != 5001323.3 &
                Trawl_ID != 5001345.1 &
                Trawl_ID != 5001378.1 &
                Trawl_ID != 5001378.3 &
                Trawl_ID != 5001407.1 &
                Trawl_ID != 5001407.3 &
                Trawl_ID != 5001418.1 &
                Trawl_ID != 5001418.3 &
                Trawl_ID != 5001458.7 &
                Trawl_ID != 5000406.2 &
                Trawl_ID != 5000422.2 &
                Trawl_ID != 5000422.3 &
                Trawl_ID != 5000422.4 &
                Trawl_ID != 5000454.2 &
                Trawl_ID != 5000480.2 &
                Trawl_ID != 5000500.2 &
                Trawl_ID != 5000500.3 &
                Trawl_ID != 5000500.4 &
                Trawl_ID != 5000522.2 &
                Trawl_ID != 5000522.3 &
                Trawl_ID != 5000522.4 &
                Trawl_ID != 5000589.2 &
                Trawl_ID != 5000589.4 &
                Trawl_ID != 5000605.2 &
                Trawl_ID != 5000605.4 &
                Trawl_ID != 5000691.2 &
                Trawl_ID != 5000691.3 &
                Trawl_ID != 5000691.4 &
                Trawl_ID != 5000691.5 &
                Trawl_ID != 5000708.8 &
                Trawl_ID != 5000736.2 &
                Trawl_ID != 5000736.3 &
                Trawl_ID != 5000736.4 &
                Trawl_ID != 5000736.5 &
                Trawl_ID != 5000798.2 &
                Trawl_ID != 5000798.3 &
                Trawl_ID != 5000819.13 &
                Trawl_ID != 5000874.2 &
                Trawl_ID != 5000874.4 &
                Trawl_ID != 5000893.7 &
                Trawl_ID != 5001035.2 &
                Trawl_ID != 5001035.4 &
                Trawl_ID != 5001056.2 &
                Trawl_ID != 5001070.2 &
                Trawl_ID != 5001070.4 &
                Trawl_ID != 5001153.3 &
                Trawl_ID != 5001153.4 &
                Trawl_ID != 5001254.2 &
                Trawl_ID != 5001254.4 &
                Trawl_ID != 5001266.11 &
                Trawl_ID != 5001345.2 &
                Trawl_ID != 5001345.3 &
                Trawl_ID != 5001378.2 &
                Trawl_ID != 5001378.4 &
                Trawl_ID != 5001407.2 &
                Trawl_ID != 5001418.2 &
                Trawl_ID != 5000893.4 &
                Trawl_ID != 5000893.3 &
                Trawl_ID != 5000819.6 &
                Trawl_ID != 3063214.4 &
                Trawl_ID != 3057631.1 &
                Trawl_ID != 3095495.8 &
                Trawl_ID != 3090017.2 &
                Trawl_ID != 3090119.1 &
                Trawl_ID != 3108263.6 &
                Trawl_ID != 3106726.8 &
                Trawl_ID != 3095375.2 &
                Trawl_ID != 3095371.4 &
                Trawl_ID != 3093785.11 &
                Trawl_ID != 3093749.13 &
                Trawl_ID != 3095427.9 &
                Trawl_ID != 3104127.3 &
                Trawl_ID != 3142212.5)
save(logbooks_final, file = "../data/ODFW_data/logbooks_corrected")



########################################################################################################################
###########Fish Tickets#################################################################################################
# Match fish ticket species numbers to species codes in logbooks_corrected
combined <- logbook_lon_filtered
match_id <- match(fish_tickets$TICKET,
                  combined$TICKET)
fish_tickets$lon <- combined$lon[match_id]
fish_tickets$lat <- combined$lat[match_id]
fish_tickets$depth <- combined$depth.new[match_id]
fish_tickets$lat[is.na(fish_tickets$lat)] <- 0
tickets_corrected <- fish_tickets[!(fish_tickets$lat == 0), ]

species_index <-
        read.csv("/Users/howar/Documents/Oregon State/Thesis/Logbook Data/Logbook/species-code_index.csv",
                header = T)
match_id <- match(tickets_corrected$SP_CODE,
                  species_index$?..ODFW)
tickets_corrected$PACFIN <- species_index$PACFIN[match_id]
tickets_corrected$scientific_name <- species_index$scientific_name[match_id]
tickets_corrected$common_name <- species_index$common_name[match_id]
tickets_corrected$code <- species_index$code[match_id]

# Filter out non-FMP species
FMP_tickets <-
        tickets_corrected[!(tickets_corrected$code == "non-FMP"), ]
save(FMP_tickets, file  = "fish_tickets_final")
# Match vessels to their license numbers



