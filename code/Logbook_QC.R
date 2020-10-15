library(dplyr)
library(tidyr)
library(reshape)
library(reshape2)
library(MASS)
library(marmap)

setwd("C:/Users/howar/Documents/Oregon State/ORnearshore_groundfish/code/")
load("logbooks")
load("fish_tickets")
load("vessel_data")

log_location <- logbooks[-c(4, 5, 10, 17:20, 22, 23, 25:28, 30, 31, 33, 34, 38)]

# log_location$lat <- ifelse(is.na(log_location$SET_LAT),
#                           log_location$UP_LAT,log_location$SET_LAT)
# log_location$lon <- ifelse(is.na(log_location$SET_LONG),
                           # log_location$UP_LONG,log_location$SET_LONG)

log_location$lat <- log_location$SET_LAT
log_location$lon <- -abs((log_location$SET_LONG))

# Filter logbooks to depths shallower than 110 fathoms and within OR nearshore
logbook_depth <-
        filter(log_location,
               lat >= 42.0000,
               lat <= 47.0000,
               lon <= -123.9,
               lon >= -125)
logbook_depth$year <- format(as.Date(logbook_depth$TOWDATE, format = "%Y-%m-%d"), '%Y')
logbook_depth$year <- as.numeric(as.character((logbook_depth$year)))
logbook_depth$year[is.na(logbook_depth$year)] <- 0
logbook_depth <- logbook_depth[!(logbook_depth$year == 0), ]

#Build LOESS model: play with SPAN and DEGREE of the LOESS function to find the best fit
bathy.dat <- read.table('etopo1.xyz', sep = '')
names(bathy.dat) <- c('lon', 'lat', 'depth')
bathy.dat$depth[bathy.dat$depth > 0] <- NA
depth.loess <- loess(depth ~ lon * lat,
                     span = 0.01,
                     degree = 2,
                     data = bathy.dat)
summary(lm(depth.loess$fitted ~ bathy.dat$depth))

#Predict depth on the new grid
logbook_depth$depth.pred <- predict(depth.loess, newdata = logbook_depth)
#remove points on land and in extremely unreasonably shallow water
logbook_depth <- logbook_depth[logbook_depth$depth.pred <= -10, ]
#remove depths above 200 meters
logbook_final <- logbook_depth[logbook_depth$depth.pred >= -200, ]

# Remove the large footrope gear
logbook_final <- logbook_final[!(logbook_final$GEAR == 391), ]
logbook_final <- logbook_final[!is.na(logbook_final$GEAR), ]

#create unique ID for each haul
logbook_final$Trawl_ID <-
        paste(logbook_final$TICKET, logbook_final$NTOW, sep = ".")
#logbook_final$Trawl_ID <- as.numeric(as.factor(logbook_final$Trawl_ID))
logbook_final$Trawl_ID[is.na(logbook_final$Trawl_ID)] <- 0
logbook_final <- logbook_final[!(logbook_final$Trawl_ID == 0), ]

# Reshape to get a species column for each haul
logbook_select <-
        logbook_final %>% dplyr::select(Trawl_ID, c(23:131))
logbook_hauls <-
        logbook_final %>% dplyr::select(Trawl_ID, c(1:22), c(132:141))
logbook_relevant <-
        logbook_hauls %>% dplyr::select(Trawl_ID,
                                        TICKET,
                                        lon,
                                        lat,
                                        depth.pred,
                                        TOWDATE,
                                        year,
                                        DURATION,
                                        GEAR,
                                        DOCNUM)
logbook_species <- melt(logbook_select, id = "Trawl_ID")

# put the two datasets back together
match_id <- match(logbook_species$trawl_ID, logbook_relevant$trawl_ID)
logbook_species$TripID <- logbook_relevant$TripID[match_id]
combined <-
        merge(logbook_relevant,
              logbook_species,
              by = "Trawl_ID",
              all = TRUE)
colnames(combined)[5] <- "depth"
colnames(combined)[11] <- "species"
colnames(combined)[12] <- "species_weight"

# Filter out trawls still on land (identified through maps of each year)
combined_filtered <-
        filter(combined,
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
combined <- combined_filtered
save(combined, file = "logbooks_corrected")
###########Fish Tickets############################################################################################################################33

# Match fish ticket species numbers to species codes in logbooks_corrected
combined <- logbook_lon_filtered

match_id <- match(fish_tickets$TICKET, combined$TICKET)
fish_tickets$lon <- combined$lon[match_id]
fish_tickets$lat <- combined$lat[match_id]
fish_tickets$depth <- combined$depth.new[match_id]
fish_tickets$lat[is.na(fish_tickets$lat)] <- 0
tickets_corrected <- fish_tickets[!(fish_tickets$lat == 0), ]

species_index <-
        read.csv("/Users/howar/Documents/Oregon State/Thesis/Logbook Data/Logbook/species-code_index.csv",
                header = T)
match_id <- match(tickets_corrected$SP_CODE, species_index$?..ODFW)
tickets_corrected$PACFIN <- species_index$PACFIN[match_id]
tickets_corrected$scientific_name <-
        species_index$scientific_name[match_id]
tickets_corrected$common_name <- species_index$common_name[match_id]
tickets_corrected$code <- species_index$code[match_id]

# Filter out non-FMP species
FMP_tickets <-
        tickets_corrected[!(tickets_corrected$code == "non-FMP"), ]
save(FMP_tickets, file  = "fish_tickets_final")
# Match vessels to their license numbers




#### Raw data map
logbook_final <- logbook_final[logbook_final$depth.pred <= -15, ]

windows(width = 28, height = 18)
par(mfrow = c(1, 4))
plot(1, 1,
        xlim = range(logbook_final$lon, na.rm = TRUE) + c(-.5, .2),
        ylim = range(logbook_final$lat, na.rm = TRUE) + c(-0.2, .2),
        ylab = expression(paste("latitude (" ^ 0, 'N)')),
        xlab = expression(paste("longitude (" ^ 0, 'E)')),
        main = paste('Petrale Sole 1980s'))
map("worldHires",
    fill = T,
    col = "grey",
    add = T)
points(logbook_final$lon[logbook_final$year >= 1981 & logbook_final$year <= 1989],
       logbook_final$lat[logbook_final$year >= 1981 & logbook_final$year <= 1989],
       pch = ".",
       col = 'purple')
#contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=c(-10,-200),
#        labcex=0.7,add=T,col='black', labels = NULL, lwd = 2)

savePol = locator(40, type = "o")
base_polygon = data.frame(x = savePol$x, y = savePol$y)

plot(1, 1, xlim = range(logbook_final$lon, na.rm = TRUE) + c(-.5, .2),
        ylim = range(logbook_final$lat, na.rm = TRUE) + c(-.2, .2),
        ylab = expression(paste("latitude (" ^ 0, 'N)')),
        xlab = expression(paste("longitude (" ^ 0, 'E)')),
        main = paste('Petrale Sole 1990s'))
map("worldHires",
    fill = T,
    col = "grey",
    add = T)
points(logbook_final$lon[logbook_final$year >= 1990 & logbook_final$year <= 1999],
       logbook_final$lat[logbook_final$year >= 1990 & logbook_final$year <= 1999],
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

plot(1, 1, xlim = range(logbook_final$lon, na.rm = TRUE) + c(-.5, .2),
        ylim = range(logbook_final$lat, na.rm = TRUE) + c(-.2, .2),
        ylab = expression(paste("latitude (" ^ 0, 'N)')),
        xlab = expression(paste("longitude (" ^ 0, 'E)')),
        main = paste('Petrale Sole 2000s'))
map("worldHires",
    fill = T,
    col = "grey",
    add = T)
points(logbook_final$lon[logbook_final$year >= 2000 & logbook_final$year <= 2009],
       logbook_final$lat[logbook_final$year >= 2000 & logbook_final$year <= 2009],
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

plot(1, 1, xlim = range(logbook_final$lon, na.rm = TRUE) + c(-.5, .2),
        ylim = range(logbook_final$lat, na.rm = TRUE) + c(-.2, .2),
        ylab = expression(paste("latitude (" ^ 0, 'N)')),
        xlab = expression(paste("longitude (" ^ 0, 'E)')),
        main = paste('Petrale Sole 2010s'))
map("worldHires",
    fill = T,
    col = "grey",
    add = T)
points(logbook_final$lon[logbook_final$year >= 2010 & logbook_final$year <= 2017],
       logbook_final$lat[logbook_final$year >= 2010 & logbook_final$year <= 2017],
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
