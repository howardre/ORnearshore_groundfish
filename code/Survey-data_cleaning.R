library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(reshape2)
library(data.table)
library(zoo)
library(lubridate)

setwd('/Users/howar/Documents/Oregon State/Thesis/Data visualization')
load('PDO')
NPGO <- read.csv("NPGO.csv")

#renames blank spaces to NA
spp_depth_1 <- read.csv("~/Oregon State/Thesis/Data visualization/spp.total.csv", na.strings = c("", NA))

#limited to less than 200m depth
spp_depth <- filter(spp_depth_1, depth_m <= 200.0)

#removed certain columns, add classification (if wanted)
spp_select <- select(spp_depth, common_name, date_yyyymmdd, depth_m, latitude_dd, longitude_dd, performance, scientific_name, total_catch_numbers, total_catch_wt_kg, cpue_numbers_per_ha_der, cpue_kg_per_ha_der, year, trawl_id )

#removes unidentified spp
spp_named <- spp_depth[!is.na(spp_select$scientific_name),]

#narrow to Oregon coast
#narrow to about ~200m isobath
spp_OR1 <- filter(spp_named, latitude_dd >= 41.0000, latitude_dd <= 48.0000)

#remove hypoxia, slope, and video survey and hake (shelf survey = 2001, shelf & slope combo = post-2003 annual)
unique(spp_OR1$project)
spp_OR <- filter(spp_OR1,  project != "Groundfish Slope Survey",
                 project != "Hypoxia Study",
                 project != "Groundfish Shelf Survey",
                 project != "Triennial Shelf Groundfish Survey: Canada",
                 project != "Fishing Power Comparison Study",
                 project != "Opportunistic Bottom Sample",
                 project != "Fishing Gear Experiment",
                 project != "Santa Barbara Basin Study",
                 project != "Video Study",
                 project != "Nonstandard Sampling",
                 project != "Opportunistic Bottom Sampling",
                 project != "AFSC/RACE Slope Survey",
                 project != "Fishing Power Comparative Study",
                 program != "AFSC/RACE Triennial Shelf Survey (by NWFSC/FRAM)",
                 scientific_name != "Merluccius productus",
                 performance != "Unsatisfactory")

#log+1 of cpue weight and number
spp_OR$lncpue<-log(spp_OR$cpue_kg_per_ha_der + 1)
spp_OR$lncpue_n <- log(spp_OR$cpue_numbers_per_ha_der + 1)

#create separate fish and invertebrate datasets, keep only fish
Scientific_id <- read.csv("~/Oregon State/Thesis/Data visualization/Scientific_id.csv", header = T)
combined <- sort(union(levels(spp_OR$scientific_name), levels(Scientific_id$scientific_name)))
OR_labeled <- left_join(mutate(spp_OR, scientific_name=factor(scientific_name, levels=combined)),
                        mutate(Scientific_id, scientific_name=factor
                               (scientific_name, levels = combined)))
OR_fish1 <- filter(OR_labeled, ID=='fish')
OR_fish1$ID <- NULL
fish_pacfin <- read.csv("~/Oregon State/Thesis/Data visualization/Fish_class.csv", header = T)
combined1 <- sort(union(levels(OR_fish1$scientific_name), levels(fish_pacfin$scientific_name)))
OR_fish <- left_join(mutate(OR_fish1, scientific_name=factor(scientific_name, levels=combined1)),
                        mutate(fish_pacfin, scientific_name=factor
                               (scientific_name, levels = combined1)))
OR_fish$cpue_kg_per_ha_der[is.na(OR_fish$cpue_kg_per_ha_der)]<-0

#rename columns
names(OR_fish)[names(OR_fish)=="cpue_kg_per_ha_der"] <- "cpue_kg"
names(OR_fish)[names(OR_fish)=="cpue_numbers_per_ha_der"] <- "cpue_numbers"
names(OR_fish)[names(OR_fish)=="latitude_dd"] <- "latitude"
names(OR_fish)[names(OR_fish)=="longitude_dd"] <- "longitude"
names(OR_fish)[names(OR_fish)=="pacfin_spid"] <- "pacfin"

# remove unnecessary columns
OR_fish <- OR_fish[-c(1,9,10,11,14,15,17,18,19,22,23,25,28,31)]

#save to text file
save(OR_fish, file = "OR_fish")
write.table(OR_fish,"OR_fish.txt",sep="\t",row.names=FALSE)

#Clean up trawl dataset
trawl_dat <- read.csv("~/Oregon State/Thesis/Data visualization/trawl_characteristics.csv")
trawl_data_lat <- filter(trawl_dat, latitude_dd >= 41.0000, latitude_dd <= 48.0000, 
                         project != "Groundfish Slope Survey",
                         project != "Hypoxia Study",
                         project != "Groundfish Shelf Survey",
                         project != "Triennial Shelf Groundfish Survey: Canada",
                         project != "Fishing Power Comparison Study",
                         project != "Opportunistic Bottom Sample",
                         project != "Fishing Gear Experiment",
                         project != "Santa Barbara Basin Study",
                         project != "Video Study",
                         project != "Nonstandard Sampling",
                         project != "Opportunistic Bottom Sampling",
                         project != "AFSC/RACE Slope Survey",
                         project != "Fishing Power Comparative Study",
                         performance != 'Unsatisfactory')
trawl_data <- filter(trawl_data_lat, depth_hi_prec_m <= 200.0)
names(trawl_data)[names(trawl_data)=="area_swept_ha_der"] <- "area_swept"
names(trawl_data)[names(trawl_data)=="depth_hi_prec_m"] <- "depth_m"
names(trawl_data)[names(trawl_data)=="latitude_dd"] <- "latitude"
names(trawl_data)[names(trawl_data)=="longitude_dd"] <- "longitude"
names(trawl_data)[names(trawl_data)=="date_dim.year"] <- "year"
names(trawl_data)[names(trawl_data)=="temperature_at_gear_c_der"] <- "bottom_temp"

match_id<-match(trawl_data$trawl_id,OR_fish$trawl_id)
trawl_data$program<-OR_fish$program[match_id]
trawl_data <- filter(trawl_data, program != "AFSC/RACE Triennial Shelf Survey (by NWFSC/FRAM)")

# add year-day by creating fake dates to set numbers for seasonality
trawl_data$date <- ymd(trawl_data$date_yyyymmdd)
trawl_data$month_day <- format(trawl_data$date, '%m-%d')
trawl_data$year_julian <- rep(1970,nrow(trawl_data))
trawl_data$date_fake <- paste(trawl_data$year_julian, trawl_data$month_day, sep = "-")
trawl_data$day <- format(trawl_data$date, '%d')
trawl_data$date_fake <- as.Date(trawl_data$date_fake)
trawl_data$julian <- julian((trawl_data$date_fake), format = "%Y%m%d")
# remove extra unnecessary columns
trawl_data <- trawl_data[-c(3,5,6,7,8,9,10,11,13, 15, 16, 17, 18,19,20,22,23,24,25,27,29,31,32,33,34,35,38,39,40)]

# add climate indices
trawl_data$year_month <- format(trawl_data$date, '%Y%m')
match_id<-match(trawl_data$year_month,PDO$date)
trawl_data$PDO<-PDO$value[match_id]
trawl_data$PDO[is.na(trawl_data$PDO)]<-0

NPGO <- as.data.frame(NPGO)
colnames(NPGO)[1] <- "NPGO"
match_id<-match(trawl_data$year_month,NPGO$Date)
trawl_data$NPGO<-NPGO$NPGO[match_id]
save(trawl_data, file = "trawl_data")
write.table(trawl_data,"trawl_data.txt",sep="\t",row.names=FALSE)

--------------------------------------#######################------------------------------------------------------

# Exploratory plots and scientific name list

#remove duplicates
sp.cpue <- tapply(spp_OR$cpue_numbers_per_ha_der, spp_OR$scientific_name, mean)
length(sp.cpue)
barplot(sort(sp.cpue, decreasing = T)[1:10])
sc.name <- unique(spp_OR$scientific_name)
length(sc.name)
write.table(sc.name,'sc_name.txt')
#df = 

match.id <- match(id, spp_name)
classification <- id[match.id]
match.id

id.match<-match(spp.sub$scientific_name,scientific_name)
length(id.match)
dim(spp.sub)
range(id.match,na.rm=T)
is.na(id.match)

match(classification, spp_OR, nomatch = NA_integer_, incomparables = NULL)
#classes_OR <- merge(catch_OR, inverts, all=TRUE)
#inverts_OR <- filter(classes_OR, class == 1)
#blank_spp <- filter(classes_OR, is.na(common_name))
#unique(blank_spp$scientific_name)
#unique()
sc.name<-unique(spp.sub$scientific_name)
head(sc.name)
length(sc.name)
write.table(sc.name,'sc_name.txt')

#match.id use to link files, code fish vs. inverts

#example of using tapply for calculating average species cpue
sp.cpue <- tapply(inverts_OR$cpue_numbers_per_ha_der, inverts_OR$scientific_name, mean)
length(sp.cpue)
barplot(sort(sp.cpue, decreasing = T)[1:10])
#use ggplot to do error bars
#check out marmap