# Title: NMFS Groundfish Survey Data Cleaning
# Purpose: Filter and prepare survey data for use in analyses
# Date created: 10/22/2020

# Load libraries ----
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(reshape2)
library(data.table)
library(zoo)
library(lubridate)
library(RANN)

# Load data
setwd('C:/Users/howar/Documents/Oregon State/ORnearshore_groundfish/code')
load('../data/Environmental_data/PDO')
NPGO <- read.csv("../data/Environmental_data/NPGO.csv")
species_data <- read.csv("../data/NMFS_data/species_data_raw.csv", na.strings = c("", NA))
scientific_name_id <- read.csv("../data/NMFS_data/Scientific_id.csv", header = T)
species_pacfin <- read.csv("../data/NMFS_data/Fish_class.csv", header = T)
trawl_data <- read.csv("../data/NMFS_data/trawl_characteristics.csv", header = T)
load('../data/Environmental_data/lith_df.Rdata')
load('../data/Environmental_data/grain_df.Rdata')

# Clean survey sample and species data ----
# Limit to less than 200 m depth
species_depth <- filter(species_data, depth_m <= 200.0)

#removed certain columns, add classification (if wanted)
species_limited <- dplyr::select(species_depth,
                          common_name,
                          date_yyyymmdd,
                          depth_m,
                          latitude_dd,
                          longitude_dd,
                          performance,
                          project,
                          scientific_name,
                          total_catch_numbers,
                          total_catch_wt_kg,
                          cpue_numbers_per_ha_der,
                          cpue_kg_per_ha_der,
                          year,
                          trawl_id )

# Remove unidentified species
species_identified <- species_limited[!is.na(species_limited$scientific_name), ]

# Remove extraneous surveys
# Limit to Oregon fishing grounds
unique(species_identified$project)
species_oregon <- filter(species_identified,
                         latitude_dd >= 41.0000,
                         latitude_dd <= 48.0000,
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
                         scientific_name != "Merluccius productus",
                         performance != "Unsatisfactory")

# Log + 1 of CPUE weight and number
species_oregon$lncpue <- log(species_oregon$cpue_kg_per_ha_der + 1)
species_oregon$lncpue_n <- log(species_oregon$cpue_numbers_per_ha_der + 1)

# Create separate fish and invertebrate datasets, keep only fish
species_oregon$ID <- scientific_name_id$ID[match(species_oregon$scientific_name,
                                                 scientific_name_id$scientific_name)]
OR_fish <- filter(species_oregon, ID == 'fish')
OR_fish$ID <- NULL
OR_fish$Class <- species_pacfin$Class[match(OR_fish$scientific_name,
                                                 species_pacfin$scientific_name)]
OR_fish$cpue_kg_per_ha_der[is.na(OR_fish$cpue_kg_per_ha_der)] <- 0

#rename columns
names(OR_fish)[names(OR_fish)=="cpue_kg_per_ha_der"] <- "cpue_kg"
names(OR_fish)[names(OR_fish)=="cpue_numbers_per_ha_der"] <- "cpue_numbers"
names(OR_fish)[names(OR_fish)=="latitude_dd"] <- "latitude"
names(OR_fish)[names(OR_fish)=="longitude_dd"] <- "longitude"
names(OR_fish)[names(OR_fish)=="pacfin_spid"] <- "pacfin"

# Save to text file
save(OR_fish, file = "../data/NMFS_data/OR_fish")
write.table(OR_fish, "../data/NMFS_data/OR_fish.txt",sep="\t", row.names=FALSE)

#Clean up trawl dataset
trawl_data_lat <- filter(trawl_data, latitude_dd >= 41.0000, latitude_dd <= 48.0000,
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

# add year-day by creating fake dates to set numbers for seasonality
trawl_data$date <- ymd(trawl_data$date_yyyymmdd)
trawl_data$month_day <- format(trawl_data$date, '%m-%d')
trawl_data$year_julian <- rep(1970,nrow(trawl_data))
trawl_data$date_fake <- paste(trawl_data$year_julian, trawl_data$month_day, sep = "-")
trawl_data$day <- format(trawl_data$date, '%d')
trawl_data$date_fake <- as.Date(trawl_data$date_fake)
trawl_data$julian <- julian((trawl_data$date_fake), format = "%Y%m%d")
# remove extra unnecessary columns
trawl_data <- trawl_data[-c(5:11, 13, 37:40)]

# add climate indices
# PDO
setnames(PDO, old=c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"),
         new = c("01","02","03","04","05","06", "07", "08", "09", "10", "11", "12"))
PDO <- as.data.frame(PDO)
rownames(PDO) <- PDO$YEAR
PDO$YEAR <- NULL
PDO <- PDO[-c(120),]
pdo <- as.matrix(PDO)

months <- format(seq.Date(as.Date("2018-01-01"), as.Date("2018-12-12"), by = "month"), format = "%b")
PDO$Year <- seq(1900, 2018)
PDO <- melt(PDO, id.vars = "Year")
PDO$Date <- as.Date(paste(PDO$Year, PDO$variable, "01", sep = "-"),
                    format = ("%Y-%b"))
PDO <- PDO[order(PDO$Date),]
PDO$Date <- as.yearmon(paste(PDO$Year, PDO$variable), "%Y%m")
PDO$date <- as.Date(PDO$Date)
PDO$date <- format(PDO$date, '%Y%m')

PDO <- PDO[-c(1)]
PDO <- PDO[-c(1)]
PDO <- PDO[-c(2)]

save(PDO, file = '../data/Environmental_data/PDO_column')

trawl_data$year_month <- format(trawl_data$date, '%Y%m')
match_id<-match(trawl_data$year_month, PDO$date)
trawl_data$PDO<-PDO$value[match_id]
trawl_data$PDO[is.na(trawl_data$PDO)]<-0

NPGO <- as.data.frame(NPGO)
colnames(NPGO)[1] <- "NPGO"
match_id<-match(trawl_data$year_month,NPGO$Date)
trawl_data$NPGO<-NPGO$NPGO[match_id]

# Add the grain size to trawl data
trawl_data[, c(33, 34)] <- as.data.frame(RANN::nn2(grain_df[, c(3, 2)], trawl_data[, c(5, 6)], k = 1))
trawl_data[, 35] <- grain_df[trawl_data[, 33], 1]
trawl_data <- trawl_data[- c(33, 34)]
colnames(trawl_data)[33] <- "grain_size"

# Add the lithology size to trawl data
trawl_data[, c(34, 35)] <- as.data.frame(RANN::nn2(lith_df[, c(3, 2)], trawl_data[, c(5, 6)], k = 1))
trawl_data[, 36] <- lith_df[trawl_data[, 34], 1]
trawl_data <- trawl_data[- c(34, 35)]
colnames(trawl_data)[34] <- "lithology"

save(trawl_data, file = "../data/NMFS_data/trawl_data")
write.table(trawl_data,"../data/NMFS_data/trawl_data.txt",sep="\t",row.names=FALSE)


# Create separate annual and triennial datasets ----
# Filter by FMP species
commercial_filter <- filter(OR_fish, Class != "non-FMP")
match_id<-match(commercial_filter$trawl_id,trawl_data$trawl_id)
commercial_filter$catch_kg<-trawl_data$vertebrate_weight_kg[match_id]
commercial_species <- commercial_filter

# Create separate survey datasets (get rid of 2004 triennial survey)
annual <- filter(commercial_species, project != "Groundfish Triennial Shelf Survey", year > "2002")
annual_trawl <- filter(trawl_data, project != "Groundfish Triennial Shelf Survey", year > "2002")
triennial <- filter(commercial_species, project != "Groundfish Slope and Shelf Combination Survey")
triennial_trawl <- filter(trawl_data, project != "Groundfish Slope and Shelf Combination Survey")
save(annual_trawl, file = "../data/NMFS_data/annual_tows")
save(triennial_trawl, file = "../data/NMFS_data/triennial_tows")
save(annual, file = "../data/NMFS_data/annual_samples")
save(triennial, file = "../data/NMFS_data/triennial_samples")

# Annual: Remove species that appear in < 1% of trawls (47 species extracted, excluded Raja and Rajiformes)
annual <- as.data.frame(annual[order(annual$trawl_id),])
annual_matrix <- create.matrix(annual, tax.name = "scientific_name",
                               locality = "trawl_id", time.col = NULL, time = NULL,
                               abund = F, abund.col = "lncpue")
annual_matrix <- as.data.frame(annual_matrix)
annual_matrix$sum <- rowSums(annual_matrix)
annual_matrix$percentage <- ((annual_matrix$sum)/2405)*100
annual_filtered <- annual_matrix[annual_matrix$percentage >= 1,]
annual_filtered <- rownames_to_column(annual_filtered, "scientific_name")
annual_vec <- as.vector(annual_filtered['scientific_name'])
filter <- c("Anoplopoma fimbria", "Atheresthes stomias",
            "Bathyraja kincaidii", "Citharichthys sordidus","Eopsetta jordani", "Gadus macrocephalus",
            "Glyptocephalus zachirus", "Hexagrammos decagrammus", "Hippoglossoides elassodon", "Hydrolagus colliei",
            "Isopsetta isolepis", "Lepidopsetta bilineata", "Microstomus pacificus","Ophiodon elongatus","Parophrys vetulus",
            "Platichthys stellatus", "Pleuronichthys decurrens", "Psettichthys melanostictus", "Raja binoculata", "Raja inornata",
            "Raja rhina","Raja stellulata","Sebastes alutus","Sebastes babcocki","Sebastes brevispinis", "Sebastes chlorostictus",
            "Sebastes crameri","Sebastes diploproa","Sebastes elongatus","Sebastes entomelas","Sebastes flavidus","Sebastes goodei",
            "Sebastes helvomaculatus","Sebastes jordani","Sebastes maliger","Sebastes paucispinis","Sebastes pinniger",
            "Sebastes proriger","Sebastes ruberrimus","Sebastes saxicola","Sebastes sp. (aleutianus / melanostictus)",
            "Sebastes wilsoni","Sebastes zacentrus","Sebastolobus alascanus","Squalus suckleyi")
annual_filter <- filter(annual, scientific_name %in% filter)
save(annual_filter, file = "../data/NMFS_data/annual_filtered")
write.table(annual_filter,"../data/NMFS_data/annual_filtered.txt",sep="\t",row.names=T)


# Triennial: Remove species that appear in < 1% of trawls (44 species extracted, excluded Sebastes sp.)
triennial <- triennial[order(triennial$trawl_id),]
triennial_matrix <- create.matrix(triennial, tax.name = "scientific_name",
                                  locality = "trawl_id", time.col = NULL, time = NULL,
                                  abund = F, abund.col = "lncpue")
triennial_matrix <- as.data.frame(triennial_matrix)
triennial_matrix$sum <- rowSums(triennial_matrix)
triennial_matrix$percentage <- ((triennial_matrix$sum)/1942)*100
triennial_filtered <- triennial_matrix[triennial_matrix$percentage >= 1,]
triennial_filtered <- rownames_to_column(triennial_filtered, "scientific_name")
triennial_vec <- as.vector(triennial_filtered['scientific_name'])
filter <- c("Anoplopoma fimbria", "Atheresthes stomias",
            "Bathyraja kincaidii", "Citharichthys sordidus","Eopsetta jordani", "Gadus macrocephalus",
            "Glyptocephalus zachirus", "Hexagrammos decagrammus", "Hippoglossoides elassodon", "Hydrolagus colliei",
            "Isopsetta isolepis", "Lepidopsetta bilineata", "Microstomus pacificus","Ophiodon elongatus","Parophrys vetulus",
            "Pleuronichthys decurrens", "Psettichthys melanostictus", "Raja binoculata",
            "Raja rhina","Sebastes alutus","Sebastes babcocki","Sebastes brevispinis", "Sebastes chlorostictus",
            "Sebastes crameri","Sebastes diploproa","Sebastes elongatus","Sebastes entomelas","Sebastes flavidus","Sebastes goodei",
            "Sebastes helvomaculatus","Sebastes jordani","Sebastes melanops","Sebastes paucispinis","Sebastes pinniger",
            "Sebastes proriger","Sebastes ruberrimus","Sebastes saxicola","Sebastes sp. (aleutianus / melanostictus)",
            "Sebastes wilsoni","Sebastes zacentrus","Sebastolobus alascanus","Squalus suckleyi")
triennial_filter <- filter(triennial, scientific_name %in% filter)
save(triennial_filter, file = "../data/NMFS_data/triennial_filtered")
write.table(triennial_filter,"../data/NMFS_data/triennial_filtered.txt",sep="\t",row.names=T)


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
