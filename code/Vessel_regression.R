library(dplyr)
library(tidyr)
library(reshape)
library(reshape2)
library(MASS)
library(marmap)

setwd("C:/Users/howar/Documents/Oregon State/Thesis/Logbook Data/Logbook")
load("vessel_data")
load("logbooks_corrected")
load("tickets_corrected")

# Remove vessel records missing length and horsepower
vessel_data <- vessel_data[!vessel_data$Length==0,]
vessel_data <- vessel_data[!vessel_data$Horsepower==0,]

# Create starry sole subset
subset_starry<-combined[combined$species=='STRY_ADJ',]

# Match vessel data (length, horsepower) to logbooks using document number
subset_starry <- subset_starry[!is.na(subset_starry$DOCNUM),]
colnames(subset_starry)[10] <- "Docnum"
match_id<-match(subset_starry$Docnum,vessel_data$Docnum)
subset_starry$vessel_length<-vessel_data$Length[match_id]
subset_starry$vessel_hp<-vessel_data$Horsepower[match_id]

# Filter to just survey months
subset_starry$month_day <- as.numeric(format(subset_starry$TOWDATE, '%m%d'))
subset_starry <- subset_starry[subset_starry$month_day>=517&subset_starry$month_day<= 929,]

### CPUE kg/hr
subset_starry$kg_caught <- subset_starry$species_weight*0.4535924
subset_starry <- subset_starry[!subset_starry$DURATION==0,] # remove tows with no trawl duration
subset_starry$CPUE <- subset_starry$kg_caught/subset_starry$DURATION
subset_starry <- subset_starry[!is.na(subset_starry$CPUE),]
subset_starry$lncpue <- log(subset_starry$CPUE+1)
subset_starry <- subset_starry[subset_starry$depth<=-5,] # remove unreasonably shallow tows

# linear regression of length for each decade
subset_starry_eighties <- subset_starry[subset_starry$year<1990,]
subset_starry_nineties <- subset_starry[subset_starry$year>1989&subset_starry$year<2000,]
subset_starry_thousands <- subset_starry[subset_starry$year>1999&subset_starry$year<2010,]
subset_starry_teens <- subset_starry[subset_starry$year>2009&subset_starry$year<2018,]
eighties_lm <- lm(formula = lncpue ~ vessel_length, data = subset_starry_eighties)
nineties_lm <- lm(formula = lncpue ~ vessel_length, data = subset_starry_nineties)
thousands_lm <- lm(formula = lncpue ~ vessel_length, data = subset_starry_thousands)
teens_lm <- lm(formula = lncpue ~ vessel_length, data = subset_starry_teens)

windows()
par(mfrow=c(2,2))
plot(x= subset_starry_eighties$vessel_length, y=subset_starry_eighties$lncpue, xlab = "Length", ylab = "lnCPUE", main = "1980s")
abline(eighties_lm, col = "red")
plot(x= subset_starry_nineties$vessel_length, y=subset_starry_nineties$lncpue, xlab = "Length", ylab = "lnCPUE", main = "1990s")
abline(nineties_lm, col = "red")
plot(x= subset_starry_thousands$vessel_length, y=subset_starry_thousands$lncpue, xlab = "Length", ylab = "lnCPUE", main = "2000s")
abline(thousands_lm, col = "red")
plot(x= subset_starry_teens$vessel_length, y=subset_starry_teens$lncpue, xlab = "Length", ylab = "lnCPUE", main = "2010s")
abline(teens_lm, col = "red")