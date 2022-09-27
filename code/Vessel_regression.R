library(dplyr)
library(here)

load(here("data/ODFW_data", "vessel_data"))
load(here("data/ODFW_data", "logbooks_corrected"))

# Remove vessel records missing length and horsepower
vessel_data <- vessel_data[!vessel_data$Length == 0, ]
vessel_data <- vessel_data[!vessel_data$Horsepower == 0, ]

# Create starry sole subset
subset_starry <- logbooks_final[logbooks_final$species == 'STRY_ADJ', ]

# Match vessel data (length, horsepower) to logbooks using document number
subset_starry <- subset_starry[!is.na(subset_starry$DOCNUM), ]
colnames(subset_starry)[10] <- "Docnum"
match_id <- match(subset_starry$Docnum, vessel_data$Docnum)
subset_starry$vessel_length <- vessel_data$Length[match_id]
subset_starry$vessel_hp <- vessel_data$Horsepower[match_id]

# Filter to just survey months
subset_starry$month_day <- as.numeric(format(subset_starry$TOWDATE, '%m%d'))
subset_starry <- subset_starry[subset_starry$month_day >= 517 &
                  subset_starry$month_day <= 929, ]

### CPUE kg/hr
subset_starry$kg_caught <- subset_starry$species_weight * 0.4535924
subset_starry <- subset_starry[!subset_starry$DURATION == 0, ] # remove tows with no trawl duration
subset_starry$CPUE <- subset_starry$kg_caught / subset_starry$DURATION
subset_starry <- subset_starry[!is.na(subset_starry$CPUE), ]
subset_starry$lncpue <- log(subset_starry$CPUE + 1)
subset_starry <- subset_starry[subset_starry$depth <= -5, ] # remove unreasonably shallow tows

# Decade subset
subset_starry_eighties <- filter(subset_starry,
                                 !is.na(vessel_length),
                                 !is.na(CPUE),
                                 year <= 1989)
subset_starry_nineties <- filter(subset_starry,
                                 !is.na(vessel_length),
                                 !is.na(CPUE),
                                 year >= 1990 & year <= 2001)
subset_starry_thousands <- filter(subset_starry,
                                 !is.na(vessel_length),
                                 !is.na(CPUE),
                                 year >= 2002 & year <= 2009)
subset_starry_teens <- filter(subset_starry,
                                  !is.na(vessel_length),
                                  !is.na(CPUE),
                                  year >= 2010 & year <= 2018)

par(mfrow = c(2, 2))
hist(subset_starry$lncpue,
     main = "log(CPUE + 1)")
hist(log(subset_starry$CPUE + mean(subset_starry$CPUE) * 0.1),
     main = "log(CPUE + 10% of mean CPUE)")
hist(sqrt(subset_starry$CPUE),
     main = "sqrt(CPUE)")
hist((subset_starry$CPUE ^ (1 / 4)),
     main = "4th rt(CPUE)")

# Sensitivity analysis
eighties_test <- nlsLM(log(x + CPUE) ~ vessel_length,
                       start = list(x = 1),
                       data = subset_starry_eighties) # not working for response variable

# Linear regression
regression_sensitivity <- function(subset){
  eighties_lm <- lm(formula = lncpue ~ vessel_length, data = subset)
  eighties_lm_log <- lm(formula = log(CPUE + (mean(CPUE) * 0.1)) ~ vessel_length, data = subset)
  eighties_lm_sqrt <- lm(formula = sqrt(CPUE) ~ vessel_length, data = subset)
  eighties_lm_frth <- lm(formula = (CPUE)^(1/4) ~ vessel_length, data = subset)
  lm_list <- list(eighties_lm, eighties_lm_log, eighties_lm_sqrt, eighties_lm_frth)
  return(lm_list)
  }

lm_plots <- function(subset, title, list){
  plot(x = subset$vessel_length,
       y = subset$lncpue,
       xlab = "Length",
       ylab = "log(CPUE + 1)",
       main = title)
  abline(list[[1]], col = "red")
  plot(x = subset$vessel_length,
       y = log(subset$CPUE + (mean(subset$CPUE) * 0.1)),
       xlab = "Length",
       ylab = "log(CPUE + 10% mean)",
       main = title)
  abline(list[[2]], col = "red")
  plot(x = subset$vessel_length,
       y = sqrt(subset$CPUE),
       xlab = "Length",
       ylab = "sqrt(CPUE)",
       main =  title)
  abline(list[[3]], col = "red")
  plot(x = subset$vessel_length,
       y = subset$CPUE ^ (1 / 4),
       xlab = "Length",
       ylab = "fourth root CPUE",
       main = title)
  abline(list[[4]], col = "red")
}

eighties_list <- regression_sensitivity(subset_starry_eighties)
nineties_list <- regression_sensitivity(subset_starry_nineties)
thousands_list <- regression_sensitivity(subset_starry_thousands)
teens_list <- regression_sensitivity(subset_starry_teens)

windows()
par(mfrow = c(2, 2))
lm_plots(subset_starry_eighties, "1980s", eighties_list)

windows()
par(mfrow = c(2, 2))
lm_plots(subset_starry_nineties, "1990s", nineties_list)

windows()
par(mfrow = c(2, 2))
lm_plots(subset_starry_thousands, "2000s", thousands_list)

windows()
par(mfrow = c(2, 2))
lm_plots(subset_starry_teens, "2010s", teens_list)
