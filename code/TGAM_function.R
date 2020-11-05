# SHOUld work for fishes with year in the gam. Couldn't figure out indexing for line
# 23 though
# Go through each year to find the threshold year for sample ----
get_tgam_yr_aic <- function(df, yr) {
  df$thr <- ifelse(df$year <= yr, 'before', 'after')
  aic_year <-  gam(pres ~ factor(year) +
                     s(longitude, latitude, by = factor(thr)) +
                     s(julian),
                   data = df,
                   family = binomial)$aic
  return_list <- list(aic_year, yr)
}

# Returns AIC difference, best year, best TGAM and Reference GAM ----
get_tgam <- function(df) {
  ref_gam <- gam(pres ~ factor(year) +
                   s(julian) +
                   s(longitude, latitude),
                 data = samp,
                 family = binomial)
  all_tgams <- data.frame(aic_year = rep(0, length(years)),
                          year = rep(0, length(years)))
  all_tgams <- map(years, ~ get_tgam_yr_aic(samp, .x))
  best_tgam <- all_tgams[[which.min(sapply(1:length(all_tgams),function(x)min(all_tgams[[x]]$aic_year)))]]
  diff <- ref_gam$aic - best_tgam$aic_year
  # now get actual tgam
  df$thr <- ifelse(df$year <= best_tgam$yr, 'before', 'after')
  final_tgam <-  gam(pres ~ factor(year) +
                          s(longitude, latitude, by = factor(thr)) +
                          s(julian),
                        data = df,
                        family = binomial)
  return_list <- list(ref_gam, final_tgam, best_tgam$yr, diff)
}
