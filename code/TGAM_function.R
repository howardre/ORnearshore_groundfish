# SHOUld work for fishes with year in the gam. Couldn't figure out indexing for line
# 23 though
# Go through each year to find the threshold year for sample ----
get_tgam_yr_aic <- function(df, yr) {
  df$thr <- ifelse(df$year <= yr, 'before', 'after')
  aic_year <-  gam(
    pres ~ factor(year) + s(longitude, latitude, by = factor(thr)) + s(julian),
    data = df,
    family = binomial)$aic
  return_list <- list(aic_yr, yr)
}

# Returns AIC difference, best year, best TGAM and Reference GAM ----
get_gams <- function(df) {
  ref_gam <- gam(
    pres ~ factor(year) + s(julian) + s(longitude, latitude),
    data = samp,
    family = binomial)
  all_tgams <- data.frame(
    aic_yr = rep(0, length(years)),
    year= rep(0, length(years)))
  all_tgams <- map(years, ~get_tgam_aic(samp, .x))
  best_tgam <- all_tgams$year[which(sort(all_tgams$aic_yr)[1]),] # can NOT figure out this indexing. want the year that matches with the lowest aic
  diff <- ref_gam$aic - best_tgam$aic_yr
  # now get actual tgam
  df$thr <- ifelse(df$year <= best_tgam$year, 'before', 'after')
  best_tgam_gam <-  gam(
    pres ~ factor(year) + s(longitude, latitude, by = factor(thr)) + s(julian),
    data = df,
    family = binomial)
  return_list <- list(ref_gam, best_tgam_gam, best_tgam$year, diff)
}
