# Go through each year to find the threshold year for sample ----
# This is used in the get_tgam function
get_yr_aic <- function(df, yr) {
  df$thr <- ifelse(df$year <= yr, 'before', 'after')
  ctrl <- list(nthreads = 6)
  aic_year <-  gam(pres ~ factor(year) +
                     s(longitude, latitude, by = factor(thr)) +
                     s(julian),
                   control = ctrl,
                   data = df,
                   family = binomial)$aic # gives the AIC values
  return_list <- data.frame(aic_year, yr)
}

# Returns AIC difference, best year, best TGAM and reference GAM ----
get_tgam <- function(df, years) {
  ctrl <- list(nthreads = 6)
  ref_gam <- gam(pres ~ factor(year) +
                   s(julian) +
                   s(longitude, latitude),
                 control = ctrl,
                 data = df,
                 family = binomial) # this is the reference gam with no threshold
  all_tgams <- data.frame(aic_year = rep(0, length(years)),
                          year = rep(0, length(years))) # contains all of the gams for each year
  all_tgams <- future_map(years, ~ get_yr_aic(df, .x))
  best_tgam <- all_tgams[[which.min(sapply(1:length(all_tgams),
                                           function(x) min(all_tgams[[x]]$aic_year)))]]
  diff <- ref_gam$aic - best_tgam$aic_year
  df$thr <- ifelse(df$year <= best_tgam$yr, 'before', 'after')
  final_tgam <-  gam(pres ~ factor(year) +
                          s(longitude, latitude, by = factor(thr)) +
                          s(julian),
                     control = ctrl,
                     data = df,
                     family = binomial)
  all_aic <- sapply(all_tgams, function(x) {as.numeric(x[1])})
  return_list <- list(ref_gam, final_tgam, best_tgam$yr, diff, all_aic)
}

# Get threshold year using formula without year ----
get_aic_woyear <- function(df, yr) {
  ctrl <- list(nthreads = 6)
  df$thr <- ifelse(df$year <= yr, 'before', 'after')
  aic_year <-  gam(pres ~ factor(thr) +
                     s(longitude, latitude, by = factor(thr)) +
                     s(julian),
                   control = ctrl,
                   data = df,
                   family = binomial)$aic
  return_list <- data.frame(aic_year, yr)
}

# Returns AIC difference, best year, best TGAM and Reference GAM ----
get_tgam_woyear <- function(df, years) {
  ctrl <- list(nthreads = 6)
  ref_gam <- gam(pres ~ s(julian) +
                   s(longitude, latitude),
                 control = ctrl,
                 data = df,
                 family = binomial)
  all_tgams <- data.frame(aic_year = rep(0, length(years)),
                          year = rep(0, length(years)))
  all_tgams <- future_map(years, ~ get_aic_woyear(df, .x))
  best_tgam <- all_tgams[[which.min(sapply(1:length(all_tgams),
                                           function(x) min(all_tgams[[x]]$aic_year)))]]
  diff <- ref_gam$aic - best_tgam$aic_year
  df$thr <- ifelse(df$year <= best_tgam$yr, 'before', 'after')
  final_tgam <-  gam(pres ~ factor(thr) +
                       s(longitude, latitude, by = factor(thr)) +
                       s(julian),
                     control = ctrl,
                     data = df,
                     family = binomial)
  all_aic <- sapply(all_tgams, function(x) {as.numeric(x[1])})
  return_list <- list(ref_gam, final_tgam, best_tgam$yr, diff, all_aic)
}

# Use to plot the final results ----
plot_AIC <- function(tgam, years) {
  plot(years,
       tgam[[5]], # need y axis to be the AIC range for all GAMs
       type = 'b',
       xlab = 'Year',
       ylab = 'AIC',
       main = deparse(substitute(tgam)),
       cex.main = 1.4,
       cex.lab = 1.2,
       cex.axis = 1.2)
  abline(v = tgam[[3]], lty = 2)
  abline(h = AIC(tgam[[1]]), lty = 2)
}
