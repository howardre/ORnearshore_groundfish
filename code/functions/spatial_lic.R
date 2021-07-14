# Fill empty cells with 0s or NAs

spatial_lic <- function(logbooks_cells, logbook_data,
                        survey_cells, survey_data,
                        year1, year2){
  nstations1 = NA * (1:nrow(grid_lon))

  for (i in 1:length(nstations1)) {
    tmp = in.chull(logbook_data$lon[logbook_data$year >= year1 &
                                logbook_data$year <= year2],
                   logbook_data$lat[logbook_data$year >= year1 &
                                logbook_data$year <= year2],
                   grid_lon[i,],
                   grid_lat[i,])
    nstations1[i] = (sum(logbook_data$CPUE[logbook_data$year >= year1 &
                                      logbook_data$year <= year2] * tmp) / sum(1 * tmp)) /
      (sum(logbook_data$CPUE[logbook_data$year >= year1 & logbook_data$year <= year2]))
    points(logbook_data$longitude[tmp],
           logbook_data$latitude[tmp],
           col = i,
           pch = 16)
    polygon(grid_lon[i,], grid_lat[i,])
  }

  nstations2 = NA * (1:nrow(grid_lon))

  for (i in 1:length(nstations2)) {
    tmp = in.chull(survey_data$lon[survey_data$year >= year1 &
                                survey_data$year <= year2],
                   survey_data$lat[survey_data$year >= year1 &
                                survey_data$year <= year2],
                   grid_lon[i,],
                   grid_lat[i,])
    nstations2[i] = (sum(survey_data$CPUE[survey_data$year >= year1 &
                                      survey_data$year <= year2] * tmp) / sum(1 * tmp)) /
      (sum(survey_data$CPUE[survey_data$year >= year1 & survey_data$year <= year2]))
    points(survey_data$longitude[tmp],
           survey_data$latitude[tmp],
           col = i,
           pch = 16)
    polygon(grid_lon[i,], grid_lat[i,])
  }

  decade_grid <- as.matrix(logbooks_cells * survey_cells, na.rm = T) /
    (sqrt(sum(logbooks_cells ^ 2, na.rm = T) * sum(survey_cells ^ 2, na.rm = T)))

  logbook_cells <-  matrix(nstations1,
                           ncol = length(zlat),
                           nrow = length(zlon),
                           byrow = F)
  survey_cells <- matrix(nstations2,
                         ncol = length(zlat),
                         nrow = length(zlon),
                         byrow = F)
  decade_grid[is.nan(survey_cells)] <- 0
  decade_grid[is.nan(logbook_cells)] <- NA
  return(decade_grid)
  }
