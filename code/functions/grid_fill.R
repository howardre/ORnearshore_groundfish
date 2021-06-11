# Manually populates grid cells with data
grid_fill <- function(lat_res, lon_res, year1, year2, data) {
  nlat = lat_res # determine resolution of grid
  nlon = lon_res
  latd = seq(42, 47, length.out = nlat)
  lond = seq(-125, -123.9, length.out = nlon)

  # make dataframe of just lon
  grid_lon = data.frame(
    lon1 = rep(lond[-length(lond)], (nlat - 1)),
    lon2 = rep(lond[-1], (nlat - 1)),
    lon3 = rep(lond[-1], (nlat - 1)),
    lon4 = rep(lond[-length(lond)], (nlat - 1)))

  # lat dataframe
  grid_lat = data.frame(
    lat1 = sort(rep(latd[-length(latd)], (nlon - 1))),
    lat2 = sort(rep(latd[-length(latd)], (nlon - 1))),
    lat3 = sort(rep(latd[-1], (nlon - 1))),
    lat4 = sort(rep(latd[-1], (nlon - 1))))

  # fill grid cells
  nstations = NA * (1:nrow(grid_lon))
  for (i in 1:length(nstations)) {
    tmp = in.chull(data$lon[data$year > year1 &
                                      data$year < year2],
                   data$lat[data$year > year1 &
                                      data$year < year2], grid_lon[i, ], grid_lat[i, ])
    nstations[i] = sum(data$hauls[data$year > year1 &
                                             data$year < year2] * tmp) # This decides what goes into each grid pixel
    points(data$lon[tmp],
           data$lat[tmp],
           col = i,
           pch = 16)
    polygon(grid_lon[i, ], grid_lat[i, ])
  }
  }
