## Fill proportion of biomass in grid cells

biomass_grid <- function(subset, year1, year2) {
  nstations = NA * (1:nrow(grid_lon))

  for (i in 1:length(nstations)) {
    tmp = in.chull(subset$lon[subset$year >= year1 &
                                subset$year <= year2],
                   subset$lat[subset$year >= year1 &
                                subset$year <= year2],
                   grid_lon[i,],
                   grid_lat[i,])
    nstations[i] = (sum(subset$CPUE[subset$year >= year1 &
                                      subset$year <= year2] * tmp) / sum(1 * tmp)) /
      (sum(subset$CPUE[subset$year >= year1 & subset$year <= year2]))
    points(subset$longitude[tmp],
           subset$latitude[tmp],
           col = i,
           pch = 16)
    polygon(grid_lon[i,], grid_lat[i,])
  }
  matrix(nstations,
         ncol = length(zlat),
         nrow = length(zlon),
         byrow = F)
}
