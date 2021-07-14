## Fill grid cells with corresponding points

biomass_fillpts <- function(subset, year1, year2){
plot(subset$lon[subset$year >= year1 &
                  subset$year <= year2],
     subset$lat[subset$year >= year1 &
                  subset$year <= year2],
     pch = '.',
     ylim = c(42, 47),
     xlim = c(-125,-123.9))
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
  # This gives proportion of biomass per grid cell (divide avg. weight for cell by biomass for that decade)
  }

  points(subset$longitude[tmp],
         subset$latitude[tmp],
         col = i,
         pch = 16)
  polygon(grid_lon[i,], grid_lat[i,])
  }
