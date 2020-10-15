# Single Species
species_plot <- function(lower_yr, upper_yr) {
  plot(1, 1, xlim = range(species_subset$lon, na.rm = TRUE) + c(-.5, .2),
       ylim = range(species_subset$lat, na.rm = TRUE) + c(-.2, .2),
       ylab = "latitude °N",
       xlab = "longitude °W",
       main = paste(lower_yr, 's', sep = ""))
  map("worldHires",
      fill = T,
      col = "grey",
      add = T)
  points(species_subset$lon[species_subset$year >= lower_yr &
                              species_subset$year <= upper_yr &
                              species_subset$pres == 1],
         species_subset$lat[species_subset$year >= lower_yr &
                              species_subset$year <= upper_yr &
                              species_subset$pres == 1],
         pch = ".",
         col = 'purple')
  contour(unique(bathy.dat$lon),
          sort(unique(bathy.dat$lat)),
          bathy.mat,
          levels = c(-10, -200),
          labcex = 0.7,
          add = T,
          col = 'black',
          labels = NULL,
          lwd = 2)
}

# Plot the maps to check
#windows(width = 28, height = 18)
#par(mfrow = c(1, 4))
#species_plot(1980, 1989)
#species_plot(1990, 1999)
#species_plot(2000, 2009)
#species_plot(2010, 2017)

#rm(list = ls())
