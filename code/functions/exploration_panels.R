## Logbook exploration 4-panel map
exploration_panels <- function(subset, bathy_dat, bathy_mat, year1, year2, title){
plot(1,
     1,
     xlim = range(subset$lon, na.rm = TRUE) + c(-.5, .2),
     ylim = range(subset$lat, na.rm = TRUE) + c(-0.2, .2),
     ylab = expression(paste("latitude (" ^ 0, 'N)')),
     xlab = expression(paste("longitude (" ^ 0, 'W)')),
     main = title)
  maps::map("worldHires",
            fill = T,
            col = "grey",
            add = T)
points(subset$lon[subset$year >= year1 &
                            subset$year <= year2 & subset$pres == 1],
       subset$lat[subset$year >= year1 &
                            subset$year <= year2 & subset$pres == 1],
       pch = ".",
       col = 'purple')
contour(unique(bathy_dat$lon),
        sort(unique(bathy_dat$lat)),
        bathy_mat,
        levels = c(-10,-200),
        labcex = 0.7,
        add = T,
        col = 'black',
        labels = NULL,
        lwd = 2)
}
