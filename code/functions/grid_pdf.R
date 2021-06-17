grid_pdf <- function(lat_res, lon_res, year1, year2, zmax,
                     data, z_matrix, title, bathy_dat, bathy_mat) {
  nlat = lat_res # determine resolution of grid
  nlon = lon_res
  latd = seq(42, 47, length.out = nlat)
  lond = seq(-125, -123.9, length.out = nlon)
  zlat <- (latd[1:(length(latd) - 1)] + latd[2:length(latd)]) / 2
  zlon <- (lond[1:(length(lond) - 1)] + lond[2:length(lond)]) / 2
  data_subset <- length(unique(data$year[data$year > year1 &
                                           data$year < year2]))
  image(zlon,
        zlat,
        z_matrix,
        col = viridis(100, option = "A", direction = -1),
        xlim = c(-125,-123.6),
        ylim = c(42, 47),
        main = title,
        ylab = " ",
        xlab = " ",
        cex.axis = 1.7,
        cex.main = 2,
        zlim = c(0, zmax))
  rect(par("usr")[1], par("usr")[3],
       par("usr")[2], par("usr")[4],
       col = "aliceblue") # Color
  par(new = TRUE)
  image(zlon,
        zlat,
        z_matrix,
        col = viridis(100, option = "A", direction = -1),
        xlim = c(-125,-123.6),
        ylim = c(42, 47),
        main = title,
        ylab = " ",
        xlab = " ",
        cex.axis = 1.7,
        cex.main = 2,
        zlim = c(0, zmax))
  map("worldHires",
      fill = T,
      col = "grey",
      add = T)
  contour(unique(bathy_dat$lon),
          sort(unique(bathy_dat$lat)),
          bathy_mat,
          lwd = 1,
          levels = -100,
          col = "gray19",
          add = T)
  contour(unique(bathy_dat$lon),
          sort(unique(bathy_dat$lat)),
          bathy_mat,
          lwd = 1,
          levels = -200,
          col = "gray19",
          add = T)
  points(-124.0535, 44.6368,
         pch = 20,
         cex = 1.5)
  text(-124.0535, 44.6368, "Newport",
       adj = c(0, 1.7),
       cex = 1.6)
  points(-123.8313, 46.1879,
         pch = 20,
         cex = 1.6)
  text(-123.88028, 46.13361, "Astoria",
       adj = c(0.15, 1.2),
       cex = 1.6)
  points(-124.3, 43.3,
         pch = 20,
         cex = 1.6)
  text(-124.3, 43.3, "Charleston",
       adj = c(-0.1, 0.5),
       cex = 1.6)
}
