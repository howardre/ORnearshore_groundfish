## Create manuscript four panel figure

lic_sup <- function(matrix, matrix2, title,
                    bathy_dat, bathy_mat) {
  image(zlon,
        zlat,
        matrix,
        col = hcl.colors(100, "YlOrRd", rev = T),
        xlim = c(-125, -123.6),
        ylim = c(42, 47),
        main = title,
        ylab = " ",
        xlab = " ",
        cex.axis = 2,
        cex.main = 2,
        zlim = c(0, max(matrix2, na.rm = T)))
  rect(par("usr")[1], par("usr")[3],
       par("usr")[2], par("usr")[4],
       col = "aliceblue") # Color
  par(new = TRUE)
  image(zlon,
        zlat,
        matrix,
        col = hcl.colors(100, "YlOrRd", rev = T),
        xlim = c(-125, -123.6),
        ylim = c(42, 47),
        main = title,
        ylab = " ",
        xlab = " ",
        cex.axis = 2,
        cex.main = 2,
        zlim = c(0, max(matrix2, na.rm = T)))
  map("worldHires",
      fill = T,
      col = "grey",
      add = T)
  contour(unique(bathy_dat$lon),
          sort(unique(bathy_dat$lat)),
          bathy_mat,
          lwd = 1,
          levels = -55,
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
         cex = 1.6)
  text(-124.0535, 44.6368, "Newport",
       adj = c(0, 1.7),
       cex = 2)
  points(-123.8313, 46.1879,
         pch = 20,
         cex = 1.6)
  text(-123.88028, 46.13361, "Astoria",
       adj = c(0.15, 1.2),
       cex = 2)
  points(-124.3, 43.3,
         pch = 20,
         cex =  1.6)
  text(-124.3, 43.3, "Charleston",
       adj = c(-0.1, 0.5),
       cex = 2)
  mtext(expression(paste("Latitude ("^0,'N)')),
        side = 2,
        line = -2.2,
        cex = 1.5,
        outer = TRUE)
  mtext(expression(paste("Longitude ("^0, 'W)')),
        side = 1,
        line = -1.5,
        cex = 1.5,
        outer = TRUE)
}
