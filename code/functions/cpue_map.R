## Create maps of cpue for each species

cpue_map <- function(matrix, matrix2, color, title, bathy_dat, bathy_mat){
  image(zlon,
        zlat,
        matrix,
        col = hcl.colors(40, color, rev = F),
        xlim = c(-125, -123.6),
        ylim = c(42, 47),
        main = title,
        ylab = expression(paste("Latitude (" ^ 0, 'N)')),
        xlab = expression(paste("Longitude (" ^ 0, 'W)')),
        zlim = c(0, max(matrix2, na.rm = T)))
  map("worldHires",
      fill = T,
      col = "grey",
      add = T)
  image.plot(legend.only = T,
             col = hcl.colors(40, color, rev = F),
             legend.shrink = 0.2,
             smallplot = c(.76, .81, .09, .25),
             legend.cex = 0.7,
             axis.args = list(cex.axis = 0.9),
             legend.width = 1,
             zlim = c(0, max(matrix2, na.rm = T)),
             legend.lab = "Scaled Catch")
  contour(unique(bathy_dat$lon),
          sort(unique(bathy_dat$lat)),
          bathy_mat,
          lwd = 1,
          levels = -100,
          col = "gray40",
          add = T)
  contour(unique(bathy_dat$lon),
          sort(unique(bathy_dat$lat)),
          bathy_mat,
          lwd = 1,
          levels = -200,
          col = "gray40",
          add = T)
  points(-124.0535, 44.6368, pch = 20)
  text(-124.0535, 44.6368, "Newport", adj = c(0, 1.2))
  points(-123.8313, 46.1879, pch = 20)
  text(-123.88028, 46.13361, "Astoria", adj = c(0, 1.2))
  points(-124.3, 43.3, pch = 20)
  text(-124.3, 43.3, "Charleston", adj = c(0, 1.2))
}
