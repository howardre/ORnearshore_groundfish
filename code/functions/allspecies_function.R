# generalizable to all data frames
raw_plot <- function(lower_yr, upper_yr, df, df_name){
  plot(1,1,xlim = range(df$lon, na.rm=TRUE) + c(-.5,.2),
       ylim = range(df$lat, na.rm=TRUE) + c(-.2,.2),
       ylab = expression(paste("latitude ("^0,'N)')), # change to a copy and pasted degree sign and get rid of "expression()" portion
       xlab = expression(paste("longitude ("^0,'W)')),
       main = paste(df_name, " ", lower_yr, 's'))
  map("worldHires",fill=T,col="grey",add=T)
  points(df$lon[df$year>=lower_yr & df$year<=upper_yr & df$pres==1],
         df$lat[df$year>=lower_yr & df$year<=upper_yr & df$pres==1],
         pch=".",col='purple')
  contour(
    unique(bathy.dat$lon),
    sort(unique(bathy.dat$lat)),
    bathy.mat,
    levels = c(-10, -200),
    labcex = 0.7,
    add = T,
    col = 'black',
    labels = NULL,
    lwd = 2
  )
}

windows(width=28,height=18)
par(mfrow=c(1,4))

raw_plot(lower_yr = 1980, upper_yr =  1989, df = subset_starry, df_name = "Starry Flounder")
raw_plot(1990, 1999, suset_starry, df_name = "Starry Flounder")
raw_plot(2000, 2009, suset_starry, df_name = "Starry Flounder")
raw_plot(2010, 2017, suset_starry, df_name = "Starry Flounder")