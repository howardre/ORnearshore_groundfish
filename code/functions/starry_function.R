# for starry flounder only
grid_plot <- function(lower_yr, upper_yr){
  plot(1,1,xlim=range(subset_starry$lon,na.rm=TRUE)+c(-.5,.2),
       ylim=range(subset_starry$lat,na.rm=TRUE)+c(-.2,.2),
       ylab=expression(paste("latitude ("^0,'N)')),
       xlab=expression(paste("longitude ("^0,'W)')),
       main=paste('Starry Flounder ', lower_yr, 's'))
  map("worldHires",fill=T,col="grey",add=T)
  points(subset_starry$lon[subset_starry$year>=lower_yr & subset_starry$year<=upper_yr & subset_starry$pres==1],
         subset_starry$lat[subset_starry$year>=lower_yr & subset_starry$year<=upper_yr & subset_starry$pres==1],
         pch=".",col='purple')
  contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=c(-10,-200),
          labcex=0.7,add=T,col='black', labels = NULL, lwd = 2)
}



windows(width=28,height=18)
par(mfrow=c(1,4))

grid_plot(1980, 1989)
grid_plot(1990, 1999)
