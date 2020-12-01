setwd('/Users/howar/Documents/Oregon State/ORnearshore_groundfish/code')

source("functions/vis_gam_COLORS.R")
jet.colors <- colorRampPalette(c("#F7FCFD", "#E5F5F9", "#CCECE6", "#99D8C9",
                                 "#66C2A4", "#41AE76", "#238B45", "#006D2C", "#00441B"))

#Libraries
library(maps)
library(mapdata)
library(fields)
library(marmap)
library(colorRamps)
library(itsadug)
library(RColorBrewer)

#Load data
trawl_data <- read.delim("../data/NMFS_data/trawl_data.txt", header = T)
load("../data/bathy.dat")
load("../data/bathy.mat")

blues <- c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
greys <- c(grey(0.6), grey(0.93), grey(0.99))

OR_bathy <- getNOAA.bathy(lon1= -127, lon2= -121, lat1= 49, lat2= 39, resolution=1)

# Function for map
# Survey plot
survey_map <- function(title, bathymetry, survey, years, letter){
plot.bathy(
  bathymetry,
  image = T,
  axes = T,
  lwd = 0.03,
  land = T,
  n = 0,
  bpal = list(c(0, max(bathymetry), greys), c(min(bathymetry), 0, blues)),
  ylim = c(40.3, 48.2),
  xlim = c(-125.5,-123.5),
  xlab = ,
  ylab = ,
  main = title,
  cex.lab = 1.5,
  cex.main = 2,
  cex.axis = 1.1
)
plot(
  bathymetry,
  deep = 0,
  shallow = 0,
  lwd = 1,
  add = T
)
plot(
  bathymetry,
  deep = -50,
  shallow = -50,
  lwd = 0.4,
  drawlabels = T,
  add = T,
  col = "slategrey"
)
plot(
  bathymetry,
  deep = -200,
  shallow = -200,
  lwd = 0.4,
  drawlabels = T,
  add = T,
  col = "slategrey"
)
map.scale(-126., 40.3, cex = 1)
points(
  survey$longitude[survey$year == years],
  survey$latitude[survey$year == years],
  pch = 18,
  col = 'black',
  cex = 1.3
)
title(
  outer = outer,
  adj = .025,
  main = letter,
  cex.main = 2,
  col = "black",
  font = 2,
  line = -2
)
}
tiff("../final_figs/Figure_1.tiff", width = 11, height = 12, units = "in", res = 300)
par(mfrow = c(1, 2),
    family = 'serif',
    mar = c(4, 5, 3, .2) + .15)
survey_map("AFSC Triennial Survey", OR_bathy, trawl_data, 1995, "(A)")
survey_map("NWFSC Annual Survey", OR_bathy, trawl_data, 2018, "(B)")
dev.off()

fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {

  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright",
                          "left", "center", "right",
                          "bottomleft", "bottom", "bottomright"))

  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")

    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    }
  }

  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }

  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100

  x1 <- switch(pos,
               topleft     =x[1] + sw,
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)

  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)

  old.par <- par(xpd=NA)
  on.exit(par(old.par))

  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}

###########################################################################################################################################################
# Make temperature GAM ----
trawl_data$bottom_temp[is.na(trawl_data$bottom_temp)] <- 0
trawl_data <- trawl_data[trawl_data$bottom_temp < 10, ]
trawl_data <- trawl_data[trawl_data$bottom_temp > 4, ]

temp_gam <- gam(bottom_temp ~ s(year) +
                  s(longitude, latitude) +
                  s(julian), data = trawl_data)
summary(temp_gam)
windows()
par(mfrow=c(2,2))
gam.check(temp_gam)

pdf("temp_plot.pdf", width = 5, height = 13)
par(family = 'serif', mar=c(4,5,3,.2)+.15)
myvis.gam(temp_gam,view=c('longitude','latitude'),too.far=0.03,plot.type='contour', contour.col = "gray48",
        color="jet" ,type='response', main = paste("Nearshore Temperature"),
        xlim=range(trawl_data$longitude,na.rm=TRUE)+c(-.1,0.8),
        ylim=range(trawl_data$latitude,na.rm=TRUE)+c(-.4,.5),
        xlab=expression(paste("Longitude ("^o,'W)')),
        ylab=expression(paste("Latitude ("^o,'N)')), family = "serif",cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.1)
contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-seq(200,200,by=200),
        labcex=0.7,add=T,col='black', labels = NULL, lwd = 1.7)
map('worldHires',add=T,col='antiquewhite4',fill=T)
points(cities, col="black", pch=16)
text(cities, labels=row.names(cities), cex=1.2, font=2, offset=0.5, adj=c(-0.2,2), family = "serif")
image.plot(legend.only=T, zlim=c(6.5,8.5), col=jet.colors(100), legend.shrink = 0.2, smallplot = c(.26,.3,.08,.23),
           legend.cex = 1., legend.lab = "Temperature (C)", axis.args=list(cex.axis = .8), legend.width = 1)
dev.off()

# Temperature TGAM
years<-sort(unique(trawl_data$year))[4:22]

aic.year<-years*NA
trawl_data$thr<-NA
for(i in 1:length(years)){
  trawl_data$thr<-ifelse(trawl_data$year<=years[i],'before','after')
  gam.obj <- gam(bottom_temp~factor(year)+s(julian)+s(longitude,latitude,by=factor(thr)),
                 data=trawl_data)
  aic.year[i]<-gam.obj$aic}
best.year<-years[order(aic.year)[1]]
best.aic<-sort(aic.year)[1]

#Null model
gam.ref<-gam(bottom_temp~factor(year)+s(julian)+s(longitude,latitude),data=trawl_data)
AIC(gam.ref)

windows()
par(family = "serif")
plot(years,aic.year,type='b',xlab='Year',ylab='AIC',ylim=range(c(aic.year,AIC(gam.ref))), main = "Bottom Temperature")
abline(v=best.year,lty=2)
text(1990,3000,'Before')
text(1994,3000,'After')
abline(h=AIC(gam.ref),lty=2)

trawl_data$thr<-ifelse(trawl_data$year<=best.year,'before','after')
gam.real<- gam(bottom_temp~factor(year)+s(julian)+s(longitude,latitude,by=factor(thr)),data=trawl_data)

#TGAM Maps
pdf("temp_TGAM.pdf", width = 11, height = 12)
par(mfrow=c(1,2), family = 'serif', mar=c(4,5,3,.2)+.1)
myvis.gam(gam.real,view=c('longitude','latitude'), too.far = 0.025, plot.type='contour',color='jet',
          type='response', cond=list(thr='before'),main = 'Pre-1992', contour.col = "gray35",
          xlim = range(trawl_data$longitude,na.rm=TRUE)+c(-0.5,1), ylim=range(trawl_data$latitude,na.rm=TRUE)+c(-0.5,0.5),
          cex.lab = 1.5,cex.axis = 1.5,cex.main = 2,xlab=expression(paste("Longitude ("^0,'W)')), ylab=expression(paste("Latitude ("^0,'N)')))
contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-seq(200,200,by=200),
        labcex=0.7,add=T,col='black', labels = NULL, lwd = 1.7)
map('worldHires',add=T,col='peachpuff3',fill=T)
points(cities, col="black", pch=16)
text(cities, labels=row.names(cities), cex=1.2, font =2, offset=0.2, adj=c(-0.3,2))
image.plot(legend.only=T, zlim=c(6.4,8.7), col=jet.colors(100), legend.shrink = 0.2, smallplot = c(.26,.3,.08,.23),
           legend.cex = 1., legend.lab = "Temperature (C)", axis.args=list(cex.axis = .8), legend.width = 1)
myvis.gam(gam.real,view=c('longitude','latitude'),too.far=0.025,plot.type='contour',color='jet',
          type='response', cond=list(thr='after'),main = 'Post-1992', contour.col = "gray35",
          xlim = range(trawl_data$longitude,na.rm=TRUE)+c(-0.5,1), ylim=range(trawl_data$latitude,na.rm=TRUE)+c(-0.5,0.5),
          cex.lab = 1.5, cex.axis = 1.5, cex.main = 2,xlab=expression(paste("Longitude ("^0,'W)')), ylab=expression(paste("Latitude ("^0,'N)')))
contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-seq(200,200,by=200),
        labcex=0.7,add=T,col='black', labels = NULL, lwd = 1.7)
map('worldHires',add=T,col='peachpuff3',fill=T)
points(cities, col="black", pch=16)
text(cities, labels=row.names(cities), cex=1.2, font =2, offset=0.2, adj=c(-0.3,2))
image.plot(legend.only=T, zlim=c(6.4,8.7), col=jet.colors(100), legend.shrink = 0.2, smallplot = c(.26,.3,.08,.23),
           legend.cex = 1., legend.lab = "Temperature (C)", axis.args=list(cex.axis = .8), legend.width = 1)
dev.off()

# make map
windows(width=15, height=15)
par(mfrow=c(1,2))
plot(1,1,xlim=range(trawl_data$longitude,na.rm=TRUE),
     ylim=range(trawl_data$latitude,na.rm=TRUE),
     ylab=expression(paste("Latitude ("^0,'N)')),
     xlab=expression(paste("Longitude ("^0,'E)')),
     main=paste('2004'))
points(trawl_data$longitude[trawl_data$year==2004],
       trawl_data$latitude[trawl_data$year==2004],
       pch=18,col='black', cex = 1)
plot(1,1,xlim=range(trawl_data$longitude,na.rm=TRUE),
     ylim=range(trawl_data$latitude,na.rm=TRUE),
     ylab=expression(paste("Latitude ("^0,'N)')),
     xlab=expression(paste("Longitude ("^0,'E)')),
     main=paste('2005'))
points(trawl_data$longitude[trawl_data$year==2005],
       trawl_data$latitude[trawl_data$year==2005],
       pch=18,col='black', cex = 1)

## Presentation Location map
browns <- c("navajowhite4", "navajowhite3", "navajowhite2", "navajowhite1", "navajowhite", "moccasin")
blues <- c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")

westcoast_bathy <- getNOAA.bathy(lon1= -129, lon2= -117, lat1= 49, lat2= 30, resolution=1)

windows(width = 8, height = 12.8)
par(family = "serif")
plot.bathy(westcoast_bathy, image = T, axes = T, lwd=0.03, land = T, n = 0,
           bpal = list(c(0, max(westcoast_bathy), browns), c(min(westcoast_bathy),0,blues)),
           ylim=c(30,48.5), xlim = c(-128, -117),
           xlab=expression(paste("Longitude ("^o,'W)')),
           ylab=expression(paste("Latitude ("^o,'N)')),
           main = "Northern California Current", cex.lab = 1, cex.main = 1, cex.axis = 1.1)
points(-118.2437,34.0522, pch=20)
text(-120.2,35.3,"Los Angeles", adj=c(0,1.2))
points(-122.6750,45.5051, pch=20, cex = 0.8)
text(-122.6750,45.5051,"Portland", adj=c(0,1.2))
points(-122.3321,47.6062, pch=20)
text(-122.32,47.6062,"Seattle", adj=c(0,1.2))
points(-122.4194,37.7749, pch=20)
text(-121.4,37.7749,"San Francisco", adj=c(0,1.2))
points(-121.4944,38.5816, pch=20)
text(-121.4944,39.4,"Sacramento", adj=c(0,1.2))
points(-123.2620,44.5646, pch=20)
text(-123.2620,44.5646,"Corvallis", adj=c(0,1.2))
