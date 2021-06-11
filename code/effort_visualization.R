##Part 1: make a regular grid, identifies points that fall within each grid cell, 
#and plot results on a map

  setwd("C:/Users/howar/Documents/Oregon State/Thesis/Logbook Data/Logbook")
  load('logbooks_corrected')
  library(sgeostat)
  library(maps)
  library(mapdata)
  library(fields)
  library(marmap)
  library(dplyr)
  
  combined$month_day <- as.numeric(format(combined$TOWDATE, '%m%d'))
  combined <- combined[combined$month_day>=517&combined$month_day<= 929,]
  combined <- combined[combined$depth<=-5,]
  # generate a column to count each trawl
  combined$hauls <- combined$GEAR[combined$GEAR>0]<-1
  combined <- combined[combined$species=="PTRL_ADJ",]

#DEPTH: import data and show contour on a map
#http://maps.ngdc.noaa.gov/viewers/wcs-client/
#Choose ETOPO1 (ice)
  bathy.dat<-read.table('etopo1.xyz',sep='')
  names(bathy.dat)<-c('lon','lat','depth')
  bathy.dat$depth[bathy.dat$depth>0]<-NA#Avoid points above water
  head(bathy.dat)
  bathy.mat<-matrix(bathy.dat$depth,nrow=length(unique(bathy.dat$lon)),ncol=length(unique(bathy.dat$lat)))[,order(unique(bathy.dat$lat))]
  
#Part 1: make a regular grid and count stations within each grid cell. Plot results
  
    nlat=20#determine resolution of grid
    nlon=15
    latd=seq(42,47,length.out=nlat)
    lond=seq(-125,-123.9,length.out=nlon)
    
    
  grid.lon=data.frame(
    lon1=rep(lond[-length(lond)],(nlat-1)),
    lon2=rep(lond[-1],(nlat-1)),
    lon3=rep(lond[-1],(nlat-1)),
    lon4=rep(lond[-length(lond)],(nlat-1)))#make dataframe of just lon
  
  grid.lat=data.frame(
    lat1=sort(rep(latd[-length(latd)],(nlon-1))),
    lat2=sort(rep(latd[-length(latd)],(nlon-1))),
    lat3=sort(rep(latd[-1],(nlon-1))),
    lat4=sort(rep(latd[-1],(nlon-1))))#lat dataframe
  
  dev.new(width=4,height=10)
  plot(combined$lon,combined$lat,pch='.',ylim=c(42,47),xlim=c(-125,-123.9))
  n.stations=NA*(1:nrow(grid.lon))
  

  #for(i in 1:length(n.stations)){
    #tmp=in.chull(combined$lon,combined$lat,grid.lon[i,],grid.lat[i,])
    #n.stations[i]=sum(combined$PTRL_ADJ*tmp)#This decides what goes into each grid pixel
    #points(combined$lon[tmp],combined$lat[tmp],col=i,pch=16)
    #polygon(grid.lon[i,],grid.lat[i,])
  #}
  #map("worldHires",fill=T,col="grey",add=T)
  
  
  ##ALL YEARS
  years.data<-length(unique(combined$year))
  z.lat<-(latd[1:(length(latd)-1)]+latd[2:length(latd)])/2
  z.lon<-(lond[1:(length(lond)-1)]+lond[2:length(lond)])/2
  z.matrix<-matrix(n.stations,ncol=length(z.lat),nrow=length(z.lon),byrow=F)/years.data #n of tows/year
  
##80s
  n.stations1=NA*(1:nrow(grid.lon))
  
  for(i in 1:length(n.stations1)){
    tmp=in.chull(combined$lon[combined$year>1980&combined$year<1990],combined$lat[combined$year>1980&combined$year<1990],grid.lon[i,],grid.lat[i,])
    n.stations1[i]=sum(combined$hauls[combined$year>1980&combined$year<1990]*tmp)#This decides what goes into each grid pixel
    points(combined$lon[tmp],combined$lat[tmp],col=i,pch=16)
    polygon(grid.lon[i,],grid.lat[i,])
  }


##90s
   n.stations2=NA*(1:nrow(grid.lon))
   
   for(i in 1:length(n.stations2)){
     tmp=in.chull(combined$lon[combined$year>1989&combined$year<2000],combined$lat[combined$year>1989&combined$year<2000],grid.lon[i,],grid.lat[i,])
     n.stations2[i]=sum(combined$hauls[combined$year>1989&combined$year<2000]*tmp)#This decides what goes into each grid pixel
     points(combined$lon[tmp],combined$lat[tmp],col=i,pch=16)
     polygon(grid.lon[i,],grid.lat[i,])
   }

   
##00s
   n.stations3=NA*(1:nrow(grid.lon))
   
   for(i in 1:length(n.stations3)){
     tmp=in.chull(combined$lon[combined$year>=2000&combined$year<=2009],combined$lat[combined$year>=2000&combined$year<=2009],grid.lon[i,],grid.lat[i,])
     n.stations3[i]=sum(combined$hauls[combined$year>=2000&combined$year<=2009]*tmp)#This decides what goes into each grid pixel
     points(combined$lon[tmp],combined$lat[tmp],col=i,pch=16)
     polygon(grid.lon[i,],grid.lat[i,])
   }
   
##10s
   n.stations4=NA*(1:nrow(grid.lon))
   
   for(i in 1:length(n.stations4)){
     tmp=in.chull(combined$lon[combined$year>=2010&combined$year<=2017],combined$lat[combined$year>=2010&combined$year<=2017],grid.lon[i,],grid.lat[i,])
     n.stations4[i]=sum(combined$hauls[combined$year>=2010&combined$year<=2017]*tmp)#This decides what goes into each grid pixel
     points(combined$lon[tmp],combined$lat[tmp],col=i,pch=16)
     polygon(grid.lon[i,],grid.lat[i,])
   }
   
   ## decades
   z.lat<-(latd[1:(length(latd)-1)]+latd[2:length(latd)])/2
   z.lon<-(lond[1:(length(lond)-1)]+lond[2:length(lond)])/2
   eighties.data<-length(unique(combined$year[combined$year>1980&combined$year<1990]))
   z.matrix1<-matrix(n.stations1,ncol=length(z.lat),nrow=length(z.lon),byrow=F)/eighties.data
   z.matrix1 <- ifelse(z.matrix1==0, NA, z.matrix1)
   nineties.data<-length(unique(combined$year[combined$year>1989&combined$year<2000]))
   z.matrix2<-matrix(n.stations2,ncol=length(z.lat),nrow=length(z.lon),byrow=F)/nineties.data
   z.matrix2 <- ifelse(z.matrix2==0, NA, z.matrix2)
   thousands.data<-length(unique(combined$year[combined$year>1999&combined$year<2010]))
   z.matrix3<-matrix(n.stations3,ncol=length(z.lat),nrow=length(z.lon),byrow=F)/thousands.data
   z.matrix3 <- ifelse(z.matrix3==0, NA, z.matrix3)
   tens.data<-length(unique(combined$year[combined$year>2009&combined$year<2018]))
   z.matrix4<-matrix(n.stations4,ncol=length(z.lat),nrow=length(z.lon),byrow=F)/tens.data
   z.matrix4 <- ifelse(z.matrix4==0, NA, z.matrix4)
   
###MAKE MAPS

 #windows(width=11,height=5)
  #par(mfrow=c(1,2))
   windows(width=15,height=9)
   par(mfrow=c(1,4),family = 'serif', mar=c(4,5,3,.3)+.1)
 
   image(z.lon,z.lat,z.matrix1,col=hcl.colors(40, "RdYlBu", rev=T),xlim=c(-125,-123.6),ylim=c(42,47), main="Fishing Effort 1980s",
         ylab=expression(paste("Latitude ("^0,'N)')), xlab=expression(paste("Longitude ("^0,'W)')), zlim = c(0,230))
   map("worldHires",fill=T,col="grey",add=T)
   image.plot(legend.only=T, col=hcl.colors(40, "RdYlBu", rev=T), legend.shrink = 0.2, smallplot = c(.76,.81,.09,.25), 
              legend.cex = 0.7, axis.args=list(cex.axis = 0.9), legend.width = 1, zlim=c(0,230), legend.lab = "avg. number of tows")
   contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=1, 
           levels=-100,col= "gray19",add=T)
   contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=1, 
           levels=-200,col= "gray19",add=T)
   points(-124.0535,44.6368, pch=20)
   text(-124.0535,44.6368,"Newport", adj=c(0,1.2))
   points(-123.8313,46.1879, pch=20)
   text(-123.88028,46.13361,"Astoria", adj=c(0,1.2))
   points(-124.3,43.3, pch=20)
   text(-124.3,43.3,"Charleston", adj=c(0,1.2))

    image(z.lon,z.lat,z.matrix2,col=hcl.colors(40, "RdYlBu", rev=T),xlim=c(-125,-123.6),ylim=c(42,47), main="Fishing Effort 1990s",
          ylab=expression(paste("Latitude ("^0,'N)')), xlab=expression(paste("Longitude ("^0,'W)')), zlim=c(0,230))
    map("worldHires",fill=T,col="grey",add=T)
    image.plot(legend.only=T, col=hcl.colors(40, "RdYlBu", rev=T), legend.shrink = 0.2, smallplot = c(.76,.81,.09,.25), 
               legend.cex = 0.7, axis.args=list(cex.axis = 0.9), legend.width = 1, zlim=c(0,230), legend.lab = "avg. number of tows")
    contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=1, 
            levels=-100,col= "gray19",add=T)
    contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=1, 
            levels=-200,col= "gray19",add=T)
    points(-124.0535,44.6368, pch=20)
    text(-124.0535,44.6368,"Newport", adj=c(0,1.2))
    points(-123.8313,46.1879, pch=20)
    text(-123.88028,46.13361,"Astoria", adj=c(0,1.2))
    points(-124.3,43.3, pch=20)
    text(-124.3,43.3,"Charleston", adj=c(0,1.2))
  
    image(z.lon,z.lat,z.matrix3,col=hcl.colors(40, "RdYlBu", rev=T),xlim=c(-125,-123.6),ylim=c(42,47), main="Fishing Effort 2000s",
          ylab=expression(paste("Latitude ("^0,'N)')), xlab=expression(paste("Longitude ("^0,'W)')), zlim=c(0,230))
    map("worldHires",fill=T,col="grey",add=T)
    image.plot(legend.only=T, col=hcl.colors(40, "RdYlBu", rev=T), legend.shrink = 0.2, smallplot = c(.76,.81,.09,.25), 
               legend.cex = 0.7, axis.args=list(cex.axis = 0.9), legend.width = 1, zlim=c(0,230), legend.lab = "avg. number of tows")
    contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=1, 
            levels=-100,col= "gray19",add=T)
    contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=1, 
            levels=-200,col= "gray19",add=T)
    points(-124.0535,44.6368, pch=20)
    text(-124.0535,44.6368,"Newport", adj=c(0,1.2))
    points(-123.8313,46.1879, pch=20)
    text(-123.88028,46.13361,"Astoria", adj=c(0,1.2))
    points(-124.3,43.3, pch=20)
    text(-124.3,43.3,"Charleston", adj=c(0,1.2))
  
    image(z.lon,z.lat,z.matrix4,col=hcl.colors(40, "RdYlBu", rev=T),xlim=c(-125,-123.6),ylim=c(42,47), main="Fishing Effort 2010s",
          ylab=expression(paste("Latitude ("^0,'N)')), xlab=expression(paste("Longitude ("^0,'W)')), zlim=c(0,230))
    map("worldHires",fill=T,col="grey",add=T)
    image.plot(legend.only=T, col=hcl.colors(40, "RdYlBu", rev=T), legend.shrink = 0.2, smallplot = c(.76,.81,.09,.25), 
               legend.cex = 0.7, axis.args=list(cex.axis = 0.9), legend.width = 1, zlim=c(0,230), legend.lab = "avg. number of tows")
    contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=1, 
            levels=-100,col= "gray19",add=T)
    contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=1, 
            levels=-200,col= "gray19",add=T)
    points(-124.0535,44.6368, pch=20)
    text(-124.0535,44.6368,"Newport", adj=c(0,1.2))
    points(-123.8313,46.1879, pch=20)
    text(-123.88028,46.13361,"Astoria", adj=c(0,1.2))
    points(-124.3,43.3, pch=20)
    text(-124.3,43.3,"Charleston", adj=c(0,1.2))
    
  ## Survey
    load('/Users/howar/Documents/Oregon State/Thesis/Data Visualization/trawl_data')
    trawl_data$hauls <- trawl_data$month_day[trawl_data$year_month>0]<-1
    
    nlat=20#determine resolution of grid
    nlon=15
    latd=seq(42,47,length.out=nlat)
    lond=seq(-125,-123.9,length.out=nlon)
    
    
    grid.lon=data.frame(
      lon1=rep(lond[-length(lond)],(nlat-1)),
      lon2=rep(lond[-1],(nlat-1)),
      lon3=rep(lond[-1],(nlat-1)),
      lon4=rep(lond[-length(lond)],(nlat-1)))#make dataframe of just longitude
    
    grid.lat=data.frame(
      lat1=sort(rep(latd[-length(latd)],(nlon-1))),
      lat2=sort(rep(latd[-length(latd)],(nlon-1))),
      lat3=sort(rep(latd[-1],(nlon-1))),
      lat4=sort(rep(latd[-1],(nlon-1))))#lat dataframe
  
    #80s
    dev.new(width=4,height=10)
    plot(trawl_data$longitude[trawl_data$year>1979&trawl_data$year<1990],trawl_data$latitude[trawl_data$year>1979&trawl_data$year<1990],pch='.',ylim=c(42,47),xlim=c(-125,-123.9))
    n.stations2=NA*(1:nrow(grid.lon))
    
    for(i in 1:length(n.stations1)){
      tmp=in.chull(trawl_data$longitude[trawl_data$year>1979&trawl_data$year<1990],trawl_data$latitude[trawl_data$year>1979&trawl_data$year<1990],grid.lon[i,],grid.lat[i,])
      n.stations1[i]=sum(trawl_data$hauls[trawl_data$year>1979&trawl_data$year<1990]*tmp)#This decides what goes into each grid pixel
      points(trawl_data$longitude[tmp],trawl_data$latitude[tmp],col=i,pch=16)
      polygon(grid.lon[i,],grid.lat[i,])
    }
    
    ##90s
    dev.new(width=4,height=10)
    plot(trawl_data$longitude[trawl_data$year>1989&trawl_data$year<2000],trawl_data$latitude[trawl_data$year>1989&trawl_data$year<2000],pch='.',ylim=c(42,47),xlim=c(-125,-123.9))
    n.stations2=NA*(1:nrow(grid.lon))
    
    for(i in 1:length(n.stations2)){
      tmp=in.chull(trawl_data$longitude[trawl_data$year>1989&trawl_data$year<2000],trawl_data$latitude[trawl_data$year>1989&trawl_data$year<2000],grid.lon[i,],grid.lat[i,])
      n.stations2[i]=sum(trawl_data$hauls[trawl_data$year>1989&trawl_data$year<2000]*tmp)#This decides what goes into each grid pixel
      points(trawl_data$longitude[tmp],trawl_data$latitude[tmp],col=i,pch=16)
      polygon(grid.lon[i,],grid.lat[i,])
    }
    
    
    ##00s
    dev.new(width=4,height=10)
    plot(trawl_data$longitude[trawl_data$year>=2000&trawl_data$year<=2009],trawl_data$latitude[trawl_data$year>=2000&trawl_data$year<=2009],pch='.',ylim=c(42,47),xlim=c(-125,-123.9))
    n.stations3=NA*(1:nrow(grid.lon))
    
    for(i in 1:length(n.stations3)){
      tmp=in.chull(trawl_data$longitude[trawl_data$year>=2000&trawl_data$year<=2009],trawl_data$latitude[trawl_data$year>=2000&trawl_data$year<=2009],grid.lon[i,],grid.lat[i,])
      n.stations3[i]=sum(trawl_data$hauls[trawl_data$year>=2000&trawl_data$year<=2009]*tmp)#This decides what goes into each grid pixel
      points(trawl_data$longitude[tmp],trawl_data$latitude[tmp],col=i,pch=16)
      polygon(grid.lon[i,],grid.lat[i,])
    }
    
    ##10s
    dev.new(width=4,height=10)
    plot(trawl_data$longitude[trawl_data$year>=2010&trawl_data$year<=2017],trawl_data$latitude[trawl_data$year>=2010&trawl_data$year<=2017],pch='.',ylim=c(42,47),xlim=c(-125,-123.9))
    n.stations4=NA*(1:nrow(grid.lon))
    
    for(i in 1:length(n.stations4)){
      tmp=in.chull(trawl_data$longitude[trawl_data$year>=2010&trawl_data$year<=2017],trawl_data$latitude[trawl_data$year>=2010&trawl_data$year<=2017],grid.lon[i,],grid.lat[i,])
      n.stations4[i]=sum(trawl_data$hauls[trawl_data$year>=2010&trawl_data$year<=2017]*tmp)#This decides what goes into each grid pixel
      points(trawl_data$longitude[tmp],trawl_data$latitude[tmp],col=i,pch=16)
      polygon(grid.lon[i,],grid.lat[i,])
    }
    
    ##ALL YEARS
    years.data<-length(unique(trawl_data$year))
    z.lat<-(latd[1:(length(latd)-1)]+latd[2:length(latd)])/2
    z.lon<-(lond[1:(length(lond)-1)]+lond[2:length(lond)])/2
    
    ## decades
    eighties.data<-length(unique(trawl_data$year[trawl_data$year>=1980&trawl_data$year<=1989]))
    z.matrix1<-matrix(n.stations1,ncol=length(z.lat),nrow=length(z.lon),byrow=F)/eighties.data
    z.matrix1 <- ifelse(z.matrix1 == 0, NA, z.matrix1)
    nineties.data<-length(unique(trawl_data$year[trawl_data$year>1989&trawl_data$year<2000]))
    z.matrix2<-matrix(n.stations2,ncol=length(z.lat),nrow=length(z.lon),byrow=F)/nineties.data
    z.matrix2 <- ifelse(z.matrix2==0, NA, z.matrix2)
    thousands.data<-length(unique(trawl_data$year[trawl_data$year>1999&trawl_data$year<2010]))
    z.matrix3<-matrix(n.stations3,ncol=length(z.lat),nrow=length(z.lon),byrow=F)/thousands.data
    z.matrix3 <- ifelse(z.matrix3==0, NA, z.matrix3)
    tens.data<-length(unique(trawl_data$year[trawl_data$year>2009&trawl_data$year<2018]))
    z.matrix4<-matrix(n.stations4,ncol=length(z.lat),nrow=length(z.lon),byrow=F)/tens.data
    z.matrix4 <- ifelse(z.matrix4==0, NA, z.matrix4)
    
    ###MAKE MAPS
    
    #windows(width=11,height=5)
    #par(mfrow=c(1,2))
    windows(width=15,height=9)
    par(mfrow=c(1,4),family = 'serif', mar=c(4,5,3,.3)+.1)
    
    image(z.lon,z.lat,z.matrix1,col=hcl.colors(40, "RdYlBu", rev=T),xlim=c(-125,-123.6),ylim=c(42,47), main="Survey Effort 1980s",
          ylab=expression(paste("Latitude ("^0,'N)')), xlab=expression(paste("Longitude ("^0,'W)')), zlim=c(0,4.5))
    map("worldHires",fill=T,col="grey",add=T)
    image.plot(legend.only=T, col=hcl.colors(40, "RdYlBu", rev=T), legend.shrink = 0.2, smallplot = c(.76,.81,.09,.25), 
               legend.cex = 0.7, axis.args=list(cex.axis = 0.9), legend.width = 1, zlim=c(0,4.5), legend.lab = "avg. number of tows")
    contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=1, 
            levels=-100,col= "gray19",add=T)
    contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=1, 
            levels=-200,col= "gray19",add=T)
    points(-124.0535,44.6368, pch=20)
    text(-124.0535,44.6368,"Newport", adj=c(0,1.2))
    points(-123.8313,46.1879, pch=20)
    text(-123.88028,46.13361,"Astoria", adj=c(0,1.2))
    points(-124.3,43.3, pch=20)
    text(-124.3,43.3,"Charleston", adj=c(0,1.2))
    
    image(z.lon,z.lat,z.matrix2,col=hcl.colors(40, "RdYlBu", rev=T),xlim=c(-125,-123.6), zlim=c(0,4.5), main="Survey Effort 1990s",
          ylab=expression(paste("Latitude ("^0,'N)')), xlab=expression(paste("Longitude ("^0,'W)')))
    map("worldHires",fill=T,col="grey",add=T)
    image.plot(legend.only=T, col=hcl.colors(40, "RdYlBu", rev=T), legend.shrink = 0.2, smallplot = c(.76,.81,.09,.25), 
               legend.cex = 0.7, axis.args=list(cex.axis = 0.9), legend.width = 1, zlim=c(0,4.5), legend.lab = "avg. number of tows")
    contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=1, 
            levels=-100,col= "gray19",add=T)
    contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=1, 
            levels=-200,col= "gray19",add=T)
    points(-124.0535,44.6368, pch=20)
    text(-124.0535,44.6368,"Newport", adj=c(0,1.2))
    points(-123.8313,46.1879, pch=20)
    text(-123.88028,46.13361,"Astoria", adj=c(0,1.2))
    points(-124.3,43.3, pch=20)
    text(-124.3,43.3,"Charleston", adj=c(0,1.2))
    
    image(z.lon,z.lat,z.matrix3,col=hcl.colors(40, "RdYlBu", rev=T),xlim=c(-125,-123.6),ylim=c(42,47), main="Survey Effort 2000s",
          ylab=expression(paste("Latitude ("^0,'N)')), xlab=expression(paste("Longitude ("^0,'W)')), zlim=c(0,4.5))
    map("worldHires",fill=T,col="grey",add=T)
    image.plot(legend.only=T, col=hcl.colors(40, "RdYlBu", rev=T), legend.shrink = 0.2, smallplot = c(.76,.81,.09,.25), 
               legend.cex = 0.7, axis.args=list(cex.axis = 0.9), legend.width = 1, zlim=c(0,4.5), legend.lab = "avg. number of tows")
    contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=1, 
            levels=-100,col= "gray19",add=T)
    contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=1, 
            levels=-200,col= "gray19",add=T)
    points(-124.0535,44.6368, pch=20)
    text(-124.0535,44.6368,"Newport", adj=c(0,1.2))
    points(-123.8313,46.1879, pch=20)
    text(-123.88028,46.13361,"Astoria", adj=c(0,1.2))
    points(-124.3,43.3, pch=20)
    text(-124.3,43.3,"Charleston", adj=c(0,1.2))
    
    image(z.lon,z.lat,z.matrix4,col=hcl.colors(40, "RdYlBu", rev=T),xlim=c(-125,-123.6),ylim=c(42,47), main="Survey Effort 2010s",
          ylab=expression(paste("Latitude ("^0,'N)')), xlab=expression(paste("Longitude ("^0,'W)')), zlim=c(0,4.5))
    map("worldHires",fill=T,col="grey",add=T)
    image.plot(legend.only=T, col=hcl.colors(40, "RdYlBu", rev=T), legend.shrink = 0.2, smallplot = c(.76,.81,.09,.25), 
               legend.cex = 0.7, axis.args=list(cex.axis = 0.9), legend.width = 1, zlim=c(0,4.5), legend.lab = "avg. number of tows")
    contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=1, 
            levels=-100,col= "gray19",add=T)
    contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=1, 
            levels=-200,col= "gray19",add=T)
    points(-124.0535,44.6368, pch=20)
    text(-124.0535,44.6368,"Newport", adj=c(0,1.2))
    points(-123.8313,46.1879, pch=20)
    text(-123.88028,46.13361,"Astoria", adj=c(0,1.2))
    points(-124.3,43.3, pch=20)
    text(-124.3,43.3,"Charleston", adj=c(0,1.2))