library(dplyr)
library(tidyr)
library(reshape)
library(reshape2)
library(maps)
library(mapdata)
library(sgeostat)
library(fields)

setwd("C:/Users/howar/Documents/Oregon State/Thesis/Logbook Data/Logbook")
load("logbooks_corrected")

bathy.dat<-read.table('etopo1.xyz',sep='')
names(bathy.dat)<-c('lon','lat','depth')
bathy.dat$depth[bathy.dat$depth>0]<-NA
head(bathy.dat)
bathy.mat<-matrix(bathy.dat$depth,nrow=length(unique(bathy.dat$lon)),ncol=length(unique(bathy.dat$lat)))[,order(unique(bathy.dat$lat))]

#### Logbooks

#Filter to just starry and create p/a and year columns
subset_starry<-combined[combined$species=='STRY_ADJ',]

# Filter to just survey months
subset_starry$month_day <- as.numeric(format(subset_starry$TOWDATE, '%m%d'))
subset_starry <- subset_starry[subset_starry$month_day>=517&subset_starry$month_day<= 929,]
subset_starry <- subset_starry[subset_starry$depth<=-5,]
subset_starry$kg_caught <- subset_starry$species_weight*0.4535924
subset_starry <- subset_starry[!subset_starry$DURATION==0,]
subset_starry$CPUE <- subset_starry$kg_caught/subset_starry$DURATION

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

## Decades
# 80s
dev.new(width=4,height=10)
plot(subset_starry$lon[subset_starry$year>1980&subset_starry$year<1990],subset_starry$lat[subset_starry$year>1980&subset_starry$year<1990],
     pch='.',ylim=c(42,47),xlim=c(-125,-123.9))
n.stations1=NA*(1:nrow(grid.lon))

for(i in 1:length(n.stations1)){
  tmp=in.chull(subset_starry$lon[subset_starry$year>1980&subset_starry$year<1990],
               subset_starry$lat[subset_starry$year>1980&subset_starry$year<1990],grid.lon[i,],grid.lat[i,])
  n.stations1[i]=(sum(subset_starry$CPUE[subset_starry$year>1980&subset_starry$year<1990]*tmp)/sum(1*tmp))/
    (sum(subset_starry$CPUE[subset_starry$year>1980&subset_starry$year<1990])) #This gives proportion of biomass per grid cell (divide avg. weight for cell by biomass for that decade)
  points(subset_starry$longitude[tmp],subset_starry$latitude[tmp],col=i,pch=16)
  polygon(grid.lon[i,],grid.lat[i,])
}

# 90s
dev.new(width=4,height=10)
plot(subset_starry$lon[subset_starry$year>1989&subset_starry$year<2000],subset_starry$lat[subset_starry$year>1989&subset_starry$year<2000],
     pch='.',ylim=c(42,47),xlim=c(-125,-123.9))
n.stations2=NA*(1:nrow(grid.lon))

for(i in 1:length(n.stations2)){
  tmp=in.chull(subset_starry$lon[subset_starry$year>1989&subset_starry$year<2000],
               subset_starry$lat[subset_starry$year>1989&subset_starry$year<2000],grid.lon[i,],grid.lat[i,])
  n.stations2[i]=(sum(subset_starry$CPUE[subset_starry$year>1989&subset_starry$year<2000]*tmp)/sum(1*tmp))/
    (sum(subset_starry$CPUE[subset_starry$year>1989&subset_starry$year<2000])) #This gives proportion of biomass per grid cell
  points(subset_starry$lon[tmp],subset_starry$lat[tmp],col=i,pch=16)
  polygon(grid.lon[i,],grid.lat[i,])
}

# 00s
dev.new(width=4,height=10)
plot(subset_starry$lon[subset_starry$year>1999&subset_starry$year<2010],subset_starry$lat[subset_starry$year>1999&subset_starry$year<2010],
     pch='.',ylim=c(42,47),xlim=c(-125,-123.9))
n.stations3=NA*(1:nrow(grid.lon))

for(i in 1:length(n.stations3)){
  tmp=in.chull(subset_starry$lon[subset_starry$year>1999&subset_starry$year<2010],
               subset_starry$lat[subset_starry$year>1999&subset_starry$year<2010],grid.lon[i,],grid.lat[i,])
  n.stations3[i]=(sum(subset_starry$CPUE[subset_starry$year>1999&subset_starry$year<2010]*tmp)/sum(1*tmp))/
    (sum(subset_starry$CPUE[subset_starry$year>1999&subset_starry$year<2010])) #This gives proportion of biomass per grid cell
  points(subset_starry$lon[tmp],subset_starry$lat[tmp],col=i,pch=16)
  polygon(grid.lon[i,],grid.lat[i,])
}

# 10s
dev.new(width=4,height=10)
plot(subset_starry$lon[subset_starry$year>2009&subset_starry$year<2018],subset_starry$lat[subset_starry$year>2009&subset_starry$year<2018],
     pch='.',ylim=c(42,47),xlim=c(-125,-123.9))
n.stations4=NA*(1:nrow(grid.lon))

for(i in 1:length(n.stations4)){
  tmp=in.chull(subset_starry$lon[subset_starry$year>2009&subset_starry$year<2018],
               subset_starry$lat[subset_starry$year>2009&subset_starry$year<2018],grid.lon[i,],grid.lat[i,])
  n.stations4[i]=(sum(subset_starry$CPUE[subset_starry$year>2009&subset_starry$year<2018]*tmp)/sum(1*tmp))/
    (sum(subset_starry$CPUE[subset_starry$year>2009&subset_starry$year<2018])) #This gives proportion of biomass per grid cell
  points(subset_starry$lon[tmp],subset_starry$lat[tmp],col=i,pch=16)
  polygon(grid.lon[i,],grid.lat[i,])
}


# decades
z.lat<-(latd[1:(length(latd)-1)]+latd[2:length(latd)])/2
z.lon<-(lond[1:(length(lond)-1)]+lond[2:length(lond)])/2
eighties.data.logbooks<-length(unique(subset_starry$year[subset_starry$year>1980&subset_starry$year<1990]))
z.matrix1.logbooks<-matrix(n.stations1,ncol=length(z.lat),nrow=length(z.lon),byrow=F)/eighties.data.logbooks
nineties.data.logbooks<-length(unique(subset_starry$year[subset_starry$year>1989&subset_starry$year<2000]))
z.matrix2.logbooks<-matrix(n.stations2,ncol=length(z.lat),nrow=length(z.lon),byrow=F)/nineties.data.logbooks
thousands.data.logbooks<-length(unique(subset_starry$year[subset_starry$year>1999&subset_starry$year<2010]))
z.matrix3.logbooks<-matrix(n.stations3,ncol=length(z.lat),nrow=length(z.lon),byrow=F)/thousands.data.logbooks
tens.data.logbooks<-length(unique(subset_starry$year[subset_starry$year>2009&subset_starry$year<2018]))
z.matrix4.logbooks<-matrix(n.stations4,ncol=length(z.lat),nrow=length(z.lon),byrow=F)/tens.data.logbooks


#### Survey
load('/Users/howar/Documents/Oregon State/Thesis/Data Visualization/trawl_data')
load('/Users/howar/Documents/Oregon State/Thesis/Data Visualization/OR_fish')

subset_starry.1<-OR_fish[OR_fish$scientific_name=='Platichthys stellatus',]
match_id<-match(trawl_data$trawl_id,subset_starry.1$trawl_id)
trawl_data$CPUE<-subset_starry.1$cpue_kg[match_id]
subset_starry.1 <- trawl_data
subset_starry.1$CPUE[is.na(subset_starry.1$CPUE)]<-0

# Grids
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

### pres
### Used log(x+1) pres for kilograms
dev.new(width=4,height=10)
plot(subset_starry.1$longitude[subset_starry.1$year>1979&subset_starry.1$year<1990],subset_starry.1$latitude[subset_starry.1$year>1979&subset_starry.1$year<1990],
     pch='.',ylim=c(42,47),xlim=c(-125,-123.9))
n.stations5=NA*(1:nrow(grid.lon))

for(i in 1:length(n.stations5)){
  tmp=in.chull(subset_starry.1$longitude[subset_starry.1$year>1979&subset_starry.1$year<1990],
               subset_starry.1$latitude[subset_starry.1$year>1979&subset_starry.1$year<1990],grid.lon[i,],grid.lat[i,])
  n.stations5[i]=(sum(subset_starry.1$CPUE[subset_starry.1$year>1979&subset_starry.1$year<1990]*tmp)/sum(1*tmp))/
    (sum(subset_starry.1$CPUE[subset_starry.1$year>1979&subset_starry.1$year<1990]))#This decides what goes into each grid pixel
  points(subset_starry.1$longitude[tmp],subset_starry.1$latitude[tmp],col=i,pch=16)
  polygon(grid.lon[i,],grid.lat[i,])
}

## 90s
dev.new(width=4,height=10)
plot(subset_starry.1$longitude[subset_starry.1$year>1989&subset_starry.1$year<2000],subset_starry.1$latitude[subset_starry.1$year>1989&subset_starry.1$year<2000],
     pch='.',ylim=c(42,47),xlim=c(-125,-123.9))
n.stations6=NA*(1:nrow(grid.lon))

for(i in 1:length(n.stations6)){
  tmp=in.chull(subset_starry.1$longitude[subset_starry.1$year>1989&subset_starry.1$year<2000],
               subset_starry.1$latitude[subset_starry.1$year>1989&subset_starry.1$year<2000],grid.lon[i,],grid.lat[i,])
  n.stations6[i]=(sum(subset_starry.1$CPUE[subset_starry.1$year>1989&subset_starry.1$year<2000]*tmp)/sum(1*tmp))/
    (sum(subset_starry.1$CPUE[subset_starry.1$year>1989&subset_starry.1$year<2000]))#This decides what goes into each grid pixel
  points(subset_starry.1$longitude[tmp],subset_starry.1$latitude[tmp],col=i,pch=16)
  polygon(grid.lon[i,],grid.lat[i,])
}

## 00s
dev.new(width=4,height=10)
plot(subset_starry.1$longitude[subset_starry.1$year>1999&subset_starry.1$year<2010],subset_starry.1$latitude[subset_starry.1$year>1999&subset_starry.1$year<2010],pch='.',ylim=c(42,47),xlim=c(-125,-123.9))
n.stations7=NA*(1:nrow(grid.lon))

for(i in 1:length(n.stations7)){
  tmp=in.chull(subset_starry.1$longitude[subset_starry.1$year>1999&subset_starry.1$year<2010],
               subset_starry.1$latitude[subset_starry.1$year>1999&subset_starry.1$year<2010],grid.lon[i,],grid.lat[i,])
  n.stations7[i]=(sum(subset_starry.1$CPUE[subset_starry.1$year>1999&subset_starry.1$year<2010]*tmp)/sum(1*tmp))/
    (sum(subset_starry.1$CPUE[subset_starry.1$year>1999&subset_starry.1$year<2010]))#This decides what goes into each grid pixel
  points(subset_starry.1$longitude[tmp],subset_starry.1$latitude[tmp],col=i,pch=16)
  polygon(grid.lon[i,],grid.lat[i,])
}

## 10s
dev.new(width=4,height=10)
plot(subset_starry.1$longitude[subset_starry.1$year>2009&subset_starry.1$year<2018],subset_starry.1$latitude[subset_starry.1$year>2009&subset_starry.1$year<2018],pch='.',ylim=c(42,47),xlim=c(-125,-123.9))
n.stations8=NA*(1:nrow(grid.lon))

for(i in 1:length(n.stations8)){
  tmp=in.chull(subset_starry.1$longitude[subset_starry.1$year>2009&subset_starry.1$year<2018],
               subset_starry.1$latitude[subset_starry.1$year>2009&subset_starry.1$year<2018],grid.lon[i,],grid.lat[i,])
  n.stations8[i]=(sum(subset_starry.1$CPUE[subset_starry.1$year>2009&subset_starry.1$year<2018]*tmp)/sum(1*tmp))/
    (sum(subset_starry.1$CPUE[subset_starry.1$year>2009&subset_starry.1$year<2018]))#This decides what goes into each grid pixel
  points(subset_starry.1$longitude[tmp],subset_starry.1$latitude[tmp],col=i,pch=16)
  polygon(grid.lon[i,],grid.lat[i,])
}

##ALL YEARS
years.data<-length(unique(subset_starry.1$year))
z.lat<-(latd[1:(length(latd)-1)]+latd[2:length(latd)])/2
z.lon<-(lond[1:(length(lond)-1)]+lond[2:length(lond)])/2

## decades
eighties.data.survey<-length(unique(subset_starry.1$year[subset_starry.1$year>1980&subset_starry.1$year<1990]))
z.matrix1.survey<-matrix(n.stations5,ncol=length(z.lat),nrow=length(z.lon),byrow=F)/eighties.data.survey
nineties.data.survey<-length(unique(subset_starry.1$year[subset_starry.1$year>1989&subset_starry.1$year<2000]))
z.matrix2.survey<-matrix(n.stations6,ncol=length(z.lat),nrow=length(z.lon),byrow=F)/nineties.data.survey
thousands.data.survey<-length(unique(subset_starry.1$year[subset_starry.1$year>1999&subset_starry.1$year<2010]))
z.matrix3.survey<-matrix(n.stations7,ncol=length(z.lat),nrow=length(z.lon),byrow=F)/thousands.data.survey
tens.data.survey<-length(unique(subset_starry.1$year[subset_starry.1$year>2009&subset_starry.1$year<2018]))
z.matrix4.survey<-matrix(n.stations8,ncol=length(z.lat),nrow=length(z.lon),byrow=F)/tens.data.survey


#### Spatial indicators
# local index of collocation function (just gives singular value)
loc_collocfn <- function(prey, pred) {
  p_prey <- prey/sum(prey, na.rm = T)
  p_pred <- pred/sum(pred, na.rm = T)
  sum(p_prey*p_pred, na.rm = T)/(sqrt(sum(p_prey^2, na.rm = T)*sum(p_pred^2, na.rm = T)))
}

# local index of collocation per decade for starry sole overall
loc_collocfn(z.matrix1.logbooks, z.matrix1.survey)
loc_collocfn(z.matrix2.logbooks, z.matrix2.survey)
loc_collocfn(z.matrix3.logbooks, z.matrix3.survey)
loc_collocfn(z.matrix4.logbooks, z.matrix4.survey)

# Map out cell values for LIC
eighties.lic <- as.matrix(z.matrix1.logbooks*z.matrix1.survey, na.rm = T)/(sqrt(sum(z.matrix1.logbooks^2, na.rm = T)*sum(z.matrix1.survey^2, na.rm = T)))
nineties.lic <- as.matrix(z.matrix2.logbooks*z.matrix2.survey, na.rm = T)/(sqrt(sum(z.matrix2.logbooks^2, na.rm = T)*sum(z.matrix2.survey^2, na.rm = T)))
thousands.lic <- as.matrix(z.matrix3.logbooks*z.matrix3.survey, na.rm = T)/(sqrt(sum(z.matrix3.logbooks^2, na.rm = T)*sum(z.matrix3.survey^2, na.rm = T)))
tens.lic <- as.matrix(z.matrix4.logbooks*z.matrix4.survey, na.rm = T)/(sqrt(sum(z.matrix4.logbooks^2, na.rm = T)*sum(z.matrix4.survey^2, na.rm = T)))

n.stations1<-matrix(n.stations1,ncol=length(z.lat),nrow=length(z.lon),byrow=F)/eighties.data.logbooks
n.stations2<-matrix(n.stations2,ncol=length(z.lat),nrow=length(z.lon),byrow=F)/nineties.data.logbooks
n.stations3<-matrix(n.stations3,ncol=length(z.lat),nrow=length(z.lon),byrow=F)/thousands.data.logbooks
n.stations4<-matrix(n.stations4,ncol=length(z.lat),nrow=length(z.lon),byrow=F)/tens.data.logbooks
n.stations5<-matrix(n.stations5,ncol=length(z.lat),nrow=length(z.lon),byrow=F)/eighties.data.survey
n.stations6<-matrix(n.stations6,ncol=length(z.lat),nrow=length(z.lon),byrow=F)/nineties.data.survey
n.stations7<-matrix(n.stations7,ncol=length(z.lat),nrow=length(z.lon),byrow=F)/thousands.data.survey
n.stations8<-matrix(n.stations8,ncol=length(z.lat),nrow=length(z.lon),byrow=F)/tens.data.survey

eighties.lic[is.nan(n.stations5)]<- 0
eighties.lic[is.nan(n.stations1)]<- NA
nineties.lic[is.nan(n.stations6)]<- 0
nineties.lic[is.nan(n.stations2)]<- NA
thousands.lic[is.nan(n.stations7)]<- 0
thousands.lic[is.nan(n.stations3)]<- NA
tens.lic[is.nan(n.stations8)]<- 0
tens.lic[is.nan(n.stations4)]<- NA


windows(width=15,height=9)
par(mfrow=c(1,4),family = 'serif', mar=c(4,5,3,.3)+.1)

image(z.lon,z.lat,eighties.lic,col=hcl.colors(40, "viridis", rev=F),xlim=c(-125,-123.6),ylim=c(42,47), main="LIC 1980s",
           ylab=expression(paste("Latitude ("^0,'N)')), xlab=expression(paste("Longitude ("^0,'W)')), zlim = c(0,0.035))
map("worldHires",fill=T,col="grey",add=T)
image.plot(legend.only=T, col=hcl.colors(40, "viridis", rev=F), legend.shrink = 0.2, smallplot = c(.76,.81,.09,.25), 
           legend.cex = 0.7, axis.args=list(cex.axis = 0.9), legend.width = 1, zlim = c(0,0.035), legend.lab = "LIC")
contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=1, 
        levels=-55,col= "gray40",add=T)
contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=1, 
        levels=-200,col= "gray40",add=T)
points(-124.0535,44.6368, pch=20)
text(-124.0535,44.6368,"Newport", adj=c(0,1.2))
points(-123.8313,46.1879, pch=20)
text(-123.88028,46.13361,"Astoria", adj=c(0,1.2))
points(-124.3,43.3, pch=20)
text(-124.3,43.3,"Charleston", adj=c(0,1.2))

image(z.lon,z.lat,nineties.lic,col=hcl.colors(40, "viridis", rev=F), ylab=expression(paste("Latitude ("^0,'N)')), xlab=expression(paste("Longitude ("^0,'W)')),
           xlim=c(-125,-123.6),ylim=c(42,47), main="LIC 1990s", zlim = c(0,0.035))
map("worldHires",fill=T,col="grey",add=T)
image.plot(legend.only=T, col=hcl.colors(40, "viridis", rev=F), legend.shrink = 0.2, smallplot = c(.76,.81,.09,.25), 
           legend.cex = 0.7, axis.args=list(cex.axis = 0.9), legend.width = 1, zlim = c(0,0.035), legend.lab = "LIC")
contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=1, 
        levels=-55,col= "gray40",add=T)
contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=1, 
        levels=-200,col= "gray40",add=T)
points(-124.0535,44.6368, pch=20)
text(-124.0535,44.6368,"Newport", adj=c(0,1.2))
points(-123.8313,46.1879, pch=20)
text(-123.88028,46.13361,"Astoria", adj=c(0,1.2))
points(-124.3,43.3, pch=20)
text(-124.3,43.3,"Charleston", adj=c(0,1.2))

image(z.lon,z.lat,thousands.lic,col=hcl.colors(40, "viridis", rev=F), ylab=expression(paste("Latitude ("^0,'N)')), xlab=expression(paste("Longitude ("^0,'W)')),
           xlim=c(-125,-123.6),ylim=c(42,47), main="LIC 2000s", zlim = c(0,0.035))
map("worldHires",fill=T,col="grey",add=T)
image.plot(legend.only=T, col=hcl.colors(40, "viridis", rev=F), legend.shrink = 0.2, smallplot = c(.76,.81,.09,.25), 
           legend.cex = 0.7, axis.args=list(cex.axis = 0.9), legend.width = 1, zlim = c(0,0.035), legend.lab = "LIC")
contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=1, 
        levels=-55,col= "gray40",add=T)
contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=1, 
        levels=-200,col= "gray40",add=T)
points(-124.0535,44.6368, pch=20)
text(-124.0535,44.6368,"Newport", adj=c(0,1.2))
points(-123.8313,46.1879, pch=20)
text(-123.88028,46.13361,"Astoria", adj=c(0,1.2))
points(-124.3,43.3, pch=20)
text(-124.3,43.3,"Charleston", adj=c(0,1.2))

image(z.lon,z.lat,tens.lic,col=hcl.colors(40, "viridis", rev=F), ylab=expression(paste("Latitude ("^0,'N)')), xlab=expression(paste("Longitude ("^0,'W)')),
           xlim=c(-125,-123.6),ylim=c(42,47), main="LIC 2010s", zlim = c(0,0.035))
map("worldHires",fill=T,col="grey",add=T)
image.plot(legend.only=T, col=hcl.colors(40, "viridis", rev=F), legend.shrink = 0.2, smallplot = c(.76,.81,.09,.25), 
           legend.cex = 0.7, axis.args=list(cex.axis = 0.9), legend.width = 1, zlim = c(0,0.035), legend.lab = "LIC")
contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=1, 
        levels=-55,col= "gray40",add=T)
contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=1, 
        levels=-200,col= "gray40",add=T)
points(-124.0535,44.6368, pch=20)
text(-124.0535,44.6368,"Newport", adj=c(0,1.2))
points(-123.8313,46.1879, pch=20)
text(-123.88028,46.13361,"Astoria", adj=c(0,1.2))
points(-124.3,43.3, pch=20)
text(-124.3,43.3,"Charleston", adj=c(0,1.2))