#Install packages
install.packages(c("dismo","terra","spatstat","cowplot","ggplot2","precrec","glmnet","maxnet","ranger"))
#We also need to load the "MLR" package for this practical
install.packages("mlr")
#Function to convert raster to images for spatstat. Takes one argument "im" that should be a raster object
install.packages(c("terra","sf","mapview"))
library(terra)
library(sf)
library(spatstat) #for point process modelling and converting between raster and pixel image objectslibrary(spatstat)


#for point process modelling and converting between raster and pixel image objects
raster.as.im = function(im) {
  #get the resolution (cell size of the raster)
  r = raster::res(im)[1]
  #get the origin (bottom left corner of the raster/image)
  orig = ext(im)[c(1,3)]
  #set the coordinates of the columns which is just a series of number from zero increasing by 100 metres (the resolution of the raster) for every cell along the rows and columns.
  xx = orig[1] + seq(from=0,to=(ncol(im) - 1)*100,by=r)
  #set the coordinates of the columns
  yy = orig[2] + seq(from=0,to=(nrow(im) - 1)*100,by=r)
  
  mat=matrix(raster::values(im), ncol = ncol(im), 
             nrow = nrow(im), byrow = TRUE)[nrow(im):1, ]
  return(spatstat.geom::im(mat, 
                           xcol = xx, yrow = yy))
}



#Read and clean species data
melesmeles<- read.csv("Melesmeles.csv")

#data cleaning to ensure there are no NA
melesmeles<-melesmeles[!is.na(melesmeles$Latitude),]

#remove all points with uncertainty > 1000m
melesmeles<-melesmeles[melesmeles$Coordinate.uncertainty_m<1001,]

#create crs object
melesmeles.latlong<-data.frame(x=melesmeles$Longitude,y=melesmeles$Latitude)

#Use coordinates object to create our spatial points object
melesmeles.sp=st_as_sf(melesmeles.latlong,coords=c("x","y"),crs="epsg:4326")

plot(melesmeles.sp)



#First set the extent to the study area
scot=st_read('scotSamp.shp')

#load in the land cover map and then clip to the polygon
LCM=rast("LCMUK.tif")
#crop to the extent of the study area plus a little more (because we will lose a small amount of data in the next step)
LCM=crop(LCM,st_buffer(scot, dist= 1000))
#aggregate LCM raster
LCM=aggregate(LCM$LCMUK_1,fact=4,fun="modal")
#project squirrel data
melesmeles.sp=st_transform(melesmeles.sp,crs(LCM))

#now crop our points to the study area
melesmelesFin=melesmeles.sp[scot,]

#finally, mask the LCM to this boundary

LCM=crop(LCM,scot,mask=TRUE)
plot(LCM)
plot(melesmelesFin$geometry,add=T)



#access levels of the raster by treating them as categorical data ('factors' in R)
LCM=as.factor(LCM$LCMUK_1)

#create an vector object called reclass
reclass = c(0,1,rep(0,20))

# combine with the LCM categories into a matrix of old and new values.
RCmatrix=cbind(levels(LCM)[[1]],reclass)

RCmatrix=RCmatrix[,2:3]

#apply function to make sure new columns are numeric (here the "2" specifies that we want to apply the as.numeric function to columns, where "1" would have specified rows)
RCmatrix=apply(RCmatrix,2,FUN=as.numeric)
#Use the classify() function to asssign new values to LCM with our reclassification matrix

broadleaf=classify(LCM, RCmatrix)

#neighbourhood weights matrix to sum all available resources for each cell

#get number of picels needed to cover the 1800 metre radius 
nPix=round(900/res(LCM)[1])

#next, you need to double this number
nPix=(nPix*2)+1

#buiild weights matrix
weightsMatrix=matrix(1:nPix^2,nrow=nPix,ncol=nPix)

#get focal cell 
x=ceiling(ncol(weightsMatrix)/2)
y=ceiling(nrow(weightsMatrix)/2)


focalCell=weightsMatrix[x,y]

indFocal=which(weightsMatrix==focalCell,arr.ind = TRUE)

#compute distances
distances=list()

for(i in 1:nPix^2){
  ind.i=which(weightsMatrix==i,arr.ind=T)
  diffX=abs(ind.i[1,1]-indFocal[1,1])*res(LCM)[1]
  diffY=abs(ind.i[1,2]-indFocal[1,2])*res(LCM)[1]
  
  dist.i=sqrt(diffX^2+diffY^2)
  distances[[i]]=dist.i
  
}

#add distance values to the weights matrix
weightsMatrix[]=unlist(distances)

#set cells outside search radius to NA
weightsMatrix[weightsMatrix>900]=NA

#plot weights matrix
plot(rast(weightsMatrix))

#normalise the weights matrix by dividing all cell values by the number of cells. 
weightsMatrixNorm=weightsMatrix
weightsMatrixNorm[!is.na(weightsMatrixNorm)]=1/length(weightsMatrixNorm[!is.na(weightsMatrixNorm)])

#test to see for yourself
sum(weightsMatrixNorm,na.rm=T)

plot(rast(weightsMatrixNorm))

#sum neighbourhood values from all surrounding cells
lcm_wood_900=focal(broadleaf,w=weightsMatrixNorm,fun="sum")

plot(lcm_wood_900)