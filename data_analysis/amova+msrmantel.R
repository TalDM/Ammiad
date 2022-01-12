### Instructions to calculate AMOVA, partial Mantel tests and MSR-Mantel tests ###
## First, it is necessary to having calculated genetic distance and DGGs following the script identity.R

##LOAD R PACKAGES ###
library(data.table)
library(tidyverse)
library(deldir)
library(data.table)
library(pinfsc50)
library(vcfR)
library(ape)
library(phangorn)
library(pegas)
library(adegenet)
library(stats)
library(ggplot2)
library(phytools)

identity_threshold=0.981

latlon2m <- function(lon1,lat1,lon2,lat2) {
#convert coordinates in terms of lat-long to meters following
#https://stackoverflow.com/questions/639695/how-to-convert-latitude-or-longitude-to-meters
    R <- 6378.137                                # radius of earth in Km
    dLat <- (lat2-lat1)*pi/180
    dLon <- (lon2-lon1)*pi/180
    a <- sin((dLat/2))^2 + cos(lat1*pi/180)*cos(lat2*pi/180)*(sin(dLon/2))^2
    c <- 2 * atan2(sqrt(a), sqrt(1-a))
    d <- R * c
    return (d * 1000)                            # distance in meters
}

mydata<-read.csv2("map.csv",header=TRUE)
mydata<-mydata[mydata$pop=="ammiad" & ! is.na(mydata$pop),]
mydata<-mydata[!is.na(mydata$longitude),]
mydata<-mydata[!is.na(mydata$latitude),]
mydata$transect<-0

mydata.t<-mydata[mydata$targetsample %in% rownames(mat),]
mydata.t$year<-as.numeric(mydata.t$year)
mat.t<-mat[rownames(mat) %in% mydata.t$targetsample,]
mat.t<-mat.t[,colnames(mat.t) %in% mydata.t$targetsample]

mat.t2<-identity.matrix[rownames(identity.matrix) %in% rownames(mat.t),rownames(identity.matrix) %in% rownames(mat.t)]

mat.dummy<-data.frame(targetsample=rownames(mat.t))
mydata.t<-merge(mat.dummy,mydata.t)
#colnames(mat.t) %>% head
#mydata.t$targetsample %>% head
mat.space<-matrix(rep(0,nrow(mat.t)*nrow(mat.t)),ncol=nrow(mat.t))
for (i in 1:nrow(mat.t)){
print(i)
for (j in 1:nrow(mat.t)){
mat.space[i,j]<-latlon2m(mydata.t$longitude[i],mydata.t$latitude[i],mydata.t$longitude[j],mydata.t$longitude[j])
}}
mat.time<-matrix(rep(0,nrow(mat.t)*nrow(mat.t)),ncol=nrow(mat.t))
for (i in 1:nrow(mat.t)){
print(i)
for (j in 1:nrow(mat.t)){
mat.time[i,j]<-abs(mydata.t$year[i]-mydata.t$year[j])
}}
mat.habitat<-matrix(rep(0,nrow(mat.t)*nrow(mat.t)),ncol=nrow(mat.t))
for (i in 1:nrow(mat.t)){
print(i)
for (j in 1:nrow(mat.t)){
mat.habitat[i,j]<-as.numeric(mydata.t$habitat[i]==mydata.t$habitat[j])
}}
for (i in 1:nrow(mat.t)){
for (j in 1:nrow(mat.t)){
mat.t[i,j]<-mat.t[j,i]
}}

mat.identity<-apply(mat.t2,MARGIN=1,FUN=function(x) as.numeric(x>=0.981))
rownames(mat.identity)<-rownames(mat.t2)
colnames(mat.identity)<-colnames(mat.t2)

### Mantel tests with time, space and habitat to explain pairwise distances or DGGs ###
multi.mantel(Y=mat.t2, X=list(mat.space,mat.time,mat.habitat), nperm=1000)
multi.mantel(Y=as.dist(1-mat.t2), X=list(as.dist(mat.space),as.dist(mat.time),as.dist(1-mat.habitat)), nperm=1000)
multi.mantel(Y=mat.identity, X=list(mat.space,mat.time,mat.habitat), nperm=1000)
multi.mantel(Y=as.dist(mat.identity), X=list(as.dist(mat.space),as.dist(mat.time),as.dist(1-mat.habitat)), nperm=1000)

### AMOVA ###
mydata.t$habitat<-as.factor(mydata.t$habitat)
d.mat<-as.dist(1-mat.t2)
amova.onlyhabitat<-pegas::amova(formula=d.mat ~ habitat,data=mydata.t)

mydata.t$year<-as.numeric(mydata.t$year)
mydata.t$transects_id<-as.factor(mydata.t$transects_id)
mydata.t$h.year<-as.factor(paste(mydata.t$year,mydata.t$habitat,sep="_"))
amova.habitatyear<-pegas::amova(formula=d.mat ~ habitat/h.year,data=mydata.t)
amova.trhabitatyear<-pegas::amova(formula=d.mat ~ transects_id/habitat/h.year,data=mydata.t)
amova.onlyhabitat<-pegas::amova(formula=d.mat ~ habitat,data=pops)
amova.habitatyear2<-pegas::amova(formula=d.mat ~ habitat/h.year,data=pops)

### MORAN SPECTRAL RANDOMIZATION ###
## load packages
library(fossil)
library(spdep)
library(ade4)
library(vegan)
library(adespatial)
library(adegenet)
## run preliminary mantel
gen.dist<-cailliez(as.dist(1-mat.t2+10^(-6)))
gen.dist<-cailliez(as.dist(1-mat.t2+10^(-6)))
A<-matrix(rep(0,832*832),ncol=832)
A[mat.t2>0.981]<-1
identity.dist<-cailliez(as.dist(1-A+10^(-6)))
space.dist<-earth.dist(data.frame(mydata.t$longitude,mydata.t$latitude))
mrt<-mantel.randtest(gen.dist, as.dist(1-mat.habitat))
mrt.id<-mantel.randtest(identity.dist, as.dist(1-mat.habitat))

## run MSR-Mantel
space.dist.x<-earth.dist(data.frame(mydata.t$longitude,rep(0,nrow(mydata.t))))
space.dist.y<-earth.dist(data.frame(rep(0,nrow(mydata.t)),mydata.t$latitude))
x.min<-which.min(mydata.t$longitude)
y.min<-which.min(mydata.t$latitude)
xcoords<-as.matrix(space.dist.x)[x.min,]
ycoords<-as.matrix(space.dist.y)[y.min,]
space.vecs<-cbind(mydata.t$longitude,mydata.t$latitude)
space.vecs.m<-cbind(xcoords,ycoords)
dist.m<-dist(space.vecs.m)
rownames(space.vecs.m)<-rownames(mat.t2)
nb1 <- graph2nb(gabrielneigh(space.vecs.m,nnmult=10), sym = T)
lw1 <- nb2listw(nb1)
msr.test <- msr(mrt, lw1, 1000)
msr.test.id <- msr(mrt.id, lw1, 1000)

