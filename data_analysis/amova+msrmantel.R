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
#https://stackoverflow.com/questions/639695/how-to-convert-Latitude-or-Longitude-to-meters
    R <- 6378.137                                # radius of earth in Km
    dLat <- (lat2-lat1)*pi/180
    dLon <- (lon2-lon1)*pi/180
    a <- sin((dLat/2))^2 + cos(lat1*pi/180)*cos(lat2*pi/180)*(sin(dLon/2))^2
    c <- 2 * atan2(sqrt(a), sqrt(1-a))
    d <- R * c
    return (d * 1000)                            # distance in meters
}

mydata<-read.csv("../data/dataset_S1.csv",header=TRUE)

mydata<-mydata[mydata$Position!="Zavitan",]
mydata<-mydata[!is.na(mydata$Longitude),]
mydata<-mydata[!is.na(mydata$Latitude),]

mat<-read.table("../data/fstmat.tsv",header=TRUE,row.names=1)
mat<-as.matrix(mat)
colnames(mat)<-gsub("^X","",colnames(mat))

mydata.t<-mydata[mydata$Sample %in% rownames(mat),]
mydata.t$Year<-as.numeric(mydata.t$Year)
mat.t<-mat[rownames(mat) %in% mydata.t$Sample,]
mat.t<-mat.t[,colnames(mat.t) %in% mydata.t$Sample]

identity.matrix<-read.table("../data/identity.matrix.tsv",sep="\t",header=T)
identity.matrix<-as.matrix(identity.matrix)
colnames(identity.matrix)<-gsub("^X","",colnames(identity.matrix))

mat.t2<-identity.matrix[rownames(identity.matrix) %in% rownames(mat.t),rownames(identity.matrix) %in% rownames(mat.t)]

mat.dummy<-data.frame(targetsample=rownames(mat.t))
mydata.t<-merge(mat.dummy,mydata.t)
#colnames(mat.t) %>% head
#mydata.t$targetsample %>% head
mat.space<-matrix(rep(0,nrow(mat.t)*nrow(mat.t)),ncol=nrow(mat.t))
for (i in 1:nrow(mat.t)){
print(i)
for (j in 1:nrow(mat.t)){
mat.space[i,j]<-latlon2m(mydata.t$Longitude[i],mydata.t$Latitude[i],mydata.t$Longitude[j],mydata.t$Longitude[j])
}}
mat.time<-matrix(rep(0,nrow(mat.t)*nrow(mat.t)),ncol=nrow(mat.t))
for (i in 1:nrow(mat.t)){
print(i)
for (j in 1:nrow(mat.t)){
mat.time[i,j]<-abs(mydata.t$Year[i]-mydata.t$Year[j])
}}
mat.Habitat<-matrix(rep(0,nrow(mat.t)*nrow(mat.t)),ncol=nrow(mat.t))
for (i in 1:nrow(mat.t)){
print(i)
for (j in 1:nrow(mat.t)){
mat.Habitat[i,j]<-as.numeric(mydata.t$Habitat[i]==mydata.t$Habitat[j])
}}
for (i in 1:nrow(mat.t)){
for (j in 1:nrow(mat.t)){
mat.t[i,j]<-mat.t[j,i]
}}

mat.identity<-apply(mat.t2,MARGIN=1,FUN=function(x) as.numeric(x>=0.981))
rownames(mat.identity)<-rownames(mat.t2)
colnames(mat.identity)<-colnames(mat.t2)

### Mantel tests with time, space and Habitat to explain pairwise distances or DGGs ###
multi.mantel(Y=mat.t2, X=list(mat.space,mat.time,mat.Habitat), nperm=1000)
multi.mantel(Y=as.dist(1-mat.t2), X=list(as.dist(mat.space),as.dist(mat.time),as.dist(1-mat.Habitat)), nperm=1000)
multi.mantel(Y=mat.identity, X=list(mat.space,mat.time,mat.Habitat), nperm=1000)
multi.mantel(Y=as.dist(mat.identity), X=list(as.dist(mat.space),as.dist(mat.time),as.dist(1-mat.Habitat)), nperm=1000)

### AMOVA ###
mydata.t$Habitat<-as.factor(mydata.t$Habitat)
d.mat<-as.dist(1-mat.t2)
amova.onlyHabitat<-pegas::amova(formula=d.mat ~ Habitat,data=mydata.t)

mydata.t$Year<-as.numeric(mydata.t$Year)
mydata.t$h.Year<-as.factor(paste(mydata.t$Year,mydata.t$Habitat,sep="_"))
amova.HabitatYear<-pegas::amova(formula=d.mat ~ Habitat/h.Year,data=mydata.t)
#amova.trHabitatYear<-pegas::amova(formula=d.mat ~ Transect_order/Habitat/h.Year,data=mydata.t)
amova.onlyHabitat<-pegas::amova(formula=d.mat ~ Habitat,data=mydata.t)
amova.HabitatYear2<-pegas::amova(formula=d.mat ~ Habitat/h.Year,data=mydata.t)
mydata.t$Year<-as.factor(mydata.t$Year)
amova.onlYear<-pegas::amova(formula=d.mat ~ Year,data=mydata.t)
"
Call: pegas::amova(formula = d.mat ~ Year, data = mydata.t)
             SSD        MSD  df
Year   0.1197359 0.01496699   8
Error 12.1002699 0.01470264 823
Total 12.2200058 0.01470518 831

Variance components:
          sigma2 P.value
Year  2.8618e-06  0.4226
Error 1.4703e-02

Phi-statistics:
Year.in.GLOBAL
  0.0001946065

Variance coefficients:
     a
92.372
"
"
        Analysis of Molecular Variance

Call: pegas::amova(formula = d.mat ~ Habitat, data = mydata.t)

              SSD        MSD  df
Habitat  3.134754 0.44782200   7
Error    9.085252 0.01102579 824
Total   12.220006 0.01470518 831

Variance components:
           sigma2 P.value
Habitat 0.0043384       0
Error   0.0110258

Phi-statistics:
Habitat.in.GLOBAL
        0.2823704

Variance coefficients:
       a
100.6817

"
"Call: pegas::amova(formula = d.mat ~ Habitat/h.Year, data = mydata.t)

               SSD        MSD  df
Habitat  3.1347540 0.44782200   7
h.Year   0.7146763 0.01116682  64
Error    8.3705755 0.01101392 760
Total   12.2200058 0.01470518 831

Variance components:
            sigma2 P.value
Habitat 4.3370e-03  0.0000
h.Year  1.3254e-05  0.7902
Error   1.1014e-02

Phi-statistics:
Habitat.in.GLOBAL (Phi_CT)  h.Year.in.GLOBAL (Phi_ST)
                0.28228083                 0.28314346
h.Year.in.Habitat (Phi_SC)
                0.00120191

Variance coefficients:
        a         b         c
 11.53661  11.33873 100.68166
"

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
space.dist<-earth.dist(data.frame(mydata.t$Longitude,mydata.t$Latitude))
mrt<-mantel.randtest(gen.dist, as.dist(1-mat.Habitat))
mrt.id<-mantel.randtest(identity.dist, as.dist(1-mat.Habitat))

## run MSR-Mantel
space.dist.x<-earth.dist(data.frame(mydata.t$Longitude,rep(0,nrow(mydata.t))))
space.dist.y<-earth.dist(data.frame(rep(0,nrow(mydata.t)),mydata.t$Latitude))
x.min<-which.min(mydata.t$Longitude)
y.min<-which.min(mydata.t$Latitude)
xcoords<-as.matrix(space.dist.x)[x.min,]
ycoords<-as.matrix(space.dist.y)[y.min,]
space.vecs<-cbind(mydata.t$Longitude,mydata.t$Latitude)
space.vecs.m<-cbind(xcoords,ycoords)
dist.m<-dist(space.vecs.m)
rownames(space.vecs.m)<-rownames(mat.t2)
nb1 <- graph2nb(gabrielneigh(space.vecs.m,nnmult=10), sym = T)
lw1 <- nb2listw(nb1)
msr.test <- msr(mrt, lw1, 1000)
msr.test.id <- msr(mrt.id, lw1, 1000)
