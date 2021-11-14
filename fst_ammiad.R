#============================================================================
#==================================Fst-Habitat===============================
#============================================================================
{
#
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
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("snpStats")
#BiocManager::install("SNPRelate")
#BiocManager::install("GWASTools")
#BiocManager::install("SeqArray")
#BiocManager::install("hierfstat")
library("snpStats")
library("SNPRelate")
library("hierfstat")
#library("GWASTools")
#library("SeqArray")

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

par(mfrow=c(1,1))
mydata<-read.table("map_file.txt",header=TRUE)
mydata<-mydata[mydata$pop=="ammiad" & ! is.na(mydata$pop),]
mydata<-mydata[!is.na(mydata$longitude),]
mydata<-mydata[!is.na(mydata$latitude),]
mydata$transect<-0


#partition sites in transects with diameter of 30m (radius 15m) within each habitat
itransect<-0
distance_radius<-30
for (xhabitat in unique(mydata$habitat)){
    mydata_t<-mydata %>% filter(habitat==xhabitat,transect==0)
    mins<-mydata %>% filter(habitat==xhabitat) %>% summarise(minlat=min(latitude),minlon=min(longitude)) %>% as.numeric
    maxs<-mydata %>% filter(habitat==xhabitat) %>% summarise(maxlat=max(latitude),maxlon=max(longitude)) %>% as.numeric
    if ((maxs[1]-mins[1])>(maxs[2]-mins[2])){ lat_is_wider<-TRUE } else { lat_is_wider<-FALSE }
        while(sum(mydata_t$transect==0)>0)
        {        
        itransect<-itransect+1
        mins<-mydata_t %>% filter(transect==0) %>% summarise(minlat=min(latitude),minlon=min(longitude)) %>% as.numeric
        maxs<-mydata_t %>% filter(transect==0) %>% summarise(maxlat=max(latitude),maxlon=max(longitude)) %>% as.numeric
        #latlon2m(mins[2],mins[1],maxs[2],maxs[1]) #113m
        #choose point1: most extreme. to decide between left and bottom, I use overall distribution to understand which axis is more relevant.
        if (lat_is_wider){point1<-mydata_t %>% filter(transect==0,latitude==mins[1]) %>% select(latitude,longitude) %>% head(1) %>% as.numeric} else 
        {point1<-mydata_t %>% filter(transect==0,longitude==mins[2]) %>% select(latitude,longitude) %>% head(1) %>% as.numeric}
        dist_t<-apply(mydata_t %>% select(latitude,longitude),MARGIN=1,FUN=function(x) latlon2m(point1[2],point1[1],x[2],x[1]))
        mydata_t$transect[dist_t<distance_radius & mydata_t$transect==0]<-itransect
        }
        mydata$transect[mydata$habitat==xhabitat]<-mydata_t$transect
}


#plot map of parcels and habitats
pdf(paste0("transects",distance_radius,"m.pdf"))
par(mfrow=c(1,2))
plot(mydata$longitude,mydata$latitude,col=as.factor(mydata$habitat),pch=19,cex=0.5,xlab="longitude",ylab="latitude",main="habitat")
plot(mydata$longitude,mydata$latitude,col=as.factor(mydata$transect),pch=19,cex=0.5,xlab="longitude",ylab="latitude",main=paste("transects ",distance_radius,"m"))
dev.off()

#read vcf and calculate pairwise-Fst
vcf<-read.vcfR("vcf/ammiad_only_filtered_hets_845.vcf")
vcf.loci<-vcfR2loci(vcf)
vcf.genind<-vcfR2genind(vcf)
vcf.labels<-VCFlabels("vcf/ammiad_only_filtered_hets_845.vcf")
pops<-vcf.labels
pops[!vcf.labels %in% mydata$targetsample]<-"other"
pops<-left_join(data.frame(targetsample=pops),mydata)
pops$transect[is.na(pops$transect)]<-max(pops$transect,na.rm=TRUE)+1
snpgdsVCF2GDS("ammiad_only_filtered_hets_845.vcf", "ammiad_only_filtered_hets_845.gds")
ammiad.gds<-snpgdsOpen("ammiad_only_filtered_hets_845.gds")
transects.fst<-snpgdsFst(ammiad.gds,as.factor(pops$transect),"W&H02",verbose=TRUE)
matrixFst<-transects.fst$Beta

vcf.genind@pop<-as.factor(pops$transect)
mathierfstat<-genind2hierfstat(vcf.genind)
matFst <- pairwise.WCfst(mathierfstat,diploid=TRUE)


#Partition in between and within habitat comparisons
mydata<-mydata %>% group_by(transect) %>% summarise(long_transect=mean(longitude),lat_transect=mean(latitude)) %>% left_join(mydata)
mat_v<-c()
samehab<-c()
dist.transects_v<-c()
for (i in 1:(nrow(matFst)-1)){
for (j in (i+1):nrow(matFst)){
mat_v<-c(mat_v,matFst[i,j])
habi<-mydata$habitat[mydata$transect==rownames(matFst)[i]][1]
habj<-mydata$habitat[mydata$transect==rownames(matFst)[j]][1]
if (habi==habj){ samehab<-c(samehab,"within") } else { samehab<-c(samehab,"between") }
long_i<-mydata$long_transect[mydata$transect==rownames(matFst)[i]][1]
lat_i<-mydata$lat_transect[mydata$transect==rownames(matFst)[i]][1]
long_j<-mydata$long_transect[mydata$transect==rownames(matFst)[j]][1]
lat_j<-mydata$lat_transect[mydata$transect==rownames(matFst)[j]][1]
dist.transects_v<-c(dist.transects_v,latlon2m(long_i,lat_i,long_j,lat_j))
}
}

Fst2Nm<-function(x) (1-x)/(4*x)

#add as second axis m estimated by Fst: note that this m should be considered a ballpark estimate and not an accurate one,
#see for example: Whitlock et al.,1999, Heredity
#FST!=1/(4Nm+1) #->#Nm=(1-Fst)/(4*Fst)  

data.plot<-data.frame(distance=dist.transects_v,habitat=as.factor(samehab),Fst=mat_v,Nm=(1-mat_v)/(4*mat_v))
p<-ggplot(data.plot,aes(x=distance,y=Fst,colour=habitat))+geom_point()+scale_x_log10()+stat_smooth(method=glm,method.args=list(family=binomial))+
scale_y_continuous("Fst",sec.axis = sec_axis(~ (1 - .)/(4*.) , name = "Nm"))+theme_bw()
ggsave(p,file="Fstdecay_30m_binomial.pdf")
p<-ggplot(data.plot,aes(x=distance,y=Fst,colour=habitat))+geom_point()+scale_x_log10()+stat_smooth(method=loess,method.args=list(degree=1))+theme_bw()
#+scale_y_continuous("Fst",sec.axis = sec_axis(~ (1 - .)/(4*.) , name = "Nm"))
ggsave(p,file="Fstdecay_30m_loess.pdf")

data.obs<-data.plot[data.plot$distance<50,]
data.obs<-data.obs %>% group_by(habitat) %>% summarise(Fst=mean(Fst)) %>% ungroup
data.obs$bs<-0

samplei<-c();samplej<-c()
for (i in 1:(nrow(matFst)-1)){
for (j in (i+1):nrow(matFst)){
samplei<-c(samplei,i);samplej<-c(samplej,j)}}

bs.n<-1000
bs.v<-c()
mydata<-mydata %>% group_by(transect) %>% summarise(long_transect=mean(longitude),lat_transect=mean(latitude)) %>% left_join(mydata)
mat_v<-c()
samehab<-c()
dist.transects_v<-c()
for (bs.iteration in 1:bs.n){
newsample<-sample(1:length(samplei),replace=TRUE)
for (i in 1:length(samplei)){
i.bs<-samplei[newsample[i]]
j.bs<-samplej[newsample[i]]
mat_v<-c(mat_v,matFst[i.bs,j.bs])
habi<-mydata$habitat[mydata$transect==rownames(matFst)[i.bs]][1]
habj<-mydata$habitat[mydata$transect==rownames(matFst)[j.bs]][1]
if (habi==habj){ samehab<-c(samehab,"within") } else { samehab<-c(samehab,"between") }
long_i<-mydata$long_transect[mydata$transect==rownames(matFst)[i.bs]][1]
lat_i<-mydata$lat_transect[mydata$transect==rownames(matFst)[i.bs]][1]
long_j<-mydata$long_transect[mydata$transect==rownames(matFst)[j.bs]][1]
lat_j<-mydata$lat_transect[mydata$transect==rownames(matFst)[j.bs]][1]
dist.transects_v<-c(dist.transects_v,latlon2m(long_i,lat_i,long_j,lat_j))
bs.v<-c(bs.v,bs.iteration)
}}
data.bs<-data.frame(distance=dist.transects_v,habitat=as.factor(samehab),Fst=mat_v,bs=bs.v)
data.bs.t<-data.bs[data.bs$distance<50,]
data.plot<-data.bs.t %>% group_by(bs,habitat) %>% summarise(Fst=mean(Fst)) %>% ungroup


p<-ggplot(data.plot,aes(x=habitat,y=Fst,fill=habitat))+geom_violin()+geom_point(data=data.obs)+theme_bw()
ggsave(p,file="Fst_bootstraps_30m.pdf")

pvalues<-(data.plot %>% filter(habitat=="between") %>% select(Fst)) -( data.plot %>% filter(habitat=="within") %>% select(Fst))
pvalue<-1-sum(pvalues>0)/100 #0.031

}
