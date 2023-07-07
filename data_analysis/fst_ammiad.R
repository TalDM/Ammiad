#============================================================================
bjobs -u all -p -q my_queue
#==================================Fst-Habitat===============================
#============================================================================
#Fst between vs within
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
#https://stackoverflow.com/questions/639695/how-to-convert-Latitude-or-Longitude-to-meters
    R <- 6378.137                                # radius of earth in Km
    dLat <- (lat2-lat1)*pi/180
    dLon <- (lon2-lon1)*pi/180
    a <- sin((dLat/2))^2 + cos(lat1*pi/180)*cos(lat2*pi/180)*(sin(dLon/2))^2
    c <- 2 * atan2(sqrt(a), sqrt(1-a))
    d <- R * c
    return (d * 1000)                            # distance in meters
}

par(mfrow=c(1,1))

mydata<-read.csv("../data/dataset_S1.csv",header=TRUE)
mydata<-mydata[mydata$Position!="Zavitan",]
mydata<-mydata[!is.na(mydata$Longitude),]
mydata<-mydata[!is.na(mydata$Latitude),]

system("mkdir indiv_pops")
for (xHabitat in unique(mydata$Habitat)){
write.table(mydata$targetsample[mydata$Habitat==xHabitat],file=paste0("indiv_pops/",xHabitat,".txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)
}

#partition sites in Transect_orders with diameter of 30m (radius 15m) within each Habitat
iTransect_order<-0
distance_radius<-30
for (xHabitat in unique(mydata$Habitat)){
    mydata_t<-mydata %>% filter(Habitat==xHabitat,Transect_order==0)
    mins<-mydata %>% filter(Habitat==xHabitat) %>% summarise(minlat=min(Latitude),minlon=min(Longitude)) %>% as.numeric
    maxs<-mydata %>% filter(Habitat==xHabitat) %>% summarise(maxlat=max(Latitude),maxlon=max(Longitude)) %>% as.numeric
    if ((maxs[1]-mins[1])>(maxs[2]-mins[2])){ lat_is_wider<-TRUE } else { lat_is_wider<-FALSE }
        while(sum(mydata_t$Transect_order==0)>0)
        {
        iTransect_order<-iTransect_order+1
        mins<-mydata_t %>% filter(Transect_order==0) %>% summarise(minlat=min(Latitude),minlon=min(Longitude)) %>% as.numeric
        maxs<-mydata_t %>% filter(Transect_order==0) %>% summarise(maxlat=max(Latitude),maxlon=max(Longitude)) %>% as.numeric
        #latlon2m(mins[2],mins[1],maxs[2],maxs[1]) #113m
        #choose point1: most extreme. to decide between left and bottom, I use overall distribution to understand which axis is more relevant.
        if (lat_is_wider){point1<-mydata_t %>% filter(Transect_order==0,Latitude==mins[1]) %>% select(Latitude,Longitude) %>% head(1) %>% as.numeric} else
        {point1<-mydata_t %>% filter(Transect_order==0,Longitude==mins[2]) %>% select(Latitude,Longitude) %>% head(1) %>% as.numeric}
        dist_t<-apply(mydata_t %>% select(Latitude,Longitude),MARGIN=1,FUN=function(x) latlon2m(point1[2],point1[1],x[2],x[1]))
        mydata_t$Transect_order[dist_t<distance_radius & mydata_t$Transect_order==0]<-iTransect_order
        }
        mydata$Transect_order[mydata$Habitat==xHabitat]<-mydata_t$Transect_order
}


#plot map of parcels and Habitats
pdf(paste0("Transect_orders",distance_radius,"m.pdf"))
par(mfrow=c(1,2))
plot(mydata$Longitude,mydata$Latitude,col=as.factor(mydata$Habitat),pch=19,cex=0.5,xlab="Longitude",ylab="Latitude",main="Habitat")
plot(mydata$Longitude,mydata$Latitude,col=as.factor(mydata$Transect_order),pch=19,cex=0.5,xlab="Longitude",ylab="Latitude",main=paste("Transect_orders ",distance_radius,"m"))
dev.off()

#read vcf and calculate pairwise-Fst
vcf<-read.vcfR("vcf/ammiad_only_filtered_hets_845.vcf")
vcf.loci<-vcfR2loci(vcf)
vcf.genind<-vcfR2genind(vcf)
vcf.labels<-VCFlabels("vcf/ammiad_only_filtered_hets_845.vcf")
pops<-vcf.labels
pops[!vcf.labels %in% mydata$targetsample]<-"other"
pops<-left_join(data.frame(targetsample=pops),mydata)
pops$Transect_order[is.na(pops$Transect_order)]<-max(pops$Transect_order,na.rm=TRUE)+1
snpgdsVCF2GDS("ammiad_only_filtered_hets_845.vcf", "ammiad_only_filtered_hets_845.gds")
ammiad.gds<-snpgdsOpen("ammiad_only_filtered_hets_845.gds")
Transect_orders.fst<-snpgdsFst(ammiad.gds,as.factor(pops$Transect_order),"W&H02",verbose=TRUE)
matrixFst<-Transect_orders.fst$Beta

vcf.genind@pop<-as.factor(pops$Transect_order)
mathierfstat<-genind2hierfstat(vcf.genind)
matFst <- pairwise.WCfst(mathierfstat,diploid=TRUE)


#Partition in between and within Habitat comparisons
mydata<-mydata %>% group_by(Transect_order) %>% summarise(long_Transect_order=mean(Longitude),lat_Transect_order=mean(Latitude)) %>% left_join(mydata)
mat_v<-c()
samehab<-c()
dist.Transect_orders_v<-c()
for (i in 1:(nrow(matFst)-1)){
for (j in (i+1):nrow(matFst)){
mat_v<-c(mat_v,matFst[i,j])
habi<-mydata$Habitat[mydata$Transect_order==rownames(matFst)[i]][1]
habj<-mydata$Habitat[mydata$Transect_order==rownames(matFst)[j]][1]
if (habi==habj){ samehab<-c(samehab,"within") } else { samehab<-c(samehab,"between") }
long_i<-mydata$long_Transect_order[mydata$Transect_order==rownames(matFst)[i]][1]
lat_i<-mydata$lat_Transect_order[mydata$Transect_order==rownames(matFst)[i]][1]
long_j<-mydata$long_Transect_order[mydata$Transect_order==rownames(matFst)[j]][1]
lat_j<-mydata$lat_Transect_order[mydata$Transect_order==rownames(matFst)[j]][1]
dist.Transect_orders_v<-c(dist.Transect_orders_v,latlon2m(long_i,lat_i,long_j,lat_j))
}
}

Fst2Nm<-function(x) (1-x)/(4*x)

#add as second axis m estimated by Fst: note that this m should be considered a ballpark estimate and not an accurate one,
#see for example: Whitlock et al.,1999, Heredity
#FST!=1/(4Nm+1) #->#Nm=(1-Fst)/(4*Fst)

data.plot<-data.frame(distance=dist.Transect_orders_v,Habitat=as.factor(samehab),Fst=mat_v,Nm=(1-mat_v)/(4*mat_v))
p<-ggplot(data.plot,aes(x=distance,y=Fst,colour=Habitat))+geom_point()+scale_x_log10()+stat_smooth(method=glm,method.args=list(family=binomial))+
scale_y_continuous("Fst",sec.axis = sec_axis(~ (1 - .)/(4*.) , name = "Nm"))+theme_bw()
ggsave(p,file="Fstdecay_30m_binomial.pdf")
p<-ggplot(data.plot,aes(x=distance,y=Fst,colour=Habitat))+geom_point()+scale_x_log10()+stat_smooth(method=loess,method.args=list(degree=1))+theme_bw()
#+scale_y_continuous("Fst",sec.axis = sec_axis(~ (1 - .)/(4*.) , name = "Nm"))
ggsave(p,file="Fstdecay_30m_loess.pdf")

data.obs<-data.plot[data.plot$distance<50,]
data.obs<-data.obs %>% group_by(Habitat) %>% summarise(Fst=mean(Fst)) %>% ungroup
data.obs$bs<-0

samplei<-c();samplej<-c()
for (i in 1:(nrow(matFst)-1)){
for (j in (i+1):nrow(matFst)){
samplei<-c(samplei,i);samplej<-c(samplej,j)}}

bs.n<-1000
bs.v<-c()
mydata<-mydata %>% group_by(Transect_order) %>% summarise(long_Transect_order=mean(Longitude),lat_Transect_order=mean(Latitude)) %>% left_join(mydata)
mat_v<-c()
samehab<-c()
dist.Transect_orders_v<-c()
for (bs.iteration in 1:bs.n){
newsample<-sample(1:length(samplei),replace=TRUE)
for (i in 1:length(samplei)){
i.bs<-samplei[newsample[i]]
j.bs<-samplej[newsample[i]]
mat_v<-c(mat_v,matFst[i.bs,j.bs])
habi<-mydata$Habitat[mydata$Transect_order==rownames(matFst)[i.bs]][1]
habj<-mydata$Habitat[mydata$Transect_order==rownames(matFst)[j.bs]][1]
if (habi==habj){ samehab<-c(samehab,"within") } else { samehab<-c(samehab,"between") }
long_i<-mydata$long_Transect_order[mydata$Transect_order==rownames(matFst)[i.bs]][1]
lat_i<-mydata$lat_Transect_order[mydata$Transect_order==rownames(matFst)[i.bs]][1]
long_j<-mydata$long_Transect_order[mydata$Transect_order==rownames(matFst)[j.bs]][1]
lat_j<-mydata$lat_Transect_order[mydata$Transect_order==rownames(matFst)[j.bs]][1]
dist.Transect_orders_v<-c(dist.Transect_orders_v,latlon2m(long_i,lat_i,long_j,lat_j))
bs.v<-c(bs.v,bs.iteration)
}}
data.bs<-data.frame(distance=dist.Transect_orders_v,Habitat=as.factor(samehab),Fst=mat_v,bs=bs.v)
data.bs.t<-data.bs[data.bs$distance<50,]
data.plot<-data.bs.t %>% group_by(bs,Habitat) %>% summarise(Fst=mean(Fst)) %>% ungroup


p<-ggplot(data.plot,aes(x=Habitat,y=Fst,fill=Habitat))+geom_violin()+geom_point(data=data.obs)+theme_bw()
ggsave(p,file="Fst_bootstraps_30m.pdf")

pvalues<-(data.plot %>% filter(Habitat=="between") %>% select(Fst)) -( data.plot %>% filter(Habitat=="within") %>% select(Fst))
pvalue<-1-sum(pvalues>0)/100 #0.031

############SUPPLEMENTARY CODE###################
#----------------Supplementary code to compute Fst from data in snpgds format-------------#

mydata<-read.csv("../data/dataset_S1.csv",header=TRUE)
mydata<-mydata[mydata$Position!="Zavitan",]
mydata<-mydata[!is.na(mydata$Longitude),]
mydata<-mydata[!is.na(mydata$Latitude),]

vcf<-read.vcfR("vcf/ammiad_only_filtered_hets_845.vcf")
vcf.loci<-vcfR2loci(vcf)
vcf.genind<-vcfR2genind(vcf)
vcf.labels<-VCFlabels("vcf/ammiad_only_filtered_hets_845.vcf")
pops<-vcf.labels
pops[!vcf.labels %in% mydata$targetsample]<-"other"
pops<-left_join(data.frame(targetsample=pops),mydata)
pops<-pops[!is.na(pops$pop),]
dist.hamming(vcf)
d<-pegas::dist.hamming(vcf.loci)
save.image("dist.hamming.RData")
d<-pegas::dist.dna(vcf.loci)
pops$Transect_order[is.na(pops$Transect_order)]<-max(pops$Transect_order,na.rm=TRUE)+1
snpgdsVCF2GDS("ammiad_only_filtered_hets_845.vcf", "ammiad_only_filtered_hets_845.gds")
ammiad.gds<-snpgdsOpen("ammiad_only_filtered_hets_845.gds")
Transect_orders.fst<-snpgdsFst(ammiad.gds,as.factor(pops$Transect_order),"W&H02",verbose=TRUE)
matrixFst<-Transect_orders.fst$Beta


#----------------BASH and R code to calculate pi with vcftools------------#

#in bash
cat vcf/ammiad_only_filtered_hets_845.vcf | grep CHROM | awk '{for (i=10;i<=NF;i++){print $i}}' > vcf/samples.txt
NLINES=$( cat vcf/samples.txt | wc -l )
for i in `seq 1 $( expr $NLINES - 1 )`;do
for j in `seq $i $NLINES`;do
cat vcf/samples.txt | head -${i} | tail -1 > temp.txt
cat vcf/samples.txt | head -${j} | tail -1 >> temp.txt
vcftools --vcf vcf/ammiad_only_filtered_hets_845.vcf --keep temp.txt --site-pi --out res
pi=$( cat res.sites.pi | grep -v CHROM awk 'BEGIN{tot=0}{tot=tot+$3}END{print tot/NR}' )
echo $i $j $pi >> results.txt
done
done

#in R                      
xsamples<-unname(unlist(read.table("vcf/samples.txt")))
mat<-matrix(rep(0,845*845),nrow=845)
for ( i in 1:844){
print(i)
xpi<-read.table(paste0("results.pi/results.",i,".txt"))
mat[xpi$V2,i]<-xpi$V3
rownames(mat)<-xsamples
colnames(mat)<-xsamples
}
write.table(mat,file="fstmat.tsv",quote=FALSE,sep='\t')

