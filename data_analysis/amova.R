#Amova

mydata<-read.table("map_file.txt",header=TRUE)
mydata<-mydata[mydata$pop=="ammiad" & ! is.na(mydata$pop),]
mydata<-mydata[!is.na(mydata$longitude),]
mydata<-mydata[!is.na(mydata$latitude),]

vcf<-read.vcfR("vcf/ammiad_only_filtered_hets_845.vcf")
vcf.loci<-vcfR2loci(vcf)
vcf.genind<-vcfR2genind(vcf)
vcf.labels<-VCFlabels("vcf/ammiad_only_filtered_hets_845.vcf")
pops<-vcf.labels
pops[!vcf.labels %in% mydata$targetsample]<-"other"
pops<-left_join(data.frame(targetsample=pops),mydata)
pops<-pops[!is.na(pops$pop),]
d<-pegas::dist.hamming(vcf.loci)
save.image("dist.hamming.RData")

d<-pegas::dist.dna(vcf.loci)


pops$transect[is.na(pops$transect)]<-max(pops$transect,na.rm=TRUE)+1
snpgdsVCF2GDS("ammiad_only_filtered_hets_845.vcf", "ammiad_only_filtered_hets_845.gds")
ammiad.gds<-snpgdsOpen("ammiad_only_filtered_hets_845.gds")
transects.fst<-snpgdsFst(ammiad.gds,as.factor(pops$transect),"W&H02",verbose=TRUE)
matrixFst<-transects.fst$Beta


#calculate pi with vcftools

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

cat vcf/ammiad_only_filtered_hets_845.vcf | grep CHROM | awk '{for (i=10;i<=NF;i++){print $i}}' > vcf/samples.txt
NLINES=$( cat vcf/samples.txt | wc -l )
for i in `seq 1 $( expr $NLINES - 1 )`;do
bsub -q new-medium -J a${i} -e ~/mylogs/a.e%J -o ~/mylogs/a.o%J -R rusage[mem=1000] ./pairwise.pi.parallel.sh $i
done

#!/bin/bash
module vcftools
i=$1
NLINES=$( cat vcf/samples.txt | wc -l )
for j in `seq $i $NLINES`;do
cat vcf/samples.txt | head -${i} | tail -1 > temp${i}.txt
cat vcf/samples.txt | head -${j} | tail -1 >> temp.${i}.txt
vcftools --vcf vcf/ammiad_only_filtered_hets_845.vcf --keep temp.${i}.txt --site-pi --out res.${i}
pi=$( cat res.${i}.sites.pi | grep -v CHROM | sed 's/-nan/0/g' | awk 'BEGIN{tot=0}{tot=tot+$3}END{print tot/NR}' )
echo $i $j $pi >> results.${i}.txt



vcftools vcf/ammiad_only_filtered_hets_845.vcf

mat<-matrix(rep(0,845*845),nrow=845)
for ( i in 1:844){
print(i)
xpi<-read.table(paste0("results.pi/results.",i,".txt"))
mat[xpi$V2,i]<-xpi$V3
rownames(mat)<-xsamples
colnames(mat)<-xsamples

pegas::amova(formula=d.mat ~ habitat,data=pops)
pops$habitat<-as.factor(pops$habitat)
row.names(mat.df)<-xsamples
col.names(mat.df)<-xsamples
mat<-mat[rownames(mat) %in% pops$targetsample,]
mat<-mat[,colnames(mat) %in% pops$targetsample]
d.mat<-as.dist(mat)
amova.onlyhabitat<-pegas::amova(formula=d.mat ~ habitat,data=pops)

pops$year<-as.numeric(pops$year)
pops$year<-as.factor(pops$year)
pops$h.year<-as.factor(pops$h.year)

pops$h.year<-paste(pops$habitat,pops$year,sep="_")
amova.onlyhabitat<-pegas::amova(formula=d.mat ~ habitat,data=pops)
amova.habitatyear2<-pegas::amova(formula=d.mat ~ habitat/h.year,data=pops)
amova.habitatyear2

#         Analysis of Molecular Variance
# 
# Call: pegas::amova(formula = d.mat ~ habitat/h.year, data = pops)
# 
#               SSD         MSD  df
# habitat 1.0733240 0.153332005   7
# h.year  0.2354624 0.003679100  64
# Error   2.8637634 0.003768110 760
# Total   4.1725498 0.005021119 831
# 
# Variance components:
#              sigma2 P.value
# habitat  1.4864e-03  0.0000
# h.year  -7.7155e-06  0.8492
# Error    3.7681e-03
# 
# Phi-statistics:
# habitat.in.GLOBAL (Phi_CT)  h.year.in.GLOBAL (Phi_ST)
#                0.283294291                0.281823778
# h.year.in.habitat (Phi_SC)
#               -0.002051767
# 
# Variance coefficients:
#         a         b         c
#  11.53661  11.33873 100.68166

# Call: pegas::amova(formula = d.mat ~ habitat/h.year, data = pops)
# 
#               SSD         MSD  df
# habitat 1.0733240 0.153332005   7
# h.year  0.2354624 0.003679100  64
# Error   2.8637634 0.003768110 760
# Total   4.1725498 0.005021119 831
# 
# Variance components:
#              sigma2 P.value
# habitat  1.4864e-03  0.0000
# h.year            0  0.8492
# Error    3.7681e-03
# 
# Phi-statistics:
# habitat.in.GLOBAL (Phi_CT)  h.year.in.GLOBAL (Phi_ST)
#                0.283294291                0.281823778
# h.year.in.habitat (Phi_SC)
#               -0.002051767
# 
}
