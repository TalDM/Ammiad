#---
#title: "principal component analysis using SNPRelate dirctly from VCF file"
#author: "Tal Dahan-Meir"
#date: "25/08/2021"
#---

### run principal component analysis on vcf file using SNPRelate library ###
library(SNPRelate)


vcf.fn="ammiad_845.vcf"
snpgdsVCF2GDS(vcf.fn, "ccm.gds",  method="biallelic.only")
genofile <- openfn.gds("ccm.gds")
ccm_pca<-snpgdsPCA(genofile)


### load output file to get DGG colors ###
data=read.csv("data/output_with_colors.csv")

### habitat colors vector ###
colors=c("#2a7fff","#f7a8b8","#ff5555","#ff9955","#5fd35f","#ffe680","#aa0000","#9955ff")

tab_identity=data.frame(sample.id=data$Sample, pop=factor(data$Habitat_North_divided), igg_col=factor(data$color), EV1=ccm_pca$eigenvec[,1],EV2=ccm_pca$eigenvec[,2], stringsAsFactors=FALSE)
tab_identity=tab_identity[!grepl("Zavitan",tab_identity$sample.id),]

### plot pca coloring by habitat ###

pdf("20210830_pc12_snprelate_by_habitat.pdf",height=8,width=8)
plot(tab_identity$EV1, tab_identity$EV2, col=colors[as.integer(tab_identity$pop)],xlab="PC1(18.14% explained var.)",ylab="PC2(8.02% explained var.)",pch=20)
legend("bottomright", legend=levels(tab_identity$pop), pch=19, col=colors)
dev.off()


### plot pca coloring by DGG ###

pdf("20210830_pc12_snprelate_by_igg.pdf",height=8,width=8)
plot(tab_identity$EV1, tab_identity$EV2, col=as.character(tab_identity$igg_col),xlab="PC1(18.14% explained var.)",ylab="PC2(8.02% explained var.)",pch=20)
dev.off()
