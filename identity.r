
---
title: "pairwise identity from VCF file"
author: "Tal Dahan-Meir"
date: "18/08/2021"
---

library(data.table)


### load and clean vcf ###

vcf=fread("ammiad_845.vcf")
vcf_tab_clean = vcf[,10:ncol(vcf)]
vcf_tab_clean_0_3 = t ( apply (vcf_tab_clean,1,function(x) { as.character(substring(x,1,3)) }) )
vcf_tab_clean_0_3 = replace(vcf_tab_clean_0_3,vcf_tab_clean_0_3 == "./.","none")
colnames(vcf_tab_clean_0_3) = colnames(vcf_tab_clean)


### get heterozygosity ###

homoz_het = vcf_tab_clean_0_3
homoz_het[,]="het"
homoz_het[grepl("0/0|1/1|2/2|3/3",vcf_tab_clean_0_3)] = "hom"
homoz_het[grepl("none",vcf_tab_clean_0_3)] = "none"
homoz_frac = apply(homoz_het,2,function(x){100*(sum(x=="het")/length(x))})
names(homoz_frac) = colnames(homoz_het)
sort(homoz_frac)


### get alternative allele fraction ###

alternative=vcf_tab_clean_0_3
alternative[,]="alt"
alternative[grepl("0/0",vcf_tab_clean_0_3)] = "ref"
alt_frac = apply(alternative,2,function(x){100*(sum(x=="alt")/length(x))})
sort(alt_frac)


### get covered positions per plant ###

covered = vcf_tab_clean_0_3
covered[,]="cov"
covered[grepl("none",vcf_tab_clean_0_3)] = "not"
covered_frac = apply(covered,2,function(x){100*(sum(x=="cov")/length(x))})
covered_sum = apply(covered,2,function(x){sum(x=="cov")})
names(covered_frac) = colnames(covered)

sort(covered_frac)
sort(covered_sum)


### calculate pairwise identity ###

vcf_tab_clean_0_3 = replace(vcf_tab_clean_0_3,vcf_tab_clean_0_3 == "none",NA)
similarity.matrix<-apply(vcf_tab_clean_0_3,2,function(x)colSums(x==vcf_tab_clean_0_3, na.rm=TRUE))
Ncomparisons_matrix=apply(vcf_tab_clean_0_3,2,function(x)colSums(!is.na(x==vcf_tab_clean_0_3)))
identity.matrix=similarity.matrix/Ncomparisons_matrix

