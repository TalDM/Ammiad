
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

#sort(covered_frac)
#sort(covered_sum)


### calculate pairwise identity ###

vcf_tab_clean_0_3 = replace(vcf_tab_clean_0_3,vcf_tab_clean_0_3 == "none",NA)
similarity.matrix<-apply(vcf_tab_clean_0_3,2,function(x)colSums(x==vcf_tab_clean_0_3, na.rm=TRUE))
Ncomparisons_matrix=apply(vcf_tab_clean_0_3,2,function(x)colSums(!is.na(x==vcf_tab_clean_0_3)))
identity.matrix=similarity.matrix/Ncomparisons_matrix

                          
### plot pairwise identity distribution of Ammiad (removing Zavitan controls) ###
                          
no_z=identity.matrix[!grepl(as.character("Zavitan"),rownames(identity.matrix)),!grepl("Zavitan",colnames(identity.matrix))]
#hist(no_z[lower.tri(no_z)],breaks=60,col="#a1d99b",border="#74C476")
   
                          
### for every plant, find all plants that has at least X (identity_threshold) identity to it ###
                          
require(igraph)
require(data.table)
                          
identity_threshold=0.981
my_list=list()
for (i in 1:ncol(identity.matrix))
{
relevant_vector = rownames(identity.matrix)[identity.matrix[,i] >= identity_threshold]
my_list[[i]]=relevant_vector
}
names(my_list) = colnames(identity.matrix)



overlap <- function(u, v) length(intersect(u, v)) / min(length(u), length(v)) > 0.8
adj <- sapply(my_list, function(u) sapply(my_list, overlap, u)) + 0
g <- graph_from_adjacency_matrix(adj)
memb <- components(g)$membership
s <- split(my_list, memb)
groups <- lapply(s, function(x) unique(unlist(x)))
cat(paste0("getting info from ",length(groups), " groups\n"))
#length(groups)

                 
### assign names to groups ###
                 
names(groups) <- paste("G",1:length(groups),sep="")
lengths=as.vector(sapply(groups, length))
#sort(lengths)

                 
### make table with groups and levels ### 
                 
m1=max(lengths(groups))
d1=as.data.frame(do.call(rbind,lapply(groups,'length<-',m1)),stringsAsFactors=FALSE)
names(d1)=paste0("Level",seq_along(d1))
csv=d1
 group_sizes = apply(csv,2,function(x){sum(!is.na(x))})

groups_vector = c()
for (i in 1:length(group_sizes)){
groups_vector = c(groups_vector,rep(i,group_sizes[i]))}
groups_list = apply(csv,2,function(x){as.character(x)[!is.na(x)]})
names(groups_vector) = unlist(groups_list,use.names=FALSE)
