
#---
#title: "pairwise identity from VCF file function"
#author: "Tal Dahan-Meir"
#date: "18/08/2021"
#---


### Ammiad summarize identity from vcf (ammiad_summarize_vcf) ###
### This function gets a vcf and a map (with the coordinates of the sampling points), and retrieves a summarizing table into output_csv_path ###
##########################################

ammiad_summarize_vcf <- function (vcf_path = NULL,
								 map_csv_path = NULL,
								 output_csv_path = NULL,
								 identity_threshold = 0.981)
{

if (is.null(vcf_path)) { stop("please enter path to samples .vcf file")}
if (is.null(map_csv_path)) { stop("please enter path map .csv file")}
if (is.null(output_csv_path)) { stop("please define path for output table (.csv)")}

### dependencies ###
require(igraph)
require(data.table)


### load and clean vcf ###
	
cat("loading files..\n")
vcf = fread(vcf_path)
ammiad_map = read.csv(map_csv_path, stringsAsFactors = FALSE, row.names="name")
vcf_tab_clean = vcf[,10:ncol(vcf)]
vcf_tab_clean_0_3 = t ( apply (vcf_tab_clean,1,function(x) { as.character(substring(x,1,3)) }) )
vcf_tab_clean_0_3 = replace(vcf_tab_clean_0_3,vcf_tab_clean_0_3 == "./.","none")
colnames(vcf_tab_clean_0_3) = colnames(vcf_tab_clean)
samples = colnames(vcf_tab_clean_0_3)

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

                 
### assign names to distinct genotype groups ###
                 
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
                 

### get year of plant collection, excluding controls ###
                 
cat("generating final table..\n")
year_vec = ifelse(grepl("^TTD",samples),"TTD",ifelse(grepl("Zavitan",samples),gsub("_Zavitan","",samples),substring(samples,1,4)))

                 
### get Position (A/B/C/D or TTD or Zavitan for controls) ###
                 
position_vec = ifelse(grepl("^TTD",samples),"TTD",ifelse(grepl("Zavitan",samples),"Zavitan",substring(samples,5)))


### get the group of the sample ###
                 
group_vec = c()
for (i in 1:length(samples)){
for (j in 1:length(groups)) {
if (sum(samples[i] %in% groups[[j]])>0) { group_vec=c(group_vec,names(groups)[j]); next;}
}}

                 
### get the size of groups (how many inidividulas in that group) ###

group_sizes_vec = unlist(lapply(groups,length))[group_vec]

                 
### get Coverage (sum of covered positions, from “covered_sum”) ###
                 
coverage_vec = covered_sum[samples]

                                 
### get Heterozygosity fraction (from “homoz_frac”) ###
                 
homozVec = homoz_frac[samples]

                 
### get Alternative allele count (from alt_frac) ###
                 
alt_vec = alt_frac[samples]

                 
### generate a unified table ###
                 
                 
samples_df = data.frame(Sample = samples,
						Year = year_vec,
					  Position = position_vec,
						IGG = group_vec,
						Group_size = group_sizes_vec,
						Coverage = coverage_vec,
						Heterozygosity = homozVec,
						Alternative = alt_vec,row.names=samples)


### get metadata from map.csv file ###
                 
ammiad_map_f = ammiad_map[rownames(ammiad_map) %in% samples_df$Position,]
missing_points = setdiff(rownames(ammiad_map),rownames(ammiad_map_f))
cat(paste0("missing samples for these ",length(missing_points)," points:\n"))
cat(paste0(missing_points,"\n"))
samples_df[,c("Longitude","Latitude","Height","Habitat","Habitat_order","Transect_order","dist")] = NA
matches = match(samples_df$Position,rownames(ammiad_map_f))
valid_matches = !is.na(match(samples_df$Position,rownames(ammiad_map_f)))
samples_df[valid_matches,c("Longitude","Latitude","Height","Habitat","Habitat_order","Transect_order","dist")] = ammiad_map_f[matches[valid_matches],c("lon","lat","ele","habitat","order_habitat","order_transect","dist")]

                 
### write into csv ###
                 
cat(paste0("\nwriting output file: ", output_csv_path,"\n"))
write.csv(samples_df, file = output_csv_path,row.names=FALSE)
}
                 
##################################################################################
##################################################################################


### assign colors to DGGs (ammiad_colorize_by_groups) ###
### This function gets a summarized table of samples and update it according to a colors table ###
################################
                 
ammiad_colorize_by_groups <- function (summarized_table_path = NULL,
																			colors_table_path = NULL,
																		  output_csv_path = NULL)
{
if (is.null(summarized_table_path)) {stop ("define path to summarized .vcf !")}
if (is.null(colors_table_path)) {stop ("define path to .vcf specifying colors !")}
if (is.null(output_csv_path)) {stop ("define path for output .csv file !")}

output_tab = read.csv(summarized_table_path)
colors_tab = read.csv(colors_table_path)
output_tab$color = colors_tab$color [match(output_tab$IGG,colors_tab$IGG)]
write.csv(output_tab, file = output_csv_path)
cat(paste0("added ",length(table(colors_tab$color))," different colors to ", nrow(output_tab), " samples\n"))
cat(paste0("saved at: ", output_csv_path))
}

##################################################################################
##################################################################################
                 
                 
                 
                 
                 
