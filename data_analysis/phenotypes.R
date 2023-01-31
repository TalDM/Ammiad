#================================================================================
#=====================ANALYSES WITH DISTINCTIVE PHENOTYPES=======================
#================================================================================
###DATA PREPARATION####
#impute data for missing phenotypes
library(tidyverse)
palette_from_tal<-c("East slope"="#2a7fff","Karst"="#f7a8b8","Middle north"="#ff5555","Narrow valley"="#ff9955","Plateau"="#5fd35f","South slope"="#ffe680","Upper north"="#aa0000","Valley"="#9955ff")
x.pheno<-read.csv("data/pheno_20_21_data.csv")
x.samples<-read.csv("data/845_output_per981.csv")
x.seeds<-read.csv("data/mean_median_var_seed_size.csv")

x.seeds<-x.seeds[grep("ZAVITAN",x.seeds$X,invert=TRUE),]
x.seeds<-x.seeds[grep("N2",x.seeds$X,invert=TRUE),]

x.seeds.str<-strsplit(as.character(x.seeds$X),split="G")
x.blocks<-sapply(x.seeds.str,function(x) x[1]) %>% gsub("BLOCK","",.) %>% as.numeric
x.gg<-sapply(x.seeds.str,function(x) x[2]) %>% paste0("G",.)
x.seeds<-data.frame(old_dgg=x.gg,block=x.blocks,MeanArea=x.seeds$Area,MeanLength=x.seeds$Length,MeanWidth=x.seeds$Width)

new_dgg_abundance<-table(x.samples$new_dgg,x.samples$Habitat_North_divided)/sum(table(x.samples$new_dgg,x.samples$Habitat_North_divided))
new_dgg_abundance<-as.data.frame(new_dgg_abundance) %>% pivot_wider(names_from="Var2",values_from="Freq")
names(new_dgg_abundance)[1]<-"new_dgg"
names(new_dgg_abundance)[2:length(names(new_dgg_abundance))]<-paste0("ab_",names(new_dgg_abundance)[2:length(names(new_dgg_abundance))])

#=============imputation taking "unique phenotypes"
mythreshold<-0.3 #0.3
maxtot<-function(x) { y<-x/sum(x,na.rm=T); y[y>mythreshold]}

res<-data.frame(new_dgg=c(),habitats=c())
for ( i in unique(x.samples$new_dgg) ){
my.temp<-x.samples %>% filter(new_dgg==i) %>% dplyr::select(Habitat_North_divided) %>% table %>% maxtot
my.lat<-x.samples %>% filter(new_dgg==i) %>% summarise(lat=median(Latitude))
my.long<-x.samples %>% filter(new_dgg==i) %>% summarise(lat=median(Longitude))
if (length(my.temp)>0)
  {
  tempres<-data.frame(new_dgg=i,names(my.temp),my.long,my.lat)
#  latx<-c();longx<-c()
#  for (ix in names(my.temp)){
#    x.samples %>% filter(new_dgg==i) %>% select(Longitude)
  res<-rbind(res,tempres)
  }
}
names(res)<-c("new_dgg","habitats","longitude","latitude")

#=============imputation taking  phenotypes using an abundance threshold to exclude singletons and rare by chance which might not be representative
x.samples.counts<-table(x.samples$new_dgg,x.samples$Habitat_North_divided) %>% as.data.frame %>% group_by(Var2) %>% summarise(Var1=Var1,Freq=Freq,tot=sum(Freq),prop=Freq/tot) %>% filter(Freq>0) %>% ungroup
x.samples.counts<-table(x.samples$new_dgg,x.samples$Habitat_North_divided) %>% as.data.frame %>% group_by(Var2) %>% summarise(Var1=Var1,Freq=Freq,tot=sum(Freq),prop=Freq/tot) %>% filter(Freq>0) %>% ungroup
habitat1xdgg<-table(x.samples$new_dgg,x.samples$Habitat_North_divided) %>% as.data.frame %>% group_by(Var2) %>% summarise(Var1=Var1,Freq=Freq,tot=sum(Freq),prop=Freq/tot) %>% ungroup %>% group_by(Var1) %>% summarise(Var1=Var1,Var2=Var2,Freq=Freq,tot=sum(Freq),prop=Freq/tot,mymax=max(prop)) %>% filter(mymax==prop)
habitat1xdgg<-habitat1xdgg[!duplicated(habitat1xdgg$Var1),1:2] #G10 and G65 cannot be assigned univocally, I remove at "random".
names(x.samples.counts)[1:2]<-c("habitats","new_dgg")
names(habitat1xdgg)<-c("habitats","new_dgg")
x.samples.counts<-x.samples.counts %>% select("new_dgg","habitats")
res<-x.samples.counts %>% inner_join(res)

x.pheno<-x.pheno[grep("trans",x.pheno$germination_date,invert=TRUE),]
x.pheno<-x.pheno %>% right_join(res)
x.pheno<-x.pheno %>% left_join(x.seeds)
x.pheno<-x.pheno %>% left_join(new_dgg_abundance)
x.ldata<-x.pheno[,(which(names(x.pheno)=="new_dgg")):length(names(x.pheno))]
allvars<-c("coleptile_color","erect_1_5","pigment_1_3","main_length","main_flag_leaf_width",
"main_number_of_spikes","number_of_tillers","spike_angle_1_3","weight_10_seeds",
"MeanLength","MeanWidth","anthesis_time")

strin2date<-function(xx) sapply(as.character(xx), function(x)
{
  if (!is.na(x))
  {
     y<-strsplit(x,split="");
     year=as.numeric(paste0(y[[1]][1:4],collapse=""));
     month=as.numeric(paste0(y[[1]][5:6],collapse=""));
     day=as.numeric(paste0(y[[1]][7:8],collapse=""));
     #c(year,month,day)
     365*year+30*month+day
  } else return(NA)
}
)

germination_date<-unname(unlist(strin2date(x.ldata$germination_date)))
anthesis_date<-unname(unlist(strin2date(x.ldata$anthesis_date)))
anthesis_time=anthesis_date-germination_date
germination_date<-germination_date-min(germination_date,na.rm=T)
x.ldata$germination_date<-germination_date
x.ldata$anthesis_time<-anthesis_time
x.ldata$habitats<-as.factor(x.ldata$habitats)
x.ldata$coleptile_color[x.ldata$coleptile_color==""]<-NA
#x.ldata$coleptile_color<-as.factor(x.ldata$coleptile_color)
x.ldata$block<-as.factor(x.ldata$block)
x.ldata$new_dgg<-as.factor(x.ldata$new_dgg)
x.ldata$latitude<-c(scale(x.ldata$latitude))
x.ldata$longitude<-c(scale(x.ldata$longitude))

currenthabitats<-c("E_facing_slope","Karst","Middle_N_facing_slope","Narrow_valley","Plateau","S_facing_slope","Upper_N_facing_slope","Valley")
x.ldata$habitats<-as.character(x.ldata$habitats)
for (i in 1:length(names(palette_from_tal))){
x.ldata$habitats[x.ldata$habitats==currenthabitats[i]]<-names(palette_from_tal)[i]
}
x.ldata$habitats<-as.factor(x.ldata$habitats)
h2.data<-x.ldata
h2.data<-h2.data[!duplicated(paste(h2.data$new_dgg,h2.data$block,sep="_")),]
x.ldata$coleptile_color<-dplyr::na_if(x.ldata$coleptile_color,"")
h2.data<-h2.data %>% mutate(coleptile_color=as.numeric(as.factor(coleptile_color))-1)
x.ldata$anthesis_date<-NULL
i_vars<-which(!names(x.ldata) %in% c("new_dgg","block","habitats","coleptile_color","longitude","latitude"))

#data for multivariate
x.mdata<-x.ldata[!is.na(x.ldata$anthesis_time) & !is.na(x.ldata$MeanWidth) & !is.na(x.ldata$main_length),]
for (i in c("germination_date","erect_1_5","pigment_1_3","main_length","main_flag_leaf_width","main_number_of_spikes",
"number_of_tillers","spike_angle_1_3","weight_10_seeds","MeanArea","MeanLength","MeanWidth","anthesis_time")){
x.mdata[[i]]<-c(scale(x.mdata[[i]]))
}

#data for regression
x.regdata<-x.ldata
x.regdata$coleptile_color[x.regdata$coleptile_color==""]<-NA
x.regdata$coleptile_color<-as.factor(x.regdata$coleptile_color)
x.regdata$coleptile_color<-as.numeric(x.regdata$coleptile_color)-1

#data for machine learning
y.ldata<-x.ldata
for (i in c("germination_date","erect_1_5","pigment_1_3","main_length","main_flag_leaf_width","main_number_of_spikes",
"number_of_tillers","spike_angle_1_3","weight_10_seeds","MeanArea","MeanLength","MeanWidth","anthesis_time")){
y.ldata[[i]]<-c(scale(y.ldata[[i]]))
}
y.ldata<-y.ldata %>% dplyr::select(!latitude & !longitude & !germination_date &! MeanArea)

x.ldata.clean<-y.ldata %>% select(c("habitats",names(y.ldata)[names(y.ldata) %in% allvars]))
x.ldata.clean<-x.ldata.clean[!is.na(x.ldata.clean$MeanLength) & !is.na(x.ldata.clean$anthesis_time) &!is.na(x.ldata.clean$main_length) & !is.na(x.ldata.clean$erect_1_5), ]
x.ldata.clean2<-x.ldata.clean %>% select(!coleptile_color)
x.ldata.clean3<-x.ldata.clean %>% filter(!is.na(coleptile_color)) # %>% mutate(coleptile_color=as.numeric(coleptile_color)-1)
x.ldata.clean4<-x.ldata.clean %>% mutate(coleptile_color=as.numeric(as.factor(coleptile_color))-1) %>% filter(!is.na(coleptile_color))

allvars.all<-allvars #c(allvars,"germination_date","MeanArea") decided to remove them from analyses, so also from table
allvars.all<-allvars.all[which(allvars.all!="coleptile_color")]
library(brms)
brmsh2<-list()
for (i in c(which(names(h2.data) %in% allvars.all)))
{
  x.temp<-h2.data[,c(1,2,i)]
  names(x.temp)[3]<-"y"
  brms_m1.1 <- brm(
    y ~ 1 + (1 | block)+(1 | new_dgg),
    data = x.temp,
    family = gaussian(),
    chains = 4, cores = 4, iter = 6000, control = list(adapt_delta = 0.97, max_treedepth = 12), seed=1
  )
  save(brms_m1.1, file = paste0("models_brms/brms_m1_",names(x.ldata)[i],".rda"))
  #plot(brms_m1.1)
  #mcmc_plot(brms_m1.1, type = "acf")
  #summary(brms_m1.1)

  v_animal <- (VarCorr(brms_m1.1, summary = FALSE)$new_dgg$sd)^2
  v_r <- (VarCorr(brms_m1.1, summary = FALSE)$residual$sd)^2
  v_block <- (VarCorr(brms_m1.1, summary = FALSE)$block$sd)^2
  h.bwt.1 <- as.mcmc(v_animal / (v_animal + v_r + v_block))
  brmsh2[[names(h2.data)[i]]]<-h.bwt.1
  print(names(h2.data)[i])
  summary(h.bwt.1)
}
save(brmsh2, file = paste0("models_brms/brmsh2.rda"))

i<-which(names(h2.data)==c("coleptile_color"))
x.temp<-h2.data[,c(1,2,i)]
names(x.temp)[3]<-"y"
x.temp$y<-as.numeric(x.temp$y)
x.temp$y <- x.temp$y -1
brms_m1.1 <- brm(
  y ~ 1 + (1 | block)+(1 | new_dgg),
  data = x.temp,
  family = bernoulli(),
  chains = 4, cores = 4, iter = 6000, control = list(adapt_delta = 0.95, max_treedepth = 12)
)
save(brms_m1.1, file = paste0("models_brms/brms_m1_",names(x.ldata)[i],".rda"))
v_animal <- (VarCorr(brms_m1.1, summary = FALSE)$new_dgg$sd)^2
v_r <- summary(brms_m1.1)$fixed$Estimate^2 #double check
v_block <- (VarCorr(brms_m1.1, summary = FALSE)$block$sd)^2
h.bwt.1 <- as.mcmc(v_animal / (v_animal + v_r + v_block))
brmsh2[[names(h2.data)[i]]]<-h.bwt.1
print(names(h2.data)[i])
summary(h.bwt.1)
save(brmsh2, file = paste0("models_brms/brmsh2.rda"))
}

for (i in 1:length(brmsh2)){
res_temp<-c(names(brmsh2[i]),summary(brmsh2[[i]])[[1]][1:2],summary(brmsh2[[i]])[[2]][1],summary(brmsh2[[i]])[[2]][2])
if (i==1){ res<-res_temp } else { res<-rbind(res,res_temp) }
}
res<-data.frame(trait=res[,1],heritability=as.numeric(res[,2]),h2_sd=as.numeric(res[,3]),CI2.5=res[,3],CI97.5=res[,4])
write.table(res, file = paste0("models_brms/res.tab"),quote=FALSE,row.names=FALSE,sep='\t')

load("models_brms/brmsh2.rda")

#density plots
for (i in 1:length(brmsh2)){
res_temp<-data.frame(phenotype=names(brmsh2[i]),heritability=unname(unlist(brmsh2[i])))
if (i==1){ res<-res_temp } else { res<-rbind(res,res_temp) }
}

p<-res%>%
  ggplot(aes(x=heritability, color=phenotype, fill=phenotype)) +
  geom_density(alpha=0.3,size=1)+
  labs(x="heritability")+facet_wrap(~phenotype)+
  theme(strip.background = element_rect(fill = "white"),axis.text.x = element_text(angle = 90, hjust=1), panel.background = element_rect(fill = 'white',color='black'),panel.grid.major = element_line(color = 'grey', linetype = 'longdash',size=0.03),legend.position = "bottom",aspect.ratio=1,strip.text.x = element_text(size = 6))
ggsave(p,file="models_brms/heritability_density.pdf")
p<-res%>% filter(phenotype!="germination_date",phenotype!="MeanArea") %>%
  ggplot(aes(x=heritability, color=phenotype, fill=phenotype)) +
  geom_density(alpha=0.3,size=1)+
  labs(x="heritability")+facet_wrap(~phenotype)+
  theme(strip.background = element_rect(fill = "white"),axis.text.x = element_text(angle = 90, hjust=1), panel.background = element_rect(fill = 'white',color='black'),panel.grid.major = element_line(color = 'grey', linetype = 'longdash',size=0.03),legend.position = "bottom",aspect.ratio=1,strip.text.x = element_text(size = 8))
ggsave(p,file="models_brms/heritability_density_filteredphenotypes.pdf")
#x-axis - genotype, ranked by mean value in block, value on y-axis with blocks as different colors
library(tidytext)
palette_from_tal<-c("East slope"="#2a7fff","Karst"="#f7a8b8","Middle north"="#ff5555","Narrow valley"="#ff9955","Plateau"="#5fd35f","South slope"="#ffe680","Upper north"="#aa0000","Valley"="#9955ff")

p<-h2.data %>% pivot_longer(cols=which(names(h2.data) %in% c(allvars,"coleptile_color")),names_to="phenotype",values_to="estimate") %>% filter(!is.na(estimate)) %>%
ggplot(aes(x=reorder_within(new_dgg, estimate, phenotype), y=estimate, color=habitats)) + geom_point() +scale_color_manual(values=palette_from_tal)+
labs(x="genotype",y="phenotype")+ facet_wrap(~phenotype, scales = "free")+
theme(strip.background = element_rect(fill = "white"),axis.text.x = element_blank(), panel.background = element_rect(fill = 'white',color='black'),panel.grid.major = element_blank(),legend.position = "bottom",aspect.ratio=1,strip.text.x = element_text(size = 8))
ggsave(p,file="models_brms/heritability_raw_talpalette.pdf")
p<-h2.data %>% pivot_longer(cols=which(names(h2.data) %in% c(allvars,"coleptile_color")),names_to="phenotype",values_to="estimate") %>% filter(!is.na(estimate)) %>%
ggplot(aes(x=reorder_within(new_dgg, estimate, phenotype), y=estimate, color=phenotype)) + geom_point() +#scale_color_manual(values=palette_from_tal)+
labs(x="genotype",y="phenotype")+ facet_wrap(~phenotype, scales = "free")+
theme(strip.background = element_rect(fill = "white"),axis.text.x = element_blank(), panel.background = element_rect(fill = 'white',color='black'),panel.grid.major = element_blank(),legend.position = "bottom",aspect.ratio=1,strip.text.x = element_text(size = 8))
ggsave(p,file="models_brms/heritability_raw.pdf")
p<-h2.data %>% pivot_longer(cols=which(names(h2.data) %in% c(allvars,"coleptile_color")),names_to="phenotype",values_to="estimate") %>% filter(!is.na(estimate)) %>%
ggplot(aes(x=reorder_within(new_dgg, estimate, phenotype), y=estimate, color=block)) + geom_point() +#scale_color_manual(values=palette_from_tal)+
labs(x="genotype",y="phenotype")+ facet_wrap(~phenotype, scales = "free")+
theme(strip.background = element_rect(fill = "white"),axis.text.x = element_blank(), panel.background = element_rect(fill = 'white',color='black'),panel.grid.major = element_blank(),legend.position = "bottom",aspect.ratio=1,strip.text.x = element_text(size = 8))
ggsave(p,file="models_brms/heritability_block.pdf")

if (FALSE){
#test habitats in space
brms.habitats<-list()
brms.space<-list()
i_predictors<-which(names(x.ldata) %in% c("new_dgg","block","habitats","longitude","latitude"))
for (i in i_vars){
x.temp<-x.ldata[,c(i,i_predictors)]
names(x.temp)[1]<-"y"
brms_space <- brm(
  y ~ 1 + (1 | block)+(1 | new_dgg)+(1 | habitats)+gp(longitude,latitude),
  data = x.temp,
  family = gaussian(),
  chains = 4, cores = 4, iter = 2000, control = list(adapt_delta = 0.99, max_treedepth = 12), seed=1,
)

v_animal <- (VarCorr(brms_space, summary = FALSE)$new_dgg$sd)^2
v_r <- (VarCorr(brms_space, summary = FALSE)$residual$sd)^2
v_block <- (VarCorr(brms_space, summary = FALSE)$block$sd)^2
v_habitats <- (VarCorr(brms_space, summary = FALSE)$habitats$sd)^2
v_space <-posterior_samples(brms_space)$sdgp_gplongitudelatitude^2
h.habitats <- as.mcmc(v_habitats / (v_habitats+ v_animal + v_r + v_block+v_space))
h.space <- as.mcmc(v_space / (v_habitats+ v_animal + v_r + v_block+v_space))
print(names(x.ldata)[i])
print(summary(brms_space))
mean(h.habitats)
brms.habitats[[names(x.ldata)[i]]]<-h.habitats
brms.space[[names(x.ldata)[i]]]<-h.space
save(brms.habitats, file = paste0("models_brms/brms.habitats.rda"))
save(brms_space, file = paste0("models_brms/brms.space.rda"))
}
}

#---------------multivariate model--------------------
library(brms)
brms_multi <- brm(
  mvbind("erect_1_5","pigment_1_3","main_length","main_flag_leaf_width","main_number_of_spikes","number_of_tillers","spike_angle_1_3","weight_10_seeds",
  "MeanLength","MeanWidth","anthesis_time") ~ (1 | habitats) + (1 | block)+(1 | new_dgg),
  data = x.mdata,
  family = gaussian(),
  control = list(adapt_delta = 0.97, max_treedepth = 12), seed=1,
  chains = 4, cores = 4, iter = 3000
)

#brms_multi <- brm(
#  mvbind("erect_1_5","pigment_1_3","main_length","main_flag_leaf_width","main_number_of_spikes","number_of_tillers","spike_angle_1_3","weight_10_seeds",
#  "MeanLength","MeanWidth","anthesis_time") ~ 1 + (1 | block)+(1 | new_dgg),
#  data = x.mdata,
#  family = gaussian(),
#  chains = 4, cores = 4, iter = 2000
#)
save(brms_multi, file = paste0("models_brms/brms_multi.rda"))

brms_multisp <- brm(
  mvbind("erect_1_5","pigment_1_3","main_length","main_flag_leaf_width","main_number_of_spikes","number_of_tillers","spike_angle_1_3","weight_10_seeds",
  "MeanLength","MeanWidth","anthesis_time") ~ (1 | habitats) + (1 | block)+(1 | new_dgg)+gp(longitude,latitude),
  data = x.mdata,
  family = gaussian(),
  control = list(adapt_delta = 0.97, max_treedepth = 12), seed=1,
  chains = 4, cores = 4, iter = 3000
)

save(brms_multisp, file = paste0("models_brms/brms_multisp.rda"))

load("models_brms/brms_multi.rda")

pdf("models_brms/brms_multi_stats.pdf")
plot(brms_multi)
dev.off()
#plot(brms_m1.1)
#mcmc_plot(brms_m1.1, type = "acf")
summary(brms_multi)
mynames<-names(posterior_samples(brms_multi))
posterior_samples.df<-posterior_samples(brms_multi)[,names(posterior_samples(brms_multi)) %in% mynames[grep("sd_",mynames)]]
posterior_samples.df<-posterior_samples.df^2
posterior_samples.df<-posterior_samples.df %>% pivot_longer(cols=1:ncol(.),names_to="var",values_to="estimate")
vars<-posterior_samples.df$var %>% gsub("new_dgg","new.dgg",.) %>%  strsplit(.,"_") %>% unlist %>% matrix(ncol=5,byrow=T) %>% .[,c(2,4)]
posterior_samples.df<-cbind(vars,posterior_samples.df)
names(posterior_samples.df)[1:2]<-c("component","trait")
posterior_samples.df2<-posterior_samples(brms_multi) %>% dplyr::select(starts_with("sigma"))
posterior_samples.df2<-posterior_samples.df2^2
names(posterior_samples.df2)<-mynames[mynames %>% grep("sigma",.)] %>% strsplit(.,"_") %>% unlist %>% matrix(ncol=2,byrow=T) %>% .[,2]
posterior_samples.df2<-posterior_samples.df2 %>% pivot_longer(cols=1:ncol(.),names_to="trait",values_to="estimate")
posterior_samples.df2$component<-"sigma"
posterior_samples.df<-rbind(posterior_samples.df,posterior_samples.df2)
posterior_samples.df$iteration<-1:4000
posterior_samples.df<-posterior_samples.df %>% filter(iteration>2000)
posterior_samples.tot<-posterior_samples.df %>% group_by(trait,iteration) %>% summarise(tot=sum(estimate)) %>% ungroup
posterior_samples.df<-posterior_samples.df %>% left_join(posterior_samples.tot) %>% mutate(prop.var=estimate/tot)
posterior_samples.df<-posterior_samples.df %>% filter(component=="habitats") %>% group_by(trait) %>% summarise(quant05=quantile(prop.var,probs=0.05),quant5=quantile(prop.var,probs=0.5),quant95=quantile(prop.var,probs=0.95))
write.table(posterior_samples.df,file="models_brms/brms_multi_var.tab",quote=FALSE,sep='\t')
# pretty much same as standard anova

#-----regression discrete-----
library(lme4)
#library(lme4)
mixedmodel=TRUE
res_phenotypes<-data.frame(phenotype=c(),pvalue=c())
for (var1 in allvars){
    x<-x.regdata[!is.na(x.regdata[[var1]]) & !is.na(x.regdata$habitats),]
    x$habitats<-as.factor(x$habitats)
    x$block<-as.factor(x$block)
    x$y<-x[[var1]]
    if (mixedmodel){
    if (var1=="coleptile_color")
      {
        mylm<-glmer(y ~ habitats+(1|block)+(1|new_dgg),data=x,family="binomial")
        mylm0<-glmer(y ~ (1|block)+(1|new_dgg),data=x,family="binomial")
      } else
      {
        #mylm<-lmer(get(var1) ~ habitats+(0+habitats|block)+(0+habitats|new_dgg),data=x) #+(1+habitats|new_dgg)
        #mylm0<-lmer(get(var1) ~ (1|block),data=x) #+(1+habitats|new_dgg)
        mylm<-lmer(y ~ habitats+(1|block)+(1|new_dgg),data=x,REML=FALSE) #+(1+habitats|new_dgg)
        mylm0<-lmer(y ~ (1|block)+(1|new_dgg),data=x,REML=FALSE) #+(1+habitats|new_dgg)
      }
      p<-anova(mylm,mylm0)[[8]][2]
    } else {
    mylm<-lm(y ~ habitats*block,data=x)
    mylm0<-lm(y ~ block,data=x)
    p<-anova(mylm,mylm0)[['Pr(>F)']][2]
  }
    print(c(var1,p))
    res_phenotypes<-rbind(res_phenotypes,c(var1,p))
  }
names(res_phenotypes)<-c("phenotypes","pvalue")
res_phenotypes$pvalue<-as.numeric(res_phenotypes$pvalue)
res_phenotypes$bonferroni<-res_phenotypes$pvalue*nrow(res_phenotypes)
res_phenotypes$bonferroni[res_phenotypes$bonferroni>1]<-1
write.table(res_phenotypes,quote=FALSE,sep="\t",file="regression_phenotypes.tab")
res_t<-read.table(file="regression_phenotypes.tab")
library(poolr)
stouffer(res_t$pvalue)
#combined p-values with:      Stouffer's method
#number of p-values combined: 12
#test statistic:              2.688 ~ N(0,1)
#adjustment:                  none
#combined p-value:            0.00359546

m.cors<-x.ldata.clean %>% mutate(coleptile_color=as.numeric(as.factor(coleptile_color))-1) %>% filter(!is.na(coleptile_color)) %>% dplyr::select(allvars) %>% cor
stouffer(res_t$pvalue, adjust = "generalized", R = mvnconv(m.cors, target = "z"))
#combined p-values with:      Stouffer's method
#number of p-values combined: 12
#test statistic:              2.265 ~ N(0,1)
#adjustment:                  Strube's method
#combined p-value:            0.0118
library(ggcorrplot)
ggsave(ggcorrplot(m.cors),file="correlation_traits.pdf")


#-----regression abundance-----
library(lme4)
x.regdata<-x.ldata
x.regdata<-x.regdata[!duplicated(paste(x.regdata$new_dgg,x.regdata$block)),]
x.regdata$coleptile_color[x.regdata$coleptile_color==""]<-NA
x.regdata$coleptile_color<-as.factor(x.regdata$coleptile_color)
x.regdata$coleptile_color<-as.numeric(x.regdata$coleptile_color)-1
#library(lme4)
mixedmodel=TRUE
res_phenotypes<-data.frame(phenotype=c(),pvalue=c())
for (var1 in allvars){
    x<-x.regdata[!is.na(x.regdata[[var1]]) & !is.na(x.regdata$habitats),]
    x$habitats<-as.factor(x$habitats)
    x$block<-as.factor(x$block)
    x$y<-x[[var1]]
    if (mixedmodel){
    if (var1=="coleptile_color")
      {
        mylm<-glmer(y ~ ab_E_facing_slope+ab_Karst+ab_Middle_N_facing_slope+ab_Narrow_valley+ab_Plateau+ab_S_facing_slope+ab_Upper_N_facing_slope+ab_Valley+(1|block)+(1|new_dgg),data=x,family="binomial")
        mylm0<-glmer(y ~ (1|block)+(1|new_dgg),data=x,family="binomial")
      } else
      {
        #mylm<-lmer(get(var1) ~ habitats+(0+habitats|block)+(0+habitats|new_dgg),data=x) #+(1+habitats|new_dgg)
        #mylm0<-lmer(get(var1) ~ (1|block),data=x) #+(1+habitats|new_dgg)
        mylm<-lmer(y ~ ab_E_facing_slope+ab_Karst+ab_Middle_N_facing_slope+ab_Narrow_valley+ab_Plateau+ab_S_facing_slope+ab_Upper_N_facing_slope+ab_Valley+(1|block)+(1|new_dgg),data=x,REML=FALSE) #+(1+habitats|new_dgg)
        mylm0<-lmer(y ~ (1|block)+(1|new_dgg),data=x,REML=FALSE) #+(1+habitats|new_dgg)
      }
      p<-anova(mylm,mylm0)[[8]][2]
    }
    print(c(var1,p))
    res_phenotypes<-rbind(res_phenotypes,c(var1,p))
  }
names(res_phenotypes)<-c("phenotypes","pvalue")
res_phenotypes$pvalue<-as.numeric(res_phenotypes$pvalue)
res_phenotypes$bonferroni<-res_phenotypes$pvalue*nrow(res_phenotypes)
res_phenotypes$bonferroni[res_phenotypes$bonferroni>1]<-1
write.table(res_phenotypes,quote=FALSE,sep="\t",file="regression_phenotypes_weighted.tab")
#write.table(res_phenotypes,quote=FALSE,sep="\t",file="regression_phenotypes_weighted_onlysingletons_out.tab")


###DATA FOR MACHINE LEARNING ANALYSES
if (FALSE){
list_env<-names(y.ldata)[!names(y.ldata) %in% c("habitats","new_dgg")]
l.ldata<-list()
l.ldata.test<-list()
l.mdata<-list()
l.mdata.test<-list()
for (i in 1:6)
{
l.mdata[[i]]<-y.ldata %>% filter(block!=i) %>% dplyr::select(!new_dgg & !block)
l.mdata.test[[i]]<-y.ldata %>% filter(block==i) %>% dplyr::select(!new_dgg & !block)
l.ldata[[i]]<-l.mdata[[i]] %>% select(!coleptile_color)
l.ldata.test[[i]]<-l.mdata.test[[i]] %>% select(!coleptile_color)
}
}

model.kmeans<-kmeans(x.ldata.clean2 %>% select(!habitats),centers=3,nstart=25)
library(caret)
model.lda2<-train(habitats ~.,method="lda2",data=x.ldata.clean2,trControl=trainControl(method="LGOCV")) #,p=0.3,number=10))
varImp(model.lda2) #less seeds
model.lda2$finalModel
model.tree<-train(habitats ~.,method="treebag",data=x.ldata.clean3,trControl=trainControl(method="LGOCV",p=0.3,number=10))
varImp(model.tree) #more seeds
model.rf<-train(habitats ~.,method="rf",data=x.ldata.clean3,trControl=trainControl(method="LGOCV"))
varImp(model.rf) #more seeds #best
model.rpart<-train(habitats ~.,method="rpart",data=x.ldata.clean3,trControl=trainControl(method="LGOCV",p=0.8,number=10))
model.bafFDA<-train(habitats ~.,method="bagFDA",data=x.ldata.clean3,trControl=trainControl(method="LGOCV",p=0.3,number=10))
model.lda2<-train(habitats ~.,method="lda2",data=x.ldata.clean4,trControl=trainControl(method="LGOCV")) #,p=0.3,number=10))
#dimen  Accuracy   Kappa
#5      0.3586325  0.2136130
#Accuracy was used to select the optimal model using the largest value.
#The final value used for the model was dimen = 5.

predicted<-predict(model.lda2$finalModel,newdata=x.ldata.clean4 %>% select(!habitats))
lda_plot <- data.frame(habitats=predicted$class, predicted$x)

p<-ggplot(lda_plot, aes(LD1, LD2)) +
  geom_point(aes(color = habitats)) +
  theme(strip.background = element_rect(fill = "white"),axis.text.x = element_text(angle = 90, hjust=1), panel.background = element_rect(fill = 'white',color='black'),panel.grid.major = element_line(color = 'white', linetype = 'longdash',size=0),aspect.ratio=1)
ggsave(p,file="lda_LD1LD2.pdf")
p<-ggplot(lda_plot, aes(LD3, LD4)) +
  geom_point(aes(color = habitats)) +
  theme(strip.background = element_rect(fill = "white"),axis.text.x = element_text(angle = 90, hjust=1), panel.background = element_rect(fill = 'white',color='black'),panel.grid.major = element_line(color = 'white', linetype = 'longdash',size=0),aspect.ratio=1)
ggsave(p,file="lda_LD3LD4.pdf")
coefficients<-as.data.frame(model.lda2$finalModel[[4]])
write.table(coefficients,quote=FALSE,sep="\t",file="lda_v1_coefficients.tab")
pdf("importance.lda.pdf")
plot(varImp(model.lda2))
dev.off()


plot(model.tree$finalModel)

#x.ldata<-#x.pheno[,((which(names(x.pheno)=="block"))+1):length(names(x.pheno))]
x.ldata<-l.ldata[[1]]
x.ldata.clean<-x.ldata[!is.na(x.ldata$MeanLength) & !is.na(x.ldata$anthesis_time) &!is.na(x.ldata$main_length) & !is.na(x.ldata$erect_1_5), ]
model.rf<-train(habitats ~.,method="rf",data=x.ldata.clean,trControl=trainControl(method="LGOCV",p=0.1,number=10))
model.rf<-train(habitats ~.,method="rf",data=x.ldata.clean,trControl=trainControl(method="LGOCV",p=0.1,number=10))

model.rf<-train(habitats ~.,method="rf",data=x.ldata.clean,trControl=trainControl(method="LGOCV",p=0.3,number=10))
model.rf<-train(habitats ~.,method="rf",data=x.ldata.clean,trControl=trainControl(method="LGOCV",p=0.3,number=10))
sum(x.ldata.clean$habitats==predict(model.rf))/length(predict(model.rf))

sum(x.ldata.clean$habitats==predict(model.rf))/length(predict(model.rf))
varImp(model.rf)[[1]] %>% apply(MARGIN=1,FUN=mean)

model.rf<-train(method="pca",data=x.ldata.clean)
ctrl <- trainControl(method = "repeatedcv",
                     repeats = 3,
                     classProbs = TRUE,
                     preProcOptions = list(thresh = 0.85))#, #or list(pcaComp = 7)
                     #summaryFunction = twoClassSummary)

multinomFit <- train(habitats~., x.ldata.clean2,
                     method = "multinom",
#                     family=binomial,
                     metric = "ROC",
                     verbose = TRUE,
                     preProcess=c("center", "scale", "pca"))
                     #,
                     #trControl = ctrl)

,trControl=trainControl(method="LGOCV",p=0.1,number=10))

#----------decision tree---------
library(MASS)
library(rpart)
library(caret)
library(rpart.plot)
model.tree<-rpart(habitats ~.,method="class",data=x.ldata)
pdf("model.tree.pdf")
rpart.plot(model.tree,tweak=1.1,snip=FALSE,extra=5)
dev.off()
x.ldata.temp<-x.ldata[!is.na(x.ldata$MeanArea) & !is.na(x.ldata$anthesis_time) &!is.na(x.ldata$main_length) & !is.na(x.ldata$erect_1_5), ]
model.tree<-rpart(habitats ~.,method="class",data=x.ldata.temp)
predicted<-predict(model.tree,type="class")
tree.table<-table(x.ldata.temp$habitats,predicted)
sum(diag(tree.table))/sum(tree.table) #0.4575937

#----------linear discriminant analysis---------

model.lda.l<-list()
performance.test.l<-c()
for (i in 1:6)
{
model<-lda(habitats~., data=l.ldata[[i]])
model.lda.l[[i]]<-model
predicted.classes<-model %>% predict(l.ldata[[i]]) %>% .[[1]]
mean(l.ldata[[i]]$habitats==predicted.classes,na.rm=T) #0.4122288
predicted.classes<-model %>% predict(l.ldata.test[[i]]) %>% .[[1]]
performance.test<-mean(l.ldata.test[[i]]$habitats==predicted.classes,na.rm=T) #0.4122288
performance.test.l<-c(performance.test.l,performance.test)
print(performance.test)
}

lda_plot <- data.frame(habitats=predict(model)$class, predict(model)$x)


p<-ggplot(lda_plot, aes(LD1, LD2)) +
  geom_point(aes(color = habitats)) + scale_color_manual(values=palette_from_tal) +
  theme(strip.background = element_rect(fill = "white"),axis.text.x = element_text(angle = 45, hjust=1), panel.background = element_rect(fill = 'white',color='black'),panel.grid.major = element_line(color = 'white', linetype = 'longdash',size=0),aspect.ratio=1)
ggsave(p,file="lda_LD1LD2_talette.pdf")
p<-ggplot(lda_plot, aes(LD3, LD4)) +
  geom_point(aes(color = habitats)) + scale_color_manual(values=palette_from_tal) +
  theme(strip.background = element_rect(fill = "white"),axis.text.x = element_text(angle = 45, hjust=1), panel.background = element_rect(fill = 'white',color='black'),panel.grid.major = element_line(color = 'white', linetype = 'longdash',size=0),aspect.ratio=1)
ggsave(p,file="lda_LD3LD4_talette.pdf")
coefficients<-as.data.frame(model[[4]])
write.table(coefficients,quote=FALSE,sep="\t",file="lda_v1_coefficients.tab")

#----------quadratic discriminant analysis---------
model.qda <- qda(habitats~., data=x.ldata)
predictions <- model.qda %>% predict(x.ldata)
# Model accuracy
mean(predictions$class == x.ldata$habitats,na.rm=T) #0.6134122
model.qda <- qda(habitats~., data=x.ldata)
x.ldata.noseed<-x.ldata[,which(! names(x.ldata) %in% c("weight_10_seeds","MeanArea","MeanLength","MeanWidth"))]
model.qda.noseed <- qda(habitats~., data=x.ldata.noseed)
predictions <- model.qda.noseed %>% predict(x.ldata.noseed)
# Model accuracy
mean(predictions$class == x.ldata.noseed$habitats,na.rm=T) #0.4554264

