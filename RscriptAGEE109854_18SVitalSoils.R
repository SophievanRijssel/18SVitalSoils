rm(list=ls())
library(ggplot2)
library(vegan)
library(reshape2)
library(Hmisc)
library(plotrix)
library(phyloseq)
library(MASS)
#library(SpiecEasi)
# library(bioDist) #package ???bioDist??? is not available (for R version 3.2.2)
library(igraph)
library(car)
library(coin)
#library(edgeR)
library(formatR)
library(gridExtra)
library(gplots)
library(indicspecies)
library(sciplot)
library(ape)
library(grid)
#library(RVAideMemoire)
library(gridBase)
library(TukeyC)
library(corrplot)
#library(userfriendlyscience)
library(caret)
library(multcompView)
library(tidyverse)
library(data.table)
library("ape")
library("car")
library(plyr)
library("dplyr")
library("scales")
library("grid")
library(ade4)
library(tidyr)
library(jsonlite)   
library(survival)           
library(stringi)
library(psych)
library(corrplot)
#library(ggpubr)
#library(ggtern)
library(RColorBrewer)
library(gridExtra)
library("permute")
library("lattice")
library("cowplot")
library("reshape2")
library("VennDiagram")
library(lmerTest)
library(geosphere)
#library(plotrix)
library(ggplot2)
library(ggpubr)
#library(NetCoMi)
library(genefilter)
library(ggnetwork)
library(intergraph)
#library(ggnet)
library(genefilter)
library(intergraph)
library(GGally)
library(microbiome)
library(network)
#library(ggnet)
library(igraph)
library(metagenomeSeq)
#library(metagMisc)
library(wesanderson)
library(nlme)
library(lme4)
library(lmerTest)
library(ggplot2)
library(MASS)

set.seed(123)


#####OPEN Field data file #####
Sampledata.Nematodes.Field<-read.csv("Sampledata.Nematodes.Field.plusFGNEW.csv",sep=";")
Sampledata.Protists.Field<-read.csv("Sampledata.Protists.Field.plusFGNEW.csv",sep=";")

DataAGEE109854<-read.csv(xxxx)


#####data transformation
Sampledata.Nematodes.Field$X<-NULL
Sampledata.Protists.Field$X<-NULL

Sampledata.Nematodes.Clay.Field<-subset(Sampledata.Nematodes.Field,Soil=="Clay")
Sampledata.Nematodes.Sand.Field<-subset(Sampledata.Nematodes.Field,Soil=="Sand")
Sampledata.Protists.Clay.Field<-subset(Sampledata.Protists.Field,Soil=="Clay")
Sampledata.Protists.Sand.Field<-subset(Sampledata.Protists.Field,Soil=="Sand")
DataAGEE109854$percSOM<-DataAGEE109854$SOM*100
DataAGEE109854$pH.std<-scale(DataAGEE109854$pH)
DataAGEE109854$SOM.std<-scale(DataAGEE109854$SOM)
DataAGEE109854$Soil_Origin<-factor(DataAGEE109854$Soil_Origin,levels=c("Flevoland","Zeeland","North","Middle","South"))


#####Phyloseq 18s#####
taxtab18s.dt<-read.csv("taxa_18S_PR2_240.csv",header=TRUE,row.names=1,text=text,sep=";",stringsAsFactors=FALSE,check.names = FALSE)
seqtab18s<-readRDS("seqtab_18S_PR2_240.rds")
samples18s<-read.csv("samples_18S_PR2_240.csv",sep=";",row.names=1)
seqtab18s<-t(as.matrix(seqtab18s))
taxtab18s<-as.matrix(taxtab18s.dt)
phyloseq18spr2 <- phyloseq(otu_table(seqtab18s,taxa_are_rows=TRUE), sample_data(samples18s), tax_table(taxtab18s))

phyloseq18spr2
phyloseq18spr2= subset_samples(phyloseq18spr2,  sample_names(phyloseq18spr2)!= "193")
phyloseq18spr2= subset_samples(phyloseq18spr2,  sample_names(phyloseq18spr2)!= "194")
phyloseq18spr2= subset_samples(phyloseq18spr2,  sample_names(phyloseq18spr2)!= "195")
phyloseq18spr2= subset_samples(phyloseq18spr2,  sample_names(phyloseq18spr2)!= "81")
phyloseq18spr2= subset_samples(phyloseq18spr2,  sample_names(phyloseq18spr2)!= "189")
phyloseq18spr2= subset_samples(phyloseq18spr2,  sample_names(phyloseq18spr2)!= "169")
phyloseq18spr2= subset_samples(phyloseq18spr2,  sample_names(phyloseq18spr2)!= "170")
phyloseq18spr2= subset_samples(phyloseq18spr2,  sample_names(phyloseq18spr2)!= "171")
phyloseq18spr2filt<-prune_samples(sample_sums(phyloseq18spr2)>=1000, phyloseq18spr2)
phyloseq18spr2filt<-prune_taxa(taxa_sums(phyloseq18spr2filt)>3, phyloseq18spr2filt)



phyloseq18spr2filt
phynem<-subset_samples(phyloseq18spr2filt,Nematodes=="Nematodes")
phyprot<-subset_samples(phyloseq18spr2filt,!Nematodes=="Nematodes")
sample_names(phynem) <- sub("Index-N","", sample_names(phynem))
sample_names(phynem) <- sub("_F_filt.fastq.gz","", sample_names(phynem))
sample_names(phyprot) <- sub("Index-P","", sample_names(phyprot))
sample_names(phyprot) <- sub("_F_filt.fastq.gz","", sample_names(phyprot))

#####Rarecurve##### 
#
# jpeg('rarecurve18Snematodeisolates.jpg')
# rarecurve(t(otu_table(phynem)), step=50, cex=0.5)
# dev.off
#
# rarecurve(t(otu_table(phyprot)), step=50, cex=0.5)
# jpeg('rarecurve18Ssoil.jpg')
# rarecurve(t(otu_table(phyprot)), step=50, cex=0.5)
# dev.off

#####Split phyloseq nematodes and protists#####

sample_names(phynem)
phynem<-subset_taxa(phynem,Class=="Nematoda")
phynem
taxnem<-tax_table(phynem)
otunem<-otu_table(phynem)
samplenem<-sample_data(DataAGEE109854)
phynemnew<-phyloseq(otunem, taxnem, samplenem)
phynem<-phynemnew

phynem<-subset_samples(phynem,!Soil_Management=="NA")
taxa_names(phynem) <- paste(seq(ntaxa(phynem)),as.matrix(tax_table(phynem)[,7]),as.matrix(tax_table(phynem)[,9]))

taxprot<-tax_table(phyprot)
otuprot<-otu_table(phyprot)
sampleprot<-sample_data(DataAGEE109854)
phyprotnew<-phyloseq(otuprot, taxprot, sampleprot)
phyprot<-phyprotnew


#phyprot<-merge_phyloseq(phyprot,sample_data(DataAGEE109854))
phyprot<-subset_taxa(phyprot,Kingdom_OLD=="Protists")
phyprot<-subset_samples(phyprot,!Full_code=="NA")
phyprot<-prune_samples(sample_sums(phyprot)>500,phyprot)

taxa_names(phyprot) <- paste(seq(ntaxa(phyprot)),as.matrix(tax_table(phyprot)[,7]),as.matrix(tax_table(phyprot)[,9]))




#####Make all subdatasets#####
filter <- phyloseq::genefilter_sample(phynem,filterfun_sample(function(x) x >= 1),A = 3)
PhyseqNematodes<- prune_taxa(filter, phynem)

filter <- phyloseq::genefilter_sample(phyprot,filterfun_sample(function(x) x >= 1),A = 3)
PhyseqProtists<- prune_taxa(filter, phyprot)

PhyseqNematodes<-prune_taxa(taxa_sums(PhyseqNematodes)>0,PhyseqNematodes)
PhyseqProtists<-prune_taxa(taxa_sums(PhyseqProtists)>0,PhyseqProtists)
library(metagMisc)
PhyseqNematodescss<-phyloseq_transform_css(PhyseqNematodes)
PhyseqProtistscss<-phyloseq_transform_css(PhyseqProtists)

PhyseqNematodesClaycss<-subset_samples(PhyseqNematodescss,Soil=="Clay")
PhyseqProtistsClaycss<-subset_samples(PhyseqProtistscss,Soil=="Clay")

PhyseqNematodesSandcss<-subset_samples(PhyseqNematodescss,Soil=="Sand")
PhyseqProtistsSandcss<-subset_samples(PhyseqProtistscss,Soil=="Sand")


PhyseqNematodesClaycss<-prune_taxa(taxa_sums(PhyseqNematodesClaycss)>0,PhyseqNematodesClaycss)
PhyseqProtistsClaycss<-prune_taxa(taxa_sums(PhyseqProtistsClaycss)>0,PhyseqProtistsClaycss)
PhyseqNematodesSandcss<-prune_taxa(taxa_sums(PhyseqNematodesSandcss)>0,PhyseqNematodesSandcss)
PhyseqProtistsSandcss<-prune_taxa(taxa_sums(PhyseqProtistsSandcss)>0,PhyseqProtistsSandcss)


Sampledata.Nematodes<-as(sample_data(PhyseqNematodescss), "data.frame")
Sampledata.Protists<-as(sample_data(PhyseqProtistscss), "data.frame")
all(Sampledata.Nematodes$Full_code %in% Sampledata.Protists$Full_code)
Sampledata.Nematodes$Full_code[!Sampledata.Nematodes$Full_code %in% Sampledata.Protists$Full_code]

subset(Sampledata.Protists)
Sampledata.Nematodes.Clay<-as(sample_data(PhyseqNematodesClaycss), "data.frame")
Sampledata.Protists.Clay<-as(sample_data(PhyseqProtistsClaycss), "data.frame")

Sampledata.Nematodes.Sand<-as(sample_data(PhyseqNematodesSandcss), "data.frame")
Sampledata.Protists.Sand<-as(sample_data(PhyseqProtistsSandcss), "data.frame")

Sampledata.Nematodes$Group<-"Nematodes"
Sampledata.Nematodes.Clay$Group<-"Nematodes"
Sampledata.Nematodes.Sand$Group<-"Nematodes"

Sampledata.Protists$Group<-"Protists"
Sampledata.Protists.Clay$Group<-"Protists"
Sampledata.Protists.Sand$Group<-"Protists"

OTU.Nematodes.css<-as.data.frame(t(otu_table(PhyseqNematodescss)))
OTU.Protists.css<-as.data.frame(t(otu_table(PhyseqProtistscss)))

OTU.Nematodes.Clay.css<-as.data.frame(t(otu_table(PhyseqNematodesClaycss)))
OTU.Protists.Clay.css<-as.data.frame(t(otu_table(PhyseqProtistsClaycss)))

OTU.Nematodes.Sand.css<-as.data.frame(t(otu_table(PhyseqNematodesSandcss)))
OTU.Protists.Sand.css<-as.data.frame(t(otu_table(PhyseqProtistsSandcss)))


PhyseqNematodes
PhyseqPhyseqNematodesSpecies<-tax_glom(PhyseqNematodes,taxrank="Species",NArm=FALSE)
PhyseqPhyseqNematodesGenus<-tax_glom(PhyseqNematodes,taxrank="Genus",NArm=FALSE)
PhyseqPhyseqNematodesFamily<-tax_glom(PhyseqNematodes,taxrank="Family",NArm=FALSE)


write.csv(as.data.frame(tax_table(PhyseqPhyseqNematodesSpecies)),"nematodespeciestaxatable.csv")

taxa_sums(PhyseqNematodes)
ntaxa
taxa_names()
View(tax_table(PhyseqNematodes)) #check for NAs

PhyseqBacteriaPhyla #check for number of phyla

PhyseqPhyseqNematodesSpecies
View(tax_table(PhyseqPhyseqNematodesSpecies))



PhyseqPhyseqNematodesGenus


PhyseqProtists
PhyseqPhyseqProtistsSpecies<-tax_glom(PhyseqProtists,taxrank="Species",NArm=FALSE)
PhyseqPhyseqProtistsGenus<-tax_glom(PhyseqProtists,taxrank="Genus",NArm=FALSE)
PhyseqPhyseqProtistsFamily<-tax_glom(PhyseqProtists,taxrank="Family",NArm=FALSE)
PhyseqPhyseqProtistsOrder<-tax_glom(PhyseqProtists,taxrank="Order",NArm=FALSE)
PhyseqPhyseqProtistsClass<-tax_glom(PhyseqProtists,taxrank="Class",NArm=FALSE)
PhyseqPhyseqProtistsDivision<-tax_glom(PhyseqProtists,taxrank="Division",NArm=FALSE)
PhyseqPhyseqProtistsSupergroups<-tax_glom(PhyseqProtists,taxrank="Supergroup",NArm=FALSE)



View(tax_table(PhyseqPhyseqProtistsSpecies))
View(tax_table(PhyseqPhyseqProtistsGenus))
View(tax_table(PhyseqPhyseqProtistsSpecies))
View(tax_table(PhyseqPhyseqProtistsOrder))
View(tax_table(PhyseqPhyseqProtistsClass))

write.csv(as.data.frame(tax_table(PhyseqPhyseqProtistsSpecies)),"protistspeciestaxatable.csv")

#####FGSubPhyseqs#####
PhyseqBacterivorousNematodes<-subset_taxa(PhyseqNematodes,Function_FeedingType=="Bacterivores")
PhyseqEprhfeederNematodes<-subset_taxa(PhyseqNematodes,Function_FeedingType=="Epidermal/root hair feeders")
PhyseqFungivorousNematodes<-subset_taxa(PhyseqNematodes,Function_FeedingType=="Fungivores")
PhyseqOmnivorousNematodes<-subset_taxa(PhyseqNematodes,Function_FeedingType=="Omnivores")
PhyseqPlantPathogenicNematodes<-subset_taxa(PhyseqNematodes,Function_FeedingType=="Plant Pathogen")
PhyseqPredatorousNematodes<-subset_taxa(PhyseqNematodes,Function_FeedingType=="Predators")

PhyseqParasiticProtists<-subset_taxa(PhyseqProtists, Function_FeedingType=="Parasite")
PhyseqPhagotrophicProtists<-subset_taxa(PhyseqProtists, Function_FeedingType=="Phagotroph")
PhyseqPhototrophicProtists<-subset_taxa(PhyseqProtists, Function_FeedingType=="Phototroph")
PhyseqPlantPathogenicProtists<-subset_taxa(PhyseqProtists, Function_FeedingType=="Plant Pathogen")
PhyseqSaprotrophicProtists<-subset_taxa(PhyseqProtists, Function_FeedingType=="Saprotroph")


#####Plant pathogenic groups#####
PhyseqPlantPathogenicNematodes

PhyseqPlantPathogenicProtists

PhyseqPlantPathogenicNematodessp<-tax_glom(PhyseqPlantPathogenicNematodes,taxrank="Species",NArm=FALSE)
PhyseqPlantPathogenicNematodesgen<-tax_glom(PhyseqPlantPathogenicNematodes,taxrank="Genus",NArm=FALSE)

PhyseqPlantPathogenicProtistsgen<-tax_glom(PhyseqPlantPathogenicProtists,taxrank="Genus",NArm=FALSE)
PhyseqPlantPathogenicNematodessp
PhyseqPlantPathogenicNematodesgen
PhyseqPlantPathogenicProtistsgen
ASV.PlPatNem<-data.frame(t(otu_table(PhyseqPlantPathogenicNematodesgen)))
ASV.PlPatProt<-data.frame(t(otu_table(PhyseqPlantPathogenicProtistsgen)))

colSums(ASV.PlPatNem==0)
colSums(!ASV.PlPatNem==0)
colSums(ASV.PlPatProt==0)
colSums(!ASV.PlPatProt==0)
colSums(ASV.PlPatNem==0)/nrow(ASV.PlPatNem)*100
colSums(ASV.PlPatProt==0)/nrow(ASV.PlPatProt)*100
# go for 40%; nematodes only Pratylenchus, protists 4 genera.
plpatnemdata<-subset(ASV.PlPatNem,select=c("X11.Chromadorea_X.Pratylenchus_thornei"))
#t(otu_table(PhyseqPlantPathogenicNematodessp)[1:4])
plpatprotdata<-data.frame(t(otu_table(PhyseqPlantPathogenicProtistsgen)[1:4]))
samplefield<-subset(DataAGEE109854,select=c("ID_sample","Fieldnr"))
samplefield$ID_sample<-as.factor(samplefield$ID_sample)
plpatprotdata$ID_sample<-as.factor(row.names(plpatprotdata))
plpatnemdata$ID_sample<-as.factor(row.names(plpatnemdata))
plpatdata<-left_join(plpatprotdata,plpatnemdata,by="ID_sample")
plpatdata<-left_join(plpatdata,samplefield,by="ID_sample")
plpatdata[is.na(plpatdata)] <- 0

plpatdatafield<-aggregate.data.frame(plpatdata,by=list(plpatdata$Fieldnr),FUN=mean)
names(plpatdatafield)
plpatdatafield$ID_sample<-NULL
plpatdatafield$Group.1<-NULL
names(plpatdatafield)<-c("Polymyxa","Pythium","Spongospora","Peronosporales","Pratylenchus","Fieldnr")

Fielddat<-subset(DataAGEE109854,Pseudorep=="1",select=c("Fieldnr","Soil","ConOrg","Yearssince.as.if"))

plpatdatafield<-left_join(plpatdatafield,Fielddat,by="Fieldnr")

plpatdatafield.hell<-plpatdatafield
plpatdatafield.hell[1:5]<- vegan::decostand(as.matrix(plpatdatafield[1:5]),method="hellinger")
plpatdatafieldclay.hell<-subset(plpatdatafield.hell,Soil=="Clay")
plpatdatafieldsand.hell<-subset(plpatdatafield.hell,Soil=="Sand")

plpatdatafieldpoisson<-plpatdatafield
plpatdatafieldpoisson[1:5]<-round(plpatdatafield[1:5], digits = 0)
plpatdatafieldclaypoisson<-subset(plpatdatafieldpoisson,Soil=="Clay")
plpatdatafieldsandpoisson<-subset(plpatdatafieldpoisson,Soil=="Sand")

plpatdatafieldclay<-subset(plpatdatafield,Soil=="Clay")
plpatdatafieldsand<-subset(plpatdatafield,Soil=="Sand")

hist(plpatdatafield.hell$Polymyxa)
hist(plpatdatafield.hell$Pythium) 
hist(plpatdatafield.hell$Spongospora)
hist(plpatdatafield.hell$Peronosporales)
hist(plpatdatafield.hell$Pratylenchus)

hist(plpatdatafield$Polymyxa)
hist(plpatdatafield$Pythium)
hist(plpatdatafield$Spongospora)
hist(plpatdatafield$Peronosporales)
hist(plpatdatafield$Pratylenchus)

hist(plpatdatafieldpoisson$Polymyxa)
hist(plpatdatafieldpoisson$Pythium)
hist(plpatdatafieldpoisson$Spongospora)
hist(plpatdatafieldpoisson$Peronosporales)
hist(plpatdatafieldpoisson$Pratylenchus)

summary(m1 <- glm(formula =Polymyxa~ConOrg*Yearssince.as.if,family="poisson",data=plpatdatafieldclaypoisson))
summary(m1 <- glm(formula =Pythium~ConOrg*Yearssince.as.if,family="poisson",data=plpatdatafieldclaypoisson))
summary(m1 <- glm(formula =Spongospora~ConOrg*Yearssince.as.if,family="poisson",data=plpatdatafieldclaypoisson))
summary(m1 <- glm(formula =Peronosporales~ConOrg*Yearssince.as.if,family="poisson",data=plpatdatafieldclaypoisson))
summary(m1 <- glm(formula =Pratylenchus~ConOrg*Yearssince.as.if,family="poisson",data=plpatdatafieldclaypoisson))

summary(m1 <- glm(formula =Polymyxa~ConOrg*Yearssince.as.if,family="poisson",data=plpatdatafieldsandpoisson))
summary(m1 <- glm(formula =Pythium~ConOrg*Yearssince.as.if,family="poisson",data=plpatdatafieldsandpoisson))
summary(m1 <- glm(formula =Spongospora~ConOrg*Yearssince.as.if,family="poisson",data=plpatdatafieldsandpoisson))
summary(m1 <- glm(formula =Peronosporales~ConOrg*Yearssince.as.if,family="poisson",data=plpatdatafieldsandpoisson))
summary(m1 <- glm(formula =Pratylenchus~ConOrg*Yearssince.as.if,family="poisson",data=plpatdatafieldsandpoisson))
#all have huge overdispersion because residual deviance / residual degrees of freedom should be around 1.
# so then we move to quasipoisson
summary(m1 <- glm(formula =Polymyxa~ConOrg*Yearssince.as.if,family="quasipoisson",data=plpatdatafieldclaypoisson))
summary(m1 <- glm(formula =Pythium~ConOrg*Yearssince.as.if,family="quasipoisson",data=plpatdatafieldclaypoisson))
summary(m1 <- glm(formula =Spongospora~ConOrg*Yearssince.as.if,family="quasipoisson",data=plpatdatafieldclaypoisson))
summary(m1 <- glm(formula =Peronosporales~ConOrg*Yearssince.as.if,family="quasipoisson",data=plpatdatafieldclaypoisson))
summary(m1 <- glm(formula =Pratylenchus~ConOrg*Yearssince.as.if,family="quasipoisson",data=plpatdatafieldclaypoisson))

summary(m1 <- glm(formula =Polymyxa~ConOrg*Yearssince.as.if,family="quasipoisson",data=plpatdatafieldsandpoisson))
summary(m1 <- glm(formula =Pythium~ConOrg*Yearssince.as.if,family="quasipoisson",data=plpatdatafieldsandpoisson))
summary(m1 <- glm(formula =Spongospora~ConOrg*Yearssince.as.if,family="quasipoisson",data=plpatdatafieldsandpoisson))
summary(m1 <- glm(formula =Peronosporales~ConOrg*Yearssince.as.if,family="quasipoisson",data=plpatdatafieldsandpoisson)) #this one is the only with a dispersion parameter of 3; larger than 1.5 which means something has to be done to correct for it in the quasipoisson; or we first just try the negative binomial
summary(m1 <- glm(formula =Pratylenchus~ConOrg*Yearssince.as.if,family="quasipoisson",data=plpatdatafieldsandpoisson))
#if dispersion parameter is larger than 15 or 20 which it is in almost all cases we should move to negative binomial (Zuur et al; Ecological modelling in R; 2016).

#clay
summary(m1 <- glm.nb(formula =Polymyxa~ConOrg*Yearssince.as.if,data=plpatdatafieldclaypoisson,link = "log"))
summary(m2 <- glm.nb(formula =Pythium~ConOrg*Yearssince.as.if,data=plpatdatafieldclaypoisson,link = "log")) #time
summary(m2 <- glm.nb(formula =Pythium~ConOrg,data=plpatdatafieldclaypoisson,link = "log")) #time
summary(m3 <- glm.nb(formula =Spongospora~ConOrg*Yearssince.as.if,data=plpatdatafieldclaypoisson,link = "log")) #man+time
ggplot(aes(x=Yearssince.as.if,y=log(Spongospora),color=ConOrg),data=plpatdatafieldclaypoisson)+geom_point()+geom_smooth(method="lm")
summary(m4 <- glm.nb(formula =Peronosporales~ConOrg*Yearssince.as.if,data=plpatdatafieldclaypoisson,link = "log"))
summary(m5 <- glm.nb(formula =Pratylenchus~ConOrg*Yearssince.as.if,data=plpatdatafieldclaypoisson,link = "log"))

#sand
summary(m6 <- glm.nb(formula =Polymyxa~ConOrg*Yearssince.as.if,data=plpatdatafieldsandpoisson,link = "log")) #time
summary(m6 <- glm.nb(formula =Polymyxa~ConOrg,data=plpatdatafieldsandpoisson,link = "log")) #time
summary(m7 <- glm.nb(formula =Pythium~ConOrg*Yearssince.as.if,data=plpatdatafieldsandpoisson,link = "log")) #management
summary(m7 <- glm.nb(formula =Pythium~ConOrg,data=plpatdatafieldsandpoisson,link = "log")) #management
summary(m8 <- glm.nb(formula =Spongospora~ConOrg*Yearssince.as.if,data=plpatdatafieldsandpoisson,link = "log"))
summary(m9 <- glm.nb(formula =Peronosporales~ConOrg*Yearssince.as.if,data=plpatdatafieldsandpoisson,link = "log"))
summary(m10 <- glm.nb(formula =Pratylenchus~ConOrg*Yearssince.as.if,data=plpatdatafieldsandpoisson,link = "log"))
#dispersion parameters are all quite OK!!!!
op <- par(mfrow = c(2,2))
plot(m1)
plot(m2)
plot(m3)
plot(m4)
plot(m5)
plot(m6)
plot(m7)
plot(m8)
plot(m9)
plot(m10)
#looks all fine I would say
par(op)
summary(m1 <- glm.nb(formula =Polymyxa~ConOrg,data=plpatdatafield,link = log))
summary(m1 <- glm.nb(formula =Pythium~ConOrg,data=plpatdatafield,link = log))
summary(m1 <- glm.nb(formula =Spongospora~ConOrg,data=plpatdatafield,link = log))
summary(m1 <- glm.nb(formula =Peronosporales~ConOrg,data=plpatdatafield,link = log))
summary(m1 <- glm.nb(formula =Pratylenchus~ConOrg,data=plpatdatafield,link = log))
#seems to do better the negative binomial test on the "raw (=css?)" data?
glm.nb(formula = daysabs ~ math + prog, data = dat, init.theta = 1.032713156,
     link = log)

anova(lm(Polymyxa~ConOrg,data=plpatdatafieldclay))
anova(lm(Polymyxa~ConOrg,data=plpatdatafieldsand))

anova(lm(Pythium~ConOrg,data=plpatdatafield))
anova(lm(Pythium~ConOrg,data=plpatdatafieldclay))
anova(lm(Pythium~ConOrg,data=plpatdatafieldsand))
aggregate(Pythium~ConOrg,data=plpatdatafieldsand,FUN=mean) #**

anova(lm(Spongospora~ConOrg,data=plpatdatafield))
anova(lm(Spongospora~ConOrg,data=plpatdatafieldclay))
anova(lm(Spongospora~ConOrg,data=plpatdatafieldsand))

anova(lm(Peronosporales~ConOrg,data=plpatdatafield))
anova(lm(Peronosporales~ConOrg,data=plpatdatafieldclay))
anova(lm(Peronosporales~ConOrg,data=plpatdatafieldsand))
aggregate(Peronosporales~ConOrg,data=plpatdatafieldsand,FUN=mean) #*

anova(lm(Pratylenchus~ConOrg,data=plpatdatafield))
anova(lm(Pratylenchus~ConOrg,data=plpatdatafieldclay))
anova(lm(Pratylenchus~ConOrg,data=plpatdatafieldsand))

#####composition protist functional groups#####
PhyseqProtistsrelab<-transform_sample_counts(PhyseqProtists, function(x) x / sum(x) )
PhyseqProtistsrelabClay<-subset_samples(PhyseqProtistsrelab,Soil=="Clay")
PhyseqProtistsrelabSand<-subset_samples(PhyseqProtistsrelab,Soil=="Sand")

ProtistsClaydf<-psmelt(PhyseqProtistsrelabClay)
ProtistsClaydfFGsum<-aggregate(Abundance~Sample*Soil_Management*Function_FeedingType*Fieldnr,data=ProtistsClaydf,FUN=sum)
ProtistsClaydfFGsum<-ProtistsClaydfFGsum[!ProtistsClaydfFGsum$Abundance==0,]
ProtistsClaydfFGsum
ProtistsClaydfdfa<-aggregate(Abundance~Fieldnr*Soil_Management*Function_FeedingType,data=ProtistsClaydfFGsum,FUN=mean)
ProtistsClaydfdfaParasite<-ProtistsClaydfFGsum[ProtistsClaydfdfa$Function_FeedingType=="Parasite",]
ProtistsClaydfdfaPhagotroph<-ProtistsClaydfFGsum[ProtistsClaydfdfa$Function_FeedingType=="Phagotroph",]
ProtistsClaydfdfaPhototroph<-ProtistsClaydfFGsum[ProtistsClaydfdfa$Function_FeedingType=="Phototroph",]
ProtistsClaydfdfaPlantPathogen<-ProtistsClaydfFGsum[ProtistsClaydfdfa$Function_FeedingType=="Plant Pathogen",]
ProtistsClaydfdfaSaprotroph<-ProtistsClaydfFGsum[ProtistsClaydfdfa$Function_FeedingType=="Saprotroph",]

t.test(Abundance~Soil_Management,ProtistsClaydfdfaParasite)
t.test(Abundance~Soil_Management,ProtistsClaydfdfaPhagotroph)
t.test(Abundance~Soil_Management,ProtistsClaydfdfaPhototroph)
t.test(Abundance~Soil_Management,ProtistsClaydfdfaPlantPathogen)
t.test(Abundance~Soil_Management,ProtistsClaydfdfaSaprotroph)

phynem
phyprot

#sum per fg
datafigS3ProtistsClay<-aggregate(Abundance~Soil_Management*Function_FeedingType,ProtistsClaydfFGsum,FUN=mean)
ggplot(data=datafigS3ProtistsClay)+
  geom_bar(aes(x=Soil_Management, y=Abundance, fill=Function_FeedingType),stat="identity") +
  labs(x="Management",y="Relative abundance per function (%)") +
  scale_fill_manual(values = phylum_colors[10:15])+
  scale_color_manual(values = phylum_colors[10:15])+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),axis.title=element_text(size=20,face="bold"),legend.title=element_text(size=15),
        legend.text=element_text(size=15),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())

ProtistsSanddf<-psmelt(PhyseqProtistsrelabSand)
ProtistsSanddfFGsum<-aggregate(Abundance~Sample*Soil_Management*Function_FeedingType*Fieldnr,data=ProtistsSanddf,FUN=sum)
ProtistsSanddfFGsum
ProtistsSanddfFGsum<-ProtistsSanddfFGsum[!ProtistsSanddfFGsum$Abundance==0,]
ProtistsSanddfdfa<-aggregate(Abundance~Fieldnr*Soil_Management*Function_FeedingType,data=ProtistsSanddfFGsum,FUN=mean)
ProtistsSanddfdfaParasite<-ProtistsSanddfFGsum[ProtistsSanddfdfa$Function_FeedingType=="Parasite",]
ProtistsSanddfdfaPhagotroph<-ProtistsSanddfFGsum[ProtistsSanddfdfa$Function_FeedingType=="Phagotroph",]
ProtistsSanddfdfaPhototroph<-ProtistsSanddfFGsum[ProtistsSanddfdfa$Function_FeedingType=="Phototroph",]
ProtistsSanddfdfaPlantPathogen<-ProtistsSanddfFGsum[ProtistsSanddfdfa$Function_FeedingType=="Plant Pathogen",]
ProtistsSanddfdfaSaprotroph<-ProtistsSanddfFGsum[ProtistsSanddfdfa$Function_FeedingType=="Saprotroph",]

t.test(Abundance~Soil_Management,ProtistsSanddfdfaParasite)
t.test(Abundance~Soil_Management,ProtistsSanddfdfaPhagotroph)
t.test(Abundance~Soil_Management,ProtistsSanddfdfaPhototroph)
t.test(Abundance~Soil_Management,ProtistsSanddfdfaPlantPathogen)
t.test(Abundance~Soil_Management,ProtistsSanddfdfaSaprotroph)


#sum per fg
datafigS3ProtistsSand<-aggregate(Abundance~Soil_Management*Function_FeedingType,ProtistsSanddfFGsum,FUN=mean)
ggplot(data=datafigS3ProtistsSand)+
  geom_bar(aes(x=Soil_Management, y=Abundance, fill=Function_FeedingType),stat="identity") +
  labs(x="Management",y="Relative abundance per function (%)") +
  scale_fill_manual(values = phylum_colors[10:15])+
  scale_color_manual(values = phylum_colors[10:15])+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),axis.title=element_text(size=20,face="bold"),legend.title=element_text(size=15),
        legend.text=element_text(size=15),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())



#####nematode count diagrams#####
PhyseqNematodesrelab<-transform_sample_counts(PhyseqNematodes, function(x) x / sum(x) )
PhyseqNematodesrelabClay<-subset_samples(PhyseqNematodesrelab,Soil=="Clay")
PhyseqNematodesrelabSand<-subset_samples(PhyseqNematodesrelab,Soil=="Sand")
TAX.Nematodes.Clay<-as.matrix(tax_table(PhyseqNematodesrelabClay))
TAX.Nematodes.Sand<-as.matrix(tax_table(PhyseqNematodesrelabSand))

OTU.Nematodes.relab.Clay<-as.data.frame(t(otu_table(PhyseqNematodesrelabClay)))
OTU.Nematodes.relab.Sand<-as.data.frame(t(otu_table(PhyseqNematodesrelabSand)))

nemcountsClay<-Sampledata.Nematodes.Clay$nema.count.gram
OTU.Nematodes.Clay.matrix<-as.matrix(OTU.Nematodes.relab.Clay)

OTU.Nematodes.Clay.Count<-nemcountsClay*OTU.Nematodes.Clay.matrix
PhyseqNematodesCountsClay<-phyloseq(otu_table(OTU.Nematodes.Clay.Count,taxa_are_rows = FALSE),tax_table=(TAX.Nematodes.Clay),sample_data(Sampledata.Nematodes.Clay))

NematodesCountsClaydf<-psmelt(PhyseqNematodesCountsClay)
NematodesCountsClaydfFGsum<-aggregate(Abundance~Sample*Soil_Management*Function_FeedingType*Fieldnr,data=NematodesCountsClaydf,FUN=sum)
NematodesCountsClaydfFGsum
NematodesCountsClaydfFGsum<-NematodesCountsClaydfFGsum[!NematodesCountsClaydfFGsum$Abundance==0,]
NematodesClaydfFGsum<-NematodesCountsClaydfFGsum
NematodesClaydfdfa<-aggregate(Abundance~Fieldnr*Soil_Management*Function_FeedingType*Fieldnr,data=NematodesClaydfFGsum,FUN=mean)
NematodesClaydfdfaBacterivores<-NematodesClaydfdfa[NematodesClaydfdfa$Function_FeedingType=="Bacterivores",]
NematodesClaydfdfaFungivores<-NematodesClaydfdfa[NematodesClaydfdfa$Function_FeedingType=="Fungivores",]
NematodesClaydfdfaPredators<-NematodesClaydfdfa[NematodesClaydfdfa$Function_FeedingType=="Predators",]
NematodesClaydfdfaPlantPathogen<-NematodesClaydfdfa[NematodesClaydfdfa$Function_FeedingType=="Plant Pathogen",]
NematodesClaydfdfaOmnivores<-NematodesClaydfdfa[NematodesClaydfdfa$Function_FeedingType=="Omnivores",]

t.test(Abundance~Soil_Management,NematodesClaydfdfaBacterivores)
t.test(Abundance~Soil_Management,NematodesClaydfdfaFungivores)
t.test(Abundance~Soil_Management,NematodesClaydfdfaPredators)
t.test(Abundance~Soil_Management,NematodesClaydfdfaPlantPathogen)
t.test(Abundance~Soil_Management,NematodesClaydfdfaOmnivores)

#sum per fg
datafigXClay<-aggregate(Abundance~Soil_Management*Function_FeedingType,NematodesCountsClaydfFGsum,FUN=mean)
ggplot(data=datafigXClay)+
  geom_bar(aes(x=Soil_Management, y=Abundance, fill=Function_FeedingType),stat="identity") +
  labs(x="Management",y="Abundance per function\n (average of samples)") +
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),axis.title=element_text(size=20,face="bold"),legend.title=element_text(size=15),
        legend.text=element_text(size=15),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())

nemcountssand<-Sampledata.Nematodes.Sand$nema.count.gram
OTU.Nematodes.Sand.matrix<-as.matrix(OTU.Nematodes.relab.Sand)
OTU.Nematodes.Sand.Count<-nemcountssand*OTU.Nematodes.Sand.matrix
PhyseqNematodesCountsSand<-phyloseq(otu_table(OTU.Nematodes.Sand.Count,taxa_are_rows = FALSE),tax_table=(TAX.Nematodes.Sand),sample_data(Sampledata.Nematodes.Sand))
NematodesCountsSanddf<-psmelt(PhyseqNematodesCountsSand)
NematodesCountsSanddfFGsum<-aggregate(Abundance~Sample*Soil_Management*Function_FeedingType*Fieldnr,data=NematodesCountsSanddf,FUN=sum)
NematodesCountsSanddfFGsum
NematodesCountsSanddfFGsum<-NematodesCountsSanddfFGsum[!NematodesCountsSanddfFGsum$Abundance==0,]
NematodesSanddfFGsum<-NematodesCountsSanddfFGsum
NematodesSanddfdfa<-aggregate(Abundance~Fieldnr*Soil_Management*Function_FeedingType,data=NematodesSanddfFGsum,FUN=mean)
NematodesSanddfdfaBacterivores<-NematodesSanddfdfa[NematodesSanddfdfa$Function_FeedingType=="Bacterivores",]
NematodesSanddfdfaFungivores<-NematodesSanddfdfa[NematodesSanddfdfa$Function_FeedingType=="Fungivores",]
NematodesSanddfdfaPredators<-NematodesSanddfdfa[NematodesSanddfdfa$Function_FeedingType=="Predators",]
NematodesSanddfdfaPlantPathogen<-NematodesSanddfdfa[NematodesSanddfdfa$Function_FeedingType=="Plant Pathogen",]
NematodesSanddfdfaOmnivores<-NematodesSanddfdfa[NematodesSanddfdfa$Function_FeedingType=="Omnivores",]

t.test(Abundance~Soil_Management,NematodesSanddfdfaBacterivores)
t.test(Abundance~Soil_Management,NematodesSanddfdfaFungivores)
t.test(Abundance~Soil_Management,NematodesSanddfdfaPredators)
t.test(Abundance~Soil_Management,NematodesSanddfdfaPlantPathogen)
t.test(Abundance~Soil_Management,NematodesSanddfdfaOmnivores)


# not sure whether this works; this is relative abundance for feeding groups or so right? (it is not anymore the total abundances per sample)
# allnemcounts<-rbind(NematodesClaydfdfa,NematodesSanddfdfa)
# allnemcountstot<-aggregate(Abundance~Fieldnr*Soil_Management,allnemcounts,FUN=mean)
# summary(aov(Abundance~Soil_Management,data=allnemcountstot))
# aggregate(Abundance~Soil_Management,data=allnemcountstot,FUN=mean)
# 
# allnemcountstot2<-left_join(allnemcountstot,Sampledata.Nematodes.Field[1:8],by="Fieldnr")
# t.test(Abundance~Soil,allnemcountstot2)
# ggplot(allnemcountstot2,aes(x=Soil,y=Abundance,fill=Management))+geom_boxplot()

t.test(Abundance~Management,allnemcountstot2)
summary(aov(Abundance~Soil*Management,allnemcountstot2))

allnemcountstot2clay<-subset(allnemcountstot2,Soil=="Clay")
t.test(Abundance~Management,allnemcountstot2clay)
allnemcountstot2sand<-subset(allnemcountstot2,Soil=="Sand")
t.test(Abundance~Management,allnemcountstot2sand)

#sum per fg
datafigXSand<-aggregate(Abundance~Soil_Management*Function_FeedingType,NematodesCountsSanddfFGsum,FUN=mean)
ggplot(data=datafigXSand)+
  geom_bar(aes(x=Soil_Management, y=Abundance, fill=Function_FeedingType),stat="identity") +
  labs(x="Management",y="Abundance per function\n (average of samples)") +
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),axis.title=element_text(size=20,face="bold"),legend.title=element_text(size=15),
        legend.text=element_text(size=15),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())

           
#####reads after filtering and rarefying#####
#reads on sample-level
Nematodes.Reads = sample_sums(PhyseqNematodes)
Protists.Reads = sample_sums(PhyseqProtists)

Sampledata.Nematodes<-cbind(Sampledata.Nematodes,Nematodes.Reads)
Sampledata.Protists<-cbind(Sampledata.Protists,Protists.Reads)

hist(Nematodes.Reads)
hist(Protists.Reads)
hist(Nematodes.Reads,breaks=100)
hist(Protists.Reads,breaks=100)
summary(Sampledata.Nematodes$Nematodes.Reads)
summary(Sampledata.Protists$Protists.Reads)

standard_error(na.omit(Sampledata.Nematodes$Nematodes.Reads))
standard_error(na.omit(Sampledata.Protists$Protists.Reads))

anova(lm(Nematodes.Reads~ConOrg*Soil,data=Sampledata.Nematodes))
anova(lm(Protists.Reads~ConOrg*Soil,data=Sampledata.Protists))

PhyseqNematodes
PhyseqProtists


# PhyseqNematodesPhyla<-tax_glom(PhyseqNematodes,taxrank="Phylum",NArm=FALSE)
# PhyseqNematodesPhyla #check for number of phyla
#
# PhyseqProtistsPhyla<-tax_glom(PhyseqProtists,taxrank="Phylum",NArm=FALSE)
# PhyseqProtistsPhyla
# #View(tax_table(PhyseqProtistsPhyla))

# #####General diversity bootstrapped  - takes some time to run...#####
# #raw diversity
# Nematodes.raw.div<-estimate_richness(PhyseqNematodes,measures=c("Observed","Shannon"))
# names(Nematodes.raw.div)[names(Nematodes.raw.div) == 'Observed'] <- 'Observed.raw'
# names(Nematodes.raw.div)[names(Nematodes.raw.div) == 'Shannon'] <- 'Shannon.raw'
# Sampledata.Nematodes<-cbind(Sampledata.Nematodes,Nematodes.raw.div)
# Protists.raw.div<-estimate_richness(PhyseqProtists,measures=c("Observed","Shannon"))
# names(Protists.raw.div)[names(Protists.raw.div) == 'Observed'] <- 'Observed.raw'
# names(Protists.raw.div)[names(Protists.raw.div) == 'Shannon'] <- 'Shannon.raw'
# Sampledata.Protists<-cbind(Sampledata.Protists,Protists.raw.div)
# 
# anova(lm(Observed.raw~ConOrg*Soil,data=Sampledata.Nematodes))
# anova(lm(Observed.raw~ConOrg*Soil,data=Sampledata.Protists))
# anova(lm(Shannon.raw~ConOrg*Soil,data=Sampledata.Nematodes))
# anova(lm(Shannon.raw~ConOrg*Soil,data=Sampledata.Protists))
# anova(glm(Observed.raw~ConOrg*Soil,data=Sampledata.Protists))
# summary(glm(Observed.raw~ConOrg*Soil*Protists.Reads,data=Sampledata.Protists))
# ggplot(Sampledata.Nematodes,aes(x=Nematodes.Reads,y=Observed.raw,color=Soil))+geom_point()+geom_smooth(method="lm")
# ggplot(Sampledata.Protists,aes(x=Protists.Reads,y=Observed.raw,color=Soil))+geom_point()+geom_smooth(method="lm")
# 
# 
# #Hiiesalu/Tedersoo approach
# Sampledata.Protists$Observed.resid<-resid(lm(Observed.raw~sqrt(Protists.Reads),data=Sampledata.Protists))
# Sampledata.Nematodes$Observed.resid<-resid(lm(Observed.raw~sqrt(Nematodes.Reads),data=Sampledata.Nematodes))
# 
# Sampledata.Protists$Shannon.resid<-resid(lm(Shannon.raw~sqrt(Protists.Reads),data=Sampledata.Protists))
# Sampledata.Nematodes$Shannon.resid<-resid(lm(Shannon.raw~sqrt(Nematodes.Reads),data=Sampledata.Nematodes))
# 
# anova(lm(Observed.resid~ConOrg*Soil,data=Sampledata.Nematodes))
# anova(lm(Observed.resid~ConOrg*Soil,data=Sampledata.Protists))
# anova(lm(Shannon.resid~ConOrg*Soil,data=Sampledata.Nematodes))
# anova(lm(Shannon.resid~ConOrg*Soil,data=Sampledata.Protists))
# 
# #1000*rarefied; then averaged per sample
# 
# #nematodes to 850; protists to 500
# richbootNematodes<-as.data.frame(matrix(nrow = 187, ncol = 1001))
# richbootNematodes$Full_code<-sample_data(PhyseqNematodes)$Full_code
# richbootNematodesdf<-as.data.frame(matrix(nrow = 187, ncol = 1))
# richbootNematodesdf$Full_code<-sample_data(PhyseqNematodes)$Full_code
# rownames(richbootNematodesdf)<-richbootNematodesdf$Full_code
# richbootNematodesdf$Full_code<-NULL
# richbootNematodesdf[1]<-NULL
# 
# #Nematodes observed
# for (i in 1:1000) {
#   PhyseqNematodesrar<-rarefy_even_depth(PhyseqNematodes,sample.size=850)
#   richbootNematodes[i]<-estimate_richness(PhyseqNematodesrar,measures=c("Observed"))
#   richbootNematodesdf<-cbind(richbootNematodesdf,richbootNematodes[i])
# }
# Sampledata.Nematodes$Observed.rarmean<-rowMeans(richbootNematodesdf)
# 
# divbootNematodes<-as.data.frame(matrix(nrow = 187, ncol = 1001))
# divbootNematodes$Full_code<-sample_data(PhyseqNematodes)$Full_code
# divbootNematodesdf<-as.data.frame(matrix(nrow = 187, ncol = 1))
# divbootNematodesdf$Full_code<-sample_data(PhyseqNematodes)$Full_code
# rownames(divbootNematodesdf)<-divbootNematodesdf$Full_code
# divbootNematodesdf$Full_code<-NULL
# divbootNematodesdf[1]<-NULL
# 
# #Nematodes Shannon
# for (i in 1:1000) {
#   PhyseqNematodesrar<-rarefy_even_depth(PhyseqNematodes,sample.size=850)
#   divbootNematodes[i]<-estimate_richness(PhyseqNematodesrar,measures=c("Shannon"))
#   divbootNematodesdf<-cbind(divbootNematodesdf,divbootNematodes[i])
# }
# Sampledata.Nematodes$Shannon.rarmean<-rowMeans(divbootNematodesdf)
# 
# 
# richbootProtists<-as.data.frame(matrix(nrow = 199, ncol = 1001))
# richbootProtists$Full_code<-sample_data(PhyseqProtists)$Full_code
# richbootProtistsdf<-as.data.frame(matrix(nrow = 199, ncol = 1))
# richbootProtistsdf$Full_code<-sample_data(PhyseqProtists)$Full_code
# rownames(richbootProtistsdf)<-richbootProtistsdf$Full_code
# richbootProtistsdf$Full_code<-NULL
# richbootProtistsdf[1]<-NULL
# 
# #Protists observed
# for (i in 1:1000) {
#   PhyseqProtistsrar<-rarefy_even_depth(PhyseqProtists,sample.size=500)
#   richbootProtists[i]<-estimate_richness(PhyseqProtistsrar,measures=c("Observed"))
#   richbootProtistsdf<-cbind(richbootProtistsdf,richbootProtists[i])
# }
# Sampledata.Protists$Observed.rarmean<-rowMeans(richbootProtistsdf)
# 
# divbootProtists<-as.data.frame(matrix(nrow = 199, ncol = 1001))
# divbootProtists$Full_code<-sample_data(PhyseqProtists)$Full_code
# divbootProtistsdf<-as.data.frame(matrix(nrow = 199, ncol = 1))
# divbootProtistsdf$Full_code<-sample_data(PhyseqProtists)$Full_code
# rownames(divbootProtistsdf)<-divbootProtistsdf$Full_code
# divbootProtistsdf$Full_code<-NULL
# divbootProtistsdf[1]<-NULL
# 
# #Protists Shannon
# for (i in 1:1000) {
#   PhyseqProtistsrar<-rarefy_even_depth(PhyseqProtists,sample.size=500)
#   divbootProtists[i]<-estimate_richness(PhyseqProtistsrar,measures=c("Shannon"))
#   divbootProtistsdf<-cbind(divbootProtistsdf,divbootProtists[i])
# }
# Sampledata.Protists$Shannon.rarmean<-rowMeans(divbootProtistsdf)
# 
# 
# 
# 
# #####General PCoA#####
PCoA.Nematodes<-pcoa(vegdist(OTU.Nematodes.css,method="bray"))
PCoA.Protists<-pcoa(vegdist(OTU.Protists.css,method="bray"))
PCoA.Nematodes$values$Relative_eig[1:5]
PCoA.Protists$values$Relative_eig[1:5]
adonis2(vegdist(OTU.Nematodes.css,method="bray")~Soil*ConOrg,data=Sampledata.Nematodes,permutations=999)
adonis2(vegdist(OTU.Protists.css,method="bray")~Soil*ConOrg,data=Sampledata.Protists,permutations=999)

adonis2(vegdist(OTU.Nematodes.css,method="bray")~Soil,data=Sampledata.Nematodes,permutations=999)
adonis2(vegdist(OTU.Protists.css,method="bray")~Soil,data=Sampledata.Protists,permutations=999)
adonis2(vegdist(OTU.Nematodes.css,method="bray")~ConOrg,data=Sampledata.Nematodes,permutations=999)
adonis2(vegdist(OTU.Protists.css,method="bray")~ConOrg,data=Sampledata.Protists,permutations=999)


anova(betadisper(vegdist(OTU.Nematodes.css,method="bray"),Sampledata.Nematodes$Soil))
anova(betadisper(vegdist(OTU.Protists.css,method="bray"),Sampledata.Protists$Soil))
anova(betadisper(vegdist(OTU.Nematodes.css,method="bray"),Sampledata.Nematodes$ConOrg))
anova(betadisper(vegdist(OTU.Protists.css,method="bray"),Sampledata.Protists$ConOrg))
 

# #####PCoAS per soil type#####
ASV.Nematodes.Clay.css<-as.data.frame(t(otu_table(PhyseqNematodesClaycss)))
ASV.Protists.Clay.css<-as.data.frame(t(otu_table(PhyseqProtistsClaycss)))
ASV.Nematodes.Sand.css<-as.data.frame(t(otu_table(PhyseqNematodesSandcss)))
ASV.Protists.Sand.css<-as.data.frame(t(otu_table(PhyseqProtistsSandcss)))

SoilSpecificPCoA.Nematodes.Clay<-pcoa(vegdist(ASV.Nematodes.Clay.css,method="bray"))
SoilSpecificPCoA.Protists.Clay<-pcoa(vegdist(ASV.Protists.Clay.css,method="bray"))

SoilSpecificPCoA.Nematodes.Sand<-pcoa(vegdist(ASV.Nematodes.Sand.css,method="bray"))
SoilSpecificPCoA.Protists.Sand<-pcoa(vegdist(ASV.Protists.Sand.css,method="bray"))
#prop explained:
SoilSpecificPCoA.Nematodes.Clay$values$Relative_eig[1:5]
SoilSpecificPCoA.Protists.Clay$values$Relative_eig[1:5]
SoilSpecificPCoA.Nematodes.Sand$values$Relative_eig[1:5]
SoilSpecificPCoA.Protists.Sand$values$Relative_eig[1:5]


#
Sampledata.Nematodes.Clay<-as(sample_data(PhyseqNematodesClaycss), "data.frame")
Sampledata.Protists.Clay<-as(sample_data(PhyseqProtistsClaycss), "data.frame")

Sampledata.Nematodes.Sand<-as(sample_data(PhyseqNematodesSandcss), "data.frame")
Sampledata.Protists.Sand<-as(sample_data(PhyseqProtistsSandcss), "data.frame")

#statistics adonis (permanova):
adonis2(vegdist(ASV.Nematodes.Clay.css,method="bray")~ConOrg,data=Sampledata.Nematodes.Clay,permutations=999)
adonis2(vegdist(ASV.Protists.Clay.css,method="bray")~ConOrg,data=Sampledata.Protists.Clay,permutations=999)
adonis2(vegdist(ASV.Nematodes.Sand.css,method="bray")~ConOrg,data=Sampledata.Nematodes.Sand,permutations=999)
adonis2(vegdist(ASV.Protists.Sand.css,method="bray")~ConOrg,data=Sampledata.Protists.Sand,permutations=999)


anova(betadisper(vegdist(ASV.Nematodes.Clay.css,method="bray"),Sampledata.Nematodes.Clay$ConOrg))
anova(betadisper(vegdist(ASV.Nematodes.Sand.css,method="bray"),Sampledata.Nematodes.Sand$ConOrg))

anova(betadisper(vegdist(ASV.Protists.Clay.css,method="bray"),Sampledata.Protists.Clay$ConOrg))
anova(betadisper(vegdist(ASV.Protists.Sand.css,method="bray"),Sampledata.Protists.Sand$ConOrg))


Sampledata.Nematodes.Clay$SoilSpecificPCoA1 <- SoilSpecificPCoA.Nematodes.Clay$vectors[,c(1)]
Sampledata.Nematodes.Clay$SoilSpecificPCoA2 <- SoilSpecificPCoA.Nematodes.Clay$vectors[,c(2)]
Sampledata.Nematodes.Sand$SoilSpecificPCoA1 <- SoilSpecificPCoA.Nematodes.Sand$vectors[,c(1)]
Sampledata.Nematodes.Sand$SoilSpecificPCoA2 <- SoilSpecificPCoA.Nematodes.Sand$vectors[,c(2)]

Sampledata.Protists.Clay$SoilSpecificPCoA1 <- SoilSpecificPCoA.Protists.Clay$vectors[,c(1)]
Sampledata.Protists.Clay$SoilSpecificPCoA2 <- SoilSpecificPCoA.Protists.Clay$vectors[,c(2)]
Sampledata.Protists.Sand$SoilSpecificPCoA1 <- SoilSpecificPCoA.Protists.Sand$vectors[,c(1)]
Sampledata.Protists.Sand$SoilSpecificPCoA2 <- SoilSpecificPCoA.Protists.Sand$vectors[,c(2)]

Sampledata.Nematodes.Clay.Selection<-subset(Sampledata.Nematodes.Clay,select=c("Fieldnr","SoilSpecificPCoA1","SoilSpecificPCoA2"))
Sampledata.Nematodes.Sand.Selection<-subset(Sampledata.Nematodes.Sand,select=c("Fieldnr","SoilSpecificPCoA1","SoilSpecificPCoA2"))

Sampledata.Protists.Clay.Selection<-subset(Sampledata.Protists.Clay,select=c("Fieldnr","SoilSpecificPCoA1","SoilSpecificPCoA2"))
Sampledata.Protists.Sand.Selection<-subset(Sampledata.Protists.Sand,select=c("Fieldnr","SoilSpecificPCoA1","SoilSpecificPCoA2"))

Sampledata.Nematodes.Clay.Field.Data<-aggregate.data.frame(Sampledata.Nematodes.Clay.Selection[2:3],by=list(Sampledata.Nematodes.Clay.Selection$Fieldnr),na.rm=TRUE, na.action=NULL,FUN=mean)#check whether this is correct
colnames(Sampledata.Nematodes.Clay.Field.Data)<-c("Fieldnr","SoilSpecificPCoA1","SoilSpecificPCoA2")
Sampledata.Nematodes.Sand.Field.Data<-aggregate.data.frame(Sampledata.Nematodes.Sand.Selection[2:3],by=list(Sampledata.Nematodes.Sand.Selection$Fieldnr),na.rm=TRUE, na.action=NULL,FUN=mean)#check whether this is correct
colnames(Sampledata.Nematodes.Sand.Field.Data)<-c("Fieldnr","SoilSpecificPCoA1","SoilSpecificPCoA2")

Sampledata.Protists.Clay.Field.Data<-aggregate.data.frame(Sampledata.Protists.Clay.Selection[2:3],by=list(Sampledata.Protists.Clay.Selection$Fieldnr),na.rm=TRUE, na.action=NULL,FUN=mean)#check whether this is correct
colnames(Sampledata.Protists.Clay.Field.Data)<-c("Fieldnr","SoilSpecificPCoA1","SoilSpecificPCoA2")
Sampledata.Protists.Sand.Field.Data<-aggregate.data.frame(Sampledata.Protists.Sand.Selection[2:3],by=list(Sampledata.Protists.Sand.Selection$Fieldnr),na.rm=TRUE, na.action=NULL,FUN=mean)#check whether this is correct
colnames(Sampledata.Protists.Sand.Field.Data)<-c("Fieldnr","SoilSpecificPCoA1","SoilSpecificPCoA2")
Nematodes.soilspec.SoilSpecificPCoAs<-rbind(Sampledata.Nematodes.Clay.Field.Data,Sampledata.Nematodes.Sand.Field.Data)
Protists.soilspec.SoilSpecificPCoAs<-rbind(Sampledata.Protists.Clay.Field.Data,Sampledata.Protists.Sand.Field.Data)

Sampledata.Nematodes.Field<-merge(Sampledata.Nematodes.Field,Nematodes.soilspec.SoilSpecificPCoAs,by="Fieldnr")
Sampledata.Protists.Field<-merge(Sampledata.Protists.Field,Protists.soilspec.SoilSpecificPCoAs,by="Fieldnr")


#####ANOVA nematode counts#####
Sampledata.Nematodes$nema.count.gram
nemacountsfield<-aggregate(nema.count.gram~Fieldnr,data=Sampledata.Nematodes,na.rm=TRUE, na.action=NULL,FUN=mean)

Sampledata.Nematodes.Field<-full_join(Sampledata.Nematodes.Field,nemacountsfield)

summary(aov(nema.count.gram~Soil*Management,data=Sampledata.Nematodes.Field))
aggregate(nema.count.gram~Soil,data=Sampledata.Nematodes.Field,function(x) c(mean = mean(x), se = standard_error(x)))



#####clay vs sand abiotics#####
t.test(pH~Soil,Sampledata.Nematodes.Field)
t.test(SOM~Soil,Sampledata.Nematodes.Field)
t.test(PLFABacteria.PLFABacteria~Soil,Sampledata.Nematodes.Field)
t.test(PLFAFungi.PLFAFungi~Soil,Sampledata.Nematodes.Field)
t.test(Reads~Soil,Sampledata.Nematodes.Field)

t.test(Reads~Soil,Sampledata.Protists.Field)
t.test(Shannon.rar~Soil,Sampledata.Protists.Field)

summary(aov(pH~Soil*Management,Sampledata.Nematodes.Field))
summary(aov(SOM~Soil*Management,Sampledata.Nematodes.Field))
summary(aov(PLFABacteria.PLFABacteria~Soil*Management,Sampledata.Nematodes.Field))
summary(aov(PLFAFungi.PLFAFungi~Soil*Management,Sampledata.Nematodes.Field))

anova(lm(pH~Soil*Management,Sampledata.Nematodes.Field))
anova(lm(SOM~Soil*Management,Sampledata.Nematodes.Field))
anova(lm(PLFABacteria.PLFABacteria~Soil*Management,Sampledata.Nematodes.Field))
anova(lm(PLFAFungi.PLFAFungi~Soil*Management,Sampledata.Nematodes.Field))
anova(lm(nema.count.gram~Soil*Management,Sampledata.Nematodes.Field))

#this in a table
aggregate(pH~Soil*Management,Sampledata.Nematodes.Field,function(x) c(mean = mean(x), se = standard_error(x)))
aggregate(SOM~Soil*Management,Sampledata.Nematodes.Field,function(x) c(mean = mean(x), se = standard_error(x)))
aggregate(PLFABacteria.PLFABacteria~Soil*Management,Sampledata.Nematodes.Field,function(x) c(mean = mean(x), se = standard_error(x)))
aggregate(PLFAFungi.PLFAFungi~Soil*Management,Sampledata.Nematodes.Field,function(x) c(mean = mean(x), se = standard_error(x)))
aggregate(nema.count.gram~Soil*Management,Sampledata.Nematodes.Field,function(x) c(mean = mean(x), se = standard_error(x)))


#analyses with time per soil type; not a 3-way interaction
# anova(lm(pH~Soil*Management*Yearssince.as.if,Sampledata.Nematodes.Field))
# anova(lm(SOM~Soil*Management*Yearssince.as.if,Sampledata.Nematodes.Field))
# anova(lm(PLFABacteria.PLFABacteria~Soil*Management*Yearssince.as.if,Sampledata.Nematodes.Field))
# anova(lm(PLFAFungi.PLFAFungi~Soil*Management*Yearssince.as.if,Sampledata.Nematodes.Field))

anova(lm(pH~Management*Yearssince.as.if,Sampledata.Nematodes.Clay.Field))
anova(lm(SOM~Management*Yearssince.as.if,Sampledata.Nematodes.Clay.Field))
anova(lm(PLFABacteria.PLFABacteria~Management*Yearssince.as.if,Sampledata.Nematodes.Clay.Field))
anova(lm(PLFAFungi.PLFAFungi~Management*Yearssince.as.if,Sampledata.Nematodes.Clay.Field))
anova(lm(nema.count.gram~Management*Yearssince.as.if,Sampledata.Nematodes.Clay.Field))


anova(lm(pH~Management*Yearssince.as.if,Sampledata.Nematodes.Sand.Field))
anova(lm(SOM~Management*Yearssince.as.if,Sampledata.Nematodes.Sand.Field))
anova(lm(PLFABacteria.PLFABacteria~Management*Yearssince.as.if,Sampledata.Nematodes.Sand.Field))
anova(lm(PLFAFungi.PLFAFungi~Management*Yearssince.as.if,Sampledata.Nematodes.Sand.Field))
anova(lm(nema.count.gram~Management*Yearssince.as.if,Sampledata.Nematodes.Sand.Field))

#####Field level Read differences between soil types and management#####
LM3.Nematodes.Reads<-lm(Reads~Soil*Management,data=Sampledata.Nematodes.Field)
LM3.Protists.Reads<-lm(Reads~Soil*Management,data=Sampledata.Protists.Field)

anova(LM3.Nematodes.Reads)
anova(LM3.Protists.Reads)

aggregate(Reads~Soil,data=Sampledata.Nematodes.Field,function(x) c(mean = mean(x), se = standard_error(x)))
aggregate(Reads~Soil,data=Sampledata.Protists.Field,function(x) c(mean = mean(x), se = standard_error(x)))
aggregate(Reads~Management,data=Sampledata.Protists.Field,function(x) c(mean = mean(x), se = standard_error(x)))

aggregate(Reads~Soil*Management,data=Sampledata.Nematodes.Field,function(x) c(mean = mean(x), se = standard_error(x)))
aggregate(Reads~Soil*Management,data=Sampledata.Protists.Field,function(x) c(mean = mean(x), se = standard_error(x)))

#####field level diversity differences between soil types#####
#raw reads
AOVObservedSoilNematodes<-aov(Observed.raw~Soil,data=Sampledata.Nematodes.Field)
AOVObservedSoilProtists<-aov(Observed.raw~Soil,data=Sampledata.Protists.Field)
AOVShannonSoilNematodes<-aov(Shannon.raw~Soil,data=Sampledata.Nematodes.Field)
AOVShannonSoilProtists<-aov(Shannon.raw~Soil,data=Sampledata.Protists.Field)

summary(AOVObservedSoilNematodes)
summary(AOVObservedSoilProtists)
summary(AOVShannonSoilNematodes)
summary(AOVShannonSoilProtists)

aggregate(Observed.raw~Soil,data=Sampledata.Nematodes.Field,function(x) c(mean = mean(x), se = standard_error(x)))
aggregate(Observed.raw~Soil,data=Sampledata.Protists.Field,function(x) c(mean = mean(x), se = standard_error(x)))
aggregate(Shannon.raw~Soil,data=Sampledata.Nematodes.Field,function(x) c(mean = mean(x), se = standard_error(x)))
aggregate(Shannon.raw~Soil,data=Sampledata.Protists.Field,function(x) c(mean = mean(x), se = standard_error(x)))

#residuals
AOVObservedSoilNematodes<-aov(Observed.resid~Soil,data=Sampledata.Nematodes.Field)
AOVObservedSoilProtists<-aov(Observed.resid~Soil,data=Sampledata.Protists.Field)
AOVShannonSoilNematodes<-aov(Shannon.resid~Soil,data=Sampledata.Nematodes.Field)
AOVShannonSoilProtists<-aov(Shannon.resid~Soil,data=Sampledata.Protists.Field)

summary(AOVObservedSoilNematodes)
summary(AOVObservedSoilProtists)
summary(AOVShannonSoilNematodes)
summary(AOVShannonSoilProtists)

aggregate(Observed.resid~Soil,data=Sampledata.Nematodes.Field,function(x) c(mean = mean(x), se = standard_error(x)))
aggregate(Observed.resid~Soil,data=Sampledata.Protists.Field,function(x) c(mean = mean(x), se = standard_error(x)))
aggregate(Shannon.resid~Soil,data=Sampledata.Nematodes.Field,function(x) c(mean = mean(x), se = standard_error(x)))
aggregate(Shannon.resid~Soil,data=Sampledata.Protists.Field,function(x) c(mean = mean(x), se = standard_error(x)))

#rarefied data
AOVObservedSoilNematodes<-aov(Observed.rarmean~Soil*Management,data=Sampledata.Nematodes.Field)
AOVObservedSoilProtists<-aov(Observed.rarmean~Soil*Management,data=Sampledata.Protists.Field)
AOVShannonSoilNematodes<-aov(Shannon.rarmean~Soil*Management,data=Sampledata.Nematodes.Field)
AOVShannonSoilProtists<-aov(Shannon.rarmean~Soil*Management,data=Sampledata.Protists.Field)

summary(AOVObservedSoilNematodes)
summary(AOVObservedSoilProtists)
summary(AOVShannonSoilNematodes)
summary(AOVShannonSoilProtists)

t.test(Observed.rarmean~Soil,data=Sampledata.Nematodes.Field)
t.test(Observed.rarmean~Soil,data=Sampledata.Protists.Field)
t.test(Shannon.rarmean~Soil,data=Sampledata.Nematodes.Field)
t.test(Shannon.rarmean~Soil,data=Sampledata.Protists.Field)


t.test(Observed.rarmean~Management,data=Sampledata.Nematodes.Field)
t.test(Observed.rarmean~Management,data=Sampledata.Protists.Field)
t.test(Shannon.rarmean~Management,data=Sampledata.Nematodes.Field)
t.test(Shannon.rarmean~Management,data=Sampledata.Protists.Field)


t.test(Observed.rarmean~Management,data=Sampledata.Protists.Clay.Field)
t.test(Observed.rarmean~Management,data=Sampledata.Protists.Sand.Field)
t.test(Shannon.rarmean~Management,data=Sampledata.Protists.Clay.Field)
t.test(Shannon.rarmean~Management,data=Sampledata.Protists.Sand.Field)


aggregate(Observed.rarmean~Soil,data=Sampledata.Nematodes.Field,function(x) c(mean = mean(x), se = standard_error(x)))
aggregate(Observed.rarmean~Soil,data=Sampledata.Protists.Field,function(x) c(mean = mean(x), se = standard_error(x)))
aggregate(Shannon.rarmean~Soil,data=Sampledata.Nematodes.Field,function(x) c(mean = mean(x), se = standard_error(x)))
aggregate(Shannon.rarmean~Soil,data=Sampledata.Protists.Field,function(x) c(mean = mean(x), se = standard_error(x)))

aggregate(Observed.rarmean~Soil*Management,data=Sampledata.Nematodes.Field,function(x) c(mean = mean(x), se = standard_error(x)))
aggregate(Observed.rarmean~Soil*Management,data=Sampledata.Protists.Field,function(x) c(mean = mean(x), se = standard_error(x)))
aggregate(Shannon.rarmean~Soil*Management,data=Sampledata.Nematodes.Field,function(x) c(mean = mean(x), se = standard_error(x)))
aggregate(Shannon.rarmean~Soil*Management,data=Sampledata.Protists.Field,function(x) c(mean = mean(x), se = standard_error(x)))

#####field level composition differences between soil types#####
AOVPCoA1SoilNematodes<-aov(General_PCoA1~Soil,data=Sampledata.Nematodes.Field)
AOVPCoA1SoilProtists<-aov(General_PCoA1~Soil,data=Sampledata.Protists.Field)
AOVPCoA2SoilNematodes<-aov(General_PCoA2~Soil,data=Sampledata.Nematodes.Field)
AOVPCoA2SoilProtists<-aov(General_PCoA2~Soil,data=Sampledata.Protists.Field)

summary(AOVPCoA1SoilNematodes)
summary(AOVPCoA1SoilProtists)
summary(AOVPCoA2SoilNematodes)
summary(AOVPCoA2SoilProtists)


#####field level general reads and composition differences between regions#####
AOVPCoA1RegionNematodes<-aov(General_PCoA1~Region,data=Sampledata.Nematodes.Field)
AOVPCoA1RegionProtists<-aov(General_PCoA1~Region,data=Sampledata.Protists.Field)
AOVPCoA2RegionNematodes<-aov(General_PCoA2~Region,data=Sampledata.Nematodes.Field)
AOVPCoA2RegionProtists<-aov(General_PCoA2~Region,data=Sampledata.Protists.Field)

summary(AOVPCoA1RegionNematodes)
summary(AOVPCoA1RegionProtists)
summary(AOVPCoA2RegionNematodes)
summary(AOVPCoA2RegionProtists)

TukeyAOVPCoA1RegionNematodes<-TukeyHSD(AOVPCoA1RegionNematodes)
TukeyAOVPCoA1RegionProtists<-TukeyHSD(AOVPCoA1RegionProtists)
TukeyAOVPCoA2RegionNematodes<-TukeyHSD(AOVPCoA2RegionNematodes)
TukeyAOVPCoA2RegionProtists<-TukeyHSD(AOVPCoA2RegionProtists)

TukeyRegioncombipval<-as.data.frame(cbind(TukeyAOVPCoA1RegionNematodes$Region[,4],
                                          TukeyAOVPCoA1RegionProtists$Region[,4],
                                          TukeyAOVPCoA2RegionNematodes$Region[,4],
                                          TukeyAOVPCoA2RegionProtists$Region[,4]))
colnames(TukeyRegioncombipval)<-c("PCoA1Nematodes","PCoA1Protists","PCoA2Nematodes","PCoA2Protists")
format(TukeyRegioncombipval,scientific = FALSE)
#within clay 1/4 differs (pcoa2nem), within sand never > only random structure needed in clay


#####lmers observed richness#####

#test whether random intercept is needed
# M1.gls<-gls(Observed.rarmean~Management*Yearssince.as.if*pH.std*SOM.std,data=Sampledata.Nematodes.Clay.Field)
# M1.lme<-lme(Observed.rarmean~Management*Yearssince.as.if*pH.std*SOM.std,random= ~1|Region,method="REML", data=Sampledata.Nematodes.Clay.Field)
# anova(M1.gls,M1.lme)#p=0.99, so random structure not needed
# E2<-resid(M1.gls,type="normalized")
# F2 <- fitted(M1.gls)
# plot(x=F2,y=E2)
# 
# M1.gls<-gls(Observed.rarmean~Management*Yearssince.as.if*pH.std*SOM.std,data=Sampledata.Protists.Clay.Field)
# M1.lme<-lme(Observed.rarmean~Management*Yearssince.as.if*pH.std*SOM.std,random= ~1|Region,method="REML", data=Sampledata.Protists.Clay.Field)
# anova(M1.gls,M1.lme)#p=0.99, so random structure not needed
# E2<-resid(M1.gls,type="normalized")
# F2 <- fitted(M1.gls)
# plot(x=F2,y=E2)
# #etc Zuur pag. 134

LM3.Nematodes.Clay.Observed<-lm(Observed.rarmean~Management*Yearssince.as.if,data=Sampledata.Nematodes.Clay.Field)
LM3.Protists.Clay.Observed<-lm(Observed.rarmean~Management*Yearssince.as.if,data=Sampledata.Protists.Clay.Field)
LM3.Nematodes.Sand.Observed<-lm(Observed.rarmean~Management*Yearssince.as.if,data=Sampledata.Nematodes.Sand.Field)
LM3.Protists.Sand.Observed<-lm(Observed.rarmean~Management*Yearssince.as.if,data=Sampledata.Protists.Sand.Field)
anova(LM3.Nematodes.Clay.Observed)
anova(LM3.Protists.Clay.Observed)
anova(LM3.Nematodes.Sand.Observed)
anova(LM3.Protists.Sand.Observed)

LM3.Observed.All<-cbind(anova(LM3.Nematodes.Clay.Observed),
                         anova(LM3.Protists.Clay.Observed),
                         anova(LM3.Nematodes.Sand.Observed),
                         anova(LM3.Protists.Sand.Observed))

# plot(LME3.Nematodes.Clay.Observed)
# plot(LME3.Nematodes.Sand.Observed)
#
# res_LME3.Nematodes.Clay.Observed<-residuals(LME3.Nematodes.Clay.Observed)
# res_LME3.Nematodes.Sand.Observed<-residuals(LME3.Nematodes.Sand.Observed)
#
# ggqqplot(res_LME3.Nematodes.Clay.Observed)
# ggqqplot(res_LME3.Nematodes.Sand.Observed)

LM4.Nematodes.Clay.Observed<-lm(Observed.rarmean~Management*Yearssince.as.if*pH.std*SOM.std,data=Sampledata.Nematodes.Clay.Field)
LM4.Protists.Clay.Observed<-lm(Observed.rarmean~Management*Yearssince.as.if*pH.std*SOM.std,data=Sampledata.Protists.Clay.Field)
LM4.Nematodes.Sand.Observed<-lm(Observed.rarmean~Management*Yearssince.as.if*pH.std*SOM.std,data=Sampledata.Nematodes.Sand.Field)
LM4.Protists.Sand.Observed<-lm(Observed.rarmean~Management*Yearssince.as.if*pH.std*SOM.std,data=Sampledata.Protists.Sand.Field)

stepAIC(LM4.Nematodes.Clay.Observed)
stepAIC(LM4.Protists.Clay.Observed)
stepAIC(LM4.Nematodes.Sand.Observed)
stepAIC(LM4.Protists.Sand.Observed)

LM5.Nematodes.Clay.Observed<-lm(Observed.rarmean~Management + Yearssince.as.if + 
                                  pH.std + SOM.std + Management:Yearssince.as.if + Management:pH.std + 
                                  Yearssince.as.if:pH.std + Management:SOM.std + Yearssince.as.if:SOM.std + 
                                  pH.std:SOM.std + Management:Yearssince.as.if:pH.std + Management:Yearssince.as.if:SOM.std + 
                                  Yearssince.as.if:pH.std:SOM.std,
                                data=Sampledata.Nematodes.Clay.Field)
LM5.Protists.Clay.Observed<-lm(Observed.rarmean~Management + Yearssince.as.if + 
                                 pH.std + SOM.std + Management:pH.std + Yearssince.as.if:pH.std + 
                                 Yearssince.as.if:SOM.std + pH.std:SOM.std,
                               data=Sampledata.Protists.Clay.Field)
LM5.Nematodes.Sand.Observed<-lm(Observed.rarmean~Management * Yearssince.as.if * pH.std * SOM.std,data=Sampledata.Nematodes.Sand.Field)
LM5.Protists.Sand.Observed<-lm(Observed.rarmean~Management + Yearssince.as.if + 
                                 pH.std + SOM.std + Management:Yearssince.as.if + Management:pH.std + 
                                 Yearssince.as.if:pH.std + Management:SOM.std + Yearssince.as.if:SOM.std + 
                                 pH.std:SOM.std + Management:Yearssince.as.if:pH.std + Management:Yearssince.as.if:SOM.std + 
                                 Management:pH.std:SOM.std + Yearssince.as.if:pH.std:SOM.std,data=Sampledata.Protists.Sand.Field)

anova(LM5.Nematodes.Clay.Observed)
anova(LM5.Protists.Clay.Observed)
anova(LM5.Nematodes.Sand.Observed) 
anova(LM5.Protists.Sand.Observed) 

#####same for Shannon, PCOA1 and PCOA2#####
####Shannon#####
# M1.gls<-gls(Shannon.rarmean~Management*Yearssince.as.if,data=Sampledata.Nematodes.Clay.Field)
# M1.lme<-lme(Shannon.rarmean~Management*Yearssince.as.if,random= ~1|Region,method="REML", data=Sampledata.Nematodes.Clay.Field)
# anova(M1.gls,M1.lme)#p=0.99, so random structure not needed
# E2<-resid(M1.gls,type="normalized")
# F2 <- fitted(M1.gls)
# plot(x=F2,y=E2)
# M1.gls<-gls(Shannon.rarmean~Management*Yearssince.as.if,data=Sampledata.Protists.Clay.Field)
# M1.lme<-lme(Shannon.rarmean~Management*Yearssince.as.if,random= ~1|Region,method="REML", data=Sampledata.Protists.Clay.Field)
# anova(M1.gls,M1.lme)#p=0.99, so random structure not needed
# E2<-resid(M1.gls,type="normalized")
# F2 <- fitted(M1.gls)
# plot(x=F2,y=E2)

LM3.Nematodes.Clay.Shannon<-lm(Shannon.rarmean~Management*Yearssince.as.if,data=Sampledata.Nematodes.Clay.Field)
LM3.Protists.Clay.Shannon<-lm(Shannon.rarmean~Management*Yearssince.as.if,data=Sampledata.Protists.Clay.Field)
LM3.Nematodes.Sand.Shannon<-lm(Shannon.rarmean~Management*Yearssince.as.if,data=Sampledata.Nematodes.Sand.Field)
LM3.Protists.Sand.Shannon<-lm(Shannon.rarmean~Management*Yearssince.as.if,data=Sampledata.Protists.Sand.Field)
anova(LM3.Nematodes.Clay.Shannon)
anova(LM3.Protists.Clay.Shannon)
anova(LM3.Nematodes.Sand.Shannon)
anova(LM3.Protists.Sand.Shannon)

aggregate(Shannon.rarmean~Management,data=Sampledata.Nematodes.Clay.Field,FUN=mean)
aggregate(Shannon.rarmean~Management,data=Sampledata.Protists.Clay.Field,FUN=mean)
aggregate(Shannon.rarmean~Management,data=Sampledata.Nematodes.Sand.Field,FUN=mean)
aggregate(Shannon.rarmean~Management,data=Sampledata.Protists.Sand.Field,FUN=mean)

LM3.Shannon.All<-cbind(anova(LM3.Nematodes.Clay.Shannon),
                        anova(LM3.Protists.Clay.Shannon),
                        anova(LM3.Nematodes.Sand.Shannon),
                        anova(LM3.Protists.Sand.Shannon))

LM4.Nematodes.Clay.Shannon<-lm(Shannon.rarmean~Management*Yearssince.as.if*pH.std*SOM.std,data=Sampledata.Nematodes.Clay.Field)
LM4.Protists.Clay.Shannon<-lm(Shannon.rarmean~Management*Yearssince.as.if*pH.std*SOM.std,data=Sampledata.Protists.Clay.Field)
LM4.Nematodes.Sand.Shannon<-lm(Shannon.rarmean~Management*Yearssince.as.if*pH.std*SOM.std,data=Sampledata.Nematodes.Sand.Field)
LM4.Protists.Sand.Shannon<-lm(Shannon.rarmean~Management*Yearssince.as.if*pH.std*SOM.std,data=Sampledata.Protists.Sand.Field)

stepAIC(LM4.Nematodes.Clay.Shannon)
stepAIC(LM4.Protists.Clay.Shannon)
stepAIC(LM4.Nematodes.Sand.Shannon)
stepAIC(LM4.Protists.Sand.Shannon)

LM5.Nematodes.Clay.Shannon<-lm(Shannon.rarmean~Management + Yearssince.as.if + 
                                 pH.std + SOM.std + Management:Yearssince.as.if + Management:pH.std + 
                                 Yearssince.as.if:pH.std + Management:SOM.std + Yearssince.as.if:SOM.std + 
                                 pH.std:SOM.std + Management:Yearssince.as.if:pH.std + Management:Yearssince.as.if:SOM.std,data=Sampledata.Nematodes.Clay.Field)
LM5.Protists.Clay.Shannon<-lm(Shannon.rarmean~Management + Yearssince.as.if + 
                                pH.std + SOM.std + Management:pH.std + Management:SOM.std + 
                                Yearssince.as.if:SOM.std + pH.std:SOM.std + Management:pH.std:SOM.std,data=Sampledata.Protists.Clay.Field)
LM5.Nematodes.Sand.Shannon<-lm(Shannon.rarmean~Management + Yearssince.as.if + 
                                 pH.std + SOM.std + Management:Yearssince.as.if + Yearssince.as.if:pH.std + 
                                 Management:SOM.std + Yearssince.as.if:SOM.std + pH.std:SOM.std + 
                                 Management:Yearssince.as.if:SOM.std + Yearssince.as.if:pH.std:SOM.std,data=Sampledata.Nematodes.Sand.Field)
LM5.Protists.Sand.Shannon<-lm(Shannon.rarmean~Management + Yearssince.as.if + 
                                pH.std + SOM.std + Management:Yearssince.as.if + Management:pH.std + 
                                Yearssince.as.if:pH.std + Management:SOM.std + Yearssince.as.if:SOM.std + 
                                pH.std:SOM.std + Management:Yearssince.as.if:pH.std + Management:Yearssince.as.if:SOM.std + 
                                Management:pH.std:SOM.std + Yearssince.as.if:pH.std:SOM.std,data=Sampledata.Protists.Sand.Field)

anova(LM5.Nematodes.Clay.Shannon)
anova(LM5.Protists.Clay.Shannon)
anova(LM5.Nematodes.Sand.Shannon)
anova(LM5.Protists.Sand.Shannon)


#####PCOA1#####
#gls + REML + anova
# M1.gls<-gls(SoilSpecificPCoA1~Management*Yearssince.as.if,data=Sampledata.Nematodes.Clay.Field)
# M1.lme<-lme(SoilSpecificPCoA1~Management*Yearssince.as.if,random= ~1|Region,method="REML", data=Sampledata.Nematodes.Clay.Field)
# anova(M1.gls,M1.lme)#p=1, so random structure IS needed for nematodes!!
# E2<-resid(M1.lme,type="normalized")
# F2 <- fitted(M1.lme)
# plot(x=F2,y=E2)
# M1.gls<-gls(SoilSpecificPCoA1~Management*Yearssince.as.if,data=Sampledata.Protists.Clay.Field)
# M1.lme<-lme(SoilSpecificPCoA1~Management*Yearssince.as.if,random= ~1|Region,method="REML", data=Sampledata.Protists.Clay.Field)
# anova(M1.gls,M1.lme)#p=0.001, so random structure is not needed for protists
# E2<-resid(M1.gls,type="normalized")
# F2 <- fitted(M1.gls)
# plot(x=F2,y=E2)

LME3.Nematodes.Clay.PCoA1<-lme(SoilSpecificPCoA1~Management*Yearssince.as.if,random= ~1|Region,method="REML", data=Sampledata.Nematodes.Clay.Field)
anova(LME3.Nematodes.Clay.PCoA1,type="sequential")
summary(LME3.Nematodes.Clay.PCoA1)
LM3.Protists.Clay.PCoA1<-lm(SoilSpecificPCoA1~Management*Yearssince.as.if,data=Sampledata.Protists.Clay.Field)
anova(LM3.Protists.Clay.PCoA1)
LM3.Nematodes.Sand.PCoA1<-lm(SoilSpecificPCoA1~Management*Yearssince.as.if,data=Sampledata.Nematodes.Sand.Field)
anova(LM3.Nematodes.Sand.PCoA1)
LM3.Protists.Sand.PCoA1<-lm(SoilSpecificPCoA1~Management*Yearssince.as.if,data=Sampledata.Protists.Sand.Field)
anova(LM3.Protists.Sand.PCoA1)

LME4.Nematodes.Clay.PCoA1<-lme(SoilSpecificPCoA1~Management*Yearssince.as.if*pH.std*SOM.std,random= ~1|Region,method="ML", data=Sampledata.Nematodes.Clay.Field)
LM4.Protists.Clay.PCoA1<-lm(SoilSpecificPCoA1~Management*Yearssince.as.if*pH.std*SOM.std,data=Sampledata.Protists.Clay.Field)
LM4.Nematodes.Sand.PCoA1<-lm(SoilSpecificPCoA1~Management*Yearssince.as.if*pH.std*SOM.std,data=Sampledata.Nematodes.Sand.Field)
LM4.Protists.Sand.PCoA1<-lm(SoilSpecificPCoA1~Management*Yearssince.as.if*pH.std*SOM.std,data=Sampledata.Protists.Sand.Field)

stepAIC(LME4.Nematodes.Clay.PCoA1)
stepAIC(LM4.Protists.Clay.PCoA1)
stepAIC(LM4.Nematodes.Sand.PCoA1)
stepAIC(LM4.Protists.Sand.PCoA1)

LME5.Nematodes.Clay.PCoA1<-lme(SoilSpecificPCoA1~Management*Yearssince.as.if*pH.std*SOM.std,random= ~1|Region,method="REML", data=Sampledata.Nematodes.Clay.Field)
LM5.Protists.Clay.PCoA1<-lm(SoilSpecificPCoA1~Management + Yearssince.as.if + 
                              pH.std + SOM.std + Management:Yearssince.as.if + Management:pH.std + 
                              Yearssince.as.if:pH.std + Management:SOM.std + Yearssince.as.if:SOM.std + 
                              pH.std:SOM.std + Management:Yearssince.as.if:pH.std + Management:Yearssince.as.if:SOM.std,data=Sampledata.Protists.Clay.Field)
LM5.Nematodes.Sand.PCoA1<-lm(SoilSpecificPCoA1~Management*Yearssince.as.if*pH.std*SOM.std,data=Sampledata.Nematodes.Sand.Field)
LM5.Protists.Sand.PCoA1<-lm(SoilSpecificPCoA1~Management*Yearssince.as.if*pH.std*SOM.std,data=Sampledata.Protists.Sand.Field)

anova(LME5.Nematodes.Clay.PCoA1,type="sequential")
anova(LM5.Protists.Clay.PCoA1)
anova(LM5.Nematodes.Sand.PCoA1)
anova(LM5.Protists.Sand.PCoA1)

#####PCOA2#####
# M1.gls<-gls(SoilSpecificPCoA2~Management*Yearssince.as.if,data=Sampledata.Nematodes.Clay.Field)
# M1.lme<-lme(SoilSpecificPCoA2~Management*Yearssince.as.if,random= ~1|Region,method="REML", data=Sampledata.Nematodes.Clay.Field)
# anova(M1.gls,M1.lme)#p=1, so random structure is not needed!!
# E2<-resid(M1.gls,type="normalized")
# F2 <- fitted(M1.gls)
# plot(x=F2,y=E2)
# M1.gls<-gls(SoilSpecificPCoA2~Management*Yearssince.as.if,data=Sampledata.Protists.Clay.Field)
# M1.lme<-lme(SoilSpecificPCoA2~Management*Yearssince.as.if,random= ~1|Region,method="REML", data=Sampledata.Protists.Clay.Field)
# anova(M1.gls,M1.lme)#p=0.001, so random structure IS needed for protists!!
# E2<-resid(M1.lme,type="normalized")
# F2 <- fitted(M1.lme)
# plot(x=F2,y=E2)


LM3.Nematodes.Clay.PCoA2<-lm(SoilSpecificPCoA2~Management*Yearssince.as.if,data=Sampledata.Nematodes.Clay.Field)
LME3.Protists.Clay.PCoA2<-lme(SoilSpecificPCoA2~Management*Yearssince.as.if,random= ~1|Region,method="REML",data=Sampledata.Protists.Clay.Field)
LM3.Nematodes.Sand.PCoA2<-lm(SoilSpecificPCoA2~Management*Yearssince.as.if,data=Sampledata.Nematodes.Sand.Field)
LM3.Protists.Sand.PCoA2<-lm(SoilSpecificPCoA2~Management*Yearssince.as.if,data=Sampledata.Protists.Sand.Field)
anova(LM3.Nematodes.Clay.PCoA2)
anova(LME3.Protists.Clay.PCoA2)
anova(LM3.Nematodes.Sand.PCoA2)
anova(LM3.Protists.Sand.PCoA2)



LM4.Nematodes.Clay.PCoA2<-lm(SoilSpecificPCoA2~Management*Yearssince.as.if*pH.std*SOM.std,data=Sampledata.Nematodes.Clay.Field)
LME4.Protists.Clay.PCoA2<-lme(SoilSpecificPCoA2~Management*Yearssince.as.if*pH.std*SOM.std,random= ~1|Region,method="ML",data=Sampledata.Protists.Clay.Field)
LM4.Nematodes.Sand.PCoA2<-lm(SoilSpecificPCoA2~Management*Yearssince.as.if*pH.std*SOM.std,data=Sampledata.Nematodes.Sand.Field)
LM4.Protists.Sand.PCoA2<-lm(SoilSpecificPCoA2~Management*Yearssince.as.if*pH.std*SOM.std,data=Sampledata.Protists.Sand.Field)
stepAIC(LM4.Nematodes.Clay.PCoA2)
stepAIC(LME4.Protists.Clay.PCoA2)
stepAIC(LM4.Nematodes.Sand.PCoA2)
stepAIC(LM4.Protists.Sand.PCoA2)

LM5.Nematodes.Clay.PCoA2<-lm(SoilSpecificPCoA2~Management + Yearssince.as.if + pH.std +      
                               SOM.std + Management:pH.std + Management:SOM.std + Yearssince.as.if:SOM.std +      
                               pH.std:SOM.std + Management:pH.std:SOM.std ,data=Sampledata.Nematodes.Clay.Field)
LME5.Protists.Clay.PCoA2<-lme(SoilSpecificPCoA2~Management + Yearssince.as.if + pH.std +      
                                SOM.std + Management:pH.std + Management:SOM.std + Yearssince.as.if:SOM.std +      
                                pH.std:SOM.std + Management:pH.std:SOM.std,random= ~1|Region,method="REML",data=Sampledata.Protists.Clay.Field)
LM5.Nematodes.Sand.PCoA2<-lm(SoilSpecificPCoA2~Management*Yearssince.as.if*pH.std*SOM.std,data=Sampledata.Nematodes.Sand.Field)
LM5.Protists.Sand.PCoA2<-lm(SoilSpecificPCoA2~Management + Yearssince.as.if + 
                              pH.std + SOM.std + Management:Yearssince.as.if + Management:pH.std + 
                              Yearssince.as.if:pH.std + Management:SOM.std + Yearssince.as.if:SOM.std + 
                              pH.std:SOM.std + Management:Yearssince.as.if:pH.std + Management:Yearssince.as.if:SOM.std + 
                              Yearssince.as.if:pH.std:SOM.std,data=Sampledata.Protists.Sand.Field)
anova(LM5.Nematodes.Clay.PCoA2)
anova(LME5.Protists.Clay.PCoA2,type="sequential")
anova(LM5.Nematodes.Sand.PCoA2)
anova(LM5.Protists.Sand.PCoA2)
#####create text output file for all final models#####
LM3richdiv<-rbind(LM3.Observed.All,LM3.Shannon.All)
write.csv(LM3richdiv,"LM3richdiv.csv")

#sink(file = "MEMbasemicroorg_outputfieldlevel.txt")
anova(LME3.Nematodes.Clay.PCoA1,type="sequential")
anova(LM3.Protists.Clay.PCoA1)
anova(LM3.Nematodes.Sand.PCoA1)
anova(LM3.Protists.Sand.PCoA1)
anova(LM3.Nematodes.Clay.PCoA2)
anova(LME3.Protists.Clay.PCoA2,type="sequential")
anova(LM3.Nematodes.Sand.PCoA2)
anova(LM3.Protists.Sand.PCoA2)
#sink(file=NULL)

#sink(file = "MEMselectedmicroorg_outputfieldlevel.txt")
anova(LM5.Nematodes.Clay.Observed)
anova(LM5.Protists.Clay.Observed)
anova(LM5.Nematodes.Sand.Observed) 
anova(LM5.Protists.Sand.Observed) 
anova(LM5.Nematodes.Clay.Shannon)
anova(LM5.Protists.Clay.Shannon)
anova(LM5.Nematodes.Sand.Shannon)
anova(LM5.Protists.Sand.Shannon)
anova(LME5.Nematodes.Clay.PCoA1,type="sequential")
anova(LM5.Protists.Clay.PCoA1)
anova(LM5.Nematodes.Sand.PCoA1)
anova(LM5.Protists.Sand.PCoA1)
anova(LM5.Nematodes.Clay.PCoA2)
anova(LME5.Protists.Clay.PCoA2,type="sequential")
anova(LM5.Nematodes.Sand.PCoA2)
anova(LM5.Protists.Sand.PCoA2)
#sink(file=NULL)

#####multimodels without interactions#####
#nematodes clay
M6.gls<-gls(Observed.rarmean~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Nematodes.Clay.Field)
M6.lme<-lme(Observed.rarmean~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),random= ~1|Region,method="REML", data=Sampledata.Nematodes.Clay.Field)
anova(M6.gls,M6.lme)
M6.Nematodes.Clay.Observed<-lm(Observed.rarmean~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Nematodes.Clay.Field)
M7.Nematodes.Clay.Observed<-stepAIC(M6.Nematodes.Clay.Observed)
anova(M7.Nematodes.Clay.Observed)

M6.gls<-gls(Shannon.rarmean~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Nematodes.Clay.Field)
M6.lme<-lme(Shannon.rarmean~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),random= ~1|Region,method="REML", data=Sampledata.Nematodes.Clay.Field)
anova(M6.gls,M6.lme)
M6.Nematodes.Clay.Shannon<-lm(Shannon.rarmean~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Nematodes.Clay.Field)
M7.Nematodes.Clay.Shannon<-stepAIC(M6.Nematodes.Clay.Shannon)
anova(M7.Nematodes.Clay.Shannon)

M6.gls<-gls(SoilSpecificPCoA1~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Nematodes.Clay.Field)
M6.lme<-lme(SoilSpecificPCoA1~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),random= ~1|Region,method="REML", data=Sampledata.Nematodes.Clay.Field)
anova(M6.gls,M6.lme)
M6.Nematodes.Clay.PCoA1<-lme(SoilSpecificPCoA1~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),random= ~1|Region,method="ML", data=Sampledata.Nematodes.Clay.Field)
stepAIC(M6.Nematodes.Clay.PCoA1)
M7.Nematodes.Clay.PCoA1<-lme(SoilSpecificPCoA1~Management + Yearssince.as.if + pH.std + SOM.std + Management:SOM.std ,random= ~1|Region,method="REML", data=Sampledata.Nematodes.Clay.Field)
anova(M7.Nematodes.Clay.PCoA1,type="sequential")

M6.gls<-gls(SoilSpecificPCoA2~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Nematodes.Clay.Field)
M6.lme<-lme(SoilSpecificPCoA2~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),random= ~1|Region,method="REML", data=Sampledata.Nematodes.Clay.Field)
anova(M6.gls,M6.lme)
M6.Nematodes.Clay.PCoA2<-lm(SoilSpecificPCoA2~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Nematodes.Clay.Field)
M7.Nematodes.Clay.PCoA2<-stepAIC(M6.Nematodes.Clay.PCoA2)
anova(M7.Nematodes.Clay.PCoA2)


#nematodes sand
M6.gls<-gls(Observed.rarmean~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Nematodes.Sand.Field)
M6.lme<-lme(Observed.rarmean~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),random= ~1|Region,method="REML", data=Sampledata.Nematodes.Sand.Field)
anova(M6.gls,M6.lme)
M6.Nematodes.Sand.Observed<-lm(Observed.rarmean~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Nematodes.Sand.Field)
M7.Nematodes.Sand.Observed<-stepAIC(M6.Nematodes.Sand.Observed)
anova(M7.Nematodes.Sand.Observed)

M6.gls<-gls(Shannon.rarmean~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Nematodes.Sand.Field)
M6.lme<-lme(Shannon.rarmean~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),random= ~1|Region,method="REML", data=Sampledata.Nematodes.Sand.Field)
anova(M6.gls,M6.lme)
M6.Nematodes.Sand.Shannon<-lm(Shannon.rarmean~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Nematodes.Sand.Field)
M7.Nematodes.Sand.Shannon<-stepAIC(M6.Nematodes.Sand.Shannon)
anova(M7.Nematodes.Sand.Shannon)

M6.gls<-gls(SoilSpecificPCoA1~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Nematodes.Sand.Field)
M6.lme<-lme(SoilSpecificPCoA1~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),random= ~1|Region,method="REML", data=Sampledata.Nematodes.Sand.Field)
anova(M6.gls,M6.lme)
M6.Nematodes.Sand.PCoA1<-lm(SoilSpecificPCoA1~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Nematodes.Sand.Field)
M7.Nematodes.Sand.PCoA1<-stepAIC(M6.Nematodes.Sand.PCoA1)
anova(M7.Nematodes.Sand.PCoA1)

M6.gls<-gls(SoilSpecificPCoA2~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Nematodes.Sand.Field)
M6.lme<-lme(SoilSpecificPCoA2~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),random= ~1|Region,method="REML", data=Sampledata.Nematodes.Sand.Field)
anova(M6.gls,M6.lme)
M6.Nematodes.Sand.PCoA2<-lm(SoilSpecificPCoA2~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Nematodes.Sand.Field)
M7.Nematodes.Sand.PCoA2<-stepAIC(M6.Nematodes.Sand.PCoA2)
anova(M7.Nematodes.Sand.PCoA2)

#Protists clay
M6.gls<-gls(Observed.rarmean~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Protists.Clay.Field)
M6.lme<-lme(Observed.rarmean~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),random= ~1|Region,method="REML", data=Sampledata.Protists.Clay.Field)
anova(M6.gls,M6.lme)
M6.Protists.Clay.Observed<-lm(Observed.rarmean~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Protists.Clay.Field)
M7.Protists.Clay.Observed<-stepAIC(M6.Protists.Clay.Observed)
anova(M7.Protists.Clay.Observed)

M6.gls<-gls(Shannon.rarmean~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Protists.Clay.Field)
M6.lme<-lme(Shannon.rarmean~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),random= ~1|Region,method="REML", data=Sampledata.Protists.Clay.Field)
anova(M6.gls,M6.lme)
M6.Protists.Clay.Shannon<-lm(Shannon.rarmean~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Protists.Clay.Field)
M7.Protists.Clay.Shannon<-stepAIC(M6.Protists.Clay.Shannon)
anova(M7.Protists.Clay.Shannon)

M6.gls<-gls(SoilSpecificPCoA1~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Protists.Clay.Field)
M6.lme<-lme(SoilSpecificPCoA1~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),random= ~1|Region,method="REML", data=Sampledata.Protists.Clay.Field)
anova(M6.gls,M6.lme)
M6.Protists.Clay.PCoA1<-lm(SoilSpecificPCoA1~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Protists.Clay.Field)
M7.Protists.Clay.PCoA1<-stepAIC(M6.Protists.Clay.PCoA1)
anova(M7.Protists.Clay.PCoA1)

M6.gls<-gls(SoilSpecificPCoA2~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Protists.Clay.Field)
M6.lme<-lme(SoilSpecificPCoA2~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),random= ~1|Region,method="REML", data=Sampledata.Protists.Clay.Field)
anova(M6.gls,M6.lme)
M6.Protists.Clay.PCoA2<-lme(SoilSpecificPCoA2~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),random= ~1|Region,method="ML", data=Sampledata.Protists.Clay.Field)
stepAIC(M6.Protists.Clay.PCoA2)
M7.Protists.Clay.PCoA2<-lme(SoilSpecificPCoA2~Management + Yearssince.as.if + PLFABacteria.std +PLFAFungi.std ,random= ~1|Region,method="REML", data=Sampledata.Protists.Clay.Field)
anova(M7.Protists.Clay.PCoA2,type="sequential")


#Protists sand
M6.gls<-gls(Observed.rarmean~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Protists.Sand.Field)
M6.lme<-lme(Observed.rarmean~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),random= ~1|Region,method="REML", data=Sampledata.Protists.Sand.Field)
anova(M6.gls,M6.lme)
M6.Protists.Sand.Observed<-lm(Observed.rarmean~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Protists.Sand.Field)
M7.Protists.Sand.Observed<-stepAIC(M6.Protists.Sand.Observed)
anova(M7.Protists.Sand.Observed)

M6.gls<-gls(Shannon.rarmean~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Protists.Sand.Field)
M6.lme<-lme(Shannon.rarmean~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),random= ~1|Region,method="REML", data=Sampledata.Protists.Sand.Field)
anova(M6.gls,M6.lme)
M6.Protists.Sand.Shannon<-lm(Shannon.rarmean~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Protists.Sand.Field)
M7.Protists.Sand.Shannon<-stepAIC(M6.Protists.Sand.Shannon)
anova(M7.Protists.Sand.Shannon)

M6.gls<-gls(SoilSpecificPCoA1~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Protists.Sand.Field)
M6.lme<-lme(SoilSpecificPCoA1~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),random= ~1|Region,method="REML", data=Sampledata.Protists.Sand.Field)
anova(M6.gls,M6.lme)
M6.Protists.Sand.PCoA1<-lm(SoilSpecificPCoA1~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Protists.Sand.Field)
M7.Protists.Sand.PCoA1<-stepAIC(M6.Protists.Sand.PCoA1)
anova(M7.Protists.Sand.PCoA1)

M6.gls<-gls(SoilSpecificPCoA2~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Protists.Sand.Field)
M6.lme<-lme(SoilSpecificPCoA2~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),random= ~1|Region,method="REML", data=Sampledata.Protists.Sand.Field)
anova(M6.gls,M6.lme)
M6.Protists.Sand.PCoA2<-lm(SoilSpecificPCoA2~Management*(Yearssince.as.if+pH.std+SOM.std+PLFABacteria.std+PLFAFungi.std),data=Sampledata.Protists.Sand.Field)
M7.Protists.Sand.PCoA2<-stepAIC(M6.Protists.Sand.PCoA2)
anova(M7.Protists.Sand.PCoA2)



#####NEW plots field level data#####
Sampledata.Nematodes.Field.sub<-subset(Sampledata.Nematodes.Field,select=c("Fieldnr","General_PCoA1","General_PCoA2","Group", "Management" ,                     "Pair","Region","Soil","Soil_Management","SoilSpecificPCoA1","SoilSpecificPCoA2" ))

Sampledata.Protists.Field.sub<-subset(Sampledata.Protists.Field,select=c("Fieldnr","General_PCoA1","General_PCoA2","Group", "Management" ,                     "Pair","Region","Soil","Soil_Management","SoilSpecificPCoA1","SoilSpecificPCoA2" ))



genpcoa<-dplyr::bind_rows(Sampledata.Nematodes.Field.sub,Sampledata.Protists.Field.sub)

genpcoa$Group<-factor(genpcoa$Group,levels=c("Nematodes","Protists"))
genpcoa$Region<-factor(genpcoa$Region,levels=c("Flevoland","Zeeland","DrentheOverijssel","Gelderland","BrabantLimburg"))

#first PCoA run on original data! not separate for clay and sand!!! (So also statistics will be different table!)
ggplot(genpcoa,aes(x=General_PCoA1,y=General_PCoA2,color=Soil,shape=Management))+
  geom_point(size=2)+
  #xlim(-1.5,1.5)+
  facet_wrap(.~Group,scales="free")+
  scale_color_manual(values = c("#336666","#CC9966"))+
  scale_shape_manual(values = c(21,16))+
  theme_bw()#+
#theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("pcoa1pcoa2.pdf",width = 9, height = 5)

ggplot(genpcoa,aes(x=General_PCoA1,y=General_PCoA2,color=Region))+
  geom_point(size=2)+
  #xlim(-1.5,1.5)+
  facet_wrap(.~Group,scales="free")+
  scale_color_manual(values = wes_palette("Cavalcanti1"))+
  theme_bw()#+
#theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("pcoa1pcoa2region.pdf",width = 9, height = 5)


#####FIGURES
datafigures<-rbind(Sampledata.Nematodes.Clay.Field,Sampledata.Nematodes.Sand.Field,
                   Sampledata.Protists.Clay.Field,Sampledata.Protists.Sand.Field)

datafigures$Group<-factor(datafigures$Group,levels=c("Nematodes","Protists"))





ggplot(datafigures,aes(x=Yearssince.as.if,y=SoilSpecificPCoA1,color=Soil,shape=Management,linetype=Management))+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  xlim(0,25)+
  facet_wrap(Soil~Group,scales="free")+
  scale_color_manual(values = c("black","black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  labs(x="Years of organic management",y='PCoA1') +
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("PCoA1time.pdf",width = 9, height = 9.3)

ggplot(datafigures,aes(x=Yearssince.as.if,y=SoilSpecificPCoA2,color=Soil,shape=Management,linetype=Management))+
  geom_point(size=2)+
  xlim(0,25)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(Soil~Group,scales="free")+
  scale_color_manual(values = c("black","black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  labs(x="Years of organic management",y='PCoA2') +
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("PCoA2time.pdf",width = 9, height = 9.3)




ggplot(datafigures,aes(x=percSOM,y=SoilSpecificPCoA1,color=Soil,shape=Management,linetype=Management))+
  geom_point(size=2)+
  xlim(0,8)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(Soil~Group,scales="free")+
  scale_color_manual(values = c("black","black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  labs(x="Soil Organic Matter (%)",y='PCoA1') +
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("PCoA1som.pdf",width = 9, height = 9.3)
ggplot(datafigures,aes(x=percSOM,y=SoilSpecificPCoA2,color=Soil,shape=Management,linetype=Management))+
  geom_point(size=2)+
  xlim(0,8)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(Soil~Group,scales="free")+
  scale_color_manual(values = c("black","black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  labs(x="Soil Organic Matter (%)",y='PCoA2') +
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("PCoA2som.pdf",width = 9, height = 9.3)

ggplot(datafigures,aes(x=pH,y=SoilSpecificPCoA1,color=Soil,shape=Management,linetype=Management))+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(Soil~Group,scales="free")+
  scale_color_manual(values = c("black","black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  labs(x="pH",y='PCoA1') +
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("PCoA1ph.pdf",width = 9, height = 9.3)
ggplot(datafigures,aes(x=pH,y=SoilSpecificPCoA2,color=Soil,shape=Management,linetype=Management))+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(Soil~Group,scales="free")+
  scale_color_manual(values = c("black","black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  labs(x="pH",y='PCoA2') +
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("PCoA2ph.pdf",width = 9, height = 9.3)



# #####new plots black&white#####



ggplot(datafigures,aes(x=Management,y=Observed.rarmean,fill=Soil_Management))+
  geom_boxplot()+
  geom_jitter(position=position_jitter(width=.1, height=0))+
  facet_wrap(Soil~Group,scales="free")+
  scale_fill_manual(values = c("white","#336666","white","#CC9966"))+
  labs(x="Management",y='Observed ASV richness') +
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("richnessboxnemprot.pdf",width = 9, height = 9.3)
ggplot(datafigures,aes(x=Management,y=Shannon.rarmean,fill=Soil_Management))+
  geom_boxplot()+
  geom_jitter(position=position_jitter(width=.1, height=0))+
  facet_wrap(Soil~Group,scales="free")+
  scale_fill_manual(values = c("white","#336666","white","#CC9966"))+
  labs(x="Management",y="Shannon's diversity index") +
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("shannonboxnemprot.pdf",width = 9, height = 9.3)

ggplot(datafigures,aes(x=Yearssince.as.if,y=Observed.rarmean,shape=Management,linetype=Management,color="black"))+
  geom_point(size=2)+
  xlim(0,25)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(Soil~Group,scales="free")+
  scale_color_manual(values = c("black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  labs(x="Years of organic management",y='Observed ASV richness') +
  theme_bw()+
  theme(legend.position = "none",
        axis.title=element_text(size=20,face="bold"),
        axis.text=element_text(size=15),
        strip.text=element_blank()
  )

ggsave("richnesstimenemprot.pdf",width = 9, height = 9.3)
ggplot(datafigures,aes(x=Yearssince.as.if,y=Shannon.rarmean,shape=Management,linetype=Management,color="black"))+
  geom_point(size=2)+
  xlim(0,25)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(Soil~Group,scales="free")+
  scale_color_manual(values = c("black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  labs(x="Years of organic management",y="Shannon's diversity index") +
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("shannontimenemprot.pdf",width = 9, height = 9.3)

