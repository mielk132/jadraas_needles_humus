### Target: Ecology Letters / New Phyt
### Prelim Title: Gadgil effect and mycorrhizal traits in pine needles and mor layer

###Authors: L A Mielke 1, J Klein 2, B Lindahl 3, R Finlay 1, A Ekblad 4, K E Clemmensen 1
### Affiliations
### 1 Dept Forest Mycology Plant Pathology - SLU Uppsala
### 2 Dept Soil & Environment - SLU Uppsala
### 3 Artdatabanken - SLU Uppsala
### 4 Örebro University

### Corresponding author
### Louis Mielke - louis.mielke@slu.se

###Jädraås Decomposition Community Data Setup

rm(list=ls()) #clear the workspace (only if you want to)

library(dplyr);library(funrar)

setwd("~/Projects/mycorrhizal removal/MeshBags/Jädraås_Gadgil_MycorrhizalTraits")


### data import
metadata <- read.delim("metadata.txt")
community <- read.table("community.txt", sep="\t", header=T, row.names = 1)
taxa <- read.table("taxa.txt", sep="\t", header=T, row.names=1) #declare NAs?


#standardize the community data by copies DNA per g substrate
community <- as.data.frame(round(make_relative(as.matrix(community))*metadata$copies_DNA_per_g_substrate))

##### MOCK ####
### subset mock community
mock <- community[metadata$treatment == "MOCK",]

##### BACKGROUND #####

######## subset background
meta_bkgrnd_humus <- metadata[metadata$treatment == "B" & metadata$substrate == "Humus",]
meta_bkgrnd_litter <- metadata[metadata$treatment == "B" & metadata$substrate == "Litter",]

bkgrnd_humus <- community[metadata$treatment == "B" & metadata$substrate == "Humus",]
bkgrnd_litter <- community[metadata$treatment == "B" & metadata$substrate == "Litter",]

##### average background copies per g substrate
h_bkgrnd_avg_copies <- as.vector(colMeans(bkgrnd_humus[-9,]))
l_bkgrnd_avg_copies <- as.vector(colMeans(bkgrnd_litter))

##### HUMUS & LITTER SUBSET ######
#### subset with samples that were incubated in experimental treatments by substrate
humusOTU <- community %>%
  filter(metadata$treatment %in% c("C","E","T","TE","DC") & metadata$substrate %in% "Humus")

humus_meta <- metadata %>%
  filter(metadata$treatment %in% c("C","E","T","TE","DC") & metadata$substrate %in% "Humus")

litterOTU <- community %>%
  filter(metadata$treatment %in% c("C","E","T","TE","DC") & metadata$substrate %in% "Litter")

litter_meta <- metadata %>%
  filter(metadata$treatment %in% c("C","E","T","TE","DC") & metadata$substrate %in% "Litter")

#subtract background from incubated samples by substrate (negative values = 0)
humusOTUnew <- humusOTU - h_bkgrnd_avg_copies
litterOTUnew <- litterOTU - l_bkgrnd_avg_copies
humusOTUnew[humusOTUnew<0] <- 0
litterOTUnew[litterOTUnew<0] <- 0

str(litterOTUnew)

litter_meta$Sample
row.names(litterOTUnew)


## pine needle litter

#hypothesis 1
litter_meta$copies_sap_whiterot <- rowSums(litterOTUnew[,taxa$guild == "saprotrophs white rot"])
litter_meta$copies_sap_other <- rowSums(litterOTUnew[,taxa$guild == "saprotrophs other"])

litter_meta$copies_fungi <- rowSums(litterOTUnew[,taxa$Kingdom == "Fungi"])

litter_meta$copies_mycena_clavicularis <- rowSums(as.matrix(litterOTUnew[,taxa$Taxon=="Mycena_clavicularis"]))
litter_meta$copies_ecto <- rowSums(as.matrix(litterOTUnew[,taxa$guild=="ectomycorrhizas"]))
litter_meta$copies_ectowr <- rowSums(as.matrix(litterOTUnew[,taxa$guild=="ectomycorrhizas white rot"]))
sum(litter_meta$copies_ectowr[litter_meta$pine == "1"])/sum(litter_meta$copies_fungi[litter_meta$pine == "1"]) #amount of white rot ectos

## humus


humus_meta$copies_fungi <- rowSums(humusOTUnew[,taxa$Kingdom == "Fungi"])

#hypothesis 1
humus_meta$copies_ecto <- rowSums(humusOTUnew[,taxa$emf == "1"])
humus_meta$copies_ericoid <- rowSums(humusOTUnew[,taxa$guild == "ericoid mycorrhizas"])

#exploration 2
humus_meta$copies_sap_whiterot <- rowSums(humusOTUnew[,taxa$guild == "saprotrophs white rot"])
humus_meta$copies_ecto_whiterot <- rowSums(humusOTUnew[,taxa$guild == "ectomycorrhizas white rot"]) 

mean(humus_meta$percent_mass_remaining[humus_meta$incubation == "5" & humus_meta$set == "B"] - humus_meta$percent_mass_remaining[humus_meta$incubation == "17" & humus_meta$set == "B"], na.rm =TRUE)

sd(humus_meta$percent_mass_remaining[humus_meta$incubation == "5" & humus_meta$set == "B"] - humus_meta$percent_mass_remaining[humus_meta$incubation == "17" & humus_meta$set == "B"], na.rm =TRUE)

mean(humus_meta$percent_mass_remaining[humus_meta$incubation == "5" & humus_meta$set == "B" & humus_meta$pine == "1" & humus_meta$shrub == "0"] - humus_meta$percent_mass_remaining[humus_meta$incubation == "17" & humus_meta$set == "B"& humus_meta$pine == "1" & humus_meta$shrub == "0"], na.rm =TRUE)

sd(humus_meta$percent_mass_remaining[humus_meta$incubation == "5" & humus_meta$set == "B" & humus_meta$pine == "1" & humus_meta$shrub == "0"] - humus_meta$percent_mass_remaining[humus_meta$incubation == "17" & humus_meta$set == "B"& humus_meta$pine == "1" & humus_meta$shrub == "0"], na.rm =TRUE)

#other interesting groups

humus_meta$copies_mycena <- rowSums(humusOTUnew[,taxa$Genus== "Mycena"])
humus_meta$copies_cortinarius <- rowSums(humusOTUnew[,taxa$Genus== "Cortinarius"])
humus_meta$copies_sap_other <- rowSums(humusOTUnew[,taxa$guild == "saprotrophs other"])
humus_meta$copies_moulds_yeasts <- rowSums(humusOTUnew[,taxa$guild == "moulds and yeasts"])
humus_meta$copies_other_root <- rowSums(humusOTUnew[,taxa$guild == "other root-associates"])
humus_meta$copies_ericoid_ecto <- rowSums(humusOTUnew[,taxa$guild == "ericoid- ectomycorrhizas"])
humus_meta$copies_piloderma <- rowSums(as.matrix(humusOTUnew[,taxa$Taxon== "Piloderma_sphaerosporum"]))
humus_meta$copies_mycena_clavicularis <- rowSums(as.matrix(humusOTUnew[,taxa$Taxon=="Mycena_clavicularis"]))
humus_meta$copies_whiterot <- rowSums(as.matrix(humusOTUnew[,taxa$Decay_type=="white rot"]))
humus_meta$copies_suillus <- rowSums(humusOTUnew[,taxa$Genus== "Suillus"])

##
write.csv(humus_meta, "humus_meta.csv")
write.csv(humusOTUnew, "humusOTU.csv")
write.csv(litter_meta, "litter_meta.csv")
write.csv(litterOTUnew, "litterOTU.csv")


