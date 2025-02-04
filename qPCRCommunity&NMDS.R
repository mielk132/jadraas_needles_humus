###Jädraås Decomposition Community Composition

rm(list=ls()) #clear the workspace (only if you want to)

library(tidyverse)
library(tibble)
library(broom)
library(funrar) #for fxn make_relative()
library(gridExtra)
library(ggpubr)
library(nlme) #?
library(betapart)
library(vegan)

setwd("~/Projects/mycorrhizal removal/MeshBags/Jädraås_Gadgil_MycorrhizalTraits")

set.seed(2022)

#### FORMAT #####

### import metadata
metadata <- read.table("metadata.txt", sep="\t", header=T, row.names=1)
metadata$set_time_trtmt_block <- paste0(metadata$set,metadata$incubation,metadata$treatment,metadata$block)
class(metadata);dim(metadata);str(metadata)
metadata$set <- factor(metadata$set, levels = c("A", "B"))

#filter out duplicate samples, negative controls, mock, and samples that were not recollected or broken
meta <- metadata[metadata$Unique_Sample != "negativecontrol",] #filter out pcr neg control
meta1 <- meta[meta$treatment != "MOCK",] #filter out mock
meta2 <- meta1[meta1$treatment != "Extra",]  # filter out samples that were sequenced twice
meta3 <- meta2[meta2$Unique_Sample != "E4Humus5B",]  #broken / lost sample
meta4 <- meta3[meta3$Unique_Sample != "E3Humus17A",]  #broken / lost sample
meta5 <- meta4[meta4$Unique_Sample != "TELitter5A",]  #broken / lost sample

meta.data =as.data.frame(meta5)
meta.data$treatment <- factor(meta.data$treatment, levels = c("B","C", "E", "T","TE","DC"))
meta.data$plot <- paste0(meta.data$treatment,meta.data$block)
meta.data$incubation <- fct_recode(as.factor(meta.data$incubation),"0" = "0", "1" = "5", "2" = "17")
str(meta.data) # 329 30

#OTU
community <- read.table("community.txt", sep="\t", header=T)
community <- data.frame(community, row.names = 1) #in one line set column 1 as row names
#str(community);community [1:5,1:5];class(community);dim(community)

community <- community[metadata$treatment != "LAB_NC",] #filter out pcr neg control
community <- community[meta$treatment != "MOCK",] #filter out mock
community <- community[meta1$treatment != "Extra",] # filter out samples that were repeated
community <- community[meta2$Unique_Sample != "E4Humus5B",] #lost sample
community <- community[meta3$Unique_Sample != "E3Humus17A",] #lost sample
community <- community[meta4$Unique_Sample != "TELitter5A",] #lost sample

OTU <- t(community)
dim(OTU) #807 329

### import taxonomy table and functional groups
taxa <- read.table("taxa.txt", sep="\t", header=T, row.names=1)
taxa$guild <- factor(taxa$guild, levels = c("ectomycorrhizas white rot","ectomycorrhizas other", "ericoid mycorrhizas","ericoid- ectomycorrhizas","other root-associates","moulds and yeasts","saprotrophs white rot","saprotrophs other","unknown","plants"))
#unique(taxa$guild);class(taxa);str(taxa);levels(taxa$guild);unique(taxa$Order);unique(taxa$Subphylum)
dim(taxa) #807 21

# add number of fungal copies from qPCR
meta.data$copies_fungi <- meta.data$copies_DNA_per_g_substrate*colSums(OTU[taxa$Kingdom == "Fungi",])/ colSums(OTU) 


##### FILTER ####
# filtering for fungi only and making relative abundance
taxa_trimd <-  taxa[rowSums(OTU) >= 1,] #prune taxa (row) with sums less than 1 read
OTU_trimd <-  OTU[rowSums(OTU) >= 1,] #prune OTU (row) with sums less than 1 read
dim(OTU_trimd) #797 329
dim(taxa_trimd) #797 21

plants <- OTU_trimd[taxa_trimd$Kingdom =="Plantae",]
sum(rowSums(plants)) / sum(rowSums(OTU))

##### PLOT ####
#graph the rank abundance curve and reads per sample
# readsumsdf = data.frame(nreads = sort(rowSums(OTU_trimd), TRUE), sorted = 1:length(OTU_trimd[,1]), type = "OTUs")
# 
# readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(OTU_trimd,TRUE),sorted = 1:length(meta.data$Unique_Sample), type = "Samples"))
#               title = "Total number of reads"
#               p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) +
#               geom_bar(stat = "identity")
#               p + ggtitle(title) + facet_wrap(~type, 1, scales = "free") #+ scale_y_log10() #to make the y axis log scale

##### FILTER 2 ####         
#to keep OTUs with at least 1% in a sample
fungi <- OTU_trimd[taxa_trimd$Kingdom =="Fungi",]
taxa_fungi <- taxa_trimd[taxa_trimd$Kingdom =="Fungi",]

sum(rowSums(fungi)) / sum(rowSums(OTU_trimd)) #98.5% fungi

taxaf_trimd <-  taxa_fungi[rowSums(fungi) >= 1,] #prune taxa (row) with sums less than 1 read
fungi_trimd <-  fungi[rowSums(fungi) >= 1,] #prune OTU (row) with sums less than 1 read

#make OTU matrix with read counts into relative abundance per sample
fungi_rbund = t(make_relative(t(fungi_trimd)))

#filter out taxa that have less than 1% in every sample but one
filtered_object <- fungi_rbund[rowSums(fungi_rbund >= 0.01) > 0, ]
filtered_taxaf <- taxaf_trimd[rowSums(fungi_rbund >= 0.01) > 0, ]
fungi_rbund = t(make_relative(t(filtered_object))) # re-make relative abundance with updated OTU table

dim(fungi_rbund) #190 329
dim(filtered_taxaf) #190 21


###### MERGE ####

meta.data <- rownames_to_column(meta.data, var = "Sample_ID") #make sure column is shared with OTU table
dim(meta.data) #329 33

#combing taxa and otu table
data <- full_join(rownames_to_column(filtered_taxaf),rownames_to_column(data.frame(fungi_rbund)))
data_long <- pivot_longer(data,cols = (length(filtered_taxaf[1,])+2):length(data[1,]), #long format
                          names_to = "Sample_ID",values_to = "abundance") %>%
            full_join(.,meta.data,by = "Sample_ID") #joining with meta.data by sample

dim(data_long) #62510   55

####### STANDARDIZE ######
#multiply averages from qPCR technical replicates by realtive abundance per sample
data_long$std_fungal_copies = data_long$copies_fungi*data_long$abundance 

######### PINE NEEDLE LITTER SUMMARIZING ########

#average qPCR abundance by plot and time point over blocks
littercopies <- meta.data %>%
  filter(.,substrate == "Litter", treatment != "B")%>% 
  group_by(treatment, set, incubation) %>% 
  summarize_at(c("copies_fungi"), funs(mean = mean, n(), se = sd(.)/sqrt(n()))) %>% 
  ungroup()

#add new variable to order x axis by set (A/B) incubation (5/17) and treatment
littercopies$time_trtmt <- paste0(littercopies$set,littercopies$incubation,littercopies$treatment)

### 

#by guild
litter.guild <- data_long %>%
  group_by(substrate,guild,treatment, block, set, incubation, ) %>% 
  filter(.,substrate == "Litter", treatment != "B")%>% 
  summarize_at(c("abundance","std_fungal_copies"), funs(sum = sum))
litterguild <- litter.guild %>% 
  group_by(guild,treatment, set, incubation) %>% 
  summarize_at(c("std_fungal_copies_sum"), funs(mean = mean)) %>%
  ungroup()

#add new variable to order x axis by set (A/B) incubation (5/17) and treatment
litterguild$time_trtmt <- paste0(litterguild$set,litterguild$incubation,litterguild$treatment)
litter.guild$time_trtmt <- paste0(litter.guild$set,litter.guild$incubation,litter.guild$treatment,litter.guild$block)

#by genera
littergenera <- data_long %>%
  group_by(substrate,Genus,treatment, block, set, incubation) %>% 
  filter(.,substrate == "Litter", treatment != "B")%>% 
  summarize_at(c("std_fungal_copies"), funs(sum = sum)) %>%
  group_by(Genus,treatment, set, incubation) %>% 
  summarize_at(c("sum"), funs(mean = mean)) %>%
  ungroup()

#add new variable to order x axis by set (A/B) incubation (5/17) and treatment
littergenera$time_trtmt <- paste0(littergenera$set,littergenera$incubation,littergenera$treatment)


######### HUMUS SUMMARIZING #########

#average plot and time over blocks
humuscopies <- meta.data %>%
  filter(.,substrate == "Humus" & treatment != "B" )%>% 
  group_by(treatment, set, incubation) %>% 
  summarize_at(c("copies_fungi"), funs(mean = mean, n(), se = sd(.)/sqrt(n()))) %>% 
  ungroup()

#add new variable to order x axis by set (A/B) incubation (5/17) and treatment
humuscopies$time_trtmt <- paste0(humuscopies$set,humuscopies$incubation,humuscopies$treatment)

##
#by guild
humus.guild <- data_long %>%
  filter(.,substrate == "Humus" & treatment != "B")%>% 
  group_by(guild,treatment, block, set,incubation) %>% 
  summarize_at(c("abundance","std_fungal_copies"), funs(sum = sum))
humusguild <- humus.guild %>%
  group_by(guild, treatment, set, incubation) %>% 
  summarize_at(c("std_fungal_copies_sum"), funs(mean = mean)) %>%
  ungroup()


#add new variable to order x axis by set (A/B) incubation (5/17) and treatment
humus.guild$time_trtmt <- paste0(humus.guild$set,humus.guild$incubation,humus.guild$treatment,humus.guild$block)
humusguild$time_trtmt <- paste0(humusguild$set,humusguild$incubation,humusguild$treatment)

##### PLOT FIGURE 2A & B #####

col_guild <- c("#88CCAA", "#44AA77","#777711","#AAAA44","gray70", "gray50", "#AA7744", "#774411","pink")

####### plot fungal guilds in needles multiplied by qPCR averaged over blocks

  guild_needles_plot_qpcr <- ggplot(data = litterguild, aes(x = time_trtmt, y = mean,  fill = guild)) +
    geom_bar(position = "stack", stat="identity", color = "black") +
    geom_errorbar(data =littercopies, aes(x = time_trtmt, ymin=mean-se, ymax = mean+se),inherit.aes = FALSE,width=.2, position=position_dodge(.9))+
    theme_classic()+
    xlim("A1C","A1E","A1T", "A1TE", "A1DC","A2C","A2E","A2T","A2TE","A2DC","B1C","B1E","B1T","B1TE","B1DC","B2C","B2E","B2T","B2TE","B2DC")+
    ylab(expression(paste("\nFungal biomass (ITS copies g ",DW^-1,")")))+
  scale_fill_manual(values = col_guild)+
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "top",
          legend.key.size = unit(1, 'cm'),
          axis.line = element_line( size = 1, linetype = "solid"), 
          axis.title=element_text(size=16, color = "black"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14, color = "black"))+
    scale_y_continuous(limits = ~ c(min(.x), max(.x) + max(.x)/10), #add a bit of space at the top for std error
                       expand = (c(0,0)),labels = function(x) format(x, scientific = TRUE))
 
  guild_needles_plot_qpcr


  ####### plot fungal guilds in humus multiplied by qPCR averaged over blocks
  
  guild_humus_plot_qpcr <- ggplot(data = humusguild, aes(x = time_trtmt, y = mean,  fill = guild)) +
    geom_bar(position = "stack", stat="identity", color = "black") +
    geom_errorbar(data =humuscopies, aes(x = time_trtmt, ymin=mean-se, ymax = mean+se),inherit.aes = FALSE,width=.2, position=position_dodge(.9))+
    theme_classic()+
    scale_fill_manual(values = col_guild)+
    ylab(expression(paste("\nFungal biomass (ITS copies g ",DW^-1,")")))+
    xlim("A1C","A1E","A1T", "A1TE", "A1DC","A2C","A2E","A2T","A2TE","A2DC","B1C","B1E","B1T","B1TE","B1DC","B2C","B2E","B2T","B2TE","B2DC")+
    theme(
          plot.title = element_text(hjust = 0.5),
          legend.position = "top",
          axis.line = element_line( size = 1, linetype = "solid"), 
          axis.title=element_text(size=16, color = "black"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 14, color = "black", angle = 90),
          axis.text.y = element_text(size = 14, color = "black"))+
    scale_y_continuous(limits = ~ c(min(.x), max(.x) + max(.x)/10), expand = (c(0,0)),labels = function(x) format(x, scientific = TRUE))

  guild_humus_plot_qpcr
  
  guild_humus_plot_qpcr_replicates <- ggplot(data = humus.guild, aes(x = time_trtmt, y = std_fungal_copies_sum,  fill = guild)) +
    geom_bar(position = "stack", stat="identity", color = "black") +
    theme_classic()+
    #scale_fill_manual(values = col_guild)+
    ylab(expression(paste("\nFungal biomass (ITS copies g ",DW^-1,")")))+
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "top",
      axis.line = element_line( size = 1, linetype = "solid"), 
      axis.title=element_text(size=16, color = "black"),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 14, color = "black", angle = 270),
      axis.text.y = element_text(size = 14, color = "black"))+
    scale_y_continuous(limits = ~ c(min(.x), max(.x) + max(.x)/10), expand = (c(0,0)),labels = function(x) format(x, scientific = TRUE))
  guild_humus_plot_qpcr_replicates

##### COMBINED PLOTS 


pdf("figures/Rplot_Figure3_Humus_Needle_Composition_qPCR.pdf", width = 10, height = 10)
comboplot <- ggarrange(nrow = 2,guild_needles_plot_qpcr,guild_humus_plot_qpcr, labels = c("a)","b)"),font.label = list(size = 18, color = "black"),common.legend = TRUE,vjust = 1.2)
comboplot
dev.off()


########## PLOTS FIGURE 4A, B, C, D, E, F

##### PLOT FIGURE 4A, B, C, D #####
######### saps in pine needles ########

  litterwhiterot_sum <- data_long %>%
    filter(.,substrate == "Litter" & guild == "saprotrophs white rot"  & treatment != "B" & set == "B" & incubation == 2)%>% 
    group_by(Genus,treatment, block) %>% 
    summarize_at(c("std_fungal_copies"), funs(sum = sum))
  
  litsum_wrGenus <-  litterwhiterot_sum %>%
    group_by(Genus,treatment) %>%
    summarize_at(c("sum"), funs(mean = mean, se = sd(.)/sqrt(n()))) %>%
    ungroup()
  
  litterwr_trmt_mean_se <-   litterwhiterot_sum %>%
    group_by(treatment,block) %>%
    summarize_at(c("sum"), funs(sum = sum)) %>%
                   group_by(treatment) %>%
                   summarize_at(c("sum"), funs(mean = mean, se = sd(.)/sqrt(n()))) %>%
    ungroup()
  
### treatment averages by white rot saps in needles
  
  #plot fungal guilds in humus
  litter_wr_Genus_plotB17 <- ggplot(data =   litsum_wrGenus, aes(x = treatment, y = mean,  fill = Genus)) +
    geom_bar(position = "stack", stat="identity", color = "black") +
    geom_errorbar(data = litterwr_trmt_mean_se, aes(x = treatment, ymin = mean - se,ymax = mean+se),inherit.aes = FALSE,width=.2, position=position_dodge(.9))+
    theme_classic()+
    xlab("\n")+
    ylab(expression(paste("\nFungal biomass (ITS copies g ",DW^-1,")")))+
    theme(plot.title = element_text(hjust = 0.5),
          #legend.position = c("none"),
          axis.line = element_line( size = 1, linetype = "solid"), 
          axis.title=element_text(size=16, color = "black"),
          axis.text.x = element_text(size = 14, color = "black"),
          axis.text.y = element_text(size = 14, color = "black"))+
    scale_x_discrete(labels = c("+Shrubs\n+Pine","-Shrubs\n+Pine","+Shrubs\n-Pine","-Shrubs\n-Pine","Disturbed\ncontrol"))+
    scale_y_continuous(limits = ~ c(min(.x), max(.x) + max(.x)/10),expand = (c(0,0)),
                       labels = function(x) format(x, scientific = TRUE))+
    scale_fill_manual(values = c("#D0BCD1", "#BD5C02", "#A86E37", "#694F06", "#D1A908", "#D4D4D4", "#E4E8CF", "#F5E889"))
  litter_wr_Genus_plotB17
  
  ### now other saps
  litter_sum <- data_long %>%
    filter(.,substrate == "Litter" & guild == "saprotrophs other"  & treatment != "B" & set == "B" & incubation == 2)%>% 
    group_by(Genus,treatment, block) %>% 
    summarize_at(c("std_fungal_copies"), funs(sum = sum))
  
  litsum_Genus <-  litter_sum %>%
    group_by(Genus,treatment) %>%
    summarize_at(c("sum"), funs(mean = mean, se = sd(.)/sqrt(n()))) %>%
    ungroup()
  
  litter_trmt_mean_se <-   litter_sum %>%
    group_by(treatment,block) %>%
    summarize_at(c("sum"), funs(sum = sum)) %>%
    group_by(treatment) %>%
    summarize_at(c("sum"), funs(mean = mean, se = sd(.)/sqrt(n()))) %>%
    ungroup()

  
  pdf("figures/Rplot_SuppFig3_OtherSap_Needle_Composition_qPCR.pdf", width = 10, height = 10)
  
  #plot other sap genera in litter
  litter_othersap_Genus_plotB17 <- ggplot(data =   litsum_Genus, aes(x = treatment, y = mean,  fill = Genus)) +
    geom_bar(position = "stack", stat="identity", color = "black") +
    geom_errorbar(data = litter_trmt_mean_se, aes(x = treatment, ymin = mean - se,ymax = mean+se),inherit.aes = FALSE,width=.2, position=position_dodge(.9))+
    theme_classic()+
    xlab("\n")+
    ylab(expression(paste("\nFungal biomass (ITS copies g ",DW^-1,")")))+
    theme(plot.title = element_text(hjust = 0.5),
          #legend.position = c("none"),
          axis.line = element_line( size = 1, linetype = "solid"), 
          axis.title=element_text(size=16, color = "black"),
          axis.text.x = element_text(size = 14, color = "black"),
          axis.text.y = element_text(size = 14, color = "black"))+
    scale_x_discrete(labels = c("+Shrubs\n+Pine","-Shrubs\n+Pine","+Shrubs\n-Pine","-Shrubs\n-Pine","Disturbed\ncontrol"))+
    scale_y_continuous(limits = ~ c(min(.x), max(.x) + max(.x)/10),expand = (c(0,0)),
                       labels = function(x) format(x, scientific = TRUE))+
    scale_fill_manual(values = c("#D0BCD1", "#BD5C02", "#A86E37", "#698F06", "#D1A908", "grey70", "#E4E8CF", "#F5E889","#D1BCD1", "#BD6C02", "#A91E37", "#695F06", "#D1A918", "#D4D4D5", "#E4E8CF", "#F5E889","#D0BCD1", "#BD7C02", "#A88E37", "#696F06", "#D1A928", "#D4D4D6", "#E4E8CF", "#F5E889","#D0BCD1", "#BD8C02", "#A89E37", "#697F06", "#D1A938", "#D4D4D7", "#E4E8CF", "#F5E889","#D0BCD1", "#BD9C02", "#A90E37", "#6B5151F1", "#D1A948", "#D4D4D8"))
  litter_othersap_Genus_plotB17
  
  dev.off()
 
   {
    
    guild_needles <- ggplot(data = litter.guild, aes(x = time_trtmt, y = abundance_sum,  fill = guild)) +
      geom_bar(position = "stack", stat="identity", color = "black") +
      #geom_errorbar(data =littercopies, aes(x = time_trtmt, ymin=mean-se, ymax = mean+se),inherit.aes = FALSE,width=.2, position=position_dodge(.9))+
      theme_classic()+
      ylab("\nITS2 copies per g pine needles\n")+
      theme(panel.border= element_rect(colour = "black", size=1, fill=NA),
            panel.background = element_rect(colour = "black", size=1, fill=NA),
            plot.title = element_text(hjust = 0.5),
            legend.position = c("right"),
            axis.line = element_line( size = 1, linetype = "solid"), 
            axis.title=element_text(size=16, color = "black"),
            axis.title.x = element_blank(),
            axis.text.x = element_text(size = 6, angle = 90),
            axis.text.y = element_text(size = 14, color = "black"))+
      scale_y_continuous(limits = ~ c(min(.x), max(.x)), #add a bit of space at the top for std error
                         expand = (c(0,0)),labels = function(x) format(x, scientific = TRUE))
    #+scale_fill_manual(values = col_guild)
    guild_needles
    
    pdf("figures/Rplot_Figure2_Needles_Composition_qPCR.pdf", width = 11, height = 11)
    guild_needles
    dev.off()
    
  } ####plot just relative abundance of fungal guilds in needles for all plots
  
######### ecto in humus ########

#col_guild <- c("#331F05", "#5C4721", "#8F6D5B", "#999176", "#CCCCCC", "grey70","#F7EDED")

  #RELATIVE ABUNDANCE OF Ectomycorrhizal SPECIES PER TREATMENT B17
  
  humusEMFGenus_sum <- data_long %>%
    filter(.,substrate == "Humus" & emf == "1" & treatment != "B" & set == "B" & incubation == 2)%>% 
    group_by(Genus,treatment, block) %>% 
    summarize_at(c("std_fungal_copies"), funs(sum = sum))
  
  humusEMFGenus_mean <-  humusEMFGenus_sum  %>%
    group_by(Genus,treatment) %>%
    summarize_at(c("sum"), funs(mean = mean, se = sd(.)/sqrt(n()))) %>%
    ungroup()
  
  humusEMF_trmt_mean_se <- humusEMFGenus_sum %>%
    group_by(treatment,block) %>%
    summarize_at(c("sum"), funs(sum = sum)) %>%
    group_by(treatment) %>%
    summarize_at(c("sum"), funs(mean = mean, se = sd(.)/sqrt(n()))) %>%
    ungroup()
  
  #plot fungal guilds in humus
  humus_EMFGenera_plotB17 <- ggplot(data = humusEMFGenus_mean, aes(x = treatment, y = mean,  fill = Genus)) +
    geom_bar(position = "stack", stat="identity", color = "black") +
    geom_errorbar(data = humusEMF_trmt_mean_se, aes(x = treatment, ymin = mean - se,ymax = mean+se),inherit.aes = FALSE,width=.2, position=position_dodge(.9))+
    theme_classic()+
    xlab("\n")+
    ylab(expression(paste("\nFungal biomass (ITS copies g ",DW^-1,")")))+
    theme(
          plot.title = element_text(hjust = 0.5),
          #legend.position = c("none"),
          axis.line = element_line( size = 1, linetype = "solid"), 
          axis.title=element_text(size=16, color = "black"),
          axis.text.x = element_text(size = 14, color = "black"),
          axis.text.y = element_text(size = 14, color = "black"))+
    scale_x_discrete(labels = c("+Shrubs\n+Pine","-Shrubs\n+Pine","+Shrubs\n-Pine","-Shrubs\n-Pine","Disturbed\ncontrol"))+
    scale_y_continuous(limits = ~ c(min(.x), max(.x) + max(.x)/10),expand = (c(0,0)),
                       labels = function(x) format(x, scientific = TRUE))+
    # scale_fill_manual(values = c("#1B0854", "#941843", "#3C8A6B", "#E0C009", "#96D158", "#E69D4A", "#ED6B74", "#94DCE0", "#BE7CCF", "#7F7E85"))
    scale_fill_manual(values = c("#E1E4E6", "#244F80", "#6174B3", "#71A5C7", "#8EC4D1", "#48C2A5", "#86BAB2", "#B9C8FA"))
  humus_EMFGenera_plotB17

 
  ##################  ericoid, ericoid-ectos and other root associated fungi in humus ########
  
  humusRootGenus_sum <- data_long %>%
    filter(.,substrate == "Humus" & guild == "ericoid mycorrhizas" & treatment != "B" & set == "B" & incubation == 2) %>% 
    group_by(Genus,treatment, block) %>% 
    summarize_at(c("std_fungal_copies"), funs(sum = sum))
  
  humusRootGenus_mean <-   humusRootGenus_sum  %>%
    group_by(Genus,treatment) %>%
    summarize_at(c("sum"), funs(mean = mean, se = sd(.)/sqrt(n()))) %>%
    ungroup()
  
  humusRoot_trmt_mean_se <- humusRootGenus_sum %>%
    group_by(treatment,block) %>%
    summarize_at(c("sum"), funs(sum = sum)) %>%
    group_by(treatment) %>%
    summarize_at(c("sum"), funs(mean = mean, se = sd(.)/sqrt(n()))) %>%
    ungroup()
  
  #plot fungal guilds in humus
  humus_RootGenera_plotB17 <- ggplot(data = humusRootGenus_mean, aes(x = treatment, y = mean,  fill = Genus)) +
    geom_bar(position = "stack", stat="identity", color = "black") +
    geom_errorbar(data = humusRoot_trmt_mean_se, aes(x = treatment, ymin = mean - se,ymax = mean+se),inherit.aes = FALSE,width=.2, position=position_dodge(.9))+
    theme_classic()+
    xlab("\n")+
    ylab(expression(paste("\nFungal biomass (ITS copies g ",DW^-1,")")))+
    theme(
          plot.title = element_text(hjust = 0.5),
          #legend.position = c("none"),
          axis.line = element_line( size = 1, linetype = "solid"), 
          axis.title=element_text(size=16, color = "black"),
          axis.text.x = element_text(size = 14, color = "black"),
          axis.text.y = element_text(size = 14, color = "black"))+
    scale_x_discrete(labels = c("+Shrubs\n+Pine","-Shrubs\n+Pine","+Shrubs\n-Pine","-Shrubs\n-Pine","Disturbed\ncontrol"))+
    scale_y_continuous(limits = ~ c(min(.x), max(.x) + max(.x)/10),expand = (c(0,0)),
                       labels = function(x) format(x, scientific = TRUE))+
    scale_fill_manual(values = c("#CAD1CF", "#0A590E", "#497A00", "#54A800", "#86D179", "#88A300", "#AED900", "#C7D900", "#EBCA11", "#D1C371", "#B39E00", "#737309", "#3D4506"))
  humus_RootGenera_plotB17 

  
  ##### COMBINED PLOTS 
  
  
  pdf("figures/Rplot_Figure4ABC_Humus_Needle_GUILDComposition_qPCR.pdf", width = 12, height = 40)
  comboplot <- ggarrange(nrow = 1,ncol = 3,litter_wr_Genus_plotB17, humus_EMFGenera_plotB17,humus_RootGenera_plotB17 , labels = c("a)","b)","c)"),font.label = list(size = 18, color = "black"),vjust = 1.2)
  comboplot
  dev.off()
  
  
  
  ## COMBINE #####
  humus_RootGenera_plotB17 
  litter_wr_Genus_plotB17
  humus_EMFGenera_plotB17
  
#### PLOT by PLOT ####
  humusRootGenus_sum$plot <- paste0(humusRootGenus_sum$treatment,humusRootGenus_sum$block)
  #plot fungal guilds in humus
  humus_RootGenera_plot <- ggplot(data = humusRootGenus_sum, aes(x = plot , y = sum,  fill = Genus)) +
    geom_bar(position = "stack", stat="identity", color = "black") +
    #geom_errorbar(data = humusRoot_trmt_mean_se, aes(x = treatment, ymin = mean - se,ymax = mean+se),inherit.aes = FALSE,width=.2, position=position_dodge(.9))+
    theme_classic()+
    xlab("\n")+
    ylab(expression(paste("\nFungal biomass (ITS copies g ",DW^-1,")")))+
    theme(
      plot.title = element_text(hjust = 0.5),
      #legend.position = c("none"),
      axis.line = element_line( size = 1, linetype = "solid"), 
      axis.title=element_text(size=16, color = "black"),
      axis.text.x = element_text(size = 14, color = "black"),
      axis.text.y = element_text(size = 14, color = "black"))+
    scale_x_discrete(labels = c("+Shrubs\n+Pine","-Shrubs\n+Pine","+Shrubs\n-Pine","-Shrubs\n-Pine","Disturbed\ncontrol"))+
    scale_y_continuous(limits = ~ c(min(.x), max(.x) + max(.x)/10),expand = (c(0,0)),
                       labels = function(x) format(x, scientific = TRUE))+
    scale_fill_manual(values = c("#CAD1CF", "#0A590E", "#497A00", "#54A800", "#86D179", "#88A300", "#AED900", "#C7D900", "#EBCA11", "#D1C371", "#B39E00", "#737309", "#3D4506"))
  humus_RootGenera_plotB17 
  
  

######### BACKGROUND ######

#average
bkgrndguild_avg <- data_long %>%
  filter(.,treatment == "B")%>% 
  group_by(substrate,guild,label) %>% 
  summarize_at(c("std_fungal_copies","DNA_ug_per_g_substrate"), funs(sum = sum)) %>%
  group_by(substrate,guild) %>% 
  summarize_at(c("std_fungal_copies_sum","DNA_ug_per_g_substrate_sum"), funs(mean = mean, se = sd(.)/sqrt(n()))) %>%
  mutate(across(2:3,function(x) format(x, scientific = TRUE))) %>%
  ungroup()

write.csv(bkgrndguild_avg,"tables/bkgrndguild_avg2.csv")

#look at mean and std error of humus and litter biomass proxy
bkgrnd_tbl <- data_long %>%
  filter(.,treatment == "B")%>% 
  group_by(substrate)%>% 
  summarize_at(c("std_fungal_copies","DNA_ug_per_g_substrate"), funs(mean = mean, se = sd(.)/sqrt(n()))) %>%
  ungroup()
bkgrnd_tbl

dna <- data_long %>%
  filter(.,treatment != "B")%>% 
  group_by(substrate)%>% 
  summarize_at(c("std_fungal_copies","DNA_ug_per_g_substrate"), funs(mean = mean, se = sd(.)/sqrt(n()))) %>%
  ungroup()
dna

dna <- data_long %>%
  filter(.,treatment != "B" & incubation == "2" & set == "B")%>% 
  group_by(substrate)%>% 
  summarize_at(c("std_fungal_copies","DNA_ug_per_g_substrate"), funs(mean = mean, se = sd(.)/sqrt(n()))) %>%
  ungroup()
dna
#Humus 2788 ± 470 copies per g substrate
#Litter 2295856 ± 423299 copies per g substrate

##### PLOT #####

#by guild
bkgrndguildL <- data_long %>%
  filter(.,treatment == "B", substrate == "Litter")%>% 
  group_by(guild,label) %>% 
  summarize_at(c("std_fungal_copies"), funs(sum = sum)) %>%
  ungroup()

#plot fungal guilds in needles
p_guildbkgrndL <- ggplot(data = bkgrndguildL, aes(x = label, y = sum,  fill = guild)) +
  geom_bar(position = "stack", stat="identity") +
  theme_classic()+
  xlab("replicate")+
  ggtitle("mean = 2.30*10^6 ± 4.23*10^5 copies")+
  ylab("\nITS2 copies per g pine needles\n")+
  theme(panel.border= element_rect(colour = "black", size=1, fill=NA),
        panel.background = element_rect(colour = "black", size=1, fill=NA),
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = c(0.2,.7),
        axis.line = element_line( size = 1, linetype = "solid"), 
        axis.title=element_text(size=16, color = "black"),
        axis.text.x = element_text(size = 16, color = "black", hjust = 0.2, vjust = 0.95),
        axis.text.y = element_text(size = 14, color = "black"))+
  scale_y_continuous(limits = ~ c(min(.x), max(.x) + 2*10^8), #add a bit of space at the top for std error
                     expand = (c(0,0)),labels = function(x) format(x, scientific = TRUE))+
  scale_fill_manual(values = col_guild)
  p_guildbkgrndL

#by guild
bkgrndguildH <- data_long %>%
  filter(.,treatment == "B", substrate == "Humus")%>% 
  group_by(guild,label) %>% 
  summarize_at(c("std_fungal_copies"), funs(sum = sum)) %>%
  ungroup()

p_guildbkgrndH <- ggplot(data = bkgrndguildH, aes(x = label, y = sum,  fill = guild)) +
  geom_bar(position = "stack", stat="identity") +
  theme_classic()+
  ylab("\nITS2 copies per g humus\n")+
  xlab("replicate")+
  ggtitle("mean = 2.79*10^3 ± 4.70*10^2 copies")+
  theme(panel.border= element_rect(colour = "black", size=1, fill=NA),
        panel.background = element_rect(colour = "black", size=1, fill=NA),
        legend.position = c(0.2,.7),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.line = element_line( size = 1, linetype = "solid"), 
        axis.title=element_text(size=16, color = "black"),
        axis.text.x = element_text(size = 16, color = "black", hjust = 0.2, vjust = 0.95),
        axis.text.y = element_text(size = 14, color = "black"))+
  scale_y_continuous(limits = ~ c(min(.x), max(.x)+10^5), #add a bit of space at the top for std error
                     expand = (c(0,0)),labels = function(x) format(x, scientific = TRUE))+
scale_fill_manual(values = col_guild)
p_guildbkgrndH

#pdf("figures/Rplot_Supplementary3_Background_Composition_qPCR.pdf", width = 11, height = 8)
bplot <- ggarrange(p_guildbkgrndH,p_guildbkgrndL, common.legend = TRUE, legend = "right", labels = c("A","B"))
annotate_figure(bplot,bottom = text_grob("Background Sample Composition\n"))
dev.off()

##### STATISTICS #######

##### pine needles ######

#### SET A ####
needle_statA5 <-  filter(meta.data, substrate == "Litter", treatment != "B", treatment != "DC", incubation == "1", set == "A")

lme <- lme(sqrt(copies_fungi) ~ pine*shrubs, random =  ~ 1|block, data = needle_statA5)
needle_qPCR_anovaA5 <- anova(lme)
needle_qPCR_anovaA5$coef <- t(coef(lme))[,1]
qqnorm(lme)

####### 17 month incubation #######
needle_statA17 <-  filter(meta.data, substrate == "Litter", treatment != "B", treatment != "DC", incubation == "2", set == "A")

lme <- lme(sqrt(copies_fungi) ~ pine*shrubs, random =  ~ 1|block, data = needle_statA17)
needle_qPCR_anovaA17 <- anova(lme)
needle_qPCR_anovaA17$coef <- t(coef(lme))[,1]
qqnorm(lme)

needle.guildA17 <- data_long %>%
  filter(.,substrate == "Litter" & treatment != "B" & treatment != "DC" & incubation == "2" & set == "A")%>% 
  group_by(treatment,pine,shrubs,guild,block) %>% 
  summarize_at(c("abundance","std_fungal_copies"), funs(sum = sum)) %>%
  ungroup()

#regression
regressionListA17 <- lapply(unique(needle.guildA17$guild),function(x) anova(lme(sqrt(std_fungal_copies_sum) ~ pine*shrubs, random = ~ 1|block, needle.guildA17[needle.guildA17$guild==x,])))
regressionListA17

qqnormPlots <- lapply(unique(needle.guildA17$guild),function(x) qqnorm(lme(sqrt(std_fungal_copies_sum) ~ pine*shrubs, random = ~ 1|block, needle.guildA17[needle.guildA17$guild==x,]), main = paste0(x)))
qqnormPlots

qqnormPlots <- lapply(unique(needle.guildA17$guild),function(x) qqnorm(lme(log(std_fungal_copies_sum+1) ~ pine*shrubs, random = ~ 1|block, needle.guildA17[needle.guildA17$guild==x,]), main = paste0(x)))
qqnormPlots

coefListA17 <- lapply(unique(needle.guildA17$guild), function(x) lme(sqrt(std_fungal_copies_sum) ~ pine*shrubs, random = ~ 1|block, needle.guildA17[needle.guildA17$guild==x,])$coefficients)

#ecto groups did not satisfy homogeneity of variance / normal distribution

#mean
needle.guildA17_ecto_wr_trtmt<- needle.guildA17 %>%
  filter(.,guild == "ectomycorrhizas white rot") %>%
  group_by(pine)%>%
  summarize_at(c("abundance_sum","std_fungal_copies_sum"), funs(mean = mean)) %>%
  ungroup()

needle.guildA17_ecto_trtmt <- needle.guildA17 %>%
  filter(.,guild == "ectomycorrhizas other") %>%
  group_by(pine)%>%
  summarize_at(c("abundance_sum","std_fungal_copies_sum"), funs(mean = mean)) %>%
  ungroup()

needle.guildA17_ecto_trtmt <- needle.guildA17 %>%
  filter(.,guild == "ectomycorrhizas other") %>%
  group_by(shrubs)%>%
  summarize_at(c("abundance_sum","std_fungal_copies_sum"), funs(mean = mean)) %>%
  ungroup()

# t.tests
needle.guildA17_ecto_wr_pine<- needle.guildA17 %>%
  filter(.,guild == "ectomycorrhizas white rot" & pine == "1")
needle.guildA17_ecto_wr_nopine<- needle.guildA17 %>%
  filter(.,guild == "ectomycorrhizas white rot" & pine == "0")

t_wr_pine <- t.test(needle.guildA17_ecto_wr_pine$std_fungal_copies_sum, needle.guildA17_ecto_wr_nopine$std_fungal_copies_sum)

needle.guildA17_ecto_wr_shrub<- needle.guildA17 %>%
  filter(.,guild == "ectomycorrhizas white rot" & pine == "1" & shrubs == "1")
needle.guildA17_ecto_wr_noshrub<- needle.guildA17 %>%
  filter(.,guild == "ectomycorrhizas white rot" & pine == "1" & shrubs == "0")

t_wr_shrub<- t.test(needle.guildA17_ecto_wr_shrub$std_fungal_copies_sum, needle.guildA17_ecto_wr_noshrub$std_fungal_copies_sum)


needle.guildA17_ecto_pine<- needle.guildA17 %>%
  filter(.,guild == "ectomycorrhizas other" & pine == "1")
needle.guildA17_ecto_nopine<- needle.guildA17 %>%
  filter(.,guild == "ectomycorrhizas other" & pine == "0")

t_wr_pine <- t.test(needle.guildA17_ecto_wr_pine$std_fungal_copies_sum, needle.guildA17_ecto_nopine$std_fungal_copies_sum)

needle.guildA17_ecto_shrub<- needle.guildA17 %>%
  filter(.,guild == "ectomycorrhizas other" & pine == "1" & shrubs == "1")
needle.guildA17_ecto_noshrub<- needle.guildA17 %>%
  filter(.,guild == "ectomycorrhizas other" & pine == "1" & shrubs == "0")

t_shrub<- t.test(needle.guildA17_ecto_wr_shrub$std_fungal_copies_sum, needle.guildA17_ecto_wr_noshrub$std_fungal_copies_sum)


##
capture.output(needle_qPCR_anovaA17,unique(needle.guildA17$guild),regressionListA17,coefListA17,t_wr,t_ecto,t_ecto_pine ,t_wr_pine, file = "needle_guild_qpcr_stat_A17.txt")

#### SET B ######
needle_statB5 <-  filter(meta.data, substrate == "Litter", treatment != "B", treatment != "DC", incubation == "1", set == "B")

lme <- lme(sqrt(copies_fungi) ~ pine*shrubs, random =  ~ 1|block, data = needle_statB5)
needle_qPCR_anovaB5 <- anova(lme)
needle_qPCR_anovaB5$coef <- t(coef(lme))[,1]
qqnorm(lme)

####### 17 month incubation
needle_statB17 <-  filter(meta.data, substrate == "Litter", treatment != "B", treatment != "DC", incubation == "2", set == "B")

#lme
lme <- lme(sqrt(copies_fungi) ~ pine*shrubs, random =  ~ 1|block, data = needle_statB17)
needle_qPCR_anovaB17 <- anova(lme)
needle_qPCR_anovaB17$coef <- t(coef(lme))[,1]
qqnorm(lme)

#summarize
needle.guildB17 <- data_long %>%
  filter(.,substrate == "Litter" & treatment != "B" & treatment != "DC" & incubation == "2" & set == "B" ,incubation == "2")%>% 
  group_by(treatment,pine,shrubs,guild,block) %>% 
  summarize_at(c("abundance","std_fungal_copies"), funs(sum = sum)) %>%
  ungroup()

#sqrt
regressionListB17_sqrt <- lapply(unique(needle.guildB17$guild),function(x) anova(lme(sqrt(std_fungal_copies_sum) ~ pine*shrubs, random = ~ 1|block, needle.guildB17[needle.guildB17$guild==x,])))
regressionListB17_sqrt 

qqnormPlotssqrt <- lapply(unique(needle.guildB17$guild),function(x) qqnorm(lme(sqrt(std_fungal_copies_sum) ~ pine*shrubs, random = ~ 1|block, needle.guildB17[needle.guildB17$guild==x,]), main = paste0(x)))
qqnormPlotssqrt

coefListB17_sqrt <- lapply(unique(needle.guildB17$guild), function(x) lme(sqrt(std_fungal_copies_sum) ~ pine*shrubs, random = ~ 1|block, needle.guildB17[needle.guildB17$guild==x,])$coefficients)

#log
regressionListB17_log1 <- lapply(unique(needle.guildB17$guild),function(x) anova(lme(log(std_fungal_copies_sum+1) ~ pine*shrubs, random = ~ 1|block, needle.guildB17[needle.guildB17$guild==x,])))
regressionListB17_log1

qqnormPlotslog <- lapply(unique(needle.guildB17$guild),function(x) qqnorm(lme(log(std_fungal_copies_sum+1) ~ pine*shrubs, random = ~ 1|block, needle.guildB17[needle.guildB17$guild==x,]), main = paste0(x)))
qqnormPlotslog

coefListB17_log1 <- lapply(unique(needle.guildB17$guild), function(x) lme(log(std_fungal_copies_sum+1) ~ pine*shrubs, random = ~ 1|block, needle.guildB17[needle.guildB17$guild==x,])$coefficients)

#ecto groups did not satisfy homogeneity of variance / normal distribution
# then I try with t.tests
needle.guildB17_ecto_wr_pine<- needle.guildB17 %>%
  filter(.,guild == "ectomycorrhizas white rot" & pine == "1")
needle.guildB17_ecto_wr_nopine<- needle.guildB17 %>%
  filter(.,guild == "ectomycorrhizas white rot" & pine == "0")

t_wr_pine <- t.test(needle.guildB17_ecto_wr_pine$std_fungal_copies_sum, needle.guildB17_ecto_wr_nopine$std_fungal_copies_sum)

needle.guildB17_ecto_wr_shrub<- needle.guildB17 %>%
  filter(.,guild == "ectomycorrhizas white rot" & pine == "1" & shrubs == "1")
needle.guildB17_ecto_wr_noshrub<- needle.guildB17 %>%
  filter(.,guild == "ectomycorrhizas white rot" & pine == "1" & shrubs == "0")

t_wr_shrub<- t.test(needle.guildB17_ecto_wr_shrub$std_fungal_copies_sum, needle.guildB17_ecto_wr_noshrub$std_fungal_copies_sum)


needle.guildB17_ecto_pine<- needle.guildB17 %>%
  filter(.,guild == "ectomycorrhizas other" & pine == "1")
needle.guildB17_ecto_nopine<- needle.guildB17 %>%
  filter(.,guild == "ectomycorrhizas other" & pine == "0")

t_wr_pine <- t.test(needle.guildB17_ecto_wr_pine$std_fungal_copies_sum, needle.guildB17_ecto_nopine$std_fungal_copies_sum)

needle.guildB17_ecto_shrub<- needle.guildB17 %>%
  filter(.,guild == "ectomycorrhizas other" & pine == "1" & shrubs == "1")
needle.guildB17_ecto_noshrub<- needle.guildB17 %>%
  filter(.,guild == "ectomycorrhizas other" & pine == "1" & shrubs == "0")

t_shrub<- t.test(needle.guildB17_ecto_wr_shrub$std_fungal_copies_sum, needle.guildB17_ecto_wr_noshrub$std_fungal_copies_sum)


capture.output(needle_qPCR_anovaB17,unique(needle.guildB17$guild),regressionListB17_sqrt,coefListB17_sqrt, regressionListB17_log1,coefListB17_log1,t_ecto,t_wr,t_wr_pine,t_ecto_pine,file = "needle_guild_qpcr_stat_B17.txt")

#means
needle.guildB17_ecto_pine <- needle.guildB17 %>%
  filter(.,guild == "ectomycorrhizas other") %>%
  group_by(pine)%>%
  summarize_at(c("abundance_sum","std_fungal_copies_sum"), funs(mean = mean)) %>%
  ungroup()

needle.guildB17_ecto_trtmt <- needle.guildB17 %>%
  filter(.,guild == "ectomycorrhizas other") %>%
  group_by(treatment)%>%
  summarize_at(c("abundance_sum","std_fungal_copies_sum"), funs(mean = mean)) %>%
  ungroup()

#### HUMUS ####
#### SET A 17 ######

#not very meaningful since there is so little biomass overall in this set and incubation
humus_statA17 <-  filter(meta.data, substrate == "Humus", treatment != "B", treatment != "DC", set == "A", incubation == "2")
dim(humus_statA17)

lme <- lme(log(copies_fungi+1) ~ pine*shrubs, random = ~1|block, data = humus_statA17)
humus_qPCR_anovaA17 <- anova(lme)
humus_qPCR_anovaA17$coef <- t(coef(lme))[,1]
qqnorm(lme)

humus.guildA17 <- data_long %>%
  filter(.,substrate == "Humus" & treatment != "B" & treatment != "DC" & set == "A" & incubation == "2")%>% 
  group_by(pine,shrubs,guild, incubation, block, plot) %>% 
  summarize_at(c("abundance","std_fungal_copies"), funs(sum = sum)) %>%
  ungroup()


#sqrt
regressionListA17 <- lapply(unique(humus.guildA17$guild), function(x) anova(lme(sqrt(std_fungal_copies_sum) ~ pine*shrubs, random = ~ 1|block, humus.guildA17[humus.guildA17$guild==x,])))
regressionListA17

qqnormPlots <- lapply(unique(humus.guildA17$guild),function(x) qqnorm(lme(sqrt(std_fungal_copies_sum) ~ pine*shrubs, random = ~ 1|block, humus.guildA17[humus.guildA17$guild==x,]),main = paste0(x)))
qqnormPlots

coefListA17 <- lapply(unique(humus.guildA17$guild),
                      function(x) summary(lme(sqrt(std_fungal_copies_sum) ~ pine*shrubs, random = ~ 1|block, humus.guildA17[humus.guildA17$guild==x,]))$coefficients)

#log (y +1)
regressionListA17 <- lapply(unique(humus.guildA17$guild), function(x) anova(lme(log(std_fungal_copies_sum+1) ~ pine*shrubs, random = ~ 1|block, humus.guildA17[humus.guildA17$guild==x,])))

qqnormPlots <- lapply(unique(humus.guildA17$guild),function(x) qqnorm(lme(log(std_fungal_copies_sum+1) ~ pine*shrubs, random = ~ 1|block, humus.guildA17[humus.guildA17$guild==x,]),main = paste0(x)))
qqnormPlots

coefListA17 <- lapply(unique(humus.guildA17$guild),
                      function(x) summary(lme(log(std_fungal_copies_sum+1) ~ pine*shrubs, random = ~ 1|block, humus.guildA17[humus.guildA17$guild==x,]))$coefficients)


### t.tests

humus.guildA17_ecto_wr_pine<- humus.guildA17 %>%
  filter(.,guild == "ectomycorrhizas white rot" & pine == "1")
humus.guildA17_ecto_wr_nopine<- humus.guildA17 %>%
  filter(.,guild == "ectomycorrhizas white rot" & pine == "0")

t_wr_pine <- t.test(humus.guildA17_ecto_wr_pine$std_fungal_copies_sum, y = humus.guildA17_ecto_wr_nopine$std_fungal_copies_sum)

humus.guildA17_ecto_pine<- humus.guildA17 %>%
  filter(.,guild == "ectomycorrhizas other" & pine == "1")
humus.guildA17_ecto_nopine<- humus.guildA17 %>%
  filter(.,guild == "ectomycorrhizas other" & pine == "0")

t_pine <- t.test(humus.guildA17_ecto_pine$std_fungal_copies_sum, y = humus.guildA17_ecto_nopine$std_fungal_copies_sum)

#shrubs
humus.guildA17_ecto_wr_shrubs<- humus.guildA17 %>%
  filter(.,guild == "ectomycorrhizas white rot" & pine == "1" & shrubs == "1")
humus.guildA17_ecto_wr_noshrubs<- humus.guildA17 %>%
  filter(.,guild == "ectomycorrhizas white rot" & pine == "1" & shrubs == "0")

t_wr_shrubs <- t.test(humus.guildA17_ecto_wr_shrubs$std_fungal_copies_sum, y = humus.guildA17_ecto_wr_noshrubs$std_fungal_copies_sum)

humus.guildA17_ecto_shrubs<- humus.guildA17 %>%
  filter(.,guild == "ectomycorrhizas other" & pine == "1" & shrubs == "1")
humus.guildA17_ecto_noshrubs<- humus.guildA17 %>%
  filter(.,guild == "ectomycorrhizas other" & pine == "1" & shrubs == "0")

t_wr_shrubs <- t.test(humus.guildA17_ecto_shrubs$std_fungal_copies_sum, y = humus.guildA17_ecto_noshrubs$std_fungal_copies_sum)


capture.output(humus_qPCR_anovaA17,unique(humus.guildA17$guild),regressionListA17,coefListA17,t_wr, t_wr_pine,file = "humus_guild_qpcr_statA17.txt")


#### SET B 17 #######

humus_statB17 <-  filter(meta.data,substrate == "Humus", treatment != "B", treatment != "DC", set == "B", incubation == "2")

lme<- lme(sqrt(copies_fungi) ~ pine*shrubs, #
          random = ~ 1|block, na.action = na.omit,
          data = humus_statB17)
humus_qPCR_anovaB17 <- anova(lme)
humus_qPCR_anovaB17$coef <- t(coef(lme))[,1]
qqnorm(lme)


humus.guildB17 <- data_long %>%
  filter(.,substrate == "Humus" & treatment != "B" & treatment != "DC", set == "B", incubation == "2")%>% 
  group_by(treatment,pine,shrubs,guild,block) %>% 
  summarize_at(c("abundance","std_fungal_copies"), funs(sum = sum)) %>%
  ungroup()

#check heteroscadiscity in log + 1
qqnormPlots <- lapply(unique(humus.guildB17$guild),function(x) qqnorm(lme(log(std_fungal_copies_sum+1) ~ pine*shrubs, random = ~ 1|block, humus.guildB17[humus.guildB17$guild==x,]), main = paste0(x)))
qqnormPlots

regressionList <- lapply(unique(humus.guildB17$guild), function(x) anova(lme(log(std_fungal_copies_sum+1) ~ pine*shrubs, random = ~ 1|block, humus.guildB17[humus.guildB17$guild==x,])))
regressionList


#sqrt transformation
qqnormPlots <- lapply(unique(humus.guildB17$guild),function(x) qqnorm(lme(sqrt(std_fungal_copies_sum) ~ pine*shrubs, random = ~ 1|block, humus.guildB17[humus.guildB17$guild==x,]), main = paste0(x)))
qqnormPlots

regressionList <- lapply(unique(humus.guildB17$guild), function(x) anova(lme(sqrt(std_fungal_copies_sum) ~ pine*shrubs, random = ~ 1|block, humus.guildB17[humus.guildB17$guild==x,])))


coefList <- lapply(unique(humus.guildB17$guild),
                   function(x) summary(lme(log(std_fungal_copies_sum+1) ~ pine*shrubs, random = ~ 1|block, humus.guildB17[humus.guildB17$guild==x,]))$coefficients)

## means
humus.guildB17_ecto_wr_trtmt<- humus.guildB17 %>%
  filter(.,guild == "ectomycorrhizas white rot") %>%
  group_by(treatment)%>%
summarize_at(c("abundance_sum","std_fungal_copies_sum"), funs(mean = mean)) %>%
  ungroup()

humus.guildB17_ecto_trtmt<- humus.guildB17 %>%
  filter(.,guild == "ectomycorrhizas other") %>%
  group_by(treatment)%>%
  summarize_at(c("abundance_sum","std_fungal_copies_sum"), funs(mean = mean)) %>%
  ungroup()

## t.tests
humus.guildB17_ecto_wr_pine<- humus.guildB17 %>%
  filter(.,guild == "ectomycorrhizas white rot" & pine == "1")
humus.guildB17_ecto_wr_nopine<- humus.guildB17 %>%
  filter(.,guild == "ectomycorrhizas white rot" & pine == "0")

t_wr_pine <- t.test(humus.guildB17_ecto_wr_pine$std_fungal_copies_sum, y = humus.guildB17_ecto_wr_nopine$std_fungal_copies_sum)

humus.guildB17_ecto_pine<- humus.guildB17 %>%
  filter(.,guild == "ectomycorrhizas other" & pine == "1")
humus.guildB17_ecto_nopine<- humus.guildB17 %>%
  filter(.,guild == "ectomycorrhizas other" & pine == "0")

t_pine <- t.test(humus.guildB17_ecto_pine$std_fungal_copies_sum, y = humus.guildB17_ecto_nopine$std_fungal_copies_sum)

#shrubs
humus.guildB17_ecto_wr_shrubs<- humus.guildB17 %>%
  filter(.,guild == "ectomycorrhizas white rot" & pine == "1" & shrubs == "1")
humus.guildB17_ecto_wr_noshrubs<- humus.guildB17 %>%
  filter(.,guild == "ectomycorrhizas white rot" & pine == "1" & shrubs == "0")

t_wr_shrubs <- t.test(humus.guildB17_ecto_wr_shrubs$std_fungal_copies_sum, y = humus.guildB17_ecto_wr_noshrubs$std_fungal_copies_sum)

humus.guildB17_ecto_shrubs<- humus.guildB17 %>%
  filter(.,guild == "ectomycorrhizas other" & pine == "1" & shrubs == "1")
humus.guildB17_ecto_noshrubs<- humus.guildB17 %>%
  filter(.,guild == "ectomycorrhizas other" & pine == "1" & shrubs == "0")

t_shrubs <- t.test(humus.guildB17_ecto_shrubs$std_fungal_copies_sum, y = humus.guildB17_ecto_noshrubs$std_fungal_copies_sum)


capture.output(humus_qPCR_anovaB17,unique(humus.guildB17$guild),regressionList,coefList, file = "humus_guild_qpcr_stat_B17.txt")

################################### Vegan (NMDS) ####
#make NMDS for each substrate
# Could rarefy the community to a sequencing depth of 1254 reads

#remove background samples
fungi_rbund_noB <- fungi_rbund[,meta.data$treatment != "B"]
meta.data_noB <- meta.data[meta.data$treatment != "B",]
meta.data_noB$treatment <- factor(meta.data_noB$treatment, levels = c('C', 'E', 'T','TE','DC'))

##### OVERALL NMDS PERMANOVA ####
spe.nmds <- metaMDS(t(fungi_rbund_noB), distance = "bray", na.rm = TRUE) #sqrt transform
spe.nmds$stress

NMDS1 <- spe.nmds$points[,1] ##also found using scores(spe.nmds)
NMDS2 <- spe.nmds$points[,2]

#Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores <- as.data.frame(scores(t(fungi_rbund_noB)))
#Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores <- as.data.frame(scores(spe.nmds, "species"))  
head(species.scores)  #look at the data

#bind columns
spe.plot <-cbind(t(fungi_rbund_noB), NMDS1, NMDS2, meta.data_noB)
na.omit <- na.omit(spe.plot)

spe.plot$set_by_incubation <- paste0(spe.plot$set,"_",spe.plot$incubation)
spe.plot <- spe.plot[order(spe.plot$set_by_incubation),]


bray <- vegdist(fungi_rbund_noB)
pinebetadisp <- betadisper(bray, group = meta.data_noB$pine)
anova(pinebetadisp) #P = 0.4942
shrubbetadisp <- betadisper(bray, group =meta.data_noB$shrub)
anova(shrubbetadisp) #P = 0.4305
setbetadisp <- betadisper(bray, group = meta.data_noB$set)
anova(setbetadisp) #0.625
duration.betadisp <- betadisper(bray, group = meta.data_noB$incubation)
anova(duration.betadisp) #0.3435

permanova <- adonis2(humus ~ meta.data_noB$shrub*meta.data_noB$pine*meta.data_noB$incubation*meta.data_noB$set, permutations = 999, method = "bray")

p1 <- ggplot(spe.plot, aes(NMDS1, NMDS2,shape = set_by_incubation, color=treatment))+
   geom_point(size=5) + 
   xlim(c(min(NMDS1)-0.045,max(NMDS1) + 0.045))+
   ylim(c(min(NMDS2)-0.045,max(NMDS2) + 0.045))+
  scale_color_manual("Treatments",values = c("#D45226","#F59B25", "#7F8A53", "#9AC287", "#C2D65C"),labels = c("+Pine +Shrubs","+Pine -Shrubs","-Pine +Shrubs","-Pine -Shrubs", " Disturbed control"))+
  scale_shape_manual(name = "Incubation",labels = c("Set 1, 5 months", "Set 1, 17 months", "Set 2, 5 months","Set 2, 17 months"), values = c(0, 15,1,16,1))+
  theme_minimal() + 
  theme(legend.title=element_blank(),axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"))+
 annotate("text",x=min(NMDS1)+max(NMDS1)/3,y=max(NMDS2)+.035,label=paste0("stress = ",round(spe.nmds$stress, digits = 3)))

p2 <-ggplot(spe.plot, aes(NMDS1, NMDS2,shape = set_by_incubation, color=substrate, label = Unique_Sample))+
  geom_point(size=5) + #geom_text()+  #add the core labels
  scale_color_manual(name = "Substrate",label = c("Humus", "Pine Litter"), values = c("#472F2F", "#8C470E"))+
  scale_shape_manual(name = "Set and Incubation",labels = c("Set A, 5 months", "Set A, 17 months", "Set B, 5 months","Set B, 17 months"), values = c(0, 15,1,16))+
  theme_minimal() +
  #annotate("text",x=-1,y=2.25,label=paste0("stress = ",round(spe.nmds$stress, digits = 3)))+
  theme(legend.title=element_blank(),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"))+
  guides(shape= "none")

 ggarrange(p1 + rremove("xlab"),p2, labels = c("a)","b)"),font.label = list(size = 16, color = "black"), combine.legend = TRUE, legend = "right", nrow = 2, vjust = 1)%>%
  ggexport(filename = "figures/Rplot_NMDS_CombinedHumus&Litter_FungalCommunity.pdf")

-----------------------

###### Litter NMDS #####
litter <- fungi_rbund_noB[,meta.data_noB$substrate == "Litter"]
metadata_litter <- meta.data_noB[meta.data_noB$substrate == "Litter",]
dim(litter) #190 159
dim(metadata_litter) #159 31

litter <-  litter[rowSums(litter) != 0,] #prune OTU (row) with sums equal to 0
#166 159
#litter <- data.frame(t(make_relative(as.matrix(litter_trimd))))

spe.nmds <- metaMDS(t(litter), distance = "bray", na.rm = TRUE)
spe.nmds$stress

NMDS1 <- spe.nmds$points[,1] ##also found using scores(spe.nmds)
NMDS2 <- spe.nmds$points[,2]

#Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores <- as.data.frame(scores(litter)) 

#look at the data
head(data.scores)
str(data.scores)
species.scores <- as.data.frame(scores(spe.nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
head(species.scores)  #look at the data

#bind columns
spe.plot <- cbind(t(litter), NMDS1, NMDS2, metadata_litter)

spe.plot$set_by_incubation <- paste0(spe.plot$set,"_",spe.plot$incubation)
spe.plot <- spe.plot[order(spe.plot$set_by_incubation),]

litter.bray <- vegdist(t(litter))

#testing the assumption that dispersion is similar within groups
trtmtbetadisp <- betadisper(litter.bray, group = metadata_litter$treatment)
anova(trtmtbetadisp) # P = 0.416
setbetadisp <- betadisper(litter.bray, group = metadata_litter$set) # does not meet homogeneity of variance
anova(setbetadisp)
durationbetadisp <- betadisper(litter.bray, group = metadata_litter$incubation) # does not meet homogeneity of variance
anova(durationbetadisp)
pinebetadisp <- betadisper(litter.bray, group = metadata_litter$pine)
anova(pinebetadisp) # P = 0.3183
shrubbetadisp <- betadisper(litter.bray, group = metadata_litter$shrub)
anova(shrubbetadisp) # P = 0.6174

permanova <- adonis2(t(litter) ~ metadata_litter$shrub*metadata_litter$pine*metadata_litter$set*metadata_litter$incubation, strata = metadata_litter$block, permutations = 999, method = "bray")
# 
capture.output(permanova,"set",anova(setbetadisp),"duration",anova(durationbetadisp),"shrub",anova(shrubbetadisp),"pine",anova(pinebetadisp), file = "permanova_needle_shrub*pine_betadispr.txt")

# write.csv(permanova,"tables/permanova_needles_shrub*pine.csv")
# permanova <- adonis2(t(litter) ~ metadata_litter$treatment, strata = metadata_litter$block, permutations = 999, method = "bray")

litternmds <- ggplot(spe.plot, aes(NMDS1, NMDS2,shape = set_by_incubation, color=treatment))+
    geom_point(size=5) +
    xlim(c(min(NMDS1)-0.045,max(NMDS1) + 0.045))+
    ylim(c(min(NMDS2)-0.045,max(NMDS2) + 0.045))+
    scale_color_manual(name = "Treatment",label = c("+Shrubs +Pine", "-Shrubs +Pine", "+Shrubs -Pine", "-Shrubs -Pine","Disturbed control"), values = c("#D45226","#F59B25", "#7F8A53", "#9AC287", "#C2D65C"))+
    scale_shape_manual(name = "Incubation",labels = c("Set 1, 5 months", "Set 1, 17 months", "Set 2, 5 months","Set 2, 17 months"), values = c(0, 15,1,16,1))+
    theme_minimal() +
    theme(legend.title=element_blank(),axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"))+
    annotate("text",x=min(NMDS1)+max(NMDS1)/3,y=max(NMDS2)+.035,label=paste0("stress = ",round(spe.nmds$stress, digits = 3)))

###### Litter 17 month Set 2 NMDS #####
litter17B <- fungi_rbund_noB[,meta.data_noB$substrate == "Litter" & meta.data_noB$incubation == "2" & meta.data_noB$set == "B" ]
metadata_litter17B <- meta.data_noB[meta.data_noB$substrate == "Litter" & meta.data_noB$incubation == "2" & meta.data_noB$set == "B",]
dim(litter17B) #190 159
dim(metadata_litter17B) #159 31

litter17B <-  litter17B[rowSums(litter17B) != 0,] #prune OTU (row) with sums equal to 0
#166 159
#litter <- data.frame(t(make_relative(as.matrix(litter_trimd))))

spe.nmds <- metaMDS(t(litter17B), distance = "bray", na.rm = TRUE)
spe.nmds$stress

NMDS1 <- spe.nmds$points[,1] ##also found using scores(spe.nmds)
NMDS2 <- spe.nmds$points[,2]

#Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores <- as.data.frame(scores(litter17B)) 

#look at the data
head(data.scores)
str(data.scores)
species.scores <- as.data.frame(scores(spe.nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
head(species.scores)  #look at the data

#bind columns
spe.plot <- cbind(t(litter17B), NMDS1, NMDS2, metadata_litter17B)

litter17B.bray <- vegdist(t(litter17B))

#testing the assumption that dispersion is similar within groups
trtmtbetadisp <- betadisper(litter17B.bray, group = metadata_litter17B$treatment)
anova(trtmtbetadisp) # P = 0.1748
pinebetadisp <- betadisper(litter17B.bray, group = metadata_litter17B$pine)
anova(pinebetadisp) # P = 0.0239
shrubbetadisp <- betadisper(litter17B.bray, group = metadata_litter17B$shrub)
anova(shrubbetadisp) # P = 0.5895

permanova <- adonis2(t(litter17B) ~ metadata_litter17B$treatment, strata = metadata_litter17B$block, permutations = 999, method = "bray")

# capture.output()
litternmdsB17 <- ggplot(spe.plot, aes(NMDS1, NMDS2, color=treatment))+
  geom_point(size=5) +
  xlim(c(min(NMDS1)-0.045,max(NMDS1) + 0.045))+
  ylim(c(min(NMDS2)-0.045,max(NMDS2) + 0.045))+
  scale_color_manual(name = "Treatment",label = c("+Shrubs +Pine", "-Shrubs +Pine", "+Shrubs -Pine", "-Shrubs -Pine","Disturbed control"), values = c("#D45226","#F59B25", "#7F8A53", "#9AC287", "#C2D65C"))+
  theme_minimal() +
  theme(legend.title=element_blank(),axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"))+
  annotate("text",x=min(NMDS1)+max(NMDS1)/3,y=max(NMDS2)+.035,label=paste0("stress = ",round(spe.nmds$stress, digits = 3)))

###### Humus NMDS #####
humus <- fungi_rbund_noB[,meta.data_noB$substrate == "Humus" & meta.data_noB$Unique_Sample != "DC1Humus17B"] #removing
metadata_humus <- meta.data_noB[meta.data_noB$substrate == "Humus" & meta.data_noB$Unique_Sample != "DC1Humus17B",] #
dim(humus) #184 157

humus <- humus[,metadata_humus$Unique_Sample != "C3Humus17B"] #removing outlier from control
metadata_humus <- metadata_humus[metadata_humus$Unique_Sample != "C3Humus17B",] #
dim(humus) #184 157

any(rowSums(humus)== "0")

humus_trimd <-  humus[rowSums(humus) != 0,] #prune OTU (row) with sums equal to 0

humus <- data.frame(t(make_relative(as.matrix(humus_trimd))))
spe.nmds <- metaMDS(humus, distance = "bray",na.rm = TRUE)

NMDS1 <- spe.nmds$points[,1] ##also found using scores(spe.nmds)
NMDS2 <- spe.nmds$points[,2]

#Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores <- as.data.frame(scores(humus)) 

#look at the data
species.scores <- as.data.frame(scores(spe.nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
head(species.scores)  #look at the data

#bind columns
spe.plot <-cbind(humus, NMDS1, NMDS2, metadata_humus)
spe.plot$set_by_incubation <- paste0(spe.plot$set,"_",spe.plot$incubation)
spe.plot <- spe.plot[order(spe.plot$set_by_incubation),]

humus.bray <- vegdist(humus)
pinebetadisp <- betadisper(humus.bray, group = metadata_humus$pine)
anova(pinebetadisp) #P = 0.4942
shrubbetadisp <- betadisper(humus.bray, group = metadata_humus$shrub)
anova(shrubbetadisp) #P = 0.4305
setbetadisp <- betadisper(humus.bray, group = metadata_humus$set)
anova(setbetadisp) #0.625
duration.betadisp <- betadisper(humus.bray, group = metadata_humus$incubation)
anova(duration.betadisp) #0.3435

permanova <- adonis2(humus ~ metadata_humus$shrub*metadata_humus$pine*metadata_humus$incubation*metadata_humus$set, permutations = 999, method = "bray")

capture.output(permanova,anova(setbetadisp),anova(duration.betadisp),anova(shrubbetadisp),anova(pinebetadisp), file = "humus_permanova_betadispr.txt")

#
humus.bray <- vegdist(humus)
trtbetadisp <- betadisper(humus.bray, group = metadata_humus$treatment)
anova(trtbetadisp) #P = 0.756
setbetadisp <- betadisper(humus.bray, group = metadata_humus$set)
anova(setbetadisp)
duration.betadisp <- betadisper(humus.bray, group = metadata_humus$incubation)
anova(duration.betadisp)

# permanova <- adonis2(humus ~ metadata_humus$treatment*metadata_humus$incubation*metadata_humus$set, permutations = 999, method = "bray")
# write.csv(permanova,"tables/permanova_humus.csv")

humusnmds <- ggplot(spe.plot, aes(NMDS1, NMDS2,shape = set_by_incubation, color=treatment))+
  geom_point(size=5) +
  xlim(c(min(NMDS1)-0.045,max(NMDS1) + 0.045))+
  ylim(c(min(NMDS2)-0.045,max(NMDS2) + 0.045))+
  scale_x_discrete()
  scale_color_manual(name = "Treatment",label = c("+Shrubs +Pine", "-Shrubs +Pine", "+Shrubs -Pine", "-Shrubs -Pine","Disturbed control"), values = c("#D45226","#F59B25", "#7F8A53", "#9AC287", "#C2D65C"))+
  scale_shape_manual(name = "Incubation",labels = c("Set 1, 5 months", "Set 1, 17 months", "Set 2, 5 months","Set 2, 17 months"), values = c(0, 15,1,16,1))+
  theme_minimal() +
  theme(legend.title=element_blank(),axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"))+
  annotate("text",x=min(NMDS1)+max(NMDS1)/4,y=max(NMDS2)+.035,label=paste0("stress = ",round(spe.nmds$stress, digits = 3)))

########
ggarrange(litternmds,humusnmds, labels = c("a)","b)"),font.label = list(size = 18, color = "black"),  legend = "right", combined.legend = TRUE, nrow = 2, vjust = 1.5)%>%
  ggexport(filename = "figures/Rplot_NMDS_Humus&Litter_FungalCommunity_Separate.pdf")
