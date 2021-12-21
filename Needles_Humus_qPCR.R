###Jädraås Decomposition Experiment - qPCR

rm(list=ls()) #clear the workspace

#Author: Louis Mielke
#Date Created: 2021 NOV 18
#Date Updated: 2021 DEC 02

# NOTES: Still need to update qPCR data by multiplying by fraction of fungi per sample

library(phyloseq);library(dplyr);library(tidyverse);library(tibble);library(viridis);library(ggplot2);library(vegan);library(ggsci);library(tidyr);library(scales);library(gridExtra);library(grid);library(lattice);library(MASS);library(lme4);library(nlme);library(emmeans);library(gamm4);library(dplyr);library(magrittr);library(reshape2);library(mgcv);library(lmerTest);library(lubridate);library(car);library(lmerTest);library(openxlsx)

setwd("~/Projects/mycorrhizal removal/MeshBags/Jädraås_Gadgil_MycorrhizalTraits")

### import metadata
meta <- read.table("metadata.txt", sep="\t", header=T, row.names=1)
class(meta);dim(meta);str(meta) #meta$start_date <- as.Date(meta$start_date)
meta$incubation <- as.factor(meta$incubation)
meta$block <- as.factor(meta$block)



### formatting the community abundance table dataset
community <- read.csv("community2.csv", header = T)
community <- data.frame(community, row.names = 1) #in one line set column 1 as row names
str(community);community [1:5,1:5];class(community);dim(community) #View(community)
OTU <- t(community) #transpose


#add plot level factor for downstream analysis
meta$plot <- paste(meta$treatment,meta$block, sep = "")
#B = background samples of intitial substrates
#C = control, 
#DC = disturbed control, 
#E = ericaceous shrub removal, 
#T = pine root exclusion, 
#TE = combined shrub removal and pine root exclusion

meta <- subset(meta, treatment %in% c("B","C","DC","E","T","TE"))
humus0 <- subset(meta, substrate %in% c("Humus"))
litter0 <- subset(meta, substrate %in% c("Litter"))


humus <- subset(humus0, treatment %in% c("C","DC","E","T","TE")) # dataset for graphing
humus.background <- subset(humus0, treatment %in% c("B")) # dataset for graphing
litter <- subset(litter0, treatment %in% c("C","DC","E","T","TE")) #dataset for graphing
litter.background <- subset(litter0, treatment %in% c("B")) # dataset for graphing



##################### qPCR STATISTICS #####


#subset by treatment
humus2 <- subset(humus0, treatment %in% c("C","E","T","TE"))
litter2 <- subset(litter0, treatment %in% c("C","E","T","TE"))
str(humus2)

####### HUMUS STATS ######
lme.h_qPCR <- lme(log(copies_DNA_per_g_substrate) ~ pine*shrubs*incubation*set, #
              random = ~ 1|block/plot, na.action = na.omit,
              data = humus2)

h_qPCR_anova <- anova(lme.h_qPCR)
h_qPCR_coef <- coef(lme.h_qPCR)
plot(resid(lme.h_qPCR))

#subset to each set A and B
humusA <- humus2[humus2$set == "A",]
humusB <- humus2[humus2$set == "B",]

#set A
lme.hA_qPCR <- lme(log(copies_DNA_per_g_substrate) ~ pine*shrubs*incubation, #
                  random = ~ 1|block/plot, na.action = na.omit,
                  data = humusA)

hA_qPCR_anova <- anova(lme.hA_qPCR)
hA_qPCR_coef <- coef(lme.hA_qPCR)
plot(resid(lme.hA_qPCR))

          #subset to each set incubation period
          humus_A5 <- humusA[humusA$incubation == "5",]
          humus_A17 <- humusA[humusA$incubation == "17",]

          lme.hA5_qPCR <- lme(log(copies_DNA_per_g_substrate) ~ pine*shrubs, #
                              random = ~ 1|block, na.action = na.omit,
                              data = humus_A5)
          
          hA5_qPCR_anova <- anova(lme.hA5_qPCR)
          hA5_qPCR_coef <- coef(lme.hA5_qPCR)
          plot(resid(lme.hA5_qPCR))
          
          lme.hA17_qPCR <- lme(log(copies_DNA_per_g_substrate) ~ pine*shrubs, #
                               random = ~ 1|block, na.action = na.omit,
                               data = humus_A17)
          
          hA17_qPCR_anova <- anova(lme.hA17_qPCR)
          hA17_qPCR_coef <- coef(lme.hA17_qPCR)
          plot(resid(lme.hA17_qPCR))
          
          
#set B
lme.hB_qPCR <- lme(log(copies_DNA_per_g_substrate) ~ pine*shrubs*incubation, #
                  random = ~ 1|block/plot, na.action = na.omit,
                  data = humusB)

hB_qPCR_anova <- anova(lme.hB_qPCR)
hB_qPCR_coef <- coef(lme.hB_qPCR)
plot(resid(lme.hB_qPCR))

            #subset to each set incubation period
            humus_B5 <- humusB[humusB$incubation == "5",]
            humus_B17 <- humusB[humusB$incubation == "17",]

            lme.hB5_qPCR <- lme(log(copies_DNA_per_g_substrate) ~ pine*shrubs, #
                               random = ~ 1|block, na.action = na.omit,
                               data = humus_B5)
            
            hB5_qPCR_anova <- anova(lme.hB5_qPCR)
            hB5_qPCR_coef <- coef(lme.hB5_qPCR)
            plot(resid(lme.hB5_qPCR))
            
            lme.hB17_qPCR <- lme(log(copies_DNA_per_g_substrate) ~ pine*shrubs, #
                               random = ~ 1|block, na.action = na.omit,
                               data = humus_B17)
            
            hB17_qPCR_anova <- anova(lme.hB17_qPCR)
            hB17_qPCR_coef <- coef(lme.hB17_qPCR)
            plot(resid(lme.hB17_qPCR))

# Create a blank workbook
humus_qpcr_anova_output <- createWorkbook()

# Add some sheets to the workbook
addWorksheet(humus_qpcr_anova_output, "Full Model")
addWorksheet(humus_qpcr_anova_output, "Set A")
addWorksheet(humus_qpcr_anova_output, "Set B")
addWorksheet(humus_qpcr_anova_output, "Set A 5 Months")
addWorksheet(humus_qpcr_anova_output, "Set A 17 Months")
addWorksheet(humus_qpcr_anova_output, "Set B 5 Months")
addWorksheet(humus_qpcr_anova_output, "Set B 17 Months")
addWorksheet(humus_qpcr_anova_output, "Coef Full Model")
addWorksheet(humus_qpcr_anova_output, "Coef Set A")
addWorksheet(humus_qpcr_anova_output, "Coef Set B")
addWorksheet(humus_qpcr_anova_output, "Coef Set A 5 Months")
addWorksheet(humus_qpcr_anova_output, "Coef Set A 17 Months")
addWorksheet(humus_qpcr_anova_output, "Coef Set B 5 Months")
addWorksheet(humus_qpcr_anova_output, "Coef Set B 17 Months")

# Write the data to the sheets
writeData(humus_qpcr_anova_output, sheet = "Full Model", x = h_qPCR_anova, rowNames = TRUE)
writeData(humus_qpcr_anova_output, sheet = "Set A", x = hA_qPCR_anova,rowNames = TRUE)
writeData(humus_qpcr_anova_output, sheet = "Set B", x = hB_qPCR_anova,rowNames = TRUE)
writeData(humus_qpcr_anova_output, sheet = "Set A 5 Months", x = hA5_qPCR_anova,rowNames = TRUE)
writeData(humus_qpcr_anova_output, sheet = "Set B 5 Months", x = hB5_qPCR_anova,rowNames = TRUE)
writeData(humus_qpcr_anova_output, sheet = "Set A 17 Months", x = hA17_qPCR_anova,rowNames = TRUE)
writeData(humus_qpcr_anova_output, sheet = "Set B 17 Months", x = hB17_qPCR_anova,rowNames = TRUE)
writeData(humus_qpcr_anova_output, sheet = "Coef Full Model", x = h_qPCR_coef, rowNames = TRUE)
writeData(humus_qpcr_anova_output, sheet = "Coef Set A", x = hA_qPCR_coef,rowNames = TRUE)
writeData(humus_qpcr_anova_output, sheet = "Coef Set B", x = hB_qPCR_coef,rowNames = TRUE)
writeData(humus_qpcr_anova_output, sheet = "Coef Set A 5 Months", x = hA5_qPCR_coef,rowNames = TRUE)
writeData(humus_qpcr_anova_output, sheet = "Coef Set B 5 Months", x = hB5_qPCR_coef,rowNames = TRUE)
writeData(humus_qpcr_anova_output, sheet = "Coef Set A 17 Months", x = hA17_qPCR_coef,rowNames = TRUE)
writeData(humus_qpcr_anova_output, sheet = "Coef Set B 17 Months", x = hB17_qPCR_coef,rowNames = TRUE)


# Export the file
saveWorkbook(humus_qpcr_anova_output, "humus_qpcr_anova_output2.xlsx")


####### PINE NEEDLE STATS ######
lme.l_qPCR <-lme(log(copies_DNA_per_g_substrate) ~ pine*shrubs*incubation*set, #
                 random = ~ 1|block/plot,
                 data = litter2, na.action = na.omit)

l_qPCR_anova <- anova(lme.l_qPCR)
l_qPCR_coef <- coef(lme.l_qPCR)
plot(resid(lme.l_qPCR))

#subset to each set A and B
litterA <- litter2[litter2$set == "A",]
litterB <- litter2[litter2$set == "B",]

#set A
lme.lA_qPCR <- lme(log(copies_DNA_per_g_substrate) ~ pine*shrubs*incubation, #
                   random = ~ 1|block/plot, na.action = na.omit,
                   data = litterA)

lA_qPCR_anova <- anova(lme.lA_qPCR)
lA_qPCR_coef <- coef(lme.lA_qPCR)
plot(resid(lme.lA_qPCR))

#subset to each set incubation period
litter_A5 <- litterA[litterA$incubation == "5",]
litter_A17 <- litterA[litterA$incubation == "17",]

lme.lA5_qPCR <- lme(log(copies_DNA_per_g_substrate) ~ pine*shrubs, #
                    random = ~ 1|block, na.action = na.omit,
                    data = litter_A5)

lA5_qPCR_anova <- anova(lme.lA5_qPCR)
lA5_qPCR_coef <- coef(lme.lA5_qPCR)
plot(resid(lme.lA5_qPCR))

lme.lA17_qPCR <- lme(log(copies_DNA_per_g_substrate) ~ pine*shrubs, #
                     random = ~ 1|block, na.action = na.omit,
                     data = litter_A17)

lA17_qPCR_anova <- anova(lme.lA17_qPCR)
lA17_qPCR_coef <- coef(lme.lA17_qPCR)
plot(resid(lme.lA17_qPCR))


#set B
lme.lB_qPCR <- lme(log(copies_DNA_per_g_substrate) ~ pine*shrubs*incubation, #
                   random = ~ 1|block/plot, na.action = na.omit,
                   data = litterB)

lB_qPCR_anova <- anova(lme.lB_qPCR)
lB_qPCR_coef <- coef(lme.lB_qPCR)
plot(resid(lme.lB_qPCR))

#subset to each set incubation period
litter_B5 <- litterB[litterB$incubation == "5",]
litter_B17 <- litterB[litterB$incubation == "17",]

lme.lB5_qPCR <- lme(log(copies_DNA_per_g_substrate) ~ pine*shrubs, #
                    random = ~ 1|block, na.action = na.omit,
                    data = litter_B5)

lB5_qPCR_anova <- anova(lme.lB5_qPCR)
lB5_qPCR_coef <- coef(lme.lB5_qPCR)
plot(resid(lme.lB5_qPCR))

lme.lB17_qPCR <- lme(log(copies_DNA_per_g_substrate) ~ pine*shrubs, #
                     random = ~ 1|block, na.action = na.omit,
                     data = litter_B17)

lB17_qPCR_anova <- anova(lme.lB17_qPCR)
lB17_qPCR_coef <- coef(lme.lB17_qPCR)
plot(resid(lme.lB17_qPCR))

# Create a blank workbook
litter_qpcr_anova_output <- createWorkbook()

# Add some sheets to the workbook
addWorksheet(litter_qpcr_anova_output, "Full Model")
addWorksheet(litter_qpcr_anova_output, "Set A")
addWorksheet(litter_qpcr_anova_output, "Set B")
addWorksheet(litter_qpcr_anova_output, "Set A 5 Months")
addWorksheet(litter_qpcr_anova_output, "Set A 17 Months")
addWorksheet(litter_qpcr_anova_output, "Set B 5 Months")
addWorksheet(litter_qpcr_anova_output, "Set B 17 Months")
addWorksheet(litter_qpcr_anova_output, "Coef Full Model")
addWorksheet(litter_qpcr_anova_output, "Coef Set A")
addWorksheet(litter_qpcr_anova_output, "Coef Set B")
addWorksheet(litter_qpcr_anova_output, "Coef Set A 5 Months")
addWorksheet(litter_qpcr_anova_output, "Coef Set A 17 Months")
addWorksheet(litter_qpcr_anova_output, "Coef Set B 5 Months")
addWorksheet(litter_qpcr_anova_output, "Coef Set B 17 Months")

# Write the data to the sheets
writeData(litter_qpcr_anova_output, sheet = "Full Model", x =l_qPCR_anova, rowNames = TRUE)
writeData(litter_qpcr_anova_output, sheet = "Set A", x = lA_qPCR_anova,rowNames = TRUE)
writeData(litter_qpcr_anova_output, sheet = "Set B", x = lB_qPCR_anova,rowNames = TRUE)
writeData(litter_qpcr_anova_output, sheet = "Set A 5 Months", x = lA5_qPCR_anova,rowNames = TRUE)
writeData(litter_qpcr_anova_output, sheet = "Set B 5 Months", x = lB5_qPCR_anova,rowNames = TRUE)
writeData(litter_qpcr_anova_output, sheet = "Set A 17 Months", x = lA17_qPCR_anova,rowNames = TRUE)
writeData(litter_qpcr_anova_output, sheet = "Set B 17 Months", x = lB17_qPCR_anova,rowNames = TRUE)
writeData(litter_qpcr_anova_output, sheet = "Coef Full Model", x = l_qPCR_coef, rowNames = TRUE)
writeData(litter_qpcr_anova_output, sheet = "Coef Set A", x = lA_qPCR_coef,rowNames = TRUE)
writeData(litter_qpcr_anova_output, sheet = "Coef Set B", x = lB_qPCR_coef,rowNames = TRUE)
writeData(litter_qpcr_anova_output, sheet = "Coef Set A 5 Months", x = lA5_qPCR_coef,rowNames = TRUE)
writeData(litter_qpcr_anova_output, sheet = "Coef Set B 5 Months", x = lB5_qPCR_coef,rowNames = TRUE)
writeData(litter_qpcr_anova_output, sheet = "Coef Set A 17 Months", x = lA17_qPCR_coef,rowNames = TRUE)
writeData(litter_qpcr_anova_output, sheet = "Coef Set B 17 Months", x = lB17_qPCR_coef,rowNames = TRUE)


# Export the file
saveWorkbook(litter_qpcr_anova_output, "litter_qpcr_anova_output.xlsx")


##################### qPCR PLOTTING ########

######## HUMUS GRAPHS ########

#make dummy variables for plotting
dummy <- tibble(set = c("A"),
                treatment = c("C","DC","E","T","TE"), 
                incubation = 0, 
                mean = mean(na.omit(humus.background$copies_DNA_per_g_substrate)),
                sd = sd(na.omit(humus.background$copies_DNA_per_g_substrate)),
                min =0,max=0,
                n=length(na.omit(humus.background$copies_DNA_per_g_substrate)),
                se=sd/sqrt(n))

dummy$treatment <- as.factor(dummy$treatment)
dummy$incubation <- as.factor(dummy$incubation)
dummy$n <- as.integer(dummy$n)


# averaging humus across treatment, incubation and set
data_grouped <- humus %>%
  group_by(set, treatment,incubation) %>% 
  na.omit() %>%
  summarize_at(c("copies_DNA_per_g_substrate"), 
               funs(mean, sd, n(), se=sd(.)/sqrt(n())))

humus_qPCR <- ungroup(data_grouped)

#subset to each set A and B
humusA_qPCR <- humus_qPCR[humus_qPCR$set == "A",]
humusB_qPCR <- humus_qPCR[humus_qPCR$set == "B",]

#combine data
humus_summaryA <- bind_rows(humusA_qPCR,dummy)
humus_summaryA$incubation <- as.numeric(as.character(humus_summaryA$incubation))

humus_summaryB <- bind_rows(humusB_qPCR,dummy)
humus_summaryB$incubation <- as.numeric(as.character(humus_summaryB$incubation))

#### Set A

ggplot(data = humus_summaryA, aes(x = incubation, y = mean, color = treatment, linetype = treatment))+ geom_line(size = 2)+
  #geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 1, position = "identity", linetype = "solid")+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 2, position = position_dodge(width=0.75), linetype = "solid")+
  scale_colour_manual(values = c("mediumseagreen", "grey", "black", "mediumseagreen","black")) +
  xlab("Months")+
  ylim(c(0,5*10^9))+
  ylab("copies per g humus")+
  theme_classic()+
  theme(legend.position="top")+
  scale_linetype_manual(values = c("solid", "dotted", "solid", "dotted", "dotted"))

#### Set B

ggplot(data = humus_summaryB, aes(x = incubation, y = mean, color = treatment, linetype = treatment))+ geom_line(size = 2)+
  #geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 1, position = "identity", linetype = "solid")+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 2, position = position_dodge(width=0.75), linetype = "solid")+
  scale_colour_manual(values = c("mediumseagreen", "grey", "black", "mediumseagreen","black")) +
  xlab("Months")+
  ylab("copies per g humus")+
  theme_classic()+
  ylim(c(0,5*10^9))+
  theme(legend.position="top")+
  scale_linetype_manual(values = c("solid", "dotted", "solid", "dotted", "dotted"))


######## PINE NEEDLES GRAPHS ######

# averaging growing season integrated results for each plot across all three years
data_grouped <- litter %>%
  group_by(set, treatment,incubation) %>% 
  na.omit() %>%
  summarize_at(c("copies_DNA_per_g_substrate"), funs(mean, sd, n(), se=sd(.)/sqrt(n())))
litter_qPCR <- ungroup(data_grouped)

#make dummy variables for plotting
dummy.litter <- tibble(set = c("A"),
                treatment = c("C","DC","E","T","TE"), 
                incubation = 0, 
                mean = mean(na.omit(litter.background$copies_DNA_per_g_substrate)),
                sd = sd(na.omit(litter.background$copies_DNA_per_g_substrate)),
                min =0,max=0,
                n=length(na.omit(litter.background$copies_DNA_per_g_substrate)),
                se=sd/sqrt(n))

dummy.litter$treatment <- as.factor(dummy.litter$treatment)
dummy.litter$incubation <- as.factor(dummy.litter$incubation)
dummy.litter$n <- as.integer(dummy.litter$n)


#subset to each set A and B
litterA_qPCR <- litter_qPCR[litter_qPCR$set == "A",]
litterB_qPCR <- litter_qPCR[litter_qPCR$set == "B",]

#combine data
litter_summaryA <- bind_rows(litterA_qPCR,dummy.litter)
litter_summaryA$incubation <- as.numeric(as.character(litter_summaryA$incubation))

litter_summaryB <- bind_rows(litterB_qPCR,dummy.litter)
litter_summaryB$incubation <- as.numeric(as.character(litter_summaryB$incubation))

#### Set A

ggplot(data = litter_summaryA, aes(x = incubation, y = mean, color = treatment, linetype = treatment))+ geom_line(size = 2)+
  #geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 1, position = "identity", linetype = "solid")+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 2, position = position_dodge(width=0.75), linetype = "solid")+
  scale_colour_manual(values = c("mediumseagreen", "grey", "black", "mediumseagreen","black")) +
  xlab("Months")+
  ylab("copies per g litter")+
  ylim(c(0,2*10^10))+
  theme_classic()+
  theme(legend.position="top")+
  scale_linetype_manual(values = c("solid", "dotted", "solid", "dotted", "dotted"))

#### Set B

ggplot(data = litter_summaryB, aes(x = incubation, y = mean, color = treatment, linetype = treatment))+ geom_line(size = 2)+
  #geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 1, position = "identity", linetype = "solid")+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 2, position = position_dodge(width=0.75), linetype = "solid")+
  scale_colour_manual(values = c("mediumseagreen", "grey", "black", "mediumseagreen","black")) +
  xlab("Months")+
  ylab("copies per g litter")+
  ylim(c(0,2*10^10))+
  theme_classic()+
  theme(legend.position="top")+
  scale_linetype_manual(values = c("solid", "dotted", "solid", "dotted", "dotted"))











####### EDits Below





######## QPCR standardized communities #######
#multiplying number of copies per g substrate by the relative abundance
copies <- meta$copies_DNA_per_g_substrate

community_qpcr <- mapply("*",as.data.frame(community[1:317,1:798]),meta$copies_DNA_per_g_substrate)

sapply(1:ncol(community),function(x) community[,x] * v[x] )
qPCRxOTU <- t(community_qpcr) #transpose
class(OTU)


######### PINE NEEDLE LITTER ########

#subset OTU table
litterOTU <- OTU[,meta$substrate == "Litter"]

#subset metadata by substrate and start date
lmeta<- meta[meta$substrate == "Litter",]
lmeta1 <- lmeta[lmeta$start_date == "2017-06-16",] #round 1 
lmeta2 <- lmeta[lmeta$start_date == "2018-06-14",] #round 2

#subset by incubation length in round 1
lmeta1.5 <- lmeta1[lmeta1$incubation_months == "5",] #5 months
lmeta1.17 <- lmeta1[lmeta1$incubation_months == "17",] #17 months

#round 2
lmeta2.5 <- lmeta2[lmeta2$incubation_months == "5",] #5 months
lmeta2.17 <- lmeta2[lmeta2$incubation_months == "17",] #17 months


#subset OTU table with metadata
litterOTU1 <- litterOTU[,lmeta$start_date == "2017-06-16"] #round 1
litterOTU2 <- litterOTU[,lmeta$start_date == "2018-06-14"] #round 2
str(litterOTU1)
#round 1
lOTU1.5 <- litterOTU1[lmeta1$incubation_months == "5",] #5 months
lOTU1.17 <- litterOTU1[lmeta1$incubation_months == "17",] #17 months

#round 2
lOTU2.5 <- litterOTU2[lmeta2$incubation_months == "5",] #5 months
lOTU2.17 <- litterOTU2[lmeta2$incubation_months == "17",] #17 months


######  format OTU table and sample metadata for phyloseq #####

lOTU1.5=otu_table(lOTU1.5,taxa_are_rows=T)
dim(lOTU1.5);class(lOTU1.5);str(lOTU1.5)

#format samples for phyloseq
lmeta1.5=sample_data(lmeta1.5)
dim(lmeta1.5);class(lmeta1.5);str(lmeta1.5)
#View(lmeta1.5)

#make phyloseq object
lp1.5 <- phyloseq(lOTU1.5,lmeta1.5,taxa)
View(lp1.5)

#format OTU table for phyloseq
lOTU1.17=otu_table(lOTU1.17,taxa_are_rows=T)
#dim(lOTU1.17);class(lOTU1.17);str(lOTU1.17)

#format samples for phyloseq
lmeta1.17=sample_data(lmeta1.17)
dim(lmeta1.17);class(lmeta1.17);str(lmeta1.17)
#View(lmeta1.17)

#make phyloseq object
lp1.17 <- phyloseq(lOTU1.17,lmeta1.17,taxa)
View(lp1.17)


lOTU2.5=otu_table(lOTU2.5,taxa_are_rows=T)
#dim(lOTU2.5);class(lOTU2.5);str(lOTU2.5)

#format samples for phyloseq
lmeta2.5=sample_data(lmeta2.5)
dim(lmeta2.5);class(lmeta2.5);str(lmeta2.5)
View(lmeta2.5)

#make phyloseq object
lp2.5 <- phyloseq(lOTU2.5,lmeta2.5,taxa)
View(lp2.5)


#format OTU table for phyloseq
lOTU2.17=otu_table(lOTU2.17,taxa_are_rows=T)
#dim(lOTU2.17);class(lOTU2.17);str(lOTU2.17)

#format samples for phyloseq
lmeta2.17=sample_data(lmeta2.17)
#dim(lmeta2.17);class(lmeta2.17);str(lmeta2.17)
View(lmeta2.17)

#make phyloseq object
lp2.17 <- phyloseq(lOTU2.17,lmeta2.17,taxa)
View(lp2.17)


########### Litter Round 1, 5 month incubation #######

#any empty samples or taxa?
any(sample_sums(lp1.5) == 0)
any(taxa_sums(lp1.5) == 0)

lp1.50 = lp1.5
lp1.5 = prune_taxa(taxa_sums(lp1.5) > 0, lp1.5)
View(lp1.5)

#graph the rank abundance curve and reads per sample
readsumsdf = data.frame(nreads = sort(taxa_sums(lp1.5), TRUE), sorted = 1:ntaxa(lp1.5), type = "OTUs")

readsumsdf = rbind(readsumsdf, 
                   data.frame(nreads = sort(sample_sums(lp1.5), TRUE),
                              sorted = 1:nsamples(lp1.5),
                              type = "Samples"))

title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + facet_wrap(~type, 1, scales = "free") # add + scale_y_log10() to make the y axis log scale

lp1.5_filtered <- filter_taxa(lp1.5, function(x) sum(x > 0.010) > 0, TRUE)

lp1.5_fungi <- subset_taxa(lp1.5, Kingdom=="Fungi")


#Phylum
psL1.5 <- tax_glom(lp1.5_fungi, "Subphylum")
psL1.5a <- transform_sample_counts(psL1.5, function(x) x / sum(x))
psL1.5b <- merge_samples(psL1.5a, "treatment")
psL1.5c <- transform_sample_counts(psL1.5b, function(x) x / sum(x))

plot_bar(psL1.5c,fill = "Subphylum")+
  geom_bar(stat="identity")+
  theme(axis.title=element_text(size=16,face="bold"),
        axis.text.y = element_text(size = 16))

#Class
psL1.5 <- tax_glom(lp1.5_fungi, "Class")
psL1.5a <- transform_sample_counts(psL1.5, function(x) x / sum(x))
psL1.5b <- merge_samples(psL1.5a, "treatment")
psL1.5c <- transform_sample_counts(psL1.5b, function(x) x / sum(x))

plot_bar(psL1.5c,fill = "Class")+
  geom_bar(stat="identity")+
  theme(axis.title=element_text(size=16,face="bold"),
        axis.text.y = element_text(size = 16))

#Order
psL1.5 <- tax_glom(lp1.5_fungi, "Order")
psL1.5a <- transform_sample_counts(psL1.5, function(x) x / sum(x))
psL1.5b <- merge_samples(psL1.5a, "treatment")
psL1.5c <- transform_sample_counts(psL1.5b, function(x) x / sum(x))


plot_bar(psL1.5c,fill = "Order") +
  geom_bar(stat="identity")+
  theme(axis.title=element_text(size=16,face="bold"),
        axis.text.y = element_text(size = 16))


# plot bar graph with standard deviation as error bars
ggplot(avgs, aes(fill=treatment, x=Family, y=mean)) + 
  geom_bar(position='dodge') +
  
  plot_bar(psL2,fill = "Genus") +
  facet_wrap(~Genus) +
  geom_bar(stat="identity")+
  theme(axis.title=element_text(size=16,face="bold"),
        axis.text.y = element_text(size = 16),
        legend.position = "none")+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd))


# Genus

psL1.5 <- tax_glom(lp1.5_fungi, "Genus")
psL1.5a <- transform_sample_counts(psL1.5, function(x) x / sum(x))
psL1.5b <- merge_samples(psL1.5a, "treatment")
psL1.5c <- transform_sample_counts(psL1.5b, function(x) x / sum(x))

plot_bar(psL1.5c,fill = "Genus") +
  facet_wrap(~Genus) +
  geom_bar(stat="identity")+
  ylim(c(0,.35))+
  theme(legend.position = "none")

#Species

#subset to top most abundant OTUs
TopNOTUs <- names(sort(taxa_sums(lp1.5_fungi), TRUE)[1:13]) #adds any same names together 
lp1.5_top12fungi <- prune_taxa(TopNOTUs, lp1.5_fungi)

psL1.5 <- tax_glom(lp1.5_top12fungi, "Taxon")
psL1.5a <- transform_sample_counts(psL1.5, function(x) x / sum(x))
psL1.5b <- merge_samples(psL1.5a, "treatment")
psL1.5c <- transform_sample_counts(psL1.5b, function(x) x / sum(x))

plot_bar(psL1.5c,fill = "Taxon") +
  facet_wrap(~Taxon) +
  geom_bar(stat="identity")+
  ylim(c(0,.35))+
  theme(legend.position = "none")      

########### Litter Round 1, 17 month incubation #######

#any empty samples or taxa?
any(sample_sums(lp1.17) == 0)
any(taxa_sums(lp1.17) == 0)

lp1.170 = lp1.17
lp1.17 = prune_taxa(taxa_sums(lp1.17) > 0, lp1.17)
View(lp1.17)

#graph the rank abundance curve and reads per sample
readsumsdf = data.frame(nreads = sort(taxa_sums(lp1.17), TRUE), sorted = 1:ntaxa(lp1.17), type = "OTUs")

#lp1.17_filtered <- filter_taxa(lp1.17, function(x) sum(x > 0.010) > 0, TRUE)
lp1.17_fungi <- subset_taxa(lp1.17, Kingdom=="Fungi")


#Phylum
psL1.17 <- tax_glom(lp1.17_fungi, "Phylum")
psL1.17a <- transform_sample_counts(psL1.17, function(x) x / sum(x))
psL1.17b <- merge_samples(psL1.17a, "treatment")
psL1.17c <- transform_sample_counts(psL1.17b, function(x) x / sum(x))


plot_bar(psL1.17c,fill = "Phylum")+
  geom_bar(stat="identity")+
  theme(axis.title=element_text(size=16,face="bold"),
        axis.text.y = element_text(size = 16))




#Order

psL1.17 <- tax_glom(lp1.17_fungi, "Order")
psL1.17a <- transform_sample_counts(psL1.17, function(x) x / sum(x))
psL1.17b <- merge_samples(psL1.17a, "treatment")
psL1.17c <- transform_sample_counts(psL1.17b, function(x) x / sum(x))


plot_bar(psL1.17c,fill = "Order") +
  geom_bar(stat="identity")+
  theme(axis.title=element_text(size=16,face="bold"),
        axis.text.y = element_text(size = 16))


# Genus

psL1.17 <- tax_glom(lp1.17_fungi, "Genus")
psL1.17a <- transform_sample_counts(psL1.17, function(x) x / sum(x))
psL1.17b <- merge_samples(psL1.17a, "treatment")
psL1.17c <- transform_sample_counts(psL1.17b, function(x) x / sum(x))

#par(mfrow=c(3,5))


plot_bar(psL1.17c,fill = "Genus") +
  facet_wrap(~Genus) +
  geom_bar(stat="identity")+
  ylim(c(0,.35))+
  theme(legend.position = "none")

#Species


#subset to top most abundant genera
TopNOTUs <- names(sort(taxa_sums(lp1.17_fungi), TRUE)[1:12]) 
lp1.17_top12fungi <- prune_taxa(TopNOTUs, lp1.17_fungi)

psL1.17 <- tax_glom(lp1.17_top12fungi, "Taxon")
psL1.17a <- transform_sample_counts(psL1.17, function(x) x / sum(x))
psL1.17b <- merge_samples(psL1.17a, "treatment")
psL1.17c <- transform_sample_counts(psL1.17b, function(x) x / sum(x))

#par(mfrow=c(3,5))


plot_bar(psL1.17c,fill = "Taxon") +
  facet_wrap(~Taxon) +
  geom_bar(stat="identity")+
  ylim(c(0,.35))+
  theme(legend.position = "none")

########### Litter Round 2, 5 month incubation #######

#any empty samples or taxa?
any(sample_sums(lp2.5) == 0)
any(taxa_sums(lp2.5) == 0)

lp2.50 = lp2.5
lp2.5 = prune_taxa(taxa_sums(lp2.5) > 0, lp2.5)
View(lp2.5)

#graph the rank abundance curve and reads per sample
readsumsdf = data.frame(nreads = sort(taxa_sums(lp2.5), TRUE), sorted = 1:ntaxa(lp2.5), type = "OTUs")

readsumsdf = rbind(readsumsdf, 
                   data.frame(nreads = sort(sample_sums(lp2.5), TRUE),
                              sorted = 1:nsamples(lp2.5),
                              type = "Samples"))

title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + facet_wrap(~type, 1, scales = "free") # add + scale_y_log10() to make the y axis log scale

#lp2.5_filtered <- filter_taxa(lp2.5, function(x) sum(x > 0.010) > 0, TRUE)
lp2.5_fungi <- subset_taxa(lp2.5, Kingdom=="Fungi")


#Phylum
psL2.5 <- tax_glom(lp2.5_fungi, "Phylum")
psL2.5a <- transform_sample_counts(psL2.5, function(x) x / sum(x))
psL2.5b <- merge_samples(psL2.5a, "treatment")
psL2.5c <- transform_sample_counts(psL2.5b, function(x) x / sum(x))

plot_bar(psL2.5c,fill = "Phylum")+
  geom_bar(stat="identity")+
  theme(axis.title=element_text(size=16,face="bold"),
        axis.text.y = element_text(size = 16))

#Subphylum
psL2.5 <- tax_glom(lp2.5_fungi, "Subphylum")
psL2.5a <- transform_sample_counts(psL2.5, function(x) x / sum(x))
psL2.5b <- merge_samples(psL2.5a, "treatment")
psL2.5c <- transform_sample_counts(psL2.5b, function(x) x / sum(x))

plot_bar(psL2.5c,fill = "Subphylum")+
  geom_bar(stat="identity")+
  theme(axis.title=element_text(size=16,face="bold"),
        axis.text.y = element_text(size = 16))

#Order
psL2.5 <- tax_glom(lp2.5_fungi, "Order")
psL2.5a <- transform_sample_counts(psL2.5, function(x) x / sum(x))
psL2.5b <- merge_samples(psL2.5a, "treatment")
psL2.5c <- transform_sample_counts(psL2.5b, function(x) x / sum(x))


plot_bar(psL2.5c,fill = "Order") +
  geom_bar(stat="identity")+
  theme(axis.title=element_text(size=16,face="bold"),
        axis.text.y = element_text(size = 16))


# Genus

psL2.5 <- tax_glom(lp2.5_fungi, "Genus")
psL2.5a <- transform_sample_counts(psL2.5, function(x) x / sum(x))
psL2.5b <- merge_samples(psL2.5a, "treatment")
psL2.5c <- transform_sample_counts(psL2.5b, function(x) x / sum(x))

plot_bar(psL2.5c,fill = "Genus") +
  facet_wrap(~Genus) +
  geom_bar(stat="identity")+
  ylim(c(0,.35))+
  theme(legend.position = "none")

#Species

#subset to top most abundant OTUs
TopNOTUs <- names(sort(taxa_sums(lp2.5_fungi), TRUE)[1:12]) 
lp2.5_top12fungi <- prune_taxa(TopNOTUs, lp2.5_fungi)

psL2.5 <- tax_glom(lp2.5_top12fungi, "Taxon")
psL2.5a <- transform_sample_counts(psL2.5, function(x) x / sum(x))
psL2.5b <- merge_samples(psL2.5a, "treatment")
psL2.5c <- transform_sample_counts(psL2.5b, function(x) x / sum(x))

plot_bar(psL2.5c,fill = "Taxon") +
  facet_wrap(~Taxon) +
  geom_bar(stat="identity")+
  ylim(c(0,.35))+
  theme(legend.position = "none")


########### Litter Round 2, 17 month incubation #######

#any empty samples or taxa?
any(sample_sums(lp2.17) == 0)
any(taxa_sums(lp2.17) == 0)

lp2.170 = lp2.17
lp2.17 = prune_taxa(taxa_sums(lp2.17) > 0, lp2.17)
View(lp2.17)

#graph the rank abundance curve and reads per sample
readsumsdf = data.frame(nreads = sort(taxa_sums(lp2.17), TRUE), sorted = 1:ntaxa(lp2.17), type = "OTUs")

readsumsdf = rbind(readsumsdf, 
                   data.frame(nreads = sort(sample_sums(lp2.17), TRUE),
                              sorted = 1:nsamples(lp2.17),
                              type = "Samples"))

title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + facet_wrap(~type, 1, scales = "free") # add + scale_y_log10() to make the y axis log scale

#lp2.17_filtered <- filter_taxa(lp2.17, function(x) sum(x > 0.010) > 0, TRUE)
lp2.17_fungi <- subset_taxa(lp2.17, Kingdom=="Fungi")


#Phylum
psL2.17 <- tax_glom(lp2.17_fungi, "Phylum")
psL2.17a <- transform_sample_counts(psL2.17, function(x) x / sum(x))
psL2.17b <- merge_samples(psL2.17a, "treatment")
psL2.17c <- transform_sample_counts(psL2.17b, function(x) x / sum(x))

plot_bar(psL2.17c,fill = "Phylum")+
  geom_bar(stat="identity")+
  theme(axis.title=element_text(size=16,face="bold"),
        axis.text.y = element_text(size = 16))


#Subhylum
psL2.17 <- tax_glom(lp2.17_fungi, "Subhylum")
psL2.17a <- transform_sample_counts(psL2.17, function(x) x / sum(x))
psL2.17b <- merge_samples(psL2.17a, "treatment")
psL2.17c <- transform_sample_counts(psL2.17b, function(x) x / sum(x))

plot_bar(psL2.17c,fill = "Subphylum")+
  geom_bar(stat="identity")+
  theme(axis.title=element_text(size=16,face="bold"),
        axis.text.y = element_text(size = 16))


#Order
psL2.17 <- tax_glom(lp2.17_fungi, "Order")
psL2.17a <- transform_sample_counts(psL2.17, function(x) x / sum(x))
psL2.17b <- merge_samples(psL2.17a, "treatment")
psL2.17c <- transform_sample_counts(psL2.17b, function(x) x / sum(x))

plot_bar(psL2.17c,fill = "Order") +
  geom_bar(stat="identity")+
  theme(axis.title=element_text(size=16,face="bold"),
        axis.text.y = element_text(size = 16))

# Genus

psL2.17 <- tax_glom(lp2.17_fungi, "Genus")
psL2.17a <- transform_sample_counts(psL2.17, function(x) x / sum(x))
psL2.17b <- merge_samples(psL2.17a, "treatment")
psL2.17c <- transform_sample_counts(psL2.17b, function(x) x / sum(x))


plot_bar(psL2.17c,fill = "Genus") +
  facet_wrap(~Genus) +
  geom_bar(stat="identity")+
  ylim(c(0,.35))+
  theme(legend.position = "none")


#species

#subset to top most abundant OTUs
TopNOTUs <- names(sort(taxa_sums(lp2.17_fungi), TRUE)[1:12]) 
lp2.17_top12fungi <- prune_taxa(TopNOTUs, lp2.17_fungi)


psL2.17 <- tax_glom(lp2.17_top12fungi, "Taxon")
psL2.17a <- transform_sample_counts(psL2.17, function(x) x / sum(x))
psL2.17b <- merge_samples(psL2.17a, "treatment")
psL2.17c <- transform_sample_counts(psL2.17b, function(x) x / sum(x))


plot_bar(psL2.17c,fill = "Taxon") +
  facet_wrap(~Taxon) +
  geom_bar(stat="identity")+
  ylim(c(0,.35))+
  theme(legend.position = "none")





















######### HUMUS #########



#subset OTU table and format the OTUs for phyloseq
HumusOTU <- OTU[,meta$substrate == "Humus"]
HumusOTU_ps <- otu_table(HumusOTU,taxa_are_rows=T)

#subset metadata by substrate and start date and format for phyloseq
hmeta<- meta[meta$substrate == "Humus",]
hmeta_ps=sample_data(hmeta)

#make phyloseq dataset for just humus
humus_ps <- phyloseq(HumusOTU_ps,hmeta_ps,taxa)
View(humus_ps)

#any empty samples or taxa?
any(sample_sums(humus_ps) == 0)
any(taxa_sums(humus_ps) == 0)

humus_phyloseq0 = humus_ps
humus_phyloseq = prune_taxa(taxa_sums(humus_ps) > 0, humus_ps)

# hellinger transformation sqrt of relative abundance
humus_phyloseq <- transform_sample_counts(humus_ps, function(x) sqrt(x / sum(x)))

humus_phyloseq = prune_taxa(taxa_sums(humus_phyloseq > 0), humus_phyloseq)

# Now perform a unconstrained correspondence analysis
ca  <- ordinate(humus_phyloseq, "CCA")
# Scree plot
plot_scree(dca, "Scree Plot for Global Patterns Correspondence Analysis")

dca12 <- plot_ordination(humus_phyloseq, ca, "Unique_Sample", color="treatment",) + geom_point(size=5)

#Taxon per plot
psH <- tax_glom(hp2.17_fungi, "Taxon")
psHa <- transform_sample_counts(psH, function(x) x / sum(x))
psHb <- merge_samples(psHa, "Unique_Sample")
psHc <- transform_sample_counts(psHb, function(x) x / sum(x))

plot_bar(psHc,fill = "Taxon")+
  geom_bar(stat="identity")+
  theme(axis.title=element_text(size=16,face="bold"),
        axis.text.y = element_text(size = 16))






#subset by round
hmeta1 <- hmeta[hmeta$start_date == "2017-06-16",] #round 1 
hmeta2 <- hmeta[hmeta$start_date == "2018-06-14",] #round 2

#subset by incubation length in round 1
hmeta1.5 <- hmeta1[hmeta1$incubation_months == "5",] #5 months
hmeta1.17 <- hmeta1[hmeta1$incubation_months == "17",] #17 months

#round 2
hmeta2.5 <- hmeta2[hmeta2$incubation_months == "5",] #5 months
hmeta2.17 <- hmeta2[hmeta2$incubation_months == "17",] #17 months


#subset OTU table with metadata
HumusOTU1 <- HumusOTU[,hmeta$start_date == "2017-06-16"] #round 1
HumusOTU2 <- HumusOTU[,hmeta$start_date == "2018-06-14"] #round 2
str(HumusOTU1)
#round 1
hOTU1.5 <- HumusOTU1[hmeta1$incubation_months == "5",] #5 months
hOTU1.17 <- HumusOTU1[hmeta1$incubation_months == "17",] #17 months

#round 2
hOTU2.5 <- HumusOTU2[hmeta2$incubation_months == "5",] #5 months
hOTU2.17 <- HumusOTU2[hmeta2$incubation_months == "17",] #17 months


###### format OTU table and sample metadata for phyloseq #####

hOTU1.5=otu_table(hOTU1.5,taxa_are_rows=T)
dim(hOTU1.5);class(hOTU1.5);str(hOTU1.5)

#format samples for phyloseq
hmeta1.5=sample_data(hmeta1.5)
dim(hmeta1.5);class(hmeta1.5);str(hmeta1.5)
#View(lmeta1.5)

#make phyloseq object
hp1.5 <- phyloseq(hOTU1.5,hmeta1.5,taxa)
View(hp1.5)

#format OTU table for phyloseq
hOTU1.17=otu_table(hOTU1.17,taxa_are_rows=T)
#dim(hOTU1.17);class(hOTU1.17);str(hOTU1.17)

#format samples for phyloseq
hmeta1.17=sample_data(hmeta1.17)
dim(hmeta1.17);class(hmeta1.17);str(hmeta1.17)
#View(hmeta1.17)

#make phyloseq object
hp1.17 <- phyloseq(hOTU1.17,hmeta1.17,taxa)
View(hp1.17)


hOTU2.5=otu_table(hOTU2.5,taxa_are_rows=T)
#dim(hOTU2.5);class(hOTU2.5);str(hOTU2.5)

#format samples for phyloseq
hmeta2.5=sample_data(hmeta2.5)
dim(hmeta2.5);class(hmeta2.5);str(hmeta2.5)
View(hmeta2.5)

#make phyloseq object
hp2.5 <- phyloseq(hOTU2.5,hmeta2.5,taxa)
View(hp2.5)


#format OTU tabhe for phyloseq
hOTU2.17=otu_table(hOTU2.17,taxa_are_rows=T)
#dim(hOTU2.17);class(hOTU2.17);str(hOTU2.17)

#format samples for phyloseq
hmeta2.17=sample_data(hmeta2.17)
#dim(hmeta2.17);class(hmeta2.17);str(hmeta2.17)
View(hmeta2.17)

#make phyhoseq object
hp2.17 <- phyloseq(hOTU2.17,hmeta2.17,taxa)
View(hp2.17)

####### Humus Round 1, 5 month incubation ######

#any empty samples or taxa?
any(sample_sums(hp1.5) == 0)
any(taxa_sums(hp1.5) == 0)

hp1.50 = hp1.5
hp1.5 = prune_taxa(taxa_sums(hp1.5) > 0, hp1.5)
View(hp1.5)

#graph the rank abundance curve and reads per sample
readsumsdf = data.frame(nreads = sort(taxa_sums(hp1.5), TRUE), sorted = 1:ntaxa(hp1.5), type = "OTUs")

readsumsdf = rbind(readsumsdf, 
                   data.frame(nreads = sort(sample_sums(hp1.5), TRUE),
                              sorted = 1:nsamples(hp1.5),
                              type = "Samples"))

title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + facet_wrap(~type, 1, scales = "free") # add + scale_y_log10() to make the y axis log scale

#hp1.5_filtered <- filter_taxa(hp1.5, function(x) sum(x > 0.010) > 0, TRUE)
hp1.5_fungi <- subset_taxa(hp1.5, Kingdom=="Fungi")


#Phylum
psH1.5 <- tax_glom(hp1.5_fungi, "Subphylum")
psH1.5a <- transform_sample_counts(psH1.5, function(x) x / sum(x))
psH1.5b <- merge_samples(psH1.5a, "treatment")
psH1.5c <- transform_sample_counts(psH1.5b, function(x) x / sum(x))

plot_bar(psH1.5c,fill = "Subphylum")+
  geom_bar(stat="identity")+
  theme(axis.title=element_text(size=16,face="bold"),
        axis.text.y = element_text(size = 16))

#Class
psH1.5 <- tax_glom(hp1.5_fungi, "Class")
psH1.5a <- transform_sample_counts(psH1.5, function(x) x / sum(x))
psH1.5b <- merge_samples(psH1.5a, "treatment")
psH1.5c <- transform_sample_counts(psH1.5b, function(x) x / sum(x))

plot_bar(psH1.5c,fill = "Class")+
  geom_bar(stat="identity")+
  theme(axis.title=element_text(size=16,face="bold"),
        axis.text.y = element_text(size = 16))

#Order
psH1.5 <- tax_glom(hp1.5_fungi, "Order")
psH1.5a <- transform_sample_counts(psH1.5, function(x) x / sum(x))
psH1.5b <- merge_samples(psH1.5a, "treatment")
psH1.5c <- transform_sample_counts(psH1.5b, function(x) x / sum(x))


plot_bar(psH1.5c,fill = "Order") +
  geom_bar(stat="identity")+
  theme(axis.title=element_text(size=16,face="bold"),
        axis.text.y = element_text(size = 16))


# plot bar graph with standard deviation as error bars
ggplot(avgs, aes(fill=treatment, x=Family, y=mean)) + 
  geom_bar(position='dodge') +
  
  plot_bar(psH2,fill = "Genus") +
  facet_wrap(~Genus) +
  geom_bar(stat="identity")+
  theme(axis.title=element_text(size=16,face="bold"),
        axis.text.y = element_text(size = 16),
        legend.position = "none")+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd))


# Genus

psH1.5 <- tax_glom(hp1.5_fungi, "Genus")
psH1.5a <- transform_sample_counts(psH1.5, function(x) x / sum(x))
psH1.5b <- merge_samples(psH1.5a, "treatment")
psH1.5c <- transform_sample_counts(psH1.5b, function(x) x / sum(x))

plot_bar(psH1.5c,fill = "Genus") +
  facet_wrap(~Genus) +
  geom_bar(stat="identity")+
  ylim(c(0,.35))+
  theme(legend.position = "none")

#Species

#subset to top most abundant OTUs
TopNOTUs <- names(sort(taxa_sums(hp1.5_fungi), TRUE)[1:13]) #adds any same names together 
hp1.5_top12fungi <- prune_taxa(TopNOTUs, hp1.5_fungi)

psH1.5 <- tax_glom(hp1.5_top12fungi, "Taxon")
psH1.5a <- transform_sample_counts(psH1.5, function(x) x / sum(x))
psH1.5b <- merge_samples(psH1.5a, "treatment")
psH1.5c <- transform_sample_counts(psH1.5b, function(x) x / sum(x))

plot_bar(psH1.5c,fill = "Taxon") +
  facet_wrap(~Taxon) +
  geom_bar(stat="identity")+
  ylim(c(0,.35))+
  theme(legend.position = "none")  

# Guild

psH1.5 <- tax_glom(hp1.5_fungi, "Guild")
psH1.5a <- transform_sample_counts(psH1.5, function(x) x / sum(x))
psH1.5b <- merge_samples(psH1.5a, "treatment")
psH1.5c <- transform_sample_counts(psH1.5b, function(x) x / sum(x))

guild1.5 <- plot_bar(psH1.5c,fill = "Guild") +
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("#56B4E9","#009E73","#E69F00", "#D55E00","grey"))+theme(legend.position="none",text = element_text(size=30))

####### Humus Round 1, 17 month incubation ######

#any empty samples or taxa?
any(sample_sums(hp1.17) == 0)
any(taxa_sums(hp1.17) == 0)

hp1.170 = hp1.17
hp1.17 = prune_taxa(taxa_sums(hp1.17) > 0, hp1.17)
#View(hp1.17)

#graph the rank abundance curve and reads per sample
readsumsdf = data.frame(nreads = sort(taxa_sums(hp1.17), TRUE), sorted = 1:ntaxa(hp1.17), type = "OTUs")

readsumsdf = rbind(readsumsdf, 
                   data.frame(nreads = sort(sample_sums(hp1.17), TRUE),
                              sorted = 1:nsamples(hp1.17),
                              type = "Samples"))

title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + facet_wrap(~type, 1, scales = "free") # add + scale_y_log10() to make the y axis log scale

#hp1.17_filtered <- filter_taxa(hp1.17, function(x) sum(x > 0.010) > 0, TRUE)
hp1.17_fungi <- subset_taxa(hp1.17, Kingdom=="Fungi")


#Phylum
psH1.17 <- tax_glom(hp1.17_fungi, "Subphylum")
psH1.17a <- transform_sample_counts(psH1.17, function(x) x / sum(x))
psH1.17b <- merge_samples(psH1.17a, "treatment")
psH1.17c <- transform_sample_counts(psH1.17b, function(x) x / sum(x))

plot_bar(psH1.17c,fill = "Subphylum")+
  geom_bar(stat="identity")+
  theme(axis.title=element_text(size=16,face="bold"),
        axis.text.y = element_text(size = 16))

#Class
psH1.17 <- tax_glom(hp1.17_fungi, "Class")
psH1.17a <- transform_sample_counts(psH1.17, function(x) x / sum(x))
psH1.17b <- merge_samples(psH1.17a, "treatment")
psH1.17c <- transform_sample_counts(psH1.17b, function(x) x / sum(x))

plot_bar(psH1.17c,fill = "Class")+
  geom_bar(stat="identity")+
  theme(axis.title=element_text(size=16,face="bold"),
        axis.text.y = element_text(size = 16))


#Order
psH1.17 <- tax_glom(hp1.17_fungi, "Order")
psH1.17a <- transform_sample_counts(psH1.17, function(x) x / sum(x))
psH1.17b <- merge_samples(psH1.17a, "treatment")
psH1.17c <- transform_sample_counts(psH1.17b, function(x) x / sum(x))


plot_bar(psH1.17c,fill = "Order") +
  geom_bar(stat="identity")+
  theme(axis.title=element_text(size=16,face="bold"),
        axis.text.y = element_text(size = 16))


# Genus

psH1.17 <- tax_glom(hp1.17_fungi, "Genus")
psH1.17a <- transform_sample_counts(psH1.17, function(x) x / sum(x))
psH1.17b <- merge_samples(psH1.17a, "treatment")
psH1.17c <- transform_sample_counts(psH1.17b, function(x) x / sum(x))

plot_bar(psH1.17c,fill = "Genus") +
  facet_wrap(~Genus) +
  geom_bar(stat="identity")+
  ylim(c(0,1))+
  theme(legend.position = "none")

#Species

#subset to top most abundant OTUs
TopNOTUs <- names(sort(taxa_sums(hp1.17_fungi), TRUE)[1:13]) #adds any same names together 
hp1.17_top12fungi <- prune_taxa(TopNOTUs, hp1.17_fungi)

psH1.17 <- tax_glom(hp1.17_top12fungi, "Taxon")
psH1.17a <- transform_sample_counts(psH1.17, function(x) x / sum(x))
psH1.17b <- merge_samples(psH1.17a, "treatment")
psH1.17c <- transform_sample_counts(psH1.17b, function(x) x / sum(x))

plot_bar(psH1.17c,fill = "Taxon") +
  facet_wrap(~Taxon) +
  geom_bar(stat="identity")+
  ylim(c(0,1))+
  theme(legend.position = "none")  

# Guild

psH1.17 <- tax_glom(hp1.17_fungi, "Guild")
psH1.17a <- transform_sample_counts(psH1.17, function(x) x / sum(x))
psH1.17b <- merge_samples(psH1.17a, "treatment")
psH1.17c <- transform_sample_counts(psH1.17b, function(x) x / sum(x))

guild1.17 <- plot_bar(psH1.17c,fill = "Guild") +
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("#56B4E9","#009E73","#F0E442", "#D55E00","grey"))+theme(legend.position = "none",text = element_text(size=30))

########### Humus Round 2, 5 month incubation #######

#any empty samples or taxa?
any(sample_sums(hp2.5) == 0)
any(taxa_sums(hp2.5) == 0)

hp2.50 = hp2.5
hp2.5 = prune_taxa(taxa_sums(hp2.5) > 0, hp2.5)
#View(hp2.5)

#graph the rank abundance curve and reads per sample
readsumsdf = data.frame(nreads = sort(taxa_sums(hp2.5), TRUE), sorted = 1:ntaxa(hp2.5), type = "OTUs")

readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(hp2.5), TRUE),
                                          sorted = 1:nsamples(hp2.5),
                                          type = "Samples", title = "Total number of reads")
                   p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
                   p + ggtitle(title) + facet_wrap(~type, 1, scales = "free") # add + scale_y_log10() to make the y axis log scale
                   
                   #hp2.5_filtered <- filter_taxa(hp2.5, function(x) sum(x > 0.010) > 0, TRUE)
                   
                   
                   hp2.5_fungi <- subset_taxa(hp2.5, Kingdom=="Fungi")
                   
                   
                   #Phylum
                   psH2.5 <- tax_glom(hp2.5_fungi, "Phylum")
                   psH2.5a <- transform_sample_counts(psH2.5, function(x) x / sum(x))
                   psH2.5b <- merge_samples(psH2.5a, "treatment")
                   psH2.5c <- transform_sample_counts(psH2.5b, function(x) x / sum(x))
                   
                   plot_bar(psH2.5c,fill = "Phylum")+
                     geom_bar(stat="identity")+
                     theme(axis.title=element_text(size=16,face="bold"),
                           axis.text.y = element_text(size = 16))
                   
                   #Subphylum
                   psH2.5 <- tax_glom(hp2.5_fungi, "Subphylum")
                   psH2.5a <- transform_sample_counts(psH2.5, function(x) x / sum(x))
                   psH2.5b <- merge_samples(psH2.5a, "treatment")
                   psH2.5c <- transform_sample_counts(psH2.5b, function(x) x / sum(x))
                   
                   plot_bar(psH2.5c,fill = "Subphylum")+
                     geom_bar(stat="identity")+
                     theme(axis.title=element_text(size=16,face="bold"),
                           axis.text.y = element_text(size = 16))
                   
                   #Class
                   psH2.5 <- tax_glom(hp2.5_fungi, "Class")
                   psH2.5a <- transform_sample_counts(psH2.5, function(x) x / sum(x))
                   psH2.5b <- merge_samples(psH2.5a, "treatment")
                   psH2.5c <- transform_sample_counts(psH2.5b, function(x) x / sum(x))
                   
                   
                   plot_bar(psH2.5c,fill = "Class") +
                     geom_bar(stat="identity")+
                     theme(axis.title=element_text(size=16,face="bold"),
                           axis.text.y = element_text(size = 16))
                   
                   
                   
                   #Order
                   psH2.5 <- tax_glom(hp2.5_fungi, "Order")
                   psH2.5a <- transform_sample_counts(psH2.5, function(x) x / sum(x))
                   psH2.5b <- merge_samples(psH2.5a, "treatment")
                   psH2.5c <- transform_sample_counts(psH2.5b, function(x) x / sum(x))
                   
                   
                   plot_bar(psH2.5c,fill = "Order") +
                     geom_bar(stat="identity")+
                     theme(axis.title=element_text(size=16,face="bold"),
                           axis.text.y = element_text(size = 16))
                   
                   
                   # Genus
                   
                   psH2.5 <- tax_glom(hp2.5_fungi, "Genus")
                   psH2.5a <- transform_sample_counts(psH2.5, function(x) x / sum(x))
                   psH2.5b <- merge_samples(psH2.5a, "treatment")
                   psH2.5c <- transform_sample_counts(psH2.5b, function(x) x / sum(x))
                   
                   plot_bar(psH2.5c,fill = "Genus") +
                     facet_wrap(~Genus) +
                     geom_bar(stat="identity")+
                     ylim(c(0,1))+
                     theme(legend.position = "none")
                   
                   #Species
                   
                   #subset to top most abundant OTUs
                   TopNOTUs <- names(sort(taxa_sums(hp2.5_fungi), TRUE)[1:12]) 
                   hp2.5_top12fungi <- prune_taxa(TopNOTUs, hp2.5_fungi)
                   
                   psH2.5 <- tax_glom(hp2.5_fungi, "Taxon")
                   psH2.5a <- transform_sample_counts(psH2.5, function(x) x / sum(x))
                   psH2.5b <- merge_samples(psH2.5a, "treatment")
                   psH2.5c <- transform_sample_counts(psH2.5b, function(x) x / sum(x))
                   
                   plot_bar(psH2.5c,fill = "Taxon") +
                     facet_wrap(~Taxon) +
                     geom_bar(stat="identity")+
                     ylim(c(0,1))+
                     theme(legend.position = "none")
                   
                   # Guild
                   
                   psH2.5 <- tax_glom(hp2.5_fungi, "Guild")
                   psH2.5a <- transform_sample_counts(psH2.5, function(x) x / sum(x))
                   psH2.5b <- merge_samples(psH2.5a, "treatment")
                   psH2.5c <- transform_sample_counts(psH2.5b, function(x) x / sum(x))
                   
                   guild2.5 <- plot_bar(psH2.5c,fill = "Guild") +
                     geom_bar(stat="identity")+
                     scale_fill_manual(values = c("#56B4E9","#009E73", "#D55E00","grey"))+
                     theme(legend.position = "none",text = element_text(size=30))
                   
                   
                   
########### Humus Round 2, 17 month incubation #######
                   
                   #any empty samples or taxa?
                   any(sample_sums(hp2.17) == 0)
                   any(taxa_sums(hp2.17) == 0)
                   
                   hp2.170 = hp2.17
                   hp2.17 = prune_taxa(taxa_sums(hp2.17) > 0, hp2.17)
                   View(hp2.17)
                   
                   #graph the rank abundance curve and reads per sample
                   readsumsdf = data.frame(nreads = sort(taxa_sums(hp2.17), TRUE), sorted = 1:ntaxa(hp2.17), type = "OTUs")
                   
                   readsumsdf = rbind(readsumsdf, 
                                      data.frame(nreads = sort(sample_sums(hp2.17), TRUE),
                                                 sorted = 1:nsamples(hp2.17),
                                                 type = "Samples"))
                   
                   title = "Total number of reads"
                   p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
                   p + ggtitle(title) + facet_wrap(~type, 1, scales = "free") # add + scale_y_log10() to make the y axis log scale
                   
                   #hp2.17_filtered <- filter_taxa(hp2.17, function(x) sum(x > 0.010) > 0, TRUE)
                   hp2.17_fungi <- subset_taxa(hp2.17, Kingdom=="Fungi")
                   
                   
                   #Phylum
                   psH2.17 <- tax_glom(hp2.17_fungi, "Phylum")
                   psH2.17a <- transform_sample_counts(psH2.17, function(x) x / sum(x))
                   psH2.17b <- merge_samples(psH2.17a, "treatment")
                   psH2.17c <- transform_sample_counts(psH2.17b, function(x) x / sum(x))
                   
                   plot_bar(psH2.17c,fill = "Phylum")+
                     geom_bar(stat="identity")+
                     theme(axis.title=element_text(size=16,face="bold"),
                           axis.text.y = element_text(size = 16))
                   
                   #Subphylum
                   psH2.17 <- tax_glom(hp2.17_fungi, "Subphylum")
                   psH2.17a <- transform_sample_counts(psH2.17, function(x) x / sum(x))
                   psH2.17b <- merge_samples(psH2.17a, "treatment")
                   psH2.17c <- transform_sample_counts(psH2.17b, function(x) x / sum(x))
                   
                   plot_bar(psH2.17c,fill = "Subphylum")+
                     geom_bar(stat="identity")+
                     theme(axis.title=element_text(size=16,face="bold"),
                           axis.text.y = element_text(size = 16))
                   
                   
                   #Class
                   psH2.17 <- tax_glom(hp2.17_fungi, "Class")
                   psH2.17a <- transform_sample_counts(psH2.17, function(x) x / sum(x))
                   psH2.17b <- merge_samples(psH2.17a, "treatment")
                   psH2.17c <- transform_sample_counts(psH2.17b, function(x) x / sum(x))
                   
                   
                   plot_bar(psH2.17c,fill = "Class") +
                     geom_bar(stat="identity")+
                     theme(axis.title=element_text(size=16,face="bold"),
                           axis.text.y = element_text(size = 16))
                   
                   #Order
                   psH2.17 <- tax_glom(hp2.17_fungi, "Order")
                   psH2.17a <- transform_sample_counts(psH2.17, function(x) x / sum(x))
                   psH2.17b <- merge_samples(psH2.17a, "treatment")
                   psH2.17c <- transform_sample_counts(psH2.17b, function(x) x / sum(x))
                   
                   
                   plot_bar(psH2.17c,fill = "Order") +
                     geom_bar(stat="identity")+
                     theme(axis.title=element_text(size=16,face="bold"),
                           axis.text.y = element_text(size = 16))
                   
                   
                   # Genus
                   
                   psH2.17 <- tax_glom(hp2.17_fungi, "Genus")
                   psH2.17a <- transform_sample_counts(psH2.17, function(x) x / sum(x))
                   psH2.17b <- merge_samples(psH2.17a, "treatment")
                   psH2.17c <- transform_sample_counts(psH2.17b, function(x) x / sum(x))
                   
                   plot_bar(psH2.17c,fill = "Genus") +
                     facet_wrap(~Genus) +
                     geom_bar(stat="identity")+
                     ylim(c(0,.35))+
                     theme(legend.position = "none")
                   
                   #Species
                   
                   #subset to top most abundant OTUs
                   TopNOTUs <- names(sort(taxa_sums(hp2.17_fungi), TRUE)[1:30]) 
                   hp2.17_top12fungi <- prune_taxa(TopNOTUs, hp2.17_fungi)
                   
                   psH2.17 <- tax_glom(hp2.17_top12fungi, "Taxon")
                   psH2.17a <- transform_sample_counts(psH2.17, function(x) x / sum(x))
                   psH2.17b <- merge_samples(psH2.17a, "treatment")
                   psH2.17c <- transform_sample_counts(psH2.17b, function(x) x / sum(x))
                   
                   plot_bar(psH2.17c,fill = "Taxon") +
                     facet_wrap(~Taxon) +
                     geom_bar(stat="identity")+
                     ylim(c(0,1))+
                     theme(legend.position = "none")
                   
                   # Guild
                   
                   psH2.17 <- tax_glom(hp2.17_fungi, "Guild")
                   psH2.17a <- transform_sample_counts(psH2.17, function(x) x / sum(x))
                   psH2.17b <- merge_samples(psH2.17a, "treatment")
                   psH2.17c <- transform_sample_counts(psH2.17b, function(x) x / sum(x))
                   
                   guild2.17 <- plot_bar(psH2.17c,fill = "Guild") +
                     geom_bar(stat="identity")+
                     scale_fill_manual(values = c("#56B4E9","#009E73","#E69F00","#F0E442","#D55E00","grey"))+
                     theme(legend.position = "none", text = element_text(size=30))
                   
                   
                   
                   ### Combined Graphs ###
                   
                   #combining plots into display
                   require(gridExtra)
                   grid.arrange(guild1.5, guild1.17, guild2.5, guild2.17, nrow=2, padding = c(0,0))
                   
    

                   
                   
#################################mvabund  #####
                   
                   library(mvabund)
                   setwd("C:/Users/csca0001/Projects/Padjelanta project/Inorg_org_N_paper/Fungal community/permutate/mvabund")
                   
                   community <- read.csv(file = "community_mvabund1.csv", header = TRUE)
                   summary(community)
                   View(community)
                   Herb_spp <- mvabund(community[,6:600])
                   summary(Herb_spp)
                   
                   
                   library(vegan)
                   library(mvabund)
                   library(reshape2)
                   library(plyr)
                   library(ggplot2)
                   #########################NC, carbon and nitrogen
                   
                   View(community)
                   mod1 <- manyglm(Herb_spp ~ N_C, data=community,family="poisson", strata= Type)
                   anova(mod1)
                   
                   mod1 <- manyglm(Herb_spp ~ N_C, data=community,family="poisson")
                   anova(mod1)
                   
                   mod1 <- manyglm(Herb_spp ~ C, data=community,family="poisson")
                   anova(mod1)
                   
                   anova_mod1<-anova(mod1, nBoot=199, test="wald", p.uni = "adjusted")
                   View(anova_mod1)
                   View(anova_mod1$uni.p)
                   write.csv(anova_mod1$uni.p, file="mvabund_all_NC.csv")
                   
                   mod1 <- manyglm(Herb_spp ~ C, data=community,family="poisson")
                   anova_mod1<-anova(mod1, nBoot=199, test="wald", p.uni = "adjusted")
                   plot(mod1)
                   write.csv(anova_mod1$uni.p, file="mvabund_all_C.csv")
                   
                   mod1 <- manyglm(Herb_spp ~ N, data=community,family="poisson")
                   anova_mod1<-anova(mod1, nBoot=199, test="wald", p.uni = "adjusted")
                   plot(mod1)
                   write.csv(anova_mod1$uni.p, file="mvabund_all_N.csv")
                   
                   
                   ####################################with genera
                   community <- read.csv(file = "community_mvabund2_genera1.csv", header = TRUE)
                   View(community)
                   fungi<- mvabund(community[,8:600])
                   View(community)
                   
                   mod1 <- manyglm(fungi ~ community$N_C, data= community, family="poisson")
                   summary(mod1,nBoot = 999,test = "LR")
                   #anova_mod1<-anova(mod1, nBoot=199, test="wald", p.uni = "adjusted")
                   anova_mod1<-anova(mod1, nBoot=199, test="LR", p.uni = "adjusted")
                   plot(mod1)
                   qqnorm(resid(mod_pois))
                   
                   #not good with poisson....
                   
                   mod_binom <- manyglm(fungi ~ community$N_C, data = community, family = 'negative_binomial')
                   anova_mod1<-anova(mod_binom, nBoot=199, test="wald", p.uni = "adjusted")
                   plot(mod_binom)
                   qqnorm(resid(mod_binom))
                   
                   write.csv(mod_binom$uni.p, file="mvabund_genera_all_C.csv")
                   
                   
                   
                   ###############with total abundance
                   community <- read.csv(file = "community_mvabund2_forest.csv", header = TRUE)
                   View(community)
                   
                   Herb_spp <- mvabund(community[,8:247])
                   par(mar=c(2,10,2,2)) # adjusts the margins
                   boxplot(community[,8:247],horizontal = TRUE,las=2, main="frequency")
                   
                   mod1 <- manyglm(Herb_spp ~ N_C, data=community,family="poisson")
                   anova(mod2, p.uni="adjusted")
                   plot(mod1)
                   
                   ###############with rel ab for heath
                   community <- read.csv(file = "community_mvabund2_heath_relab.csv", header = TRUE)
                   View(community)
                   
                   Herb_spp <- mvabund(community[,8:213])
                   
                   par(mar=c(2,10,2,2)) # adjusts the margins
                   boxplot(community[,8:247],horizontal = TRUE,las=2, main="frequency")
                   
                   mod1 <- manyglm(Herb_spp ~ N_C, data=community,family="poisson")
                   anova_mod1<-anova(mod1, p.uni="adjusted")
                   plot(mod1)
                   
                   mod1 <- manyglm(Herb_spp ~ N_C, data=community,family="negative_binomial")
                   anova_mod1<-anova(mod1, p.uni="adjusted")
                   plot(mod1)
                   
                   ###mod_pois_anova<-anova(mod2, nBoot=199, test="wald", p.uni = "adjusted")
                   plot(mod2)
                   write.csv(anova_mod1$uni.p, file="mvabund_pvalues_heath_relab.csv")
                   
                   #########################now carbon and nitrogen
                   
                   mod1 <- manyglm(Herb_spp ~ C, data=community,family="poisson")
                   anova_mod1<-anova(mod1, p.uni="adjusted")
                   plot(mod1)
                   write.csv(anova_mod1$uni.p, file="mvabund_pvalues_heath_relab_C.csv")
                   
                   mod1 <- manyglm(Herb_spp ~ N, data=community,family="poisson")
                   anova_mod1<-anova(mod1, p.uni="adjusted")
                   plot(mod1)
                   write.csv(anova_mod1$uni.p, file="mvabund_pvalues_heath_relab_N.csv")
                   
                   
                   ###############with rel ab for forest
                   community <- read.csv(file = "community_mvabund2_forest_relab.csv", header = TRUE)
                   View(community)
                   
                   Herb_spp <- mvabund(community[,8:374])
                   
                   par(mar=c(2,10,2,2)) # adjusts the margins
                   boxplot(community[,8:247],horizontal = TRUE,las=2, main="frequency")
                   
                   mod1 <- manyglm(Herb_spp ~ N_C, data=community,family="poisson")
                   anova_mod1<-anova(mod1, p.uni="adjusted")
                   plot(mod1)
                   
                   mod1 <- manyglm(Herb_spp ~ N_C, data=community,family="negative_binomial")
                   anova_mod1<-anova(mod1, p.uni="adjusted")
                   plot(mod1)
                   
                   ###mod_pois_anova<-anova(mod2, nBoot=199, test="wald", p.uni = "adjusted")
                   plot(mod2)
                   write.csv(anova_mod1$uni.p, file="mvabund_pvalues_forest_relab_N_C.csv")
                   
                   #########################now carbon and nitrogen
                   
                   mod1 <- manyglm(Herb_spp ~ C, data=community,family="poisson")
                   anova_mod1<-anova(mod1, p.uni="adjusted")
                   plot(mod1)
                   write.csv(anova_mod1$uni.p, file="mvabund_pvalues_forest_relab_C.csv")
                   
                   mod1 <- manyglm(Herb_spp ~ N, data=community,family="poisson")
                   anova_mod1<-anova(mod1, nBoot=199, test="wald")
                   anova_mod1<-anova(mod1, p.uni="adjusted")
                   plot(mod1)
                   write.csv(anova_mod1$uni.p, file="mvabund_pvalues_forest_relab_N.csv")
                   
                   
                   
                   
                   
###### CCA #######
                   
                   # Now perform a unconstrained correspondence analysis
                   gpca  <- ordinate(humus_phyloseq, "CCA")
                   # Scree plot
                   plot_scree(gpca, "Scree Plot for Global Patterns Correspondence Analysis")
                   
                   plot_ordination(humus_phyloseq, gpca, "Sample_Unique", color="treatment") + 
                     geom_line() + geom_point(size=5)
                   GP
                   data("GlobalPatterns")
                   View(GlobalPatterns)
                   