## Jädraås Mycorrhizal Removal Experiment
## 
##
#Author: Louis Mielke
#Date Created: 2021 APRIL 09
#Date Updated: 2021 DEC 01

rm(list=ls()) #clear the workspace (only if you want to)

##### load librarys ####
library(phyloseq);library(dplyr);library(tidyverse);library(tibble);library(viridis)
library(ggplot2);library(vegan);library(ggsci);library(tidyr);library(scales)
library(gridExtra);library(grid);library(lattice);library(MASS);library(lme4)
library(nlme);library(emmeans);library(gamm4);library(dplyr);library(magrittr)
library(reshape2);library(mgcv);library(lmerTest);library(lubridate);library(car)
library(lmerTest);library(openxlsx)


##### set wd #####
setwd("~/Projects/mycorrhizal removal/MeshBags/Jädraås_Gadgil_MycorrhizalTraits")

#### import metadata ####
meta <- read.table("metadata.txt", sep="\t", header=T, row.names=1)
class(meta);dim(meta);str(meta) #meta$start_date <- as.Date(meta$start_date)
meta$incubation <- as.factor(meta$incubation) #incubation months  5 or 17
meta$block <- as.factor(meta$block) 

meta=as.data.frame(meta)
class(meta);str(meta);dim(meta)

### subsetting 
meta$plot <- paste(meta$treatment,meta$block, sep = "")

#subset by substrate (Humus or Litter) for plotting
humus0 <- subset(meta, substrate %in% c("Humus"))
litter0 <- subset(meta, substrate %in% c("Litter"))

#subset by treatment for plotting: C = Control, DC = Disturbed Control, E = Shrub Removal, T = Pine Root Exclusion, TE = Pine Root Exclusion & Shrub Removal

humus <- subset(humus0, treatment %in% c("C","DC","E","T","TE")) # dataset for graphing
litter <- subset(litter0, treatment %in% c("C","DC","E","T","TE")) #dataset for graphing

#subset by treatment for statistics
meta2 <- subset(meta, treatment %in% c("C","E","T","TE")) # not including DC
humus2 <- subset(meta2, substrate %in% c("Humus")) # dataset for statistics
litter2 <- subset(meta2, substrate %in% c("Litter")) # dataset for statistics


####### humus mass loss graphs ########

## decomposition curves for humus data combined
## Humus Round 1 2017-2018 blue and green bags deployed in june 2017 - Nov 2018
## Humus Round 2 2018-2019 blue and green bags deployed in june 2018 - Nov 2019


#summary for percent mass remaining for each treatment
humus_summary <- humus %>%
  group_by(set, incubation, treatment) %>% 
  drop_na() %>% 
  summarize_at(c("percent_mass_remaining"), funs(mean, sd, min, max,n()))

#standard error
humus_summary$se <- humus_summary$sd/sqrt(humus_summary$n)
str(humus_summary)


#plotting

#make dummy variables for plotting
dummy <- tibble(incubation = 0, 
                treatment = c("C","DC","E","T","TE"), 
                mean = 100,sd = 0,
                min =0,max=0,n=0,se=0)

dummy$treatment <- as.factor(dummy$treatment)
dummy$incubation<- as.factor(dummy$incubation)
dummy$n <- as.integer(dummy$n)


#BY SET: summary for percent mass remaining for each treatment by set
humus_summary_byset <- humus %>%
  group_by(set, incubation, treatment) %>% 
  drop_na() %>% 
  summarize_at(c("percent_mass_remaining"), funs(mean, sd, min, max,n()))

#standard error
humus_summary_byset$se <- humus_summary_byset$sd/sqrt(humus_summary_byset$n)

str(humus_summary_byset)

humus_1 <- humus_summary_byset[humus_summary_byset$set == "A",]
humus_1 <- bind_rows(humus_1,dummy)
humus_1$incubation <- as.numeric(as.character(humus_1$incubation))

#Set 1
hA <- ggplot(data = humus_1, aes(x = incubation, y = mean, color = treatment, linetype = treatment))+
  geom_line(size = 2)+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 2, position = position_dodge(width=0.75), linetype = "solid")+
  scale_colour_manual(values = c("mediumseagreen", "grey", "black", "mediumseagreen","black")) +
  xlab("Months")+
  ylab("% Mass Remaining")+
  ylim(c(65,100))+
  theme_classic()+
  theme(legend.position="none", axis.title.x = element_blank(), axis.text.x=element_blank())+
  scale_linetype_manual(values = c("solid", "dotted", "solid", "dotted", "dotted"))

humus_2 <- humus_summary_byset[humus_summary_byset$set == "B",]
humus_2 <- bind_rows(humus_2,dummy)
humus_2$incubation <- as.numeric(as.character(humus_2$incubation))

#set B
hB <- ggplot(data = humus_2, aes(x = incubation, y = mean, color = treatment, linetype = treatment))+
  geom_line(size = 1.5)+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 2, position = position_dodge(width=0.75), linetype = "solid")+
  scale_colour_manual(values = c("mediumseagreen", "grey", "black", "mediumseagreen","black")) +
  xlab("Months")+
  ylab("% Mass Remaining")+
  ylim(c(65,100))+
  theme_classic()+
  theme(legend.position="none",axis.title.x=element_blank(), axis.text.x=element_blank(),axis.title.y=element_blank(), axis.text.y=element_blank())+
  scale_linetype_manual(values = c("solid", "dotted", "solid", "dotted", "dotted"))


####### pine needle litter mass loss graphs ########


#summary for percent mass remaining for each treatment
litter_summary <- litter %>%
  group_by(incubation, treatment) %>% 
  drop_na() %>% 
  summarize_at(c("percent_mass_remaining"), funs(mean, sd, min, max,n()))

#standard error
litter_summary$se <- litter_summary$sd/sqrt(litter_summary$n)

str(litter_summary)

#combine data
litter_summary <- bind_rows(litter_summary,dummy)
litter_summary$incubation <- as.numeric(as.character(litter_summary$incubation))

### PLOTTING 

#average of both rounds
ggplot(data = litter_summary, aes(x = incubation, y = mean, color = treatment, linetype = treatment))+
  geom_line(size = 2)+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 2.5, position =position_dodge(width=0.75), linetype = "solid")+
  scale_colour_manual(values = c("mediumseagreen", "grey", "black", "mediumseagreen","black")) +
  xlab("Months")+
  ylab("% Mass Remaining")+
  ylim(c(35,100))+
  theme_classic()+
  theme(legend.position="top")+
  scale_linetype_manual(values = c("solid", "dotted", "solid", "dotted", "dotted"))


#by Set: summary for percent mass remaining for each treatment AND set
litter_summary_byset <- litter %>%
  group_by(set, incubation, treatment) %>% 
  drop_na() %>% 
  summarize_at(c("percent_mass_remaining"), funs(mean, sd, min, max,n()))

#standard error
litter_summary_byset$se <- litter_summary_byset$sd/sqrt(litter_summary_byset$n)


#Set A
litter1 <- litter_summary_byset[litter_summary_byset$set == "A",]

#combine data
litter1 <- bind_rows(litter1,dummy)
litter1$incubation <- as.numeric(as.character(litter1$incubation))

lA <- ggplot(data = litter1, aes(x = incubation, y = mean, color = treatment, linetype = treatment))+
  geom_line(size = 2)+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 2, position = position_dodge(width=0.75), linetype = "solid")+
  scale_colour_manual(values = c("mediumseagreen", "grey", "black", "mediumseagreen","black"))+
  xlab("Months")+
  ylab("% Mass Remaining")+
  theme_classic()+
  ylim(c(30,100))+
  theme(legend.position="none")+
  scale_linetype_manual(values = c("solid", "dotted", "solid", "dotted", "dotted"))

#Set B
litter2 <- litter_summary_byset[litter_summary_byset$set == "B",]

#combine data
litter2 <- bind_rows(litter2,dummy)
litter2$incubation <- as.numeric(as.character(litter2$incubation))

#plot
lB <- ggplot(data = litter2, aes(x = incubation, y = mean, color = treatment, linetype = treatment))+
  geom_line(size = 2)+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 2, position = position_dodge(width=0.75), linetype = "solid")+
  scale_colour_manual(values = c("mediumseagreen", "grey", "black", "mediumseagreen","black")) +
  xlab("Months")+
  #ylab("% Mass Remaining")+
  ylim(c(30,100))+
  theme_classic()+
  theme(legend.position="none", axis.title.y=element_blank(), axis.text.y=element_blank())+
  scale_linetype_manual(values = c("solid", "dotted", "solid", "dotted", "dotted"))






#######
#combining plots into display
require(gridExtra)
grid.arrange(hA, hB, lA, lB, padding = c(0,0))


########### humus mass loss statistics #######

lme.h_mass <-lme(percent_loss ~ pine*shrubs*incubation*set, #
                 random = ~ 1|block/plot,
                 data = humus2, na.action = na.omit)

h_mass_anova <- anova(lme.h_mass)
coef_h <- coef(lme.h_mass)
plot(resid(lme.h_mass))

#subset to each set A and B
humusA <- humus2[humus2$set == "A",]
humusB <- humus2[humus2$set == "B",]

#set A
lme.h_massA <-lme(percent_loss ~ pine*shrubs*incubation, #
                 random = ~ 1|block/plot,
                 data = humusA, na.action = na.omit)

h_mass_anovaA <- anova(lme.h_massA)
coef_hA <- coef(lme.h_massA)
plot(resid(lme.h_massA))


        #subset to each set incubation period
        humus_A5 <- humusA[humusA$incubation == "5",]
        humus_A17 <- humusA[humusA$incubation == "17",]

        #Humus Set A 5 month incubation
        lme.h_massA5 <-lme(percent_loss ~ pine*shrubs, #
                          random = ~ 1|block,
                          data = humus_A5, na.action = na.omit)
        
        h_mass_anovaA5 <- anova(lme.h_massA5)
        coef_hA5 <- coef(lme.h_massA5)
        plot(resid(lme.h_massA5))
        
        #Humus Set A 17 month incubation
        lme.h_massA17 <-lme(percent_loss ~ pine*shrubs, #
                          random = ~ 1|block,
                          data = humus_A17, na.action = na.omit)
        
        h_mass_anovaA17 <- anova(lme.h_massA17)
        coef_hA17 <- coef(lme.h_massA17)
        plot(resid(lme.h_massA17))
        
        
        
        
#set B
lme.h_massB <-lme(percent_loss ~ pine*shrubs*incubation, #
                  random = ~ 1|block/plot,
                  data = humusB, na.action = na.omit)

h_mass_anovaB <- anova(lme.h_massB)
coef_hB <- coef(lme.h_massB)
plot(resid(lme.h_massB))


                #subset to each set incubation period
                humus_B5 <- humusA[humusB$incubation == "5",]
                humus_B17 <- humusA[humusB$incubation == "17",]


              #Humus Set B 5 month incubation
              lme.h_massB5 <-lme(percent_loss ~ pine*shrubs, #
                  random = ~ 1|block,
                  data = humus_B5, na.action = na.omit)

              h_mass_anovaB5 <- anova(lme.h_massB5)
              coef_hB5 <- coef(lme.h_massB5)
              plot(resid(lme.h_massB5))

              #Humus Set B 17 month incubation
              lme.h_massB17 <-lme(percent_loss ~ pine*shrubs, #
                                random = ~ 1|block,
                                data = humus_B17, na.action = na.omit)
              
              h_mass_anovaB17 <- anova(lme.h_massB17)
              coef_hB17 <- coef(lme.h_massB17)
              plot(resid(lme.h_massB17))
              
              
# Create a blank workbook
humus_mass_loss_anova_output <- createWorkbook()

# Add some sheets to the workbook
addWorksheet(humus_mass_loss_anova_output, "Full Model")
addWorksheet(humus_mass_loss_anova_output, "Set A")
addWorksheet(humus_mass_loss_anova_output, "Set B")
addWorksheet(humus_mass_loss_anova_output, "Coef Full Model")
addWorksheet(humus_mass_loss_anova_output, "Coef Set A")
addWorksheet(humus_mass_loss_anova_output, "Coef Set B")
addWorksheet(humus_mass_loss_anova_output, "Set A 5 Month")
addWorksheet(humus_mass_loss_anova_output, "Set A 17 Month")
addWorksheet(humus_mass_loss_anova_output, "Set B 5 Month")
addWorksheet(humus_mass_loss_anova_output, "Set B 17 Month")
addWorksheet(humus_mass_loss_anova_output, "Coef Set A 5 Month")
addWorksheet(humus_mass_loss_anova_output, "Coef Set A 17 Month")
addWorksheet(humus_mass_loss_anova_output, "Coef Set B 5 Month")
addWorksheet(humus_mass_loss_anova_output, "Coef Set B 17 Month")

# Write the data to the sheets
writeData(humus_mass_loss_anova_output, sheet = "Full Model", x = h_mass_anova, rowNames = TRUE)
writeData(humus_mass_loss_anova_output, sheet = "Set A", x = h_mass_anovaA,rowNames = TRUE)
writeData(humus_mass_loss_anova_output, sheet = "Set B", x = h_mass_anovaB,rowNames = TRUE)
writeData(humus_mass_loss_anova_output, sheet = "Coef Full Model", x = coef_h, rowNames = TRUE)
writeData(humus_mass_loss_anova_output, sheet = "Coef Set A", x = coef_hA,rowNames = TRUE)
writeData(humus_mass_loss_anova_output, sheet = "Coef Set B", x = coef_hB,rowNames = TRUE)
writeData(humus_mass_loss_anova_output, sheet = "Set A 5 Month", x = h_mass_anovaA5, rowNames = TRUE)
writeData(humus_mass_loss_anova_output, sheet = "Set A 17 Month", x = h_mass_anovaA17,rowNames = TRUE)
writeData(humus_mass_loss_anova_output, sheet = "Set B 5 Month", x = h_mass_anovaB5,rowNames = TRUE)
writeData(humus_mass_loss_anova_output, sheet = "Set B 17 Month", x = h_mass_anovaB17,rowNames = TRUE)
writeData(humus_mass_loss_anova_output, sheet = "Coef Set A 5 Month", x = coef_hA5, rowNames = TRUE)
writeData(humus_mass_loss_anova_output, sheet = "Coef Set A 17 Month", x = coef_hA17,rowNames = TRUE)
writeData(humus_mass_loss_anova_output, sheet = "Coef Set B 5 Month", x = coef_hB5,rowNames = TRUE)
writeData(humus_mass_loss_anova_output, sheet = "Coef Set B 17 Month", x = coef_hB17,rowNames = TRUE)

# Export the file
saveWorkbook(humus_mass_loss_anova_output, "humus_mass_loss_anova_output.xlsx")



############ pine needle litter mass loss statistics #########

lme.l_mass <-lme(percent_loss ~ pine*shrubs*incubation*set, #
                 random = ~ 1|block/plot,
                 data = litter2, na.action = na.omit)

l_mass_anova <- anova(lme.l_mass)
coef_l <- coef(lme.l_mass)
plot(resid(lme.l_mass))

#subset to each set A and B
litterA <- litter2[litter2$set == "A",]
litterB <- litter2[litter2$set == "B",]

#set A
lme.l_massA <-lme(percent_loss ~ pine*shrubs*incubation, #
                  random = ~ 1|block/plot,
                  data = litterA, na.action = na.omit)

l_mass_anovaA <- anova(lme.l_massA)
coef_lA <- coef(lme.l_massA)
plot(resid(lme.l_massA))

          #subset to each incubation in A
          litterA5 <- litterA[litterA$incubation == "5",]
          litterA17 <- litterA[litterA$incubation == "17",]

          #set A 5 Month incubation
          lme.l_massA5 <-lme(percent_loss ~ pine*shrubs, #
                            random = ~ 1|block/plot,
                            data = litterA5, na.action = na.omit)
          
          l_mass_anovaA5 <- anova(lme.l_massA5)
          coef_lA5 <- coef(lme.l_massA5)
          plot(resid(lme.l_massA5))
          
          #set A 17 Month incubation
          lme.l_massA17 <-lme(percent_loss ~ pine*shrubs, #
                             random = ~ 1|block,
                             data = litterA17, na.action = na.omit)
          
          l_mass_anovaA17 <- anova(lme.l_massA17)
          coef_lA17 <- coef(lme.l_massA17)
          plot(resid(lme.l_massA17))

#set B
lme.l_massB <-lme(percent_loss ~ pine*shrubs*incubation, #
                  random = ~ 1|block/plot,
                  data = litterB, na.action = na.omit)

l_mass_anovaB <- anova(lme.l_massB)
coef_lB <- coef(lme.l_massB)
plot(resid(lme.l_massB))


#subset to each incubation in B
litterB5 <- litterB[litterB$incubation == "5",]
litterB17 <- litterB[litterB$incubation == "17",]

#set B 5 Month incubation
lme.l_massB5 <-lme(percent_loss ~ pine*shrubs, #
                   random = ~ 1|block,
                   data = litterB5, na.action = na.omit)

l_mass_anovaB5 <- anova(lme.l_massB5)
coef_lB5 <- coef(lme.l_massB5)
plot(resid(lme.l_massB5))

#set B 17 Month incubation
lme.l_massB17 <-lme(percent_loss ~ pine*shrubs, #
                    random = ~ 1|block,
                    data = litterB17, na.action = na.omit)

l_mass_anovaB17 <- anova(lme.l_massB17)
coef_lB17 <- coef(lme.l_massB17)
plot(resid(lme.l_massB17))

# Create a blank workbook
litter_mass_loss_anova_output <- createWorkbook()

# Add some sheets to the workbook
addWorksheet(litter_mass_loss_anova_output, "Full Model")
addWorksheet(litter_mass_loss_anova_output, "Set A")
addWorksheet(litter_mass_loss_anova_output, "Set B")
addWorksheet(litter_mass_loss_anova_output, "Coef Full Model")
addWorksheet(litter_mass_loss_anova_output, "Coef Set A")
addWorksheet(litter_mass_loss_anova_output, "Coef Set B")
addWorksheet(litter_mass_loss_anova_output, "Set A 5 Month")
addWorksheet(litter_mass_loss_anova_output, "Set A 17 Month")
addWorksheet(litter_mass_loss_anova_output, "Set B 5 Month")
addWorksheet(litter_mass_loss_anova_output, "Set B 17 Month")
addWorksheet(litter_mass_loss_anova_output, "Coef Set A 5 Month")
addWorksheet(litter_mass_loss_anova_output, "Coef Set A 17 Month")
addWorksheet(litter_mass_loss_anova_output, "Coef Set B 5 Month")
addWorksheet(litter_mass_loss_anova_output, "Coef Set B 17 Month")

# Write the data to the sheets
writeData(litter_mass_loss_anova_output, sheet = "Full Model", x = l_mass_anova, rowNames = TRUE)
writeData(litter_mass_loss_anova_output, sheet = "Set A", x = l_mass_anovaA,rowNames = TRUE)
writeData(litter_mass_loss_anova_output, sheet = "Set B", x = l_mass_anovaB,rowNames = TRUE)
writeData(litter_mass_loss_anova_output, sheet = "Coef Full Model", x = coef_l, rowNames = TRUE)
writeData(litter_mass_loss_anova_output, sheet = "Coef Set A", x = coef_lA,rowNames = TRUE)
writeData(litter_mass_loss_anova_output, sheet = "Coef Set B", x = coef_lB,rowNames = TRUE)
writeData(litter_mass_loss_anova_output, sheet = "Set A 5 Month", x = l_mass_anovaA5, rowNames = TRUE)
writeData(litter_mass_loss_anova_output, sheet = "Set A 17 Month", x = l_mass_anovaA17,rowNames = TRUE)
writeData(litter_mass_loss_anova_output, sheet = "Set B 5 Month", x = l_mass_anovaB5,rowNames = TRUE)
writeData(litter_mass_loss_anova_output, sheet = "Set B 17 Month", x = l_mass_anovaB17,rowNames = TRUE)
writeData(litter_mass_loss_anova_output, sheet = "Coef Set A 5 Month", x = coef_lA5, rowNames = TRUE)
writeData(litter_mass_loss_anova_output, sheet = "Coef Set A 17 Month", x = coef_lA17,rowNames = TRUE)
writeData(litter_mass_loss_anova_output, sheet = "Coef Set B 5 Month", x = coef_lB5,rowNames = TRUE)
writeData(litter_mass_loss_anova_output, sheet = "Coef Set B 17 Month", x = coef_lB17,rowNames = TRUE)

# Export the file
saveWorkbook(litter_mass_loss_anova_output, "litter_mass_loss_anova_output.xlsx")







