## Jädraås Mycorrhizal Removal Experiment
## 
##
#Author: Louis Mielke
#Date Created: 2021 APRIL 09
#Date Updated: 2022 AUGUST 18

rm(list=ls()) #clear the workspace (only if you want to)

##### load librarys ####
library(tidyverse)
library(ggplot2)
library(openxlsx)
library(ggpubr)
library(nlme)

set.seed(2022)

##### set wd #####
setwd("~/Projects/mycorrhizal removal/MeshBags/Jädraås_Gadgil_MycorrhizalTraits")

#### import metadata ####
meta <- read.table("metadata.txt", sep="\t", header=T, row.names=1)
class(meta);dim(meta);str(meta) #meta$start_date <- as.Date(meta$start_date)
meta$incubation <- as.factor(meta$incubation) #incubation months  5 or 17
meta$block <- as.factor(meta$block) 
meta$treatment <- factor(meta$treatment, levels = c('C', 'E', 'T','TE','DC'))

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



########## PLOTTING #######

####### pine needle litter mass remaining graphs ########

#by Set: summary for percent mass remaining for each treatment AND set
litter_summary_byset <- litter %>%
  group_by(set, incubation, treatment) %>% 
  drop_na() %>% 
  summarize_at(c("percent_mass_remaining"), funs(mean, sd, min, max,n()))

#standard error
litter_summary_byset$se <- litter_summary_byset$sd/sqrt(litter_summary_byset$n)

#Set A
litterA <- litter_summary_byset[litter_summary_byset$set == "A",]

lA <- ggplot(data = litterA)+
  aes(x = incubation, y = mean, fill = treatment)+
  geom_bar(stat = "identity", position=position_dodge(), color = "black")+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = .5, position = position_dodge(.9))+
  scale_y_continuous(limits = c(0,100),expand = (c(0,0)))+
  labs(x = "Incubation length (months)\n", y = "% Mass remaining",
       title = "\nSet 1\nPine needles")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 16, color = "black"),
        plot.title = element_text(hjust = .5,size = 16,color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        #legend.justification = c("left","top"),
        legend.position = c(.85,.85),
        #panel.border= element_rect(colour = "black", size=1, fill=NA),
        #panel.background = element_rect(colour = "black", size=1, fill=NA),
        legend.key.size = unit(.5,"cm"))+
  #annotate("text", x=3.25, y=16.5, label= "(a)", size = 8)+
  scale_fill_manual("Treatments",values = c("#D45226","#F59B25", "#7F8A53", "#9AC287", "white"),labels = c("+Pine\n+Shrubs","+Pine\n-Shrubs","-Pine\n+Shrubs","-Pine\n-Shrubs", " Disturbed control"))

#old palette c("#EE9A00", "#1874CD","#A2CD5A", "#FFD700","#CD5555")

#Set B
litterB <- litter_summary_byset[litter_summary_byset$set == "B",]
#Mielke et al 2022 colors c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00")
lB <- ggplot(data = litterB)+
  aes(x = incubation, y = mean, fill = treatment)+
  geom_bar(stat = "identity", position=position_dodge(), color = "black")+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = .5, position = position_dodge(.9))+
  scale_y_continuous(limits = c(0,100),expand = (c(0,0)))+
  labs(x = "Incubation length (months)\n", y = "% Mass remaining",
       title = "\nSet 2\nPine needles")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 16, color = "black"),
        plot.title = element_text(hjust = .5,size = 16,color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        #legend.justification = c("left","top"),
        legend.position = c(.85,.85),
        #panel.border= element_rect(colour = "black", size=1, fill=NA),
        #panel.background = element_rect(colour = "black", size=1, fill=NA),
        legend.key.size = unit(.5,"cm"))+
  #annotate("text", x=3.25, y=16.5, label= "(a)", size = 8)+
  scale_fill_manual("Treatments",values = c("#D45226","#F59B25", "#7F8A53", "#9AC287", "white"),labels = c("+Pine\n+Shrubs","+Pine\n-Shrubs","-Pine\n+Shrubs","-Pine\n-Shrubs", " Disturbed control"))



####### humus mass loss graphs ########

## decomposition curves for humus data combined
## Humus Round 1 2017-2018 blue and green bags deployed in june 2017 - Nov 2018
## Humus Round 2 2018-2019 blue and green bags deployed in june 2018 - Nov 2019


#BY SET: summary for percent mass remaining for each treatment by set
humus_summary_byset <- humus %>%
  group_by(set, incubation, treatment) %>% 
  drop_na() %>% 
  summarize_at(c("percent_mass_remaining"), funs(mean, sd, min, max,n()))

#standard error
humus_summary_byset$se <- humus_summary_byset$sd/sqrt(humus_summary_byset$n)

#Humus Set A
humusA <- humus_summary_byset[humus_summary_byset$set == "A",]

hA <- ggplot(data = humusA)+
  aes(x = incubation, y = mean, fill = treatment)+
  geom_bar(stat = "identity", position=position_dodge(), color = "black")+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = .5, position = position_dodge(.9))+
  scale_y_continuous(limits = c(0,100),expand = (c(0,0)))+
  labs(x = "Incubation duration (months)\n", y = "% Mass remaining",
       title = "\nSet 1\nHumus")+
  theme_classic()+
  theme_classic()+
  theme(axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 16, color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        plot.title = element_text(hjust = .5,size = 16,color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        #legend.justification = c("left","top"),
        legend.position = c(.85,.85),
        #panel.border= element_rect(colour = "black", size=1, fill=NA),
        #panel.background = element_rect(colour = "black", size=1, fill=NA),
        legend.key.size = unit(.5,"cm"))+
  #annotate("text", x=3.25, y=16.5, label= "(a)", size = 8)+
  scale_fill_manual("Treatments",values = c("#D45226","#F59B25", "#7F8A53", "#9AC287", "white"),labels = c("+Pine\n+Shrubs","+Pine\n-Shrubs","-Pine\n+Shrubs","-Pine\n-Shrubs", " Disturbed control"))
#Humus Set B
humusB <- humus_summary_byset[humus_summary_byset$set == "B",]

hB <- ggplot(data = humusB)+
  aes(x = incubation, y = mean, fill = treatment)+
  geom_bar(stat = "identity", position=position_dodge(), color = "black")+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = .5, position = position_dodge(.9))+
  scale_y_continuous(limits = c(0,100),expand = (c(0,0)))+
  labs(x = "Incubation duration (months)\n", y = "% Mass remaining",
       title = "\nSet 2\nHumus")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 16, color = "black"),
        plot.title = element_text(hjust = .5,size = 16,color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        #legend.justification = c("left","top"),
        legend.position = c(.85,.85),
        #panel.border= element_rect(colour = "black", size=1, fill=NA),
        #panel.background = element_rect(colour = "black", size=1, fill=NA),
        legend.key.size = unit(.5,"cm"))+
  #annotate("text", x=3.25, y=16.5, label= "(a)", size = 8)+
  scale_fill_manual("Treatments",values = c("#D45226","#F59B25", "#7F8A53", "#9AC287", "white"),labels = c("+Pine\n+Shrubs","+Pine\n-Shrubs","-Pine\n+Shrubs","-Pine\n-Shrubs", " Disturbed control"))


#### combining plots
#combining plots into display

#pdf("figures/Rplot_MassLoss_Humus&Litter_Barchart.pdf", width = 10, height = 10)
combined <- ggarrange(lA + rremove("xlab"), lB + rremove("xlab") + rremove("ylab"), hA, hB + rremove("ylab"), common.legend = TRUE, legend ="top", labels = c("\na)", "\nb)","\nc)","\nd)"), font.label = list(size = 18, color = "black"), hjust = -2.5, vjust = 1)
combined
#dev.off()


####### 


########### humus mass loss statistics #######

lme.h_mass <-lme(percent_mass_remaining ~ pine*shrubs*incubation*set, #
                 random = ~ 1|block/plot,
                 data = humus2, na.action = na.omit)

h_mass_anova <- anova(lme.h_mass)
h_mass_anova$coefficients <- lme.h_mass$coefficients[["fixed"]]
plot(resid(lme.h_mass))
shapiro.test(resid(lme.h_mass))
#subset to each set A and B
humusA <- humus2[humus2$set == "A",]
humusB <- humus2[humus2$set == "B",]

#set A
lme.h_massA <-lme(percent_mass_remaining ~ pine*shrubs*incubation, #
                 random = ~ 1|block/plot,
                 data = humusA, na.action = na.omit)

h_mass_anovaA <- anova(lme.h_massA)
h_mass_anovaA$coefficients <- lme.h_massA$coefficients[["fixed"]]
plot(resid(lme.h_massA))


        #subset to each set incubation period
        humus_A5 <- humusA[humusA$incubation == "5",]
        humus_A17 <- humusA[humusA$incubation == "17",]

        #Humus Set A 5 month incubation
        lme.h_massA5 <-lme(percent_mass_remaining ~ pine*shrubs, #
                          random = ~ 1|block,
                          data = humus_A5, na.action = na.omit)
        
        h_mass_anovaA5 <- anova(lme.h_massA5)
        h_mass_anovaA5$coefficients <- lme.h_massA5$coefficients[["fixed"]]
        plot(resid(lme.h_massA5))
        
        #Humus Set A 17 month incubation
        lme.h_massA17 <-lme(percent_mass_remaining ~ pine*shrubs, #
                          random = ~ 1|block,
                          data = humus_A17, na.action = na.omit)
        
        h_mass_anovaA17 <- anova(lme.h_massA17)
        h_mass_anovaA17$coefficients <- lme.h_massA17$coefficients[["fixed"]]
        plot(resid(lme.h_massA17))
        
        
        humusDC_5A <- subset(meta, substrate %in% c("Humus") & incubation %in% c("5") & set %in% c("A") & treatment %in% c("C","DC")) # dataset for statistics with DC      
        
        t.test(percent_mass_remaining ~ treatment, data = humusDC_5A)
        #Humus Set A 5 month incubation
        lme.h_massA5_DC <-lme(percent_mass_remaining ~ treatment, #
                           random = ~ 1|block,
                           data = humusDC_5A, na.action = na.omit)
        
        h_mass_anovaDC_A5 <- anova(lme.h_massA5_DC)
        coef_hA5 <- coef(lme.h_massA5)
        plot(resid(lme.h_massA5))
        
#set B
lme.h_massB <-lme(percent_mass_remaining ~ pine*shrubs*incubation, #
                  random = ~ 1|block/plot,
                  data = humusB, na.action = na.omit)

h_mass_anovaB <- anova(lme.h_massB)
h_mass_anovaB$coefficients <- lme.h_massB$coefficients[["fixed"]]
plot(resid(lme.h_massB))


                #subset to each set incubation period
                humus_B5 <- humusA[humusB$incubation == "5",]
                humus_B17 <- humusA[humusB$incubation == "17",]


              #Humus Set B 5 month incubation
              lme.h_massB5 <-lme(percent_mass_remaining ~ pine*shrubs, #
                  random = ~ 1|block,
                  data = humus_B5, na.action = na.omit)

              h_mass_anovaB5 <- anova(lme.h_massB5)
              coef_hB5 <- coef(lme.h_massB5)
              plot(resid(lme.h_massB5))

              #Humus Set B 17 month incubation
              lme.h_massB17 <-lme(percent_mass_remaining ~ pine*shrubs, #
                                random = ~ 1|block,
                                data = humus_B17, na.action = na.omit)
              
              h_mass_anovaB17 <- anova(lme.h_massB17)
              h_mass_anovaB17$coefficients <- lme.h_massB17$coefficients[["fixed"]]
              plot(resid(lme.h_massB17))
              
              
# Create a blank workbook
humus_mass_loss_anova_output <- createWorkbook()

# Add some sheets to the workbook
addWorksheet(humus_mass_loss_anova_output, "Full Model")
addWorksheet(humus_mass_loss_anova_output, "Set A")
addWorksheet(humus_mass_loss_anova_output, "Set B")
addWorksheet(humus_mass_loss_anova_output, "Set A 5 Month")
addWorksheet(humus_mass_loss_anova_output, "Set A 17 Month")
addWorksheet(humus_mass_loss_anova_output, "Set B 5 Month")
addWorksheet(humus_mass_loss_anova_output, "Set B 17 Month")

# Write the data to the sheets
writeData(humus_mass_loss_anova_output, sheet = "Full Model", x = h_mass_anova, rowNames = TRUE)
writeData(humus_mass_loss_anova_output, sheet = "Set A", x = h_mass_anovaA,rowNames = TRUE)
writeData(humus_mass_loss_anova_output, sheet = "Set B", x = h_mass_anovaB,rowNames = TRUE)
writeData(humus_mass_loss_anova_output, sheet = "Set A 5 Month", x = h_mass_anovaA5, rowNames = TRUE)
writeData(humus_mass_loss_anova_output, sheet = "Set A 17 Month", x = h_mass_anovaA17,rowNames = TRUE)
writeData(humus_mass_loss_anova_output, sheet = "Set B 5 Month", x = h_mass_anovaB5,rowNames = TRUE)
writeData(humus_mass_loss_anova_output, sheet = "Set B 17 Month", x = h_mass_anovaB17,rowNames = TRUE)


# Export the file
saveWorkbook(humus_mass_loss_anova_output, "humus_mass_remaining_anova_output.xlsx")



############ pine needle litter mass loss statistics #########

lme.l_mass <-lme(percent_mass_remaining ~ pine*shrubs*incubation*set, #
                 random = ~ 1|block/plot,
                 data = litter2, na.action = na.omit)

l_mass_anova <- anova(lme.l_mass)
l_mass_anova$coefficients <- lme.l_mass$coefficients[["fixed"]]
plot(resid(lme.l_mass))

#subset to each set A and B
litterA <- litter2[litter2$set == "A",]
litterB <- litter2[litter2$set == "B",]

#set A
lme.l_massA <-lme(percent_mass_remaining ~ pine*shrubs*incubation, #
                  random = ~ 1|block/plot,
                  data = litterA, na.action = na.omit)

l_mass_anovaA <- anova(lme.l_massA)
l_mass_anovaA$coefficients <- lme.l_massA$coefficients[["fixed"]]
plot(resid(lme.l_massA))

          #subset to each incubation in A
          litterA5 <- litterA[litterA$incubation == "5",]
          litterA17 <- litterA[litterA$incubation == "17",]

          #set A 5 Month incubation
          lme.l_massA5 <-lme(percent_mass_remaining ~ pine*shrubs, #
                            random = ~ 1|block/plot,
                            data = litterA5, na.action = na.omit)
          
          l_mass_anovaA5 <- anova(lme.l_massA5)
          l_mass_anovaA5$coefficients <- lme.l_massA5$coefficients[["fixed"]]
          plot(resid(lme.l_massA5))
          
          #set A 17 Month incubation
          lme.l_massA17 <-lme(percent_mass_remaining ~ pine*shrubs, #
                             random = ~ 1|block,
                             data = litterA17, na.action = na.omit)
          
          l_mass_anovaA17 <- anova(lme.l_massA17)
          l_mass_anovaA17$coefficients <- lme.l_massA17$coefficients[["fixed"]]
          plot(resid(lme.l_massA17))

#set B
lme.l_massB <-lme(percent_mass_remaining ~ pine*shrubs*incubation, #
                  random = ~ 1|block/plot,
                  data = litterB, na.action = na.omit)

l_mass_anovaB <- anova(lme.l_massB)
l_mass_anovaB$coefficients <- lme.l_massB$coefficients[["fixed"]]
plot(resid(lme.l_massB))


#subset to each incubation in B
litterB5 <- litterB[litterB$incubation == "5",]
litterB17 <- litterB[litterB$incubation == "17",]

#set B 5 Month incubation
lme.l_massB5 <-lme(percent_mass_remaining ~ pine*shrubs, #
                   random = ~ 1|block,
                   data = litterB5, na.action = na.omit)

l_mass_anovaB5 <- anova(lme.l_massB5)
l_mass_anovaB5$coefficients <- lme.l_massB5$coefficients[["fixed"]]
plot(resid(lme.l_massB5))

#set B 17 Month incubation
lme.l_massB17 <-lme(percent_mass_remaining ~ pine*shrubs, #
                    random = ~ 1|block,
                    data = litterB17, na.action = na.omit)

l_mass_anovaB17 <- anova(lme.l_massB17)
l_mass_anovaB17$coefficients <- lme.l_massB17$coefficients[["fixed"]]
plot(resid(lme.l_massB17))

# Create a blank workbook
litter_mass_loss_anova_output <- createWorkbook()

# Add some sheets to the workbook
addWorksheet(litter_mass_loss_anova_output, "Full Model")
addWorksheet(litter_mass_loss_anova_output, "Set A")
addWorksheet(litter_mass_loss_anova_output, "Set B")
addWorksheet(litter_mass_loss_anova_output, "Set A 5 Month")
addWorksheet(litter_mass_loss_anova_output, "Set A 17 Month")
addWorksheet(litter_mass_loss_anova_output, "Set B 5 Month")
addWorksheet(litter_mass_loss_anova_output, "Set B 17 Month")

# Write the data to the sheets
writeData(litter_mass_loss_anova_output, sheet = "Full Model", x = l_mass_anova, rowNames = TRUE)
writeData(litter_mass_loss_anova_output, sheet = "Set A", x = l_mass_anovaA,rowNames = TRUE)
writeData(litter_mass_loss_anova_output, sheet = "Set B", x = l_mass_anovaB,rowNames = TRUE)
writeData(litter_mass_loss_anova_output, sheet = "Set A 5 Month", x = l_mass_anovaA5, rowNames = TRUE)
writeData(litter_mass_loss_anova_output, sheet = "Set A 17 Month", x = l_mass_anovaA17,rowNames = TRUE)
writeData(litter_mass_loss_anova_output, sheet = "Set B 5 Month", x = l_mass_anovaB5,rowNames = TRUE)
writeData(litter_mass_loss_anova_output, sheet = "Set B 17 Month", x = l_mass_anovaB17,rowNames = TRUE)

# Export the file
saveWorkbook(litter_mass_loss_anova_output, "litter_mass_remaining_anova_output.xlsx")

