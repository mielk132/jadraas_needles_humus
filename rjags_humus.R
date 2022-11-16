### Target: Ecology Letters
### Title:
###
### Authors: L A Mielke 1, B J Lindahl 2, J Klein 3, R Finlay 1, A Ekblad 4, K E Clemmensen 1
###
### Affiliations
### 1 Dept Forest Mycology Plant Pathology - SLU Uppsala
### 2 Dept Soil & Environment - SLU Uppsala
### 3 Artdatabanken - SLU Uppsala
### 4 Örebro University

### Corresponding author
### Louis Mielke
### louis.mielke@slu.se
###
###
###

###################################################################################

#this is the R side of the analysis. The model is written in a separate file

rm(list=ls()) #clear the workspace (only if you want to)
getwd()
#dir.create("figures")
#import libraries, data etc
#-----------------------------------------------------------#

#need these pckages installed
library(boot)
library(rjags)
library(coda)
library(scales)
library(fastDummies)
library(tidyverse)
library(corrplot)
library("runjags")

#file.choose()
#must set your working directory to wherever your JAGS models will be
setwd("~/Projects/mycorrhizal removal/MeshBags/Jädraås_Gadgil_MycorrhizalTraits")

decomp_type <- "humus"

#set <- "A"
set <- "B"

time <- "17"


#read in your data file(s)
df <- read.csv("humus_meta.csv",header=T, row.names=1) #bring in the data

#check data
names(df)


# cut & paste names in here - makes it easier when setting up your JAGS data object 
str(df) #check how R defines your variables
summary(df) #check for NA's in your explanatory variables and data range
df <- na.omit(df)

#separate by substrate, set, and incubation time
humusA17  <- subset(df, set =="A" & incubation == "17")
humusA5  <- subset(df, set =="A" & incubation == "5")
humusB17  <- subset(df, set =="B" & incubation == "17")
humusB5  <- subset(df, set =="B" & incubation == "5")
humusA  <- subset(df, set =="A")
humusB  <- subset(df, set =="B")

# 
# humusB5 <- arrange(humusB5, treatment, block)   #missing E4
# humusB17 <- arrange(humusB17, treatment, block)  
# 
# 
# #fill E4 with mean % mass remaining for that treatment
# humusB5_per_mass_remaining <- humusB5$percent_mass_remaining #make a vector
# e_avg <- mean(humusB5$percent_mass_remaining[humusB5$treatment == "E"]) 
# humusB5_per_mass_remaining <- append(humusB5_per_mass_remaining, e_avg, after = 19)
# humusB17$humusB5_per_mass_remaining <- humusB5_per_mass_remaining 

humusB17  <- filter(humusB17, treatment != "T" & treatment != "TE" & treatment != "DC")

#set up conditions for JAGS model
#-----------------------------------------------------------#

##### 1. create the data structure ####
data = NULL       #clear any old data lists that might confuse things
data = list(
  #t = rep(524, length(humusB17$start_date)), #number of growing season days 522 for 2017-2018, 524 for 2018-2019
  M0 = humusB17$original_substrate_mass_g, 
  Mt = humusB17$final_substrate_g,
  M = humusB17$percent_mass_remaining,
  nobs = nrow(humusB17),
  #plot = as.numeric(as.factor(humusB17$plot)),
  block = as.numeric(humusB17$block),
  biomass = scale(humusB17$copies_fungi),
  ecto = scale(humusB17$copies_ecto),
  ericoid = scale(humusB17$copies_ericoid),
  sap_wr = scale(humusB17$copies_sap_whiterot),
  ecto_wr = scale(humusB17$copies_ecto_whiterot),
  whiterot = scale(humusB17$copies_ecto_whiterot + humusB17$copies_sap_whiterot))

cor(data$sap_wr, data$ecto)
cor(data$ecto_wr, data$ericoid)
cor(data$ericoid, data$sap_wr)

#M <- cor(do.call(cbind.data.frame, data))
#corrplot(M, method = 'number')

### add new 
data$bm_pred <- seq(min(data$biomass), max(data$biomass), 0.1)
data$bm_ericoid_pred <- seq(min(data$ericoid), max(data$ericoid), 0.1)
data$bm_ecto_pred <- seq(min(data$ecto), max(data$ecto), 0.1)
data$bm_sap_wr_pred <- seq(min(data$sap_wr), max(data$sap_wr), 0.1)
data$bm_ecto_wr_pred <- seq(min(data$ecto_wr), max(data$ecto_wr), 0.1)
data$bm_whiterot_pred <- seq(min(data$whiterot), max(data$whiterot), 0.1)

data #look at your data file and check it looks ok
str(data)

##### What to expect: ####

plot(data$ericoid,data$M)
abline(lm(data$M ~ data$ericoid))
summary(lm(data$M ~ data$ericoid))


plot(data$ecto,data$M)
abline(lm(data$M ~ data$ecto))
summary(lm(data$M ~ data$ecto))

plot(data$sap_wr,data$M)
abline(lm(data$M ~ data$sap_wr))
summary(lm(data$M ~ data$sap_wr))

plot(data$ecto_wr,data$M)
abline(lm(data$M ~ data$ecto_wr))
summary(lm(data$M ~ data$ecto_wr))


#set up conditions for JAGS model
#-----------------------------------------------------------#

#2. create the initial values (if not provided, JAGS will do this for you)
#NOTE that it is a list within a list
#this allows you to have multiple lists that specify the initial values if you run more than 1 chain
#as with the data object you must use '=' inside a list, you cannot use '<-'

inits=list(
  list(
    sigmaM = 1,
    alpha = 0.1, ## Change if muk is logarithmized
    b_bm = 0,
    b_ecto = .1,
    b_ericoid = .1,
    b_whiterot = .1,
    b_sap_wr = .1, #white rot sap 
    b_ecto_wr = .1, #white rot ecto
    sigmablock = 1,
    k=rep(0.1,data$nobs)
  ),
  list(
    sigmaM = 0.1,
    alpha = 0.001, ## Change if muk is logarithmized
    b_bm = -0.5,
    b_ecto = -0.5,
    b_ericoid = -.5,
    b_whiterot = -.5,
    b_sap_wr = -.5, #sap white rot fungi
    b_ecto_wr = -.5,  # white rot ecto fungi
    sigmablock = 0.1,
    k=rep(0.1,data$nobs)
  ),
  list(
    sigmaM = 0.5,
    alpha = 0.01, ## Change if muk is logarithmized
    b_bm = -0.05,
    b_ecto = -0.05,
    b_ericoid = -.05,
    b_whiterot = -.05,
    b_sap_wr = -.05, #sap white rot fungi
    b_ecto_wr = -.05,  # white rot ecto fungi
    sigmablock = 0.01,
    k=rep(0.1,data$nobs)
  ))

#specify the model and compile the model
#----------------------------------------------------------------------


#3. name the JAGS model file

#model = "JAGS_humus_biomass_dnorm_model.R" #total fungal biomass
#model = "JAGS_humus_mycorrhizal_biomass_model.R" # ecto + ericoid biomass
model = "JAGS_humus_mycorrhizal_whiterot_biomass_model.R"  # ericoid + ecto_wr + sap_wr biomass

#4. compile the model
jm = jags.model(model,
                data=data, 
                n.adapt=20000, 
                inits=inits, 
                n.chains=3) 


#note you need to adjust this code if you don't provide JAGS with inits or if you increase the number of chains you run

#5. burn-in the model

burn.in=100000

update(jm, n.iter=burn.in) #this runs the number burn-in values so the model is ready for sampling

#6. generate samples from the posterior distribution CODA

samples=15000
n.thin=5

zc = coda.samples(jm,variable.names=c("alpha","sigmablock","R2","b_bm","b_ecto","b_ericoid","b_ecto_wr","b_sap_wr","b_whiterot"), n.iter=samples, thin=n.thin)

#output coda data and you can summaries these or visually inspect the chains
#if you don't know what this means read the JAGS manual or JAGS primer


#model selection, extract dic, but not reporting
ms = jags.samples(jm,
                  variable.names=c("M.new"), 
                  n.iter=samples, 
                  thin=n.thin)
PPL = sum((data$M-summary(ms$M.new,mean)$stat)^2) + sum((summary(ms$M.new,sd)$stat)^2)
PPL

#waic

zc_combined <- combine.mcmc(mcmc.objects = zc)
ecdf <-  apply(zc_combined, 2, function(x) 1-ecdf(x)(0)) #probability of positive effect

#model 1
ecdf[["b_bm"]]

#model 2
ecdf[["b_ecto"]]
ecdf[["b_ericoid"]]

#model 3
ecdf[["b_ericoid"]]
ecdf[["b_whiterot"]]

#
ecdf[["b_sap_wr"]]
ecdf[["b_ecto_wr"]]
ecdf[["b_ericoid"]]

summary(zc) #will show the mean estimate and SE and 95% CIs


# pdf("figures/MCMC_chains_humusB17_mass_remaining_linear_biomass_C_E.pdf")
# pdf("figures/MCMC_chains_humusB17_mass_remaining_linear_ericoid_ecto_biomass_C_E.pdf")
pdf("figures/MCMC_chains_humusB17_mass_remaining_linear_ericoid_ecto_sap_whiterot_biomass_C_E.pdf")
plot(zc); gelman.plot(zc) #look at the chains to see stability and mixing


zj_val <- jags.samples(jm,
                       variable.names = c("mean.y", "mean.y.new","pvalue.mean",
                                          "cv.y", "cv.y.new", "pvalue.cv",
                                          "fit", "fit.new", "pvalue"),
                       n.iter = samples,
                       thin = n.thin)

## Fit of mean:
plot(zj_val$mean.y,
     zj_val$mean.y.new,
     xlab = "mean real",
     ylab = "mean simulated",
     cex = .05)
abline(0, 1)
p <- summary(zj_val$pvalue.mean, mean)
text(x = min(zj_val$mean.y)+min(zj_val$mean.y)/10, y = min(zj_val$mean.y.new)+min(zj_val$mean.y.new)/10, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)

## Fit of variance:
plot(zj_val$cv.y,
     zj_val$cv.y.new,
     xlab = "cv real",
     ylab = "cv simulated",
     cex = .05)
abline(0,1)
p <- summary(zj_val$pvalue.cv, mean)
text(x = min(zj_val$cv.y)+min(zj_val$cv.y)/10, y = min(zj_val$cv.y.new)+min(zj_val$cv.y.new)/10, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)

## Overall fit:
plot(zj_val$fit,
     zj_val$fit.new,
     xlab = "ssq real",
     ylab = "ssq simulated",
     cex = .05)
abline(0,1)
p <- summary(zj_val$pvalue, mean)
text(x = mean(zj_val$fit), y=max(zj_val$fit.new)-max(zj_val$fit.new)/10, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)

dev.off()
#7. generate samples from the posterior distribution JAGS

#output jags.object data. For some things easier to manipulate and output from here
#if you don't know what this means, read the JAGS manual

zj = jags.samples(jm, 
                  variable.names=c("mass_bm_pred","mass_ericoid_pred","mass_ecto_pred","mass_sap_wr_pred","mass_ecto_wr_pred"), 
                  n.iter=samples, 
                  thin=n.thin)

####  BIOMASS PREDICTION #####
#plotting prediction & 95%CIs using polygon
pdf(paste0("figures/Prediction",set,time,"_model_biomass_trmts_C_E.pdf"))

pred<-summary(zj$mass_bm_pred, quantile, c(.025,.5,.975))$stat #get your prediction
x=data$bm_pred #set your x-axis relative to your x.pred prediction
y=pred
plot(x,y[2,], col="blue", xlab="total fungal biomass [scaled]", ylab="Mass remaining (% of initial mass)", cex=1.4, typ="l", tck=0.03, bty="l", ylim = c(min(y[1,]), max(y[3,]))) #plot the median prediction
polygon(c(x,rev(x)), c(y[1,], rev(y[3,])), col=alpha("blue", alpha=0.5)) #add the 95% CIs
lines(x,y[1,], lty="dashed", col="blue") #add edges to the polygon
lines(x,y[3,], lty="dashed", col="blue")
dev.off()

#### ECTO BIOMASS PREDICTION #####
pdf(paste0("figures/Prediction",set,time,"_mycorrhizal_model_ecto_biomass_trmts_C_E.pdf"))

#plotting prediction & 95%CIs using polygon
pred<-summary(zj$mass_ecto_pred, quantile, c(.025,.5,.975))$stat #get your prediction
x=data$bm_ecto_pred    #set your x-axis relative to your x.pred prediction
y=pred
plot(x,y[2,], col="blue", xlab="ectomycorrhizal biomass [scaled]", ylab="Mass remaining (% of initial mass)", cex=1.4, typ="l", tck=0.03, bty="l", ylim = c(min(y[1,]), max(y[3,]))) #plot the median prediction
polygon(c(x,rev(x)), c(y[1,], rev(y[3,])), col=alpha("grey", alpha=0.5)) #add the 95% CIs
lines(x,y[1,], lty="dashed", col="blue") #add edges to the polygon
lines(x,y[3,], lty="dashed", col="blue")
dev.off()

#### ERICOID BIOMASS PREDICTION #####
pdf(paste0("figures/Prediction",set,time,"_mycorrhizal_model_ericoid_biomass_trmts_C_E_withwhiterot.pdf"))

#plotting prediction & 95%CIs using polygon
pred<-summary(zj$mass_ericoid_pred, quantile, c(.05,.5,.95))$stat #get your prediction
pred2<-summary(zj$mass_ericoid_pred, quantile, c(.30,.5,.70))$stat #get your prediction
x=data$bm_ericoid_pred    #set your x-axis relative to your x.pred prediction
y=pred
y2 = pred2
plot(x,y[2,], col="blue", xlab="ericoid biomass [scaled]", ylab="Mass remaining (% of initial mass)", cex=1.4, typ="l", tck=0.03, bty="l", ylim = c(min(y[1,]), max(y[3,]))) #plot the median prediction
polygon(c(x,rev(x)), c(y[1,], rev(y[3,])), col=alpha("grey", alpha=0.5)) #add the 95% CIs
polygon(c(x,rev(x)), c(y2[1,], rev(y2[3,])), col=alpha("mediumseagreen", alpha=1)) #add the 70% CIs
#lines(x,y[1,], lty="dashed", col="blue") #add edges to the polygon
#lines(x,y[3,], lty="dashed", col="blue")
dev.off()

#### SAP WHITE ROT BIOMASS PREDICTION #####
pdf(paste0("figures/Prediction",set,time,"_mycorrhizal_model_sap_whiterot_biomass_trmts_C_E.pdf"))

#plotting prediction & 95%CIs using polygon
pred<-summary(zj$mass_sap_wr_pred, quantile, c(.30,.5,.70))$stat #get your prediction
x=data$bm_sap_wr_pred    #set your x-axis relative to your x.pred prediction
y=pred
plot(x,y[2,], col="blue", xlab="white rot fungal biomass [scaled]", ylab="Mass remaining (% of initial mass)", cex=1.4, typ="l", tck=0.03, bty="l", ylim = c(min(y[1,]), max(y[3,]))) #plot the median prediction
polygon(c(x,rev(x)), c(y[1,], rev(y[3,])), col=alpha("grey", alpha=0.5)) #add the 95% CIs
lines(x,y[1,], lty="dashed", col="blue") #add edges to the polygon
lines(x,y[3,], lty="dashed", col="blue")
dev.off()

#### ECTO WHITE ROT BIOMASS PREDICTION #####
pdf(paste0("figures/Prediction",set,time,"_mycorrhizal_model_ecto_whiterot_biomass_trmts_C_E.pdf"))

#plotting prediction & 95%CIs using polygon
pred<-summary(zj$mass_ecto_wr_pred, quantile, c(.025,.5,.975))$stat #get your prediction
x=data$bm_ecto_wr_pred    #set your x-axis relative to your x.pred prediction
y=pred
plot(x,y[2,], col="blue", xlab="ecto white rot biomass [scaled]", ylab="Mass remaining (% of initial mass)", cex=1.4, typ="l", tck=0.03, bty="l", ylim = c(min(y[1,]), max(y[3,]))) #plot the median prediction
polygon(c(x,rev(x)), c(y[1,], rev(y[3,])), col=alpha("blue", alpha=0.5)) #add the 95% CIs
lines(x,y[1,], lty="dashed", col="blue") #add edges to the polygon
lines(x,y[3,], lty="dashed", col="blue")
dev.off()
