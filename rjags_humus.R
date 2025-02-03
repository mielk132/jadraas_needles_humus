### Target: New Phyt
### Title: Mycorrhizal fungal guild interactions slow humus decomposition in a boreal forest
###
### Corresponding author
### Louis Mielke
### louis.mielke@stir.ac.uk
###
###
###

###################################################################################

#this is the R side of the analysis for the biomass and guild models in humus. The JAGS models are written in a separate file and are called individually

rm(list=ls()) #clear the workspace 
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
#must set your working directory to wherever the JAGS models will be
setwd("~/Projects/mycorrhizal removal/MeshBags/Jädraås_Gadgil_MycorrhizalTraits")

decomp_type <- "humus"

#Analysis only includes 17 month set B, since there was more mycelial growth (as evident by copies of ITS - Fig 1)


set <- "B" #second set of incubations (see Supp Fig 1)
#set <- "A"

time <- "17" #17 month incubation duration (see Supp Fig 1)

#read in data file(s)
df <- read.csv("humus_meta.csv",header=T, row.names=1) #bring in the data

#check data
names(df)

# cut & paste names in here - makes it easier when setting up your JAGS data object 
str(df) #check how R defines your variables
summary(df) #check for NA's in your explanatory variables and data range
df <- na.omit(df)

#separate by substrate, set, and incubation time
#humusA17  <- subset(df, set =="A" & incubation == "17")
#humusA5  <- subset(df, set =="A" & incubation == "5")
humusB17  <- subset(df, set =="B" & incubation == "17")
#humusB5  <- subset(df, set =="B" & incubation == "5")
#humusA  <- subset(df, set =="A")
#humusB  <- subset(df, set =="B")

#set up conditions for JAGS model
#-----------------------------------------------------------#

##### 1. create the data structure ####
data = NULL       #clear any old data lists that might confuse things
data = list(
  #t = rep(524, length(humusB17$start_date)), #number of growing season days 522 for 2017-2018, 524 for 2018-2019
  M0 = humusB17$original_substrate_mass_g, 
  Mt = humusB17$final_substrate_g,
  M = humusB17$percent_loss,
  nobs = nrow(humusB17),
  #plot = as.numeric(as.factor(humusB17$plot)),
  block = as.numeric(humusB17$block),

  #NOTE - only one model (Eq 2A or 2B) can be applied at a time
  
  ###### Biomass Model (Mielke et al., Equation 2A)
  #biomass = scale(humusB17$copies_fungi),
  
  ###### Guild Model (Mielke et al., Equation 2B)
  ecto = scale(humusB17$copies_ecto), #all ectos
  sap = scale(humusB17$copies_sap_other + humusB17$copies_sap_whiterot), #all filamentous saps
  molds_yeasts = scale(humusB17$copies_moulds_yeasts), #all moulds and yeasts
  ericoid = scale(humusB17$copies_ericoid)
  )

#M <- cor(do.call(cbind.data.frame, data))
#corrplot(M, method = 'number')

### MODEL 2A add new 
#data$bm_pred <- seq(min(data$biomass), max(data$biomass), 0.1)

#MODEL 2B
data$bm_ecto_pred <- seq(min(data$ecto), max(data$ecto), 0.1)
data$bm_sap_pred <- seq(min(data$sap), max(data$sap), 0.1)
data$bm_molds_yeasts_pred <- seq(min(data$molds_yeasts), max(data$molds_yeasts), 0.1)
data$bm_ericoid_pred <- seq(min(data$ericoid), max(data$ericoid), 0.1)

data #look at your data file and check it looks ok
str(data)

##### What to expect: ####

#MODEL 2B
plot(data$ericoid,data$M)
abline(lm(data$M ~ data$ericoid))
summary(lm(data$M ~ data$ericoid))

plot(data$ecto,data$M)
abline(lm(data$M ~ data$ecto))
summary(lm(data$M ~ data$ecto))

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
    b_molds_yeasts = .1,
    b_sap = 0.1,
    sigmablock = 1,
    k=rep(0.1,data$nobs)
  ),
  list(
    sigmaM = 0.1,
    alpha = 0.001, ## Change if muk is logarithmized
    b_bm = -0.5,
    b_ericoid = -0.5,
    b_ecto = -0.5,
    b_sap = -0.5,
    b_molds_yeasts = -0.5,
    sigmablock = 0.1,
    k=rep(0.1,data$nobs)
  ),
  list(
    sigmaM = 0.5,
    alpha = 0.01, ## Change if muk is logarithmized
    b_bm = -0.05,
    b_ecto = -0.05,
    b_ericoid = -0.05,
    b_sap = -0.05,
    b_molds_yeasts = -0.05,
    sigmablock = 0.01,
    k=rep(0.1,data$nobs)
  ))

#specify the model and compile the model
#----------------------------------------------------------------------


#3. name the JAGS model file

#model = "JAGS_humus_biomass_dnorm_model.R" #total fungal copies as a proxy for total biomass or 'fungal activity'
model = "JAGS_humus_guild_model.R" # ecto + sap + molds + ericoid (copies)

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

zc = coda.samples(jm,variable.names=c("alpha","sigmablock","R2","b_bm","b_ecto","b_sap","b_ericoid", "b_molds_yeasts"), n.iter=samples, thin=n.thin)

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
ecdf[["b_sap"]]
ecdf[["b_molds_yeasts"]]
ecdf[["b_ericoid"]]

summary(zc) #will show the mean estimate and SE and 95% CIs

#pdf("figures/MCMC_chains_humusB17_guild_percent_mass_loss_linear.pdf")

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
                  variable.names=c("mass_bm_pred","mass_ecto_pred", "mass_molds_yeasts_pred", "mass_sap_pred","mass_ericoid_pred"), 
                  n.iter=samples, 
                  thin=n.thin)

c("#88CCAA", "#44AA77","#777711","#AAAA44","gray70", "gray50", "#AA7744", "#774411","pink")

####  TOTAL FUNGI PREDICTION #####
#plotting prediction & 95%CIs using polygon
#pdf(paste0("figures/Prediction",set,time,"_model.pdf"))

pred<-summary(zj$mass_bm_pred, quantile, c(.025,.5,.975))$stat #get your prediction
x=data$bm_pred #set your x-axis relative to your x.pred prediction
y=pred
plot(x,y[2,], col="blue", xlab="total fungal ITS copies [scaled]", ylab="Mass Loss (% of initial mass)", cex=1.4, typ="l", tck=0.03, bty="l", ylim = c(min(y[1,]), max(y[3,]))) #plot the median prediction
polygon(c(x,rev(x)), c(y[1,], rev(y[3,])), col=alpha("blue", alpha=0.5)) #add the 95% CIs
lines(x,y[1,], lty="dashed", col="blue") #add edges to the polygon
lines(x,y[3,], lty="dashed", col="blue")
dev.off()

#### ECTOMYCORRHIZAL PREDICTION #####
pdf(paste0("figures/Prediction",set,time,"_guild_model_ectomycorrhizal.pdf"))

#plotting prediction & 95%CIs using polygon
pred<-summary(zj$mass_ecto_pred, quantile, c(.05,.5,.95))$stat #get your prediction
pred2<-summary(zj$mass_ecto_pred, quantile, c(.30,.5,.70))$stat #get your prediction
x=data$bm_ecto_pred    #set your x-axis relative to your x.pred prediction
y=pred
y2 = pred2
plot(x,y[2,], col="#44AA77", xlab="ectomycorrhizal copies [scaled]", ylab="Mass loss (% of initial mass)", cex=1.4, typ="l", tck=0.03, bty="l", ylim = c(min(y[1,]), max(y[3,]))) #plot the median prediction
polygon(c(x,rev(x)), c(y[1,], rev(y[3,])), col=alpha("grey", alpha=0.2)) #add the 95% CIs
polygon(c(x,rev(x)), c(y2[1,], rev(y2[3,])), col=alpha("#44AA77", alpha=1)) #add the 70% CIs
lines(x,y[1,], lty="dashed", col="#44AA77") #add edges to the polygon
lines(x,y[3,], lty="dashed", col="#44AA77")
dev.off()

#### ERICOID PREDICTION #####
pdf(paste0("figures/Prediction",set,time,"_guild_model_ericoid.pdf"))

#plotting prediction & 95%CIs using polygon
pred<-summary(zj$mass_ericoid_pred, quantile, c(.05,.5,.95))$stat #get your prediction
pred2<-summary(zj$mass_ericoid_pred, quantile, c(.30,.5,.70))$stat #get your prediction
x=data$bm_ericoid_pred    #set your x-axis relative to your x.pred prediction
y=pred
y2 = pred2
plot(x,y[2,], col="blue", xlab="ericoid copies [scaled]", ylab="Mass loss (% of initial mass)", cex=1.4, typ="l", tck=0.03, bty="l", ylim = c(min(y[1,]), max(y[3,]))) #plot the median prediction
polygon(c(x,rev(x)), c(y[1,], rev(y[3,])), col=alpha("grey", alpha=0.5)) #add the 95% CIs
polygon(c(x,rev(x)), c(y2[1,], rev(y2[3,])), col=alpha("#777711", alpha=1)) #add the 70% CIs
#lines(x,y[1,], lty="dashed", col="blue") #add edges to the polygon
#lines(x,y[3,], lty="dashed", col="blue")
dev.off()

#### SAPROTROPHIC PREDICTION #####

pdf(paste0("figures/Prediction",set,time,"_guild_model_saprotrophs.pdf"))

#plotting prediction & 95%CIs using polygon
pred<-summary(zj$mass_sap_pred, quantile, c(.05,.5,.95))$stat #get your prediction
pred2<-summary(zj$mass_sap_pred, quantile, c(.30,.5,.70))$stat #get your prediction
x=data$bm_sap_pred    #set your x-axis relative to your x.pred prediction
y=pred
y2 = pred2
plot(x,y[2,], col="#774411", xlab="saprotrophic ITS2 copies [scaled]", ylab="Mass Loss (% of initial mass)", cex=1.4, typ="l", tck=0.03, bty="l", ylim = c(min(y[1,]), max(y[3,]))) #plot the median prediction
polygon(c(x,rev(x)), c(y[1,], rev(y[3,])), col=alpha("grey", alpha=0.2)) #add the 95% CIs
polygon(c(x,rev(x)), c(y2[1,], rev(y2[3,])), col=alpha("#774411", alpha=1)) #add the 70% CIs
#lines(x,y[1,], lty="dashed", col="blue") #add edges to the polygon
#lines(x,y[3,], lty="dashed", col="blue")
dev.off()

#### MOLDS YEASTS PREDICTION #####

pdf(paste0("figures/Prediction",set,time,"_guild_model_molds_yeasts.pdf"))

#plotting prediction & 95%CIs using polygon
pred<-summary(zj$mass_molds_yeasts_pred, quantile, c(.05,.5,.95))$stat #get your prediction
pred2<-summary(zj$mass_molds_yeasts_pred, quantile, c(.30,.5,.70))$stat #get your prediction
x=data$bm_molds_yeasts_pred    #set your x-axis relative to your x.pred prediction
y=pred
y2 = pred2
plot(x,y[2,], col="blue", xlab="Molds and Yeasts ITS2 copies [scaled]", ylab="Mass Loss (% of initial mass)", cex=1.4, typ="l", tck=0.03, bty="l", ylim = c(min(y[1,]), max(y[3,]))) #plot the median prediction
polygon(c(x,rev(x)), c(y[1,], rev(y[3,])), col=alpha("grey", alpha=0.5)) #add the 95% CIs
polygon(c(x,rev(x)), c(y2[1,], rev(y2[3,])), col=alpha("dark grey", alpha=1)) #add the 70% CIs
#lines(x,y[1,], lty="dashed", col="blue") #add edges to the polygon
#lines(x,y[3,], lty="dashed", col="blue")
dev.off()
