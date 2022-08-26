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
library("runjags")

#file.choose()
#must set your working directory to whereever your JAGS models will be
setwd("~/Projects/mycorrhizal removal/MeshBags/Jädraås_Gadgil_MycorrhizalTraits")

decomp_type <- "litter"

set <- "A"
set <- "B"

time <- "17"


#read in your data file(s)
df <- read.csv("litter_meta.csv",header=T, row.names=1) #bring in the data

#check data
names(df)


# cut & paste names in here - makes it easier when setting up your JAGS data object 
str(df) #check how R defines your variables
summary(df) #check for NA's in your explanatory variables and data range
df <- na.omit(df)

#separate by substrate, set, and incubation time
litterA17  <- subset(df, set =="A" & incubation == "17")
litterB17  <- subset(df, set =="B" & incubation == "17")
litterA  <- subset(df, set =="A")
litterB  <- subset(df, set =="B")

#set up conditions for JAGS model
#-----------------------------------------------------------#

##### 1. create the data structure ####
data = NULL       #clear any old data lists that might confuse things
data = list(
	t = rep(524, length(litterB17$start_date)), #number of growing season days 522 for 2017-2018, 524 for 2018-2019
	M0 = litterB17$original_substrate_mass_g, 
	Mt = litterB17$final_substrate_g,
	nobs = nrow(litterB17),
	#plot = as.numeric(as.factor(litterB17$plot)),
	block = as.numeric(litterB17$block),
	biomass = scale(litterB17$copies_fungi),
	sap_whiterot = scale(litterB17$copies_sap_whiterot),
	sap_non_whiterot = scale(litterB17$copies_sap_other)
)

data$bm_pred <- seq(min(data$biomass), max(data$biomass), 0.1)
data$bm_sap_whiterot_pred <- seq(min(data$sap_whiterot), max(data$sap_whiterot), 0.1)
data$bm_sap_non_whiterot_pred <- seq(min(data$sap_non_whiterot), max(data$sap_non_whiterot), 0.1)
data #look at your data file and check it looks ok
str(data)

##### What to expect: ####

plot(data$biomass, data$Mt/data$M0)
abline(lm(data$Mt/data$M0 ~ data$biomass))
summary(lm(data$Mt/data$M0 ~ data$biomass))

plot(data$sap_whiterot, data$Mt/data$M0)
abline(lm(data$Mt/data$M0 ~ data$sap_whiterot))
summary(lm(data$Mt/data$M0 ~ data$sap_whiterot))

plot(data$sap__whiterot, data$Mt/data$M0)
abline(lm(data$Mt/data$M0 ~ data$sap_non_whiterot))
summary(lm(data$Mt/data$M0 ~ data$sap_non_whiterot))


###### 2. create the initial values (if not provided, JAGS will do this for you) ####
#NOTE that it is a list within a list
#this allows you to have multiple lists that specify the initial values if you run more than 1 chain
#as with the data object you must use '=' inside a list, you cannot use '<-'

inits=list(
list(
  sigmaM = 1,
  sigmak = .5,
  alpha = 0.1, ## Change if muk is logarithmized
  b_bm = 0,
  b_sap_whiterot = 0,
  b_sap_non_whiterot = 0,
  sigmablock = 1,
 k=rep(0.1,data$nobs) ## Change if muk is logarithmized
),
list(
 sigmaM = 0.1,
sigmak = 0.1,
alpha = 0.001, ## Change if muk is logarithmized
b_bm = 0.5,
b_sap_whiterot = 0.5,
b_sap_non_whiterot = 0.5,
sigmablock = 0.1,
k=rep(0.1,data$nobs) ## Change if muk is logarithmized
))

## If you have 3 chains you want to initialise then what would it look like?


#specify the model and compile the model
#----------------------------------------------------------------------


##### 3. name the JAGS model file #####

#NOTE: change model
model = "JAGS_pine_needle_whiterot_model.R"

## global model ---------------------------------------------------------------------------------------------                                          
load.module("dic") #adding deviance information criteria from user manual of JAGS

#4. compile the model
jm = jags.model(model,
data=data, 
n.adapt=10000, 
inits=inits, 
n.chains=2) 


#note you need to adjust this code if you don't provide JAGS with inits or if you increase the number of chains you run

#5. burn-in the model

burn.in=10000

update(jm, n.iter=burn.in) #this runs the number burn-in values so the model is ready for sampling

#6. generate samples from the posterior distribution CODA

samples=10000
n.thin=10

#checking chains
zc = coda.samples(jm,
variable.names=c("K", "sigmaM", "sigmak","alpha", "sigmablock", "b_bm","b_sap_whiterot","b_sap_non_whiterot","R2"),
n.iter=samples, 
thin=n.thin)

#output coda data and you can summaries these or visually inspect the chains
#if you don't know what this means read the JAGS manual or JAGS primer

summary(zc)

#model selection, extract dic, but not reporting
ms = jags.samples(jm,
                  variable.names=c("deviance","Mt.new"), 
                  n.iter=samples, 
                  thin=n.thin)
PPL = sum((data$Mt-summary(ms$Mt.new,mean)$stat)^2) + sum((summary(ms$Mt.new,sd)$stat)^2)
PPL

#waic

zc_combined <- combine.mcmc(zc)
ecdf <-  apply(zc_combined, 2, function(x) 1-ecdf(x)(0)) #probability of positive effect
ecdf[["b_sap_whiterot"]]
ecdf[["b_sap_non_whiterot"]]

#will show the mean estimate and SE and 95% CIs with the R2, dic and PPL
#LABEL CORRECTLY:
#write(capture.output(summary(zc),"3. ecdf",ecdf,"Posterior Predictive Loss",PPL), paste0("tables/",set,time,"_Bosetta_Ågren_pineneedle_linear_k_biomass.txt"))


zj_val <- jags.samples(jm,
                       variable.names = c("mean.y", "mean.y.new","pvalue.mean",
                                          "cv.y", "cv.y.new", "pvalue.cv",
                                          "fit", "fit.new", "pvalue","R2"),
                       n.iter = samples,
                       thin = n.thin)

pdf(paste0("figures/MCMC_chains_litter",set,time,"_Bosetta_Ågren_linear_whiterot_k.pdf"))
plot(zc); gelman.plot(zc);  #look at the chains to see stability and mixing
text(ecdf[["b_sap_whiterot"]])
## Fit of mean:
plot(zj_val$mean.y,
     zj_val$mean.y.new,
     xlab = "mean real",
     ylab = "mean simulated",
     cex = .05)
abline(0, 1)
p <- summary(zj_val$pvalue.mean, mean)
text(x = 0.5, y = 0.65, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)

## Fit of variance:
fit.variance <- plot(zj_val$cv.y,
     zj_val$cv.y.new,
     xlab = "cv real",
     ylab = "cv simulated",
     cex = .05)
abline(0,1)
p <- summary(zj_val$pvalue.cv, mean)
text(x = 0.15, y = .3, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)

## Overall fit:
overall.fit <- plot(zj_val$fit,
     zj_val$fit.new,
     xlab = "ssq real",
     ylab = "ssq simulated",
     cex = .05)
abline(0,1)
p <- summary(zj_val$pvalue, mean)
text(x = 0.2, y = 0.4, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)

dev.off()

#7. generate samples from the posterior distribution JAGS

#output jags.object data. For some things easier to manipulate and output from here
#if you don't know what this means, read the JAGS manual

zj = jags.samples(jm, 
                  variable.names=c("k_bm_pred"), 
                  n.iter=samples, 
                  thin=n.thin)

zj_wr = jags.samples(jm, 
variable.names=c("k_bm_sap_whiterot_pred"), 
n.iter=samples, 
thin=n.thin)

zj_n = jags.samples(jm, 
                  variable.names=c("k_bm_sap_non_whiterot_pred"), 
                  n.iter=samples, 
                  thin=n.thin)

##     


#### BIOMASS PREDICTION #####
pdf(paste0("figures/Prediction",set,time,"_Bosetta_Ågren_linear_fungal_biomass_k.pdf"))

#plotting prediction & 95%CIs using polygon
pred<-summary(zj$k_bm_pred, quantile, c(.025,.5,.975))$stat #get your prediction
x=data$bm_pred    #set your x-axis relative to your x.pred prediction
y=pred
plot(x,y[2,], col="blue", xlab="fungal biomass [scaled]", ylab="K [d-1]", cex=1.4, typ="l", tck=0.03, bty="l", ylim = c(min(y[1,]), max(y[3,]))) #plot the median prediction
polygon(c(x,rev(x)), c(y[1,], rev(y[3,])), col=alpha("blue", alpha=0.5)) #add the 95% CIs
lines(x,y[1,], lty="dashed", col="blue") #add edges to the polygon
lines(x,y[3,], lty="dashed", col="blue")
text(x = min(data$bm_pred)+abs(min(data$bm_pred)/2), y = max(pred)-min(pred), paste0("PPL=",round(PPL,4)), cex = 1.5)

#### SAP WHITE ROT BIOMASS PREDICTION #####
pdf(paste0("figures/Prediction_Pine_Needles",set,time,"_Bosetta_Ågren_linear_sap_whiterot_fungal_biomass_k.pdf"))

#plotting prediction & 95%CIs using polygon
pred<-summary(zj_wr$k_bm_sap_whiterot_pred, quantile, c(.025,.5,.975))$stat #get your prediction
x=data$bm_sap_whiterot_pred    #set your x-axis relative to your x.pred prediction
y=pred
plot(x,y[2,], col="blue", xlab="sap whiterot fungal biomass [scaled]", ylab="K [d-1]", cex=1.4, typ="l", tck=0.03, bty="l", ylim = c(min(y[1,]), max(y[3,]))) #plot the median prediction
polygon(c(x,rev(x)), c(y[1,], rev(y[3,])), col=alpha("blue", alpha=0.5)) #add the 95% CIs
lines(x,y[1,], lty="dashed", col="blue") #add edges to the polygon
lines(x,y[3,], lty="dashed", col="blue")
text(x = min(data$bm_sap_whiterot_pred)+abs(min(data$bm_sap_whiterot_pred)), y = max(pred)-min(pred), paste0("PPL=",round(PPL,4)), cex = 1.5)
dev.off()

#### SAP NON WHITE ROT BIOMASS PREDICTION #####
pdf(paste0("figures/Prediction_Pine_Needles",set,time,"_Bosetta_Ågren_linear_sap_non_whiterot_fungal_biomass_k.pdf"))

pred<-summary(zj_n$k_bm_sap_non_whiterot_pred, quantile, c(.025,.5,.975))$stat #get your prediction
x=data$bm_sap_non_whiterot_pred    #set your x-axis relative to your x.pred prediction
y=pred
plot(x,y[2,], col="blue", xlab="sap non white rot fungal biomass [scaled]", ylab="K [d-1]", cex=1.4, typ="l", tck=0.03, bty="l", ylim = c(min(y[1,]), max(y2[3,]))) #plot the median prediction
polygon(c(x,rev(x)), c(y[1,], rev(y[3,])), col=alpha("blue", alpha=0.5)) #add the 95% CIs
lines(x,y[1,], lty="dashed", col="blue") #add edges to the polygon
lines(x,y[3,], lty="dashed", col="blue")
text(x = min(data$bm_sap_non_whiterot_pred)+abs(min(data$bm_sap_non_whiterot_pred)/2), y = max(pred)-min(pred), paste0("PPL=",round(PPL,4)), cex = 1.5)

dev.off()
