### (Journal)
### Title:
###Authors: L A Mielke 1, B Lindahl 2, J Klein 3, K E Clemmensen 1
### Affiliations
### 1 
### 2 
### 3
### 4

###Corresponding author
###Louis Mielke
###
###louis.mielke@slu.se
###
###
###
###

###################################################################################

#this is the R side of the analysis. The model is written in a separate file

rm(list=ls()) #clear the workspace (only if you want to)
getwd()
dir.create("figures")
#import libraries, data etc
#-----------------------------------------------------------#

#need these pckages installed
library(boot)
library(rjags)
library(coda)
library(scales)
library(fastDummies)

#file.choose()
#must set your working directory to whereever your JAGS models will be
setwd("~/Projects/mycorrhizal removal/MeshBags/Jädraås_Gadgil_MycorrhizalTraits")

#read in your data file(s)
df <- read.table("metadata3.txt", sep="\t", header=T, row.names=1) #bring in the data
df <- subset(df, treatment %in% c("C","E","T","TE"))
#check data
names(df)

#making new columns with amounts or percentages of certain fungal groups
df$copies_fungi <-   df$percent_fungi*df$copies_DNA_per_g_substrate
df$copies_asco <-    df$percent_asco*df$copies_fungi  #should be positively related to mass remaining
df$copies_basid <-    df$percent_basid*df$copies_fungi #expecting to be neg related to mass remaining
df$copies_leotio <-   df$percent_leotio*df$copies_fungi  #should be positively related to mass remaining
df$copies_agaric <-   df$percent_agaric*df$copies_fungi  #expecting to be neg related to mass remaining
df$copies_mycena <-     df$percent_mycena*df$copies_fungi  #expecting to be neg related to mass remaining
df$copies_gymnopus <-   df$percent_gymnopus*df$copies_fungi  #expecting to be neg related to mass remaining 
df$copies_d_acicola <-  df$percent_d_acicola*df$copies_fungi  #expecting to be pos. related to mass remaining
df$copies_l_pinastri <- df$percent_l_pinastri*df$copies_fungi  #expecting to be neg related to mass remaining
df$copies_g_androsaceus <- df$percent_g_androsaceus*df$copies_fungi #expecting to be pos. relat. to mass remain.
df$copies_m_clavicularis <- df$percent_m_clavicularis*df$copies_fungi #expecting to be neg related to mass remaining

# cut & paste names in here - makes it easier when setting up your JAGS data object 
str(df) #check how R defines your variables
summary(df) #check for NA's in your explanatory variables and data range
df <- na.omit(df)

df$plot <- paste0(as.character(df$treatment) , as.character(df$block))
df <- na.omit(df)

#separate by substrate and set
litterA  <- subset(df, substrate=="Litter" & set =="A")
litterB  <- subset(df, substrate=="Litter" & set =="B")
humusA  <- subset(df, substrate=="Humus" & set =="A")
humusB  <- subset(df, substrate=="Humus" & set =="B")
# separate by incubation time
litterA5 <- subset(litterA, incubation == "5")
litterA17 <- subset(litterA, incubation == "17")
litterB5 <- subset(litterB, incubation == "5")
litterB17 <- subset(litterB, incubation == "17")
humusA5 <- subset(humusA, incubation == "5")
humusA17 <- subset(humusA, incubation == "17")
humusB5 <- subset(humusB, incubation == "5")
humusB17 <- subset(humusB, incubation == "17")

#PLOTTING 
#a loop for plotting percent mass remaining as y and copies or percent of the fungal group on the x axis
my_colnames <- names(litterA) 

loop.vector <- 31:53
for (i in loop.vector) {
  
  x <- litterA[,i]
  y <- litterA[,(22)]
  
  # Apply names function
  my_colname_x <- my_colnames[i]
  my_colname_y <- my_colnames[22]
  
  plot(x,y, xlab = my_colname_x, ylab = my_colname_y)+
  abline(lm(y ~ x))
}


#set up conditions for JAGS model
#-----------------------------------------------------------#

#1. create the data structure
data = NULL       #clear any old data lists that might confuse things
data = list(
	t = litterB$growing_season_days, #incubaton length (5 or 17 months)
	M0 = litterB$original_substrate_mass_g, ## Why not the same for 5 and 17 months?
	Mt = litterB$final_substrate_g,
	nobs = nrow(litterB),
	plot = as.numeric(as.factor(litterB$plot)),
	block = as.numeric(litterB$block),
	biomass = scale(litterB$copies_fungi)
) 

data$bm_pred <- seq(min(data$biomass), max(data$biomass), 0.1)
data #look at your data file and check it looks ok

## What to expect:

plot(data$biomass, data$Mt/data$M0)
abline(lm(data$Mt/data$M0 ~ data$biomass))
summary(lm(data$Mt/data$M0 ~ data$biomass))

#2. create the initial values (if not provided, JAGS will do this for you)
#NOTE that it is a list within a list
#this allows you to have multiple lists that specify the initial values if you run more than 1 chain
#as with the data object you must use '=' inside a list, you cannot use '<-'

inits=list(
list(
  sigmaM = 1,
  sigmak = .5,
  alpha = 0.1, ## Change if muk is logarithmized
  #b_sm ~ dnorm(0, 0.001)
  #b_temp ~ dnorm(0, 0.001)
  b_bm = 0,
  #b_bm2 = 0,
  sigmablock = 1,
  sigmaplot = 1,
  k=rep(0.1,data$nobs)
  #t=
),
list(
 sigmaM = 0.1,
sigmak = 0.1,
alpha = 0.001, ## Change if muk is logarithmized
#b_sm ~ dnorm(0, 0.001)
#b_temp ~ dnorm(0, 0.001)
b_bm = -0.5,
#b_bm2 = -0.5,
sigmablock = 0.1,
sigmaplot = 0.1,
 k=rep(0.1,data$nobs)
))

## If you have 3 chains you want to initialise then it would look like this


#specify the model and compile the model
#----------------------------------------------------------------------


#3. name the JAGS model file

model = "JAGS_pine_needle_model.R"

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

samples=50000
n.thin=50

zc = coda.samples(jm,
variable.names=c("K", "sigmaM", "sigmak","alpha","b_bm", "sigmablock","sigmaplot"),
n.iter=samples, 
thin=n.thin)

#output coda data and you can summaries these or visually inspect the chains
#if you don't know what this means read the JAGS manual or JAGS primer

# summary(zc) #will show the mean estimate and SE and 95% CIs

pdf("figures/MCMC_chains_litterB_Bosetta_Ågren_linearbm.pdf")
plot(zc); gelman.plot(zc) #look at the chains to see stability and mixing
dev.off()

# zj_val <- jags.samples(jm,
#                        variable.names = c("mean.y", "mean.y.new","pvalue.mean",
#                                           "cv.y", "cv.y.new", "pvalue.cv",
#                                           "fit", "fit.new", "pvalue"),
#                        n.iter = samples,
#                        thin = n.thin)
# 
# ## Fit of mean:
# plot(zj_val$mean.y,
#      zj_val$mean.y.new,
#      xlab = "mean real",
#      ylab = "mean simulated",
#      cex = .05)
# abline(0, 1)
# p <- summary(zj_val$pvalue.mean, mean)
# text(x = 0.6, y = 0.7, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)
# 
# ## Fit of variance:
# plot(zj_val$cv.y,
#      zj_val$cv.y.new,
#      xlab = "cv real",
#      ylab = "cv simulated",
#      cex = .05)
# abline(0,1)
# p <- summary(zj_val$pvalue.cv, mean)
# text(x = 0.15, y = .5, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)
# 
# ## Overall fit:
# plot(zj_val$fit,
#      zj_val$fit.new,
#      xlab = "ssq real",
#      ylab = "ssq simulated",
#      cex = .05)
# abline(0,1)
# p <- summary(zj_val$pvalue, mean)
# text(x = 0.2, y = 0.4, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)

#7. generate samples from the posterior distribution JAGS

#output jags.object data. For some things easier to manipulate and output from here
#if you don't know what this means, read the JAGS manual

zj = jags.samples(jm, 
variable.names=c("k_bm_pred"), 
n.iter=samples, 
thin=n.thin)

#plotting prediction & 95%CIs using polygon
pred<-summary(zj$k_bm_pred, quantile, c(.025,.5,.975))$stat #get your prediction
x=data$bm_pred    #set your x-axis relative to your x.pred prediction
y=pred
plot(x,y[2,], col="blue", xlab="fungal biomass [scaled]", ylab="K [d-1]", cex=1.4, typ="l", tck=0.03, bty="l", ylim = c(min(y[1,]), max(y[3,]))) #plot the median prediction
polygon(c(x,rev(x)), c(y[1,], rev(y[3,])), col=alpha("blue", alpha=0.5)) #add the 95% CIs
lines(x,y[1,], lty="dashed", col="blue") #add edges to the polygon
lines(x,y[3,], lty="dashed", col="blue")


