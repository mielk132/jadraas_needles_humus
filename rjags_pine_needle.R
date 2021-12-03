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
#setwd("~")

#read in your data file(s)
df <- read.table("metadata.txt", header=T) #bring in the data

#check your data
names(df)
# cut & paste names in here - makes it easier when setting up your JAGS data object 
str(df) #check how R defines your variables
summary(df) #check for NA's in your explanatory variables and data range

df <- dummy_cols(df, select_columns = "treatment")

#set up conditions for JAGS model
#-----------------------------------------------------------#

#1. create the data structure
data=NULL       #clear any old data lists that might confuse things
data=list(
	y=df$weight,
	xa=df$treatment_a,
	xb=df$treatment_b,
	xc=df$treatment_c,
	n.obs=length(df$weight),
	x.pred=1:30 #if creating a range of predictions
) 

data #look at your data file and check it looks ok

#NOTE you must use '=' inside a list, you cannot use '<-'

#2. create the initial values (if not provided, JAGS will do this for you)
#NOTE that it is a list within a list
#this allows you to have multiple lists that specify the initial values if you run more than 1 chain
#as with the data object you must use '=' inside a list, you cannot use '<-'

inits=list(
list(
a=0.1,
b=2
))

## If you have 3 chains you want to initialise then it would look like this

inits=list(
list(
a=0.1,
b=2
),
list(
a=5,
b=0.1
),
list(
a=-5,
b=2.5
))


#specify the model and compile the model
#----------------------------------------------------------------------


#3. name the JAGS model file

model = "name_of_JAGS_model_file.R"

#4. compile the model

jm = jags.model(model,
data=data, 
n.adapt=5000, 
inits=inits, 
n.chains=1) 

#note you need to adjust this code if you don't provide JAGS with inits or if you increase the number of chains you run


#5. burn-in the model

burn.in=10000

update(jm, n.iter=burn.in) #this runs the number burn-in values so the model is ready for sampling

#6. generate samples from the posterior distribution CODA

samples=10000
n.thin=5

zc = coda.samples(jm,
variable.names=c("a", "b"), 
n.iter=samples, 
thin=n.thin)

#output coda data and you can summaries these or visually inspect the chains
#if you don't know what this means read the JAGS manual or JAGS primer


summary(zc) #will show the mean estimate and SE and 95% CIs
plot(zc) #look at the chains to see stability and mixing


#7. generate samples from the posterior distribution JAGS

#output jags.object data. For some things easier to manipulate and output from here
#if you don't know what this means, read the JAGS manual

zj = jags.samples(jm, 
variable.names=c("output"), 
n.iter=samples, 
thin=n.thin)

#generate medians and quantiles that can be used for storing info, plotting etc.
pred<-summary(zj$output, quantile, c(.025,.5,.975))$stat



#diagnostics for convergence
#------------------------------------------------------------------------
heidel.diag(z)
gelman.diag(z) #needs at least 2 chains



#plotting prediction & 95%CIs using polygon
pred<-summary(zj$output, quantile, c(.025,.5,.975))$stat #get your prediction
x=    #set your x-axis relative to your x.pred prediction
y=pred
plot(x,y[2,], col="blue", xlab="XX", ylab="YY", cex=1.4, typ="l", tck=0.03, bty="l") #plot the median prediction
polygon(c(x,rev(x)), c(y[1,], rev(y[3,])), col=alpha("blue", alpha=0.5)) #add the 95% CIs
lines(x,y[1,], lty="dashed", col="blue") #add edges to the polygon
lines(x,y[3,], lty="dashed", col="blue")

