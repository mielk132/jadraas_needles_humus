### (Journal)
### Title 
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

 model{ 

# parameter model- uninformed priors

  sigmap~dgamma(0.001,0.001)
  sigmapk~dgamma(0.001,0.001)
  alpha
  bshrub
  bpine
  bpineshrub
  
# likelihood 
  #to do
  #t could be degree days or days. test. 

  for( i in 1:n) {
     
      Mt[i] ~ dgamma(mu[i]^2/sigmap^2,mu[i]/sigmap^2) # moment matching 
      mu[i]<- Mo[i]*exp(-k[i]*t[i])  # process model for both rounds with the second time step?
      #Mt.new[i] ~ dgamma(mu[i]^2/sigmap^2,mu[i]/sigmap^2)   
        
      k[i] ~ dgamma(muk[i]^2/sigmapk2^2,muk[i]/sigmapk^2)
      muk[i] <- alpha + bpine*pine[i]+ 
                bshrub*shrub[i]+
                bpineshrub*pine[i]*shrubs[i]+eblock[block[i]]
        
} # end of loop

  for( j in 1:nblock) {
   
     eblock[j] ~ dnorm(0,sigmablock)
    
  } # end of loop
  
  # Bayesian p values- pvalue.cv is very large (close to 1) or very small (close to 0) model fails to
   #adequately rep distribution, extreme would be <0.10 or >0.90, good representation bayesian p value #near 0.50

#cv.y <- sd(Mt[ ])/mean(Mt[ ])
#cv.y.new <- sd(Mt.new[])/mean(Mt.new[ ])
#pvalue.cv <- step(cv.y.new-cv.y)

# mean.y <-mean(Mt[])
# mean.y.new <-mean(Mt.new[])
# pvalue.mean <-step(mean.y.new - mean.y)
# 
# 
# for (j in 1:n){
#   sq[j]<-(Mt[j]-mu[j])^2
#   sq.new[j]<-(Mt.new[j]-mu[j])^2
#   
# }
# 
# fit <- sum(sq[])
# fit.new <- sum(sq.new[])
# pvalue <- step(fit-fit.new)     # bayesan p value
# 

}

 
 