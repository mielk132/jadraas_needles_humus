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

  ## Priors:

  sigmaM ~ dgamma(0.001,0.001)
  sigmak ~ dgamma(0.001,0.001)
  alpha ~ dgamma(0.001,0.001) ## Change if muk is logarithmized
  #b_sm ~ dnorm(0, 0.001)
  #b_temp ~ dnorm(0, 0.001)
  b_bm ~ dnorm(0, 0.001)
  sigmablock ~ dnorm(0, 0.001)
  
  ## likelihood: 
  
  #to do:
  #t could be degree days or days. test. 

  for(i in 1:n) {
     
      Mt[i] ~ dgamma(muM[i]^2/sigmaM^2, muM[i]/sigmaM^2) # moment matching 
      #Mt.new[i] ~ dgamma(muM[i]^2/sigmaM^2, muM[i]/sigmaM^2)   
      # muMsum <- muM + muIG
      # model for ingrowth
      muM[i] <- M0[i]*exp(-k[i]*t) ## + ingrowth in humus decomposition model # process model for both rounds with the second time step?
      ## t constant for now because only one timestep, if we add intermediate timestep this will change
        
      k[i] ~ dgamma(muk[i]^2/sigmak2^2, muk[i]/sigmak^2)
      muk[i] <- alpha + eblock[block[i]] +
                #b_sm*soilmositure[i] + ## Should be scaled!
                # b_temp*temp[i] + ## Should be scaled! ## Only needed if two timestep model
                #b_temp_2*temp[i]^2 +
                b_bm*biomass[i] ## Should be scaled!
  
  }

  for(j in 1:nblock){eblock[j] ~ dnorm(0, sigmablock)}

  ## Model validation:
  
  # Bayesian p values- pvalue.cv is very large (close to 1) or very small (close to 0) model fails to
  # adequately rep distribution, extreme would be <0.10 or >0.90, good representation bayesian p value #near 0.50

  # cv.y <- sd(Mt[])/mean(Mt[ ])
  # cv.y.new <- sd(Mt.new[])/mean(Mt.new[])
  # pvalue.cv <- step(cv.y.new-cv.y)
  
  # mean.y <-mean(Mt[])
  # mean.y.new <-mean(Mt.new[])
  # pvalue.mean <-step(mean.y.new - mean.y)
  #
  # for(j in 1:n){
  #   sq[j] <- (Mt[j]-muM[j])^2
  #   sq.new[j] <- (Mt.new[j]-muM[j])^2
  # }
  # 
  # fit <- sum(sq[])
  # fit.new <- sum(sq.new[])
  # pvalue <- step(fit-fit.new)     # bayesan p value
  
  ## Predictions:
  
  # for(k in 1:length(bm_pred)){
  #   k_bm_pred <- alpha + b_bm*bm_pred[k]
  # }
  
  
}

 
 