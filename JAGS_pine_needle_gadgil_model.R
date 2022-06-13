### Target: Ecology Letters
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

model{ 
  
  ## Priors:
  
  sigmaM ~ dt(0, pow(2.5,-2), 1)T(0,) #change to half-cauchy
  sigmak ~ dt(0, pow(2.5,-2), 1)T(0,) #change to half-cauchy
  alpha ~ dnorm(0, 0.001) ## Change if muk is logarithmized

  b_ecto ~ dnorm(0, 0.001)
  b_sap ~ dnorm(0, 0.001)
  
  sigmablock ~ dt(0, pow(2.5,-2), 1)T(0,)
  
  ## likelihood: 
  
  for(i in 1:nobs) {
    
    Mt[i] ~ dgamma(muM[i]^2/sigmaM^2, muM[i]/sigmaM^2) # moment matching
    Mt.new[i] ~ dgamma(muM[i]^2/sigmaM^2, muM[i]/sigmaM^2)
    # muMsum <- muM + muIG
    # model for ingrowth
    # muM[i] <- M0[i]*max(exp(-k[i]*t[i]), S) ## + ingrowth in humus decomposition model # old model for two time steps
    # t = growing season days and is constant if we analyze with only one time step
 
    muM[i] <- M0[i]*(1+k[i]*t[i])^(-1.19) #adapted from Bosetta & Ågren 1998, Clemmensen et al. 2013
    
    k[i] ~ dgamma(muk[i]^2/sigmak^2, muk[i]/sigmak^2)
    # muk[i] <- max(muk_real[i], 1/10E6) #why did we log
    log(muk[i]) <- alpha + eblock[block[i]] + b_sap*sap[i,1] + b_ecto*ecto[i,1] 
    }
  
  for(j in 1:max(block)){eblock[j] ~ dnorm(0, sigmablock)}
  
  ## Model validation:
  
  # Bayesian p values- pvalue.cv is very large (close to 1) or very small (close to 0) model fails to
  # adequately rep distribution, extreme would be <0.10 or >0.90, good representation bayesian p value #near 0.50

  cv.y <- sd(Mt[])/mean(Mt[ ])
  cv.y.new <- sd(Mt.new[])/mean(Mt.new[])
  pvalue.cv <- step(cv.y.new-cv.y)

  mean.y <-mean(Mt[])
  mean.y.new <-mean(Mt.new[])
  pvalue.mean <-step(mean.y.new - mean.y)

  for(j in 1:nobs){
    sq[j] <- (Mt[j]-muM[j])^2
    sq.new[j] <- (Mt.new[j]-muM[j])^2
  }

  fit <- sum(sq[])
  fit.new <- sum(sq.new[])
  pvalue <- step(fit-fit.new)     # bayesan p value

  ## Predictions:
  
  #K <- mean(k[])
  
  for(k in 1:length(bm_ecto_pred)){
   log(k_bm_ecto_pred[k]) <- alpha + b_ecto*bm_ecto_pred[k] #+  b_sap*sap[k] 
   
 }
  
  }
  
  