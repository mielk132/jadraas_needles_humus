### Target: New Phyt
### Prelim Title: Mycorrhizal fungal guild interactions slow humus decomposition in a boreal forest

### Corresponding author
### Louis Mielke - louis.mielke@stir.ac.uk

###################################################################################

model{ 
  
  ## Priors:
  
  sigmaM ~ dt(0, pow(2.5,-2), 1)T(0,) #change to half-cauchy
  sigmak ~ dt(0, pow(2.5,-2), 1)T(0,) #change to half-cauchy
  alpha ~ dnorm(0, 0.001) ## Change if muk is logarithmized

  b_bm~ dnorm(0, 0.001)
  sigmablock ~ dt(0, pow(2.5,-2), 1)T(0,)
  
  ## likelihood: 
  
  for(i in 1:nobs) {
    
    Mt[i] ~ dgamma(muM[i]^2/sigmaM^2, muM[i]/sigmaM^2) # moment matching
    Mt.new[i] ~ dgamma(muM[i]^2/sigmaM^2, muM[i]/sigmaM^2)

    muM[i] <- M0[i]*(1+k[i]*t[i])^(-1.19) #adapted from Bosetta & Ågren 1998, Clemmensen et al. 2013
    
    k[i] ~ dgamma(muk[i]^2/sigmak^2, muk[i]/sigmak^2)

    log(muk[i]) <- alpha + eblock[block[i]] + b_bm*biomass[i,1]
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
    sq_mean[j] <- (Mt[j]-mean(Mt[]))^2 # to calculate R2
  }

  fit <- sum(sq[])
  fit.new <- sum(sq.new[])
  pvalue <- step(fit-fit.new)     # bayesan p value
  R2 <- 1- sum(sq[])/sum(sq_mean[]) # coefficient of variation
  #how much of the original variation in Mt is not explained by the model / mean
  
  ## Predictions:
  
  #K <- mean(k[])

  for(k in 1:length(bm_pred)){ log(k_bm_pred[k]) <- alpha + b_bm*bm_pred[k] }

  }
  
  
