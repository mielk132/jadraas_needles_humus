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
  alpha ~ dnorm(0, 0.001) # intercept
  sigma ~ dgamma(0.001,0.001) # 
  b_bm ~ dnorm(0, 0.001) #biomass
  sigmablock ~ dt(0, pow(2.5,-2), 1)T(0,)
  
  ## likelihood: 
  
  for(l in 1:nobs) {
    
    
    M[l] ~ dnorm(mu[l],1/sigma^2)
    M.new[l] ~ dnorm(mu[l],1/sigma^2)
 
    mu[l] <- alpha + eblock[block[l]] + b_bm*biomass[l,1]
  }
  
  for(j in 1:max(block)){eblock[j] ~ dnorm(0, sigmablock)}
  
  ## Model validation:
  
  # Bayesian p values- pvalue.cv is very large (close to 1) or very small (close to 0) model fails to
  # adequately rep distribution, extreme would be <0.10 or >0.90, good representation bayesian p value #near 0.50

  cv.y <- sd(M[])/mean(M[ ])
  cv.y.new <- sd(M.new[])/mean(M.new[])
  pvalue.cv <- step(cv.y.new-cv.y)

  mean.y <-mean(M[])
  mean.y.new <-mean(M.new[])
  pvalue.mean <-step(mean.y.new - mean.y)

  for(j in 1:nobs){
    sq[j] <- (M[j]-mu[j])^2
    sq.new[j] <- (M.new[j]-mu[j])^2
    sq_mean[j] <- (M[j]-mean(M[]))^2 # to calculate R2
  }

  fit <- sum(sq[])
  fit.new <- sum(sq.new[])
  pvalue <- step(fit-fit.new)     # bayesan p value
  R2 <- 1- sum(sq[])/sum(sq_mean[]) # coefficient of variation

  ## Predictions:

  for(k in 1:length(bm_pred)){mass_bm_pred[k] <- alpha +  b_bm*bm_pred[k]
  }

}