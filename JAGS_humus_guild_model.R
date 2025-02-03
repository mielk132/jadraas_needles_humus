### Target: NewPhyt
### Prelim Title: Mycorrhizal fungal guild interactions slow humus decomposition in a boreal forest

### Corresponding author
### Louis Mielke - louis.mielke@slu.se

###################################################################################

model{
  

  ## Priors:
  alpha ~ dnorm(0, 0.001) # intercept
  beta ~ dnorm(0, 0.001) # slope
  sigma ~ dgamma(0.001,0.001) # 
  b_ecto ~ dnorm(0, 0.001)
  b_sap ~ dnorm(0, 0.001)
  b_ericoid ~ dnorm(0, 0.001)
  b_molds_yeasts ~ dnorm(0, 0.001)
  #b_archies ~ dnorm(0, 0.001)
  #b_root_fungi ~ dnorm(0, 0.001)
  sigmablock ~ dt(0, pow(2.5,-2), 1)T(0,)
  
  ## likelihood: 

  
  for(l in 1:nobs) {
   
    M[l] ~ dnorm(mu[l],1/sigma^2)
    M.new[l] ~ dnorm(mu[l],1/sigma^2)
    
    
    mu[l] <- alpha + eblock[block[l]] + b_ecto*ecto[l,1] + b_sap*sap[l,1] +  b_molds_yeasts*molds_yeasts[l,1] + b_ericoid*ericoid[l,1]  #+ b_root_fungi*root_fungi[l,1] +   b_archies*archies[l,1]  
  }
  
  for(j in 1:max(block)){eblock[j] ~ dnorm(0, sigmablock)}
       
  
  ## Model validation:
  
  # Bayesian p values- pvalue.cv is very large (close to 1) or very small (close to 0) model fails to
  # adequately rep distribution, extreme would be <0.10 or >0.90, good representation bayesian p value #near 0.50

  cv.y <- sd(M[])/mean(M[])
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
  R2 <- 1 - sum(sq[])/sum(sq_mean[]) # coefficient of variation

  # Predictions:

for(k in 1:length(bm_ecto_pred)){mass_ecto_pred[k] <- alpha + b_ecto*bm_ecto_pred[k]}
  for(k in 1:length(bm_sap_pred)){mass_sap_pred[k] <- alpha + b_sap*bm_sap_pred[k]}
  for(k in 1:length(bm_ericoid_pred)){mass_ericoid_pred[k] <- alpha + b_ericoid*bm_ericoid_pred[k]}
  for(k in 1:length(bm_molds_yeasts_pred)){mass_molds_yeasts_pred[k] <- alpha + b_molds_yeasts*bm_molds_yeasts_pred[k]}
}

  #for(k in 1:length(bm_root_fungi_pred)){mass_root_fungi_pred[k] <- alpha + b_root_fungi*bm_root_fungi_pred[k]}
  #for(k in 1:length(bm_archies_pred)){mass_archies_pred[k] <- alpha + b_archies*bm_archies_pred[k]}
