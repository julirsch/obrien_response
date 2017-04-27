# Load libraries
library(e1071)

# Function to compute MSE from regression for parameter tuning

# sample, Vector of sample gene expression values.
# sigMat, The referent profile signature matrix.
# nu_ks, Parameter for inferring proportion of support vectors; define several.

getMSE <- function(sample, sigMat, nu_ks = seq(0.1, 0.9, 0.05)){
  mses <- sapply(nu_ks, function(k){
    svr_model <- svm(sigMat, sample, sigMatscale = FALSE, kernel = "linear", type = "nu-regression", nu = k)
    predictedm <- predict(svr_model, sigMat)
    error <- sample - predictedm
    mean(error^2)
  })
  names(mses) <- as.character(nu_ks)
  mses
}

# Function to do the modeling for a given sample
# Inputs are the sample index, sample data frame and signature Matrix

modelSample <- function(idx, Y, sigMat){
  
  # Get best k
  mses <- getMSE(Y[,idx,drop = FALSE], sigMat)
  best_k <- as.numeric(names(which(min(mses) == mses)))
  
  # Run model for best k; get coefficients
  tunedModel <- svm(sigMat, Y[,idx,drop = FALSE], sigMatscale = FALSE, kernel = "linear", type = "nu-regression", nu = best_k)
  
  # Estimate coefficients; remove negative values; normalize contributions
  coefTuned <- t(tunedModel$coefs) %*% tunedModel$SV
  coefTuned[coefTuned < 0 ] <- 0
  coefTunedNormed <- round(coefTuned /sum(coefTuned), 2)
  return(coefTunedNormed)
}

