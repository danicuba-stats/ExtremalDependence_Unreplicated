### =========================================================================== ###
### Auxiliary function to compute log likelihood function for non-censored data ###
### Daniela Castro-Camilo                                                       ###
### Email: daniela.castro.camilo@gmail.com                                      ###
### =========================================================================== ###

### =========================================================================== ###
### This code was adapted from it's original version to allow for application   ###
### on un-replicated data eg. soil samples.                                     ###
### Daniela Cuba                                                                ###
### =========================================================================== ###

## Parameters
# theta = c(lambda, rho (correlation))
# data.u: data matrix in uniform scale(n columns)

# The code originally removes NA observations. In the unreplicated application, 
# NA observations are removed from the dataset prior to fitting the model. 

log.lik = function(theta, data.u){
  lbda = theta[1]
  z = apply(data.u, 1:2, F1inv, lbda = lbda)
  sigma = matrix(c(1,theta[2],
                   theta[2],1),byrow=T, ncol=2)
  n = ncol(z)
  C = t(chol(sigma)) 
  C.inv = solve(C)
  sigma.inv = t(C.inv) %*% C.inv
  
  log.f1 = function(zij, lbda){
    if(is.na(zij)) return(NA)
    log(lbda) + lbda^2/2 - lbda * zij + pnorm(zij - lbda, log.p = T)
  }
  
  log.fn = function(zi, lbda, sigma, sigma.inv, C.inv){
    ## To remove observations with NAs
    if(any(is.na(zi))){
      id.na = which(is.na(zi))
      zi = zi[!is.na(zi)]
      n = length(zi)
      sigma = sigma[-id.na, -id.na]
      sigma.inv = sigma.inv[-id.na, -id.na]
      C.inv = C.inv[-id.na, -id.na]
    }
    ##
    m1.1 = zi%*%t(C.inv)
    m1 = t(m1.1)%*%m1.1
    m1 = t(zi) %*% sigma.inv %*% zi
    m2 = t(rep(1, n)) %*% sigma.inv %*% zi
    m3 = sum(sigma.inv)
    m1.star = (m2 - lbda)/sqrt(m3)
    out = log(lbda)  - (n - 1)/2 * log(2 * pi) - sum(log(diag(C))) - (1/2) * log(m3) + ((m1.star) ^ 2 - m1)/2 + pnorm(m1.star, log.p = T)
    return(out)
  }
  
  lfn = apply(z, 1, log.fn, lbda, sigma, sigma.inv, C.inv)
  lf1 = apply(z, c(1,2), log.f1, lbda)
  value <- sum(lfn) - sum(lf1, na.rm = T)
  value
}
