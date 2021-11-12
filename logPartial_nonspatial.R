<<<<<<< HEAD
### ============================================================================== ###
### Auxiliary functions to compute log copula function for partially censored data ###
### Daniela Castro-Camilo                                                          ###
### Email: daniela.castro.camilo@gmail.com                                         ###
### ============================================================================== ###

### =========================================================================== ###
### This code was adapted from it's original version to allow for application   ###
### on un-replicated data eg. soil samples.                                     ###
### Daniela Cuba                                                                ###
### =========================================================================== ###

## Parameters
# theta = c(lambda, rho (correlation))
# data.u: data matrix in uniform scale(n columns)
# u.star: treshold in u-scale

dF1 <- function(z, theta, log = FALSE){
  lbda = theta[1]
  logdF1 <- log(lbda) -lbda*z + 0.5*lbda^2 + pnorm(z-lbda,log.p=TRUE)
  if(log){
    return(logdF1)
  } else{
    return(exp(logdF1))
  }
}

dFI <- function(z, I, theta, log = FALSE){
  lbda = theta[1]
  # browser()
  sigma = matrix(c(1,theta[2],theta[2],1),
                 byrow=T,nrow=2,ncol=2)
  D = 2
  k = length(I)
  
  sigmaII = sigma[I,I]
  C = t(chol(sigmaII))
  Cinv = solve(C)
  sigmaIIinv = t(Cinv) %*% Cinv
  CinvzI = Cinv %*% z[I]
  m1 = t(CinvzI) %*% CinvzI
  m2 = rep(1,k) %*% sigmaIIinv %*% z[I]
  m3 = sum(sigmaIIinv)
  
  
  z.I = rbind(0, z[-I] - sigma[-I,I] %*% sigmaIIinv %*% z[I])
  mu.I = as.numeric((m2 - lbda)/m3) * rbind(-1, 1 - sigma[-I,I] %*% sigmaIIinv %*% rep(1,k)) 
  sigma.I = matrix(NA, D - k + 1, D - k + 1)
  sigma.I[1,1] = m3^(-1)
  sigma.I[-1,1] = (sigma[-I,I] %*% sigmaIIinv %*% rep(1, k) - 1)/m3
  sigma.I[1,-1] = t(sigma[-I,I] %*% sigmaIIinv %*% rep(1, k) - 1)/m3
  sigma.I[-1,-1] = sigma[-I,-I] - sigma[-I,I] %*% sigmaIIinv %*% sigma[I,-I] + (sigma[-I,I] %*% sigmaIIinv %*% rep(1, k) - 1) %*% t(sigma[-I,I] %*% sigmaIIinv %*% rep(1, k) - 1)/m3
  
  if (!isSymmetric(sigma.I))
    sigma.I = as.matrix(forceSymmetric(sigma.I))
  logdFI = log(lbda) -0.5 * log(m3) -(k - 1) * 0.5 * log(2 * pi) - sum(log(diag(C))) - 0.5 * (m1 - (m2 - lbda) ^ 2/m3) + log(pmvnorm(upper = as.numeric(z.I), mean = as.numeric(mu.I), sigma = sigma.I)[1])
  
  if(log){
    return(logdFI)
  } else{
    return(exp(logdFI))
  }
}

dCI <- function(u, I, theta, log = FALSE){
  z <- qF1(u,theta)
  logdCI <- dFI(z, I, theta,log = TRUE) - sum(dF1(z[I], theta, log = TRUE))
  if(log){
    return(logdCI)
  } else{
    return(exp(logdCI))
  }
}

log.partialCn = function(theta, data.u, u.star){
  value = NULL
  if( is.null(nrow(data.u)) ) data.u = matrix(data.u, nrow = 1, ncol = length(data.u)) # if data.u has only one row
  for(i in 1:nrow(data.u)){
    u = data.u[i, ]
    if(any(is.na(u))){
      id.na = which(is.na(u))
      u.star1 = u.star[!is.na(u)]
      u1 = u[!is.na(u)]
      I = which(u1 > u.star1)
      u1[-I] = u.star1[-I]
      value[i] = dCI(u1, I, theta, log = TRUE)
    }else{
      I = which(u > u.star)
      u[-I] = u.star[-I]
      value[i] = dCI(u, I,theta, log = TRUE) 
    }
  }
  sum(value)
}
=======
### ============================================================================== ###
### Auxiliary functions to compute log copula function for partially censored data ###
### Daniela Castro-Camilo                                                          ###
### Email: daniela.castro.camilo@gmail.com                                         ###
### ============================================================================== ###

### =========================================================================== ###
### This code was adapted from it's original version to allow for application   ###
### on un-replicated data eg. soil samples.                                     ###
### Daniela Cuba                                                                ###
### =========================================================================== ###

## Parameters
# theta = c(lambda, rho (correlation))
# data.u: data matrix in uniform scale(n columns)
# u.star: treshold in u-scale

dF1 <- function(z, theta, log = FALSE){
  lbda = theta[1]
  logdF1 <- log(lbda) -lbda*z + 0.5*lbda^2 + pnorm(z-lbda,log.p=TRUE)
  if(log){
    return(logdF1)
  } else{
    return(exp(logdF1))
  }
}

dFI <- function(z, I, theta, log = FALSE){
  lbda = theta[1]
  # browser()
  sigma = matrix(c(1,theta[2],theta[2],1),
                 byrow=T,nrow=2,ncol=2)
  D = 2
  k = length(I)
  
  sigmaII = sigma[I,I]
  C = t(chol(sigmaII))
  Cinv = solve(C)
  sigmaIIinv = t(Cinv) %*% Cinv
  CinvzI = Cinv %*% z[I]
  m1 = t(CinvzI) %*% CinvzI
  m2 = rep(1,k) %*% sigmaIIinv %*% z[I]
  m3 = sum(sigmaIIinv)
  
  
  z.I = rbind(0, z[-I] - sigma[-I,I] %*% sigmaIIinv %*% z[I])
  mu.I = as.numeric((m2 - lbda)/m3) * rbind(-1, 1 - sigma[-I,I] %*% sigmaIIinv %*% rep(1,k)) 
  sigma.I = matrix(NA, D - k + 1, D - k + 1)
  sigma.I[1,1] = m3^(-1)
  sigma.I[-1,1] = (sigma[-I,I] %*% sigmaIIinv %*% rep(1, k) - 1)/m3
  sigma.I[1,-1] = t(sigma[-I,I] %*% sigmaIIinv %*% rep(1, k) - 1)/m3
  sigma.I[-1,-1] = sigma[-I,-I] - sigma[-I,I] %*% sigmaIIinv %*% sigma[I,-I] + (sigma[-I,I] %*% sigmaIIinv %*% rep(1, k) - 1) %*% t(sigma[-I,I] %*% sigmaIIinv %*% rep(1, k) - 1)/m3
  
  if (!isSymmetric(sigma.I))
    sigma.I = as.matrix(forceSymmetric(sigma.I))
  logdFI = log(lbda) -0.5 * log(m3) -(k - 1) * 0.5 * log(2 * pi) - sum(log(diag(C))) - 0.5 * (m1 - (m2 - lbda) ^ 2/m3) + log(pmvnorm(upper = as.numeric(z.I), mean = as.numeric(mu.I), sigma = sigma.I)[1])
  
  if(log){
    return(logdFI)
  } else{
    return(exp(logdFI))
  }
}

dCI <- function(u, I, theta, log = FALSE){
  z <- qF1(u,theta)
  logdCI <- dFI(z, I, theta,log = TRUE) - sum(dF1(z[I], theta, log = TRUE))
  if(log){
    return(logdCI)
  } else{
    return(exp(logdCI))
  }
}

log.partialCn = function(theta, data.u, u.star){
  value = NULL
  if( is.null(nrow(data.u)) ) data.u = matrix(data.u, nrow = 1, ncol = length(data.u)) # if data.u has only one row
  for(i in 1:nrow(data.u)){
    u = data.u[i, ]
    if(any(is.na(u))){
      id.na = which(is.na(u))
      u.star1 = u.star[!is.na(u)]
      u1 = u[!is.na(u)]
      I = which(u1 > u.star1)
      u1[-I] = u.star1[-I]
      value[i] = dCI(u1, I, theta, log = TRUE)
    }else{
      I = which(u > u.star)
      u[-I] = u.star[-I]
      value[i] = dCI(u, I,theta, log = TRUE) 
    }
  }
  sum(value)
}
>>>>>>> 8905faec071b909b4125d27c75dda9d5c3590b9c
