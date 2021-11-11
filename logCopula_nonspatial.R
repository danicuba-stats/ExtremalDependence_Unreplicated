### ========================================================================= ###
### Auxiliary function to compute log copula function for fully censored data ###
### Daniela Castro-Camilo                                                     ###
### Email: daniela.castro.camilo@gmail.com                                    ###
### ========================================================================= ###

### ========================================================================= ###
### This code has undergone changes for adaptation to a non-spatial setting   ###
### and application to unreplicated-settings (non-spatial).                   ###
### Daniela Cuba                                                              ###
### ========================================================================= ###

# Parameters
# theta = c(lambda, rho (correlation))
# u.star: threshold (vector) in uniform scale

log.Cn = function(theta, u.star){ 
  # Warning: works only for "sigma" such that diag(sigma) = 1
  lbda = theta[1]
  z = sapply(u.star, F1inv, lbda = lbda) # = z.star
  sigma=matrix(c(1,theta[2],
                 theta[2],1),byrow=T,nrow=2,ncol=2)
  n = ncol(sigma)
  mean = rep(0, n)
  probs = NULL
  condVars = lapply(1:n, mi_condMVN, mean = mean, sigma = sigma)
  unos = rep(1, (n - 1))
  sigma0 = matrix(NA, nrow = n, ncol = n)
  sigma0[n, n] = 1
  for(i in 1 : n){
    tmp = condVars[[i]]
    condVar = tmp$condVar
    C = tmp$C
    D = (unos - C); Dt = t(D)
    sigma0[1:(n - 1), 1:(n - 1)] = condVar + D %*% Dt
    sigma0[1:(n - 1), n] = -D
    sigma0[n, 1:(n - 1)] = -Dt
    upper = c(z[-i] - z[i]*unos + D * lbda, z[i] - lbda)
    set.seed(302)
    # browser()
    probs[i] = exp(lbda ^ 2/2 - lbda * z[i]) *pmvnorm(lower = -Inf, upper = upper, sigma = sigma0, algorithm = GenzBretz())[1]
  }
  set.seed(302)
  value = as.numeric(pmvnorm(lower = -Inf, upper = z, corr = sigma)[1] - sum(probs))
  log(value)
}
