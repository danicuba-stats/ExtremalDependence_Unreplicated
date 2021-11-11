### ========================================================================================================= ###
### Censored local log-likelihood and model fitting functions for the U.S. precipitation extremes application ###
### Daniela Castro-Camilo                                                                                     ###
### Email: daniela.castro.camilo@gmail.com                                                                    ###
### ========================================================================================================= ###

### ========================================================================== ###
### This code has undergone alterations from its original state to accommodate ###                              
### the possibility of unreplicated samples . All functions included here have ###
### also been altered .                                                        ###
### The Application example is a simulation study to prove the efficacy of the ###
### approach.                                                                  ###
### ========================================================================== ###

#####################################
### Libraries and auxiliary codes ###
#####################################
# install.packages(c('mvtnorm', 'R.utils', 'fields', 'numDeriv', 'condMVNorm', 'Matrix', 'matrixcalc'))
libs<-c("mvtnorm","R.utils","fields","numDeriv","condMVNorm","Matrix",
        "matrixcalc","MASS","dplyr")
lapply(libs,require,character.only=T)

source('./Tools.R')
source('./logCopula_nonspatial.R')
source('./logLikelihood_nonspatial.R')
source('./logPartial_nonspatial.R')

##############################################
### Censored local log-likelihood function ###
##############################################

## Parameters
# theta [vector]: c(lambda, range), parameter vector
# data.u [matrix 2070 x D0]: data matrix in uniform scale. D0 is the number of neighbors
# coord [matrix D0 x 2]: coordinates for all the stations
# thres [numeric]: probability for the quantile-based treshold for each location
# nu [numeric]: fixed smoothing parameter
# censorL [logical]: should we compute censorized likelihood? Default to be true

#### Alterations
# theta [vector]: c(lambda, rho (correlation)), 
# data.u [matrix n x D]: data matrix in uniform scale. n observations, D components
# thres [numeric]: probability for the quantile-based treshold for each location
# censorL [logical]: should we compute censorized likelihood? Default to be true

model.likelihood.optim = function(theta, data.u, thres, censorL){
  lbda = theta[1]
  # browser()
  if(lbda <= 0 ){
    return(1e09)
  }else{
    u.star = apply(data.u, 2, quantile, probs = thres, na.rm = T) # set threshold
    
    # Classify the data into fully censored, non-censored, or partially censored
    N1 = dim(data.u)[1]; n1 = dim(data.u)[2]
    I = matrix(rep(0, N1 * n1), N1, n1) # indicator matrix
    for(i in 1:N1){
      # cond0 = data.u[i, ] <= u.star
      cond1 = data.u[i, ] > u.star
      # I[i, which(cond0)] = 0
      I[i, which(cond1)] = 1
    }
    # browser()
    G1 = which(rowSums(I, na.rm = T) == 0) # fully censored
    G2 = which(rowSums(I, na.rm = T) == n1) # non-censored
    G3 = which(rowSums(I, na.rm = T) > 0 & rowSums(I, na.rm = T) < n1) # partially censored
    
    # Computing fully censored likelihood
    if(length(G1) > 0){
      if(censorL == TRUE){
        value1 = -length(G1) * log.Cn(theta, u.star)
      }
      else{
        value1 = -log.lik(theta, data.u[G1, ])
      }
    }
    # Computing non-censored likelihood
    if(length(G2) > 0){
      value2 = -log.lik(theta, data.u[G2, ])
    }
    # Computing partially censored likelihood
    if(length(G3) > 0){
      if(censorL == TRUE){
        value3 = -log.partialCn(theta, data.u[G3, ], u.star)
      } else {
        value3 = -log.lik(theta, data.u[G3, ])
      }
    }
    # Negative log-likelihood
    return(sum(c(value1, value2, value3))) 
  }
}

##############################
### Model fitting function ###
##############################
# Alteration - actually, this is necessary for the ease of fitting the model
# data [matrix nxD]: all the data already in uniform scale - this is due to the 
#                    complications of having modelled the marginals separately
# theta0 [vector 1x4]: c(lambda0, var1, cov12, var2), initial values
# thres [numeric]: as in model.likelihood
# censorL [logical]: as in model.likelihood
# starting: do you want a parameter grid to be fitted to find best initial values
 

fit.model.likelihood = function(data.u, thres, censorL = TRUE, starting=FALSE){
print(paste0("Start at: ",Sys.time()))
  ## data is already in uniform scale# transform to uniform margins
  
  ## Fit likelihood to a grid of potential starting values to find best initial values for optimisation
  if(starting){
    lams<-c(seq(0.1,5,0.25),seq(5,30,1))
    rs<-seq(-0.9,0.9,0.1)
    init.grid<-expand.grid(lams,rs)
    names(init.grid)<-c("lambda","rho")
    for(i in 1:nrow(init.grid)){
      init.grid$lik[i]<-tryCatch({model.likelihood.optim(theta=unlist(init.grid[i,]),
                                               data.u=data.u,
                                               thres=thres,censorL=T)},
                                 error=function(cond){
                                   return(NA)
                                 })
    }
    init.grid<-init.grid[order(init.grid$lik,decreasing=F),] # order initival values by increasing likelihood
    print("Finished exploration")
  }
  
  tmp<- list()
  
  # Fit likelihood
  for(i in 1:10) {
    if(starting){
      init.vars<-unlist(init.grid[i,c(1,2)])
    }else{
      init.vars<-c(runif(1,min=0.1,max=5),runif(1,min=-0.9,max=0.9))
    }
    # Optimise. If a mistake is found, produce NAs - typically happens if the optimisation does not converge
    tmp[[i]]<-tryCatch({optim(init.vars, model.likelihood.optim,  data.u = data.u, thres = thres,
                      censorL = censorL,
                      method="L-BFGS-B",lower=c(0.1,-0.9),upper=c(Inf,0.9))},
                      error=function(cond){
                        message(cond)
                        return(list(value=NA,par=c(NA,NA)))
                      })
    print(paste0("Optim: ",i))
  }

  # Optimise using the best initial values found in previou step
  try(tmp[[11]]<-optim(tmp[[which.min(unlist(lapply(tmp,function(x) {x$value})))]]$par, 
                   model.likelihood.optim,  data.u = data.u, thres = thres,
                   censorL = censorL,
                   method="L-BFGS-B",lower=c(0.1,-0.9),upper=c(Inf,0.9)))
  to.return<-tryCatch({tmp[[which.min(unlist(lapply(tmp, function(x) {x$value})))]]},
                      error=function(cond){
                        print("All Errors")
                        return(list(value=NA,par=c(NA,NA)))
                      })
  print(paste0("Completed at: ",Sys.time()))
  return(to.return)
}
=======
### ========================================================================================================= ###
### Censored local log-likelihood and model fitting functions for the U.S. precipitation extremes application ###
### Daniela Castro-Camilo                                                                                     ###
### Email: daniela.castro.camilo@gmail.com                                                                    ###
### ========================================================================================================= ###

### ========================================================================== ###
### This code has undergone alterations from its original state to accommodate ###                              
### the possibility of unreplicated samples . All functions included here have ###
### also been altered .                                                        ###
### The Application example is a simulation study to prove the efficacy of the ###
### approach.                                                                  ###
### ========================================================================== ###

#####################################
### Libraries and auxiliary codes ###
#####################################
# install.packages(c('mvtnorm', 'R.utils', 'fields', 'numDeriv', 'condMVNorm', 'Matrix', 'matrixcalc'))
libs<-c("mvtnorm","R.utils","fields","numDeriv","condMVNorm","Matrix",
        "matrixcalc","MASS","dplyr")
lapply(libs,require,character.only=T)

source('./Tools.R')
source('./logCopula_nonspatial.R')
source('./logLikelihood_nonspatial.R')
source('./logPartial_nonspatial.R')

##############################################
### Censored local log-likelihood function ###
##############################################

## Parameters
# theta [vector]: c(lambda, range), parameter vector
# data.u [matrix 2070 x D0]: data matrix in uniform scale. D0 is the number of neighbors
# coord [matrix D0 x 2]: coordinates for all the stations
# thres [numeric]: probability for the quantile-based treshold for each location
# nu [numeric]: fixed smoothing parameter
# censorL [logical]: should we compute censorized likelihood? Default to be true

#### Alterations
# theta [vector]: c(lambda, rho (correlation)), 
# data.u [matrix n x D]: data matrix in uniform scale. n observations, D components
# thres [numeric]: probability for the quantile-based treshold for each location
# censorL [logical]: should we compute censorized likelihood? Default to be true

model.likelihood.optim = function(theta, data.u, thres, censorL){
  lbda = theta[1]
  # browser()
  if(lbda <= 0 ){
    return(1e09)
  }else{
    u.star = apply(data.u, 2, quantile, probs = thres, na.rm = T) # set threshold
    
    # Classify the data into fully censored, non-censored, or partially censored
    N1 = dim(data.u)[1]; n1 = dim(data.u)[2]
    I = matrix(rep(0, N1 * n1), N1, n1) # indicator matrix
    for(i in 1:N1){
      # cond0 = data.u[i, ] <= u.star
      cond1 = data.u[i, ] > u.star
      # I[i, which(cond0)] = 0
      I[i, which(cond1)] = 1
    }
    # browser()
    G1 = which(rowSums(I, na.rm = T) == 0) # fully censored
    G2 = which(rowSums(I, na.rm = T) == n1) # non-censored
    G3 = which(rowSums(I, na.rm = T) > 0 & rowSums(I, na.rm = T) < n1) # partially censored
    
    # Computing fully censored likelihood
    if(length(G1) > 0){
      if(censorL == TRUE){
        value1 = -length(G1) * log.Cn(theta, u.star)
      }
      else{
        value1 = -log.lik(theta, data.u[G1, ])
      }
    }
    # Computing non-censored likelihood
    if(length(G2) > 0){
      value2 = -log.lik(theta, data.u[G2, ])
    }
    # Computing partially censored likelihood
    if(length(G3) > 0){
      if(censorL == TRUE){
        value3 = -log.partialCn(theta, data.u[G3, ], u.star)
      } else {
        value3 = -log.lik(theta, data.u[G3, ])
      }
    }
    # Negative log-likelihood
    return(sum(c(value1, value2, value3))) 
  }
}

##############################
### Model fitting function ###
##############################
# Alteration - actually, this is necessary for the ease of fitting the model
# data [matrix nxD]: all the data already in uniform scale - this is due to the 
#                    complications of having modelled the marginals separately
# theta0 [vector 1x4]: c(lambda0, var1, cov12, var2), initial values
# thres [numeric]: as in model.likelihood
# censorL [logical]: as in model.likelihood
# starting: do you want a parameter grid to be fitted to find best initial values
 

fit.model.likelihood = function(data.u, thres, censorL = TRUE, starting=FALSE){
print(paste0("Start at: ",Sys.time()))
  ## data is already in uniform scale# transform to uniform margins
  
  ## Fit likelihood to a grid of potential starting values to find best initial values for optimisation
  if(starting){
    lams<-c(seq(0.1,5,0.25),seq(5,30,1))
    rs<-seq(-0.9,0.9,0.1)
    init.grid<-expand.grid(lams,rs)
    names(init.grid)<-c("lambda","rho")
    for(i in 1:nrow(init.grid)){
      init.grid$lik[i]<-tryCatch({model.likelihood.optim(theta=unlist(init.grid[i,]),
                                               data.u=data.u,
                                               thres=thres,censorL=T)},
                                 error=function(cond){
                                   return(NA)
                                 })
    }
    init.grid<-init.grid[order(init.grid$lik,decreasing=F),] # order initival values by increasing likelihood
    print("Finished exploration")
  }
  
  tmp<- list()
  
  # Fit likelihood
  for(i in 1:10) {
    if(starting){
      init.vars<-unlist(init.grid[i,c(1,2)])
    }else{
      init.vars<-c(runif(1,min=0.1,max=5),runif(1,min=-0.9,max=0.9))
    }
    # Optimise. If a mistake is found, produce NAs - typically happens if the optimisation does not converge
    tmp[[i]]<-tryCatch({optim(init.vars, model.likelihood.optim,  data.u = data.u, thres = thres,
                      censorL = censorL,
                      method="L-BFGS-B",lower=c(0.1,-0.9),upper=c(Inf,0.9))},
                      error=function(cond){
                        message(cond)
                        return(list(value=NA,par=c(NA,NA)))
                      })
    print(paste0("Optim: ",i))
  }

  # Optimise using the best initial values found in previou step
  try(tmp[[11]]<-optim(tmp[[which.min(unlist(lapply(tmp,function(x) {x$value})))]]$par, 
                   model.likelihood.optim,  data.u = data.u, thres = thres,
                   censorL = censorL,
                   method="L-BFGS-B",lower=c(0.1,-0.9),upper=c(Inf,0.9)))
  to.return<-tryCatch({tmp[[which.min(unlist(lapply(tmp, function(x) {x$value})))]]},
                      error=function(cond){
                        print("All Errors")
                        return(list(value=NA,par=c(NA,NA)))
                      })
  print(paste0("Completed at: ",Sys.time()))
  return(to.return)
}
>>>>>>> 8905faec071b909b4125d27c75dda9d5c3590b9c
