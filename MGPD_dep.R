### ======================================================================== ###
### Fit MGPD model for extraml dependence of unreplicated data using the     ###
### Kiriliouk (2019) construction .                                          ###
### Daniela Cuba                                                             ###                                                      ###
### ======================================================================== ###

# source("./ModelDiagnosticsNewNames.r")
source("./CommonFunctions.r")
source("./code/RevExp_U_Functions.r")

# Parameters
# x_data: bivariate dataset with standard exponential margins.Set of at least one exceedance. 

MGPD.dep<-function(x_data){
  toreturn_models<-list()
  print(paste0("Start at: ",Sys.time()))
  
  ## Fit 1 - Free scale, free loc
  fit1 <- lapply(1:10,function(x){
    fit.tmp <-NULL
    k<-1
    while(is.null(fit.tmp)){
      init.vals <- c(runif(2,min=0,max=10),
                     runif(1,min=0,max=2))
      # browser()
      fit.tmp<-tryCatch(expr={fit.MGPD.RevExpU(x=x_data, 
                                               u=rep(0,2), std=T, 
                                               dep.scale.fix=F,
                                               dep.loc.fix=FALSE, 
                                               dep.start=init.vals)},
                        error=function(cond){
                          return(NULL)}
      )
      print(k)
      k<-k+1
    }
    print(x)
    return(fit.tmp)
  })
  
  initvars<- unlist(lapply(fit1, function(x) {x$nll})) %>% which.min
  
  toreturn_models[[1]] <- fit.MGPD.RevExpU(x=x_data, u=rep(0,2), std=T, 
                                           dep.scale.fix=F,
                                           dep.loc.fix=F, 
                                           dep.start=fit1[[initvars]]$mle)
  print(paste0("Fit 1 complete: ",Sys.time()))
  ## Fit 2 - Free scale, fixed loc
  fit2 <- lapply(1:10,function(x){
    fit.tmp <-NULL
    k<-1
    while(is.null(fit.tmp)){
      init.vals <- c(runif(2,min=0,max=10))
      # browser()
      fit.tmp<-tryCatch(expr={fit.MGPD.RevExpU(x=x_data, u=rep(0,2), std=T, dep.scale.fix=F,
                                               dep.loc.fix=T, dep.start=init.vals)},
                        error=function(cond){
                          return(NULL)}
      )
      print(k)
      k<-k+1
    }
    print(x)
    return(fit.tmp)
  })
  
  initvars<- unlist(lapply(fit2, function(x) {x$nll})) %>% which.min
  
  toreturn_models[[2]] <- fit.MGPD.RevExpU(x=x_data, u=rep(0,2), std=T, 
                                           dep.scale.fix=F,
                                           dep.loc.fix=T, dep.start=fit2[[initvars]]$mle)
  print(paste0("Fit 2 complete: ",Sys.time()))
  
  ## Fit 3 - Fixed scale, free location
  fit3 <- lapply(1:10,function(x){
    fit.tmp <-NULL
    k<-1
    while(is.null(fit.tmp)){
      init.vals <- c(runif(1,min=0,max=10),
                     runif(1,min=0,max=2))
      # browser()
      fit.tmp<-tryCatch(expr={fit.MGPD.RevExpU(x=x_data, u=rep(0,2), std=T, dep.scale.fix=T,
                                               dep.loc.fix=F, dep.start=init.vals)},
                        error=function(cond){
                          return(NULL)}
      )
      print(k)
      k<-k+1
    }
    print(x)
    return(fit.tmp)
  })
  
  initvars<- unlist(lapply(fit3, function(x) {x$nll})) %>% which.min
  
  toreturn_models[[3]] <- fit.MGPD.RevExpU(x=x_data, u=rep(0,2), std=T, 
                                           dep.scale.fix=T,
                                           dep.loc.fix=FALSE, 
                                           dep.start=fit3[[initvars]]$mle)
  print(paste0("Fit 3 complete: ",Sys.time()))
  
  ## Fit 4 - Fixed scale, fixed location
  fit4 <- lapply(1:10,function(x){
    fit.tmp <-NULL
    k<-1
    while(is.null(fit.tmp)){
      init.vals <- c(runif(1,min=0,max=10))
      # browser()
      fit.tmp<-tryCatch(expr={fit.MGPD.RevExpU(x=x_data, u=rep(0,2), std=T, dep.scale.fix=T,
                                               dep.loc.fix=T, dep.start=init.vals)},
                        error=function(cond){
                          return(NULL)}
      )
      print(k)
      k<-k+1
    }
    print(x)
    return(fit.tmp)
  })
  
  initvars<- unlist(lapply(fit4, function(x) {x$nll})) %>% which.min
  
  toreturn_models[[4]] <- fit.MGPD.RevExpU(x=x_data, u=rep(0,2), std=T, 
                                           dep.scale.fix=T,
                                           dep.loc.fix=T, dep.start=fit3[[initvars]]$mle)
  print(paste0("Fit 4 complete: ",Sys.time()))
  return(toreturn_models)
}