### ======================================================================== ###
### Fit Exponential Factor Copula Model to unreplicated dataset              ###
### Daniela Cuba                                                             ###                                                      ###
### ======================================================================== ###

## Load libraries and functions ##
source('./CensoredLikelihood_nonspatial.R')

## Load data ##
pair<-"PbCr"
load("./pairs_uniform.Rdata")

## Fit model to pair##
ExpFac_models<- lapply(1:length(pair), function(x){
  tmp<-tryCatch({fit.model.likelihood(pairs_uniform[[x]],thres=0.82,censorL=T,starting=T)},
                error=function(cond){
                  return(list(par=c(NA,NA),value=NA))
                })
  })

## Save output ##
# save(ExpFac_models,file="./ExpFac_models.RData")
