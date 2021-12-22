predict.mnl <- function(model, data) {
  # Function for predicting preference shares from a MNL model 
  # model: mlogit object returned by mlogit()
  # data: a data frame containing the set of designs for which you want to 
  #       predict shares.  Same format at the data used to estimate model. 
  data.model <- model.matrix(update(model$formula, 0 ~ .), data = data)[,-1]
  logitUtility <- data.model%*%model$coef
  share <- exp(logitUtility)/sum(exp(logitUtility))
  cbind(share, data)
}




sensitivity.mnl <- function(model, attrib, base.data, competitor.data) {
  # Function for creating data for a preference share-sensitivity chart
  # model: mlogit object returned by mlogit() function
  # attrib: list of vectors with attribute levels to be used in sensitivity
  # base.data: data frame containing baseline design of target product
  # competitor.data: data frame containing design of competitive set
  data <- rbind(base.data, competitor.data)
  base.share <- predict.mnl(model, data)[1,1]
  share <- NULL
  for (a in seq_along(attrib)) {
    for (i in attrib[[a]]) {
      data[1,] <- base.data
      data[1,a] <- i
      share <- c(share, predict.mnl(model, data)[1,1])
    }
  }
  data.frame(level=unlist(attrib), share=share, increase=share-base.share)
}


predict.mixed.mnl <- function(model, data, nresp=1000) {
  # Function for predicting shares from a mixed MNL model 
  # model: mlogit object returned by mlogit()
  # data: a data frame containing the set of designs for which you want to 
  #       predict shares. Same format at the data used to estimate model. 
  # Note that this code assumes all model parameters are random
  data.model <- model.matrix(update(model$formula, 0 ~ .), data = data)[,-1]
  coef.Sigma <- cov.mlogit(model)
  coef.mu <- model$coef[1:dim(coef.Sigma)[1]]
  draws <- mvrnorm(n=nresp, coef.mu, coef.Sigma)
  shares <- matrix(NA, nrow=nresp, ncol=nrow(data))
  for (i in 1:nresp) {
    utility <- data.model%*%draws[i,]
    share = exp(utility)/sum(exp(utility))
    shares[i,] <- share
  }
  cbind(colMeans(shares), data)
}


BootCI.predict.mnl <- function(model, data, nsim=500, conflevel=0.95) {
  dataModel <- model$model
  dataModel$probabilities <- NULL
  dataModel$linpred <- NULL
  idx <- dataModel$idx 
  dataModel$idx <- NULL
  dataModel <- data.frame(dataModel, idx)
  idVar <- unique(dataModel[,names(idx)[1]])
  
  bootstrapping <- function(x) {
    idbootsamp <- data.frame(sample(idVar, replace=T))
    names(idbootsamp) <- names(idx)[1]
    bootsamp <- merge(idbootsamp, dataModel, by=names(idx)[1], all.x=T)
    bootsamp[,names(idx)[1]] <- rep(1:length(table(idx[,1])), each=length(table(idx[,3])))
    bootsamp.mlogit  <- dfidx(bootsamp, idx = list(c(names(idx)[1:2]), names(idx)[3]),
                              drop.index=F)    
    bootfit <- update(model, data = bootsamp.mlogit)
    data.model <- model.matrix(update(bootfit$formula, 0 ~ .), data = data)[,-1]
    logitUtility <- data.model%*%bootfit$coef
    share <- exp(logitUtility)/sum(exp(logitUtility))
    share
  }
  
  cl <- makeCluster(detectCores())
  clusterEvalQ(cl, library(mlogit))
  clusterExport(cl, varlist=c("idVar", "dataModel", "idx", "model", "data"), 
                envir=environment())
  bootdistr <- parLapply(cl, 1:nsim, fun=bootstrapping)
  stopCluster(cl)
  
  bootdistr <- do.call(cbind, bootdistr)
  lowl <- (1-conflevel)/2
  upl <- 1-lowl  
  bootperc <- t(apply(bootdistr, 1, function(x) quantile(x, probs=c(lowl, upl))))
  pointpred <- predict.mnl(model, data)
  predictedShares <- cbind(pointpred[,1], bootperc, pointpred[,2:ncol(pointpred)])
  names(predictedShares)[1] <- "share" 
  predictedShares
}





BootCI.predict.mixed.mnl <- function(model, data, nsim=500, conflevel=0.95, nresp=1000) {
  dataModel <- model$model
  dataModel$probabilities <- NULL
  dataModel$linpred <- NULL
  idx <- dataModel$idx 
  dataModel$idx <- NULL
  dataModel <- data.frame(dataModel, idx)
  idVar <- unique(dataModel[,names(idx)[1]])
  
  bootstrapping <- function(x) {
    idbootsamp <- data.frame(sample(idVar, replace=T))
    names(idbootsamp) <- names(idx)[1]
    bootsamp <- merge(idbootsamp, dataModel, by=names(idx)[1], all.x=T)
    bootsamp[,names(idx)[1]] <- rep(1:length(table(idx[,1])), each=length(table(idx[,3])))
    bootsamp.mlogit  <- dfidx(bootsamp, idx = list(c(names(idx)[1:2]), names(idx)[3]),
                              drop.index=F)    
    bootfit <- update(model, data = bootsamp.mlogit)
    data.model <- model.matrix(update(bootfit$formula, 0 ~ .), data = data)[,-1]
    coef.Sigma <- cov.mlogit(bootfit)
    coef.mu <- bootfit$coef[1:dim(coef.Sigma)[1]]
    draws <- mvrnorm(n=nresp, coef.mu, coef.Sigma)
    shares <- matrix(NA, nrow=nresp, ncol=nrow(data))
    for (i in 1:nresp) {
      utility <- data.model%*%draws[i,]
      share <- exp(utility)/sum(exp(utility))
      shares[i,] <- share
    }
    colMeans(shares)
  }
  
  cl <- makeCluster(detectCores())
  clusterEvalQ(cl, {
    library(mlogit)
    library(MASS) })
  clusterExport(cl, varlist=c("idVar", "dataModel", "idx", "model", "data", "nresp", 
                              paste(model$call$rpar)), envir=environment())
  bootdistr <- parLapply(cl, 1:nsim, fun=bootstrapping)
  stopCluster(cl)
  
  bootdistr <- do.call(cbind, bootdistr)
  lowl <- (1-conflevel)/2
  upl <- 1-lowl  
  bootperc <- t(apply(bootdistr, 1, function(x) quantile(x, probs=c(lowl, upl))))
  pointpred <- predict.mixed.mnl(model, data, nresp)
  predictedShares <- cbind(pointpred[,1], bootperc, pointpred[,2:ncol(pointpred)])
  names(predictedShares)[1] <- "share" 
  predictedShares
}