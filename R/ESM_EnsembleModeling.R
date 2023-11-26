ESM_EnsembleModeling <- function(ESM.Mod,
                                 weighting.score,
                                 threshold = NULL,
                                 save.obj = TRUE){
  if (!weighting.score %in% c("AUC", "MaxTSS", "Boyce", 
                              "SomersD")) {
    stop("weighting score not supported! Choose one of the following: AUC, MaxTSS, Boyce, or SomersD")
  }
  if(is.null(threshold)){
    if(weighting.score == "AUC"){
      threshold = 0.5
    }else{
      threshold = 0
    }
  }
  iniwd <- getwd()
  on.exit(setwd(iniwd))
  newwd <- ESM.Mod$biva.path
  setwd(newwd)
  
  ### First generate the first level of aggregation (on algo)
  models <- ESM.Mod$models
  biva.eval <- ESM.Mod$biva.evaluations
  biva.pred <- ESM.Mod$biva.predictions
  run.names <- colnames(ESM.Mod$cv.split.table)
  resp <- ESM.Mod$resp
  cv.split.table <-ESM.Mod$cv.split.table
  
  for(i in 1:length(models)){
    ##Get the weights of the full model (which is = to mean across the cv)
    w <- as.data.frame(lapply(biva.eval,FUN = function(x,model,ws){
      x[paste0("Full.",model),ws]
    },model = models[i],ws=weighting.score))
    row.names(w) = models[i]
    w[w<threshold] = 0
    if(i == 1){
      weights.algo <- w
    }else{
      weights.algo <- rbind(weights.algo,w)
    }
    ## Then make the ensemble
    pred.EF.algo <- sapply(run.names,.bivaToEnsemble,model = models[i],
                           w = w,biva.pred = biva.pred)
    pred.EF.algo <- do.call(cbind.data.frame,pred.EF.algo)
    if(i==1){
      pred.EF <- pred.EF.algo
    }else{
      pred.EF <- cbind.data.frame(pred.EF,pred.EF.algo)
      
    }
    
  }
  
  ## Evaluate the ensembles per Run + extract w.scores
  EF.algo.eval <- .bivaEvaluation(biva = pred.EF,
                                  resp=resp, models=models,
                                  cv.split.table=cv.split.table)
  rownames(EF.algo.eval) = paste0(rownames(EF.algo.eval),".EF")
  ## Make the full ensemble
  if(length(models)>1){
    weights.EF <- EF.algo.eval[grep("Full",row.names(EF.algo.eval)),weighting.score]
    weights.EF[weights.EF<threshold] = 0
    names(weights.EF) = sub("Full.","",names(weights.EF),fixed = TRUE)
    EF <- sapply(run.names,.EFToEnsemble,
                           w = weights.EF,pred.EF = pred.EF)
    colnames(EF) = paste0(colnames(EF),".EF")
    EF.eval <- .bivaEvaluation(biva = EF,
                                    resp=resp, models="EF",
                                    cv.split.table=cv.split.table)
    EF.eval <- rbind(EF.algo.eval,EF.eval)
    
  }else{
    EF <- pred.EF
    EF.eval <- EF.algo.eval
  }
  
  obj <- list(resp = resp,
              sp.name = my.ESM$sp.name,
              cv.split.table = cv.split.table,
              models = models,
              evaluations = EF.eval,
              pred.EF.algo = pred.EF,
              weights.algo = weights.algo,
              EF = EF,
              weights.EF = weights.EF,
              modeling.id = modeling.id,
              biva.path = newwd)
  if(save.obj){
    save(obj,file= paste0("../ESM.EnsembleModeling.",modeling.id,".out"))
  }
  return(obj)
}

.bivaToEnsemble <- function(x,model,w,biva.pred){
  
    pred.algo <- as.data.frame(lapply(biva.pred,FUN = function(y,model,x){
      y[,paste0(x,".",model)]
    },model = model,x=x))
    pred.EF.run <- as.data.frame(apply(pred.algo,1,weighted.mean,
                                       w=w[colnames(pred.algo)],na.rm=T))
    colnames(pred.EF.run) = model
    return(pred.EF.run)
}
.EFToEnsemble <- function(x,pred.EF,w){
  pred.EF.algo <- pred.EF[,grep(x,colnames(pred.EF))]
  EF <- apply(pred.EF.algo,1,weighted.mean,w=w)
  return(EF)
}
