#################################################################################################################################################
## ESM_Ensemble.Modeling
#' @export
##  Description:
##    Ensemble all the bivariate models together for each modeling techniques and then
##    make an ensemble between techniques based on a weighted mean. The ensembles are
##    evaluated.
##
##  Arguments:
## @ESM.Mod: The outputs resulted from ESM_Modeling()
## @weighting.score: character. an evaluation score used to weight single models to build ensembles:"AUC","MaxTSS","Boyce" or,"SomersD"
## @threshold: numeric or NULL(Default). Threshold value of an evaluation score to select the bivariate model(s) included for building the ensemble
## @save.obj: logical. Allows or not to save each run and the final output object
##
## Values: 
##        a list containing: 
##                          data: a list with the object resp, xy, env.var and sp.name. env.var is = to the data suuplied in the argument env. 
##                                If env, was a SpatRaster, it corrsponds to the extracted values of these rasters for the species coordinates.
##                          model.info: contains the models used, their options, the combination of bivariate models (which.biva), the
##                                      modeling ID, the path to the folder where are the stored the models (biva.path), and the failed models
##                          cv.split.table: a table used to train and test models (see explanation of the argument cv.split.table)
##                          evaluations: a data.frame containing the evaluations of each run and ensemble models. The full model is evaluated
##                          using the mean values across the runs.
##                          EF.algo: a list containing: 
##                                          pred.EF.algo: a data.frame containing the predicted values for each run and the full models
##                                          for the ensemble of bivariate models across each modeling technique.
##                                          weights.algo: The weights used for each bivariate model to make the ensemble for each modeling technique
##                          EF: a list containing: (if there is only one modeling technique, EF = EF.algo)
##                                          pred.EF: a data.frame containing the predicted values for the final ESM 
##                                          weights.EF: The weights used for each modeling algorithm to make the final ESM.
##  Authors:
##          Flavien Collart based on the previous code written by Frank Breiner with the contributions of 
##          Olivier Broennimann and Flavien Collart
#################################################################################################################################################


ESM_Ensemble.Modeling <- function(ESM.Mod,
                                 weighting.score,
                                 threshold = NULL,
                                 save.obj = TRUE){
  
  ## Check some arguments
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
  
  
  ## Set the pathway
  iniwd <- getwd()
  on.exit(setwd(iniwd))
  newwd <- ESM.Mod$model.info$biva.path
  setwd(newwd)
  
  ###############
  ### First generate the first level of aggregation (on algo)
  
  models <- ESM.Mod$model.info$models
  biva.eval <- ESM.Mod$biva.evaluations
  biva.pred <- do.call(cbind,ESM.Mod$biva.predictions)
  cv.split.table <-ESM.Mod$cv.split.table
  run.names <- colnames(cv.split.table)
  resp <- ESM.Mod$data$resp
  
  for(i in 1:length(models)){
    ##Get the weights of the full model (which is = to mean across the cv)
    w <- as.data.frame(lapply(biva.eval,FUN = function(x,model,ws){
      x[paste0("Full.",model),ws]
    },model = models[i],ws=weighting.score))
    row.names(w) = models[i]
    w[w<threshold | is.na(w)] = 0
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
  #########################
  
  ### Make the ensemble across the modeling techniques. 
  if(length(models)>1){
    ## Take the weights for each modeling technique
    weights.EF <- EF.algo.eval[grep("Full",row.names(EF.algo.eval)),weighting.score]
    weights.EF[weights.EF<threshold | is.na(weights.EF)] = 0
    names(weights.EF) = sub("Full.","",names(weights.EF),fixed = TRUE)
    EF <- sapply(run.names,.EFToEnsemble,
                           w = weights.EF,pred.EF = pred.EF)
    colnames(EF) = paste0(colnames(EF),".EF")
    
    ##Evaluate the model for each run
    EF.eval <- .bivaEvaluation(biva = EF,
                                    resp=resp, models="EF",
                                    cv.split.table=cv.split.table)
    EF.eval <- rbind(EF.algo.eval,EF.eval)
    
  }else{
    EF <- pred.EF
    EF.eval <- EF.algo.eval
    weights.EF <- weights.algo
  }
  
  ### Save the outputs.
  obj <- list(data=ESM.Mod$data,
              model.info = ESM.Mod$model.info,
              cv.split.table = cv.split.table,
              evaluations = EF.eval,
              EF.algo = list(pred.EF.algo = pred.EF,
                             weights.algo = weights.algo),
              EF = list(pred.EF = EF,
                        weights.EF = weights.EF))
  if(save.obj){
    save(obj,file= paste0("../ESM.EnsembleModeling.",ESM.Mod$model.info$modeling.id,".out"))
  }
  return(obj)
}

## Allows to ensemble each bivariate model from a single algo
.bivaToEnsemble <- function(x,model,w,biva.pred){
  
    pred.algo <- biva.pred[,grep(paste0(x,".",model),colnames(biva.pred),fixed = T)]
    pred.EF.run <- as.data.frame(apply(pred.algo,1,weighted.mean,
                                       w=w[sub(paste0(".",x,".",model), "",colnames(pred.algo),fixed = T)],na.rm=T))
    ## For the w, I make sure that it is in the same order as pred.algo
    colnames(pred.EF.run) = model
    return(pred.EF.run)
}

## Allows to ensemble each algo 
.EFToEnsemble <- function(x,pred.EF,w){
  pred.EF.algo <- pred.EF[,grep(paste0(x,"."),colnames(pred.EF),fixed = T)]
  EF <- apply(pred.EF.algo,1,weighted.mean,w=w)
  return(EF)
}
