#################################################################################################################################################
#' @name ESM_Ensemble.Modeling
#' @author Flavien Collart \email{flaviencollart@hotmail.com} 
#' @title Ensemble of Small Models: Average Bivariate Models into an ESM
#' @description This function averages simple bivariate models by weighted means to Ensemble Small Models.
#' @param ESM.Mod The object returned by \code{ESM_Modeling}.
#' @param weighting.score \code{character}. An evaluation score used to weight single models to build ensembles:'AUC','MaxTSS', 'SomersD', or '
#' SBI' (if SBI = TRUE in ESM_Modeling) or 'Boyce' (if FALSE).
#' @param threshold \code{numeric} or \code{NULL}. Threshold value of an evaluation score to select the bivariate model(s) included for building 
#' the ensemble. \emph{Default}: 0.5 for AUC and 0 for the other metrics.
#' @param save.obj \code{logical}. Allows or not to save the outputs from this function.
#' @return \itemize{
#' a \code{list} containing: 
#' \item{data}: a \code{list} with the object resp, xy, env.var and sp.name. env.var is = to the data supplied in the argument env. 
#' If env, was a SpatRaster, it corresponds to the extracted values of these rasters for the species coordinates.
#' \item{model.info}: a \code{list} of the models used, their options, the combination of bivariate models (which.biva), the modeling ID, 
#' the path to the folder where are the stored the models (biva.path), and the failed models
#' \item{cv.split.table}: a \code{matrix} used to train and test models. See explanation of the argument cv.split.table
#' \item{pooling}: a \code{logical}. Does the pooling method to evaluate the models was performed?
#' \item{evaluations}: a \code{matrix} The evaluation of the ensemble based on 4 metrics: the AUC, the Somer's D (=2*AUC-1),
#' maxTSS, and the smooth Boyce Index (SBI, if SBI = TRUE in ESM_Modeling) or the regular Boyce Index (if SBI = FALSE).
#' \item{EF.algo}: a \code{list} containing \code{pred.EF.algo} which is a \code{matrix} of ensemble predictions for each run 
#' and modeling techniques; and \code{weights.algo} a \code{matrix} of weights used to generate the ensemble at the level of the algorithm.
#' \item{EF}: a \code{list} containing \code{pred.EF} which is a \code{matrix} of ESMs predictions for each run; 
#' and \code{weights.EF} a \code{matrix} of weights used to generate the ESMs. Note that if only one modelling technique is used 
#' \code{EF} will be exactly the same as \code{EF.algo}.
#' }
#' @details
#' If the pooling evaluation was selected in \code{ESM_Modeling}, the different weights will result from 
#' this evaluation. In addition, the final ensemble will also be evaluated using the pooling method (see Collart & Guisan 2023).
#' 
#' For the use of this function, please refer to the manual of ESM_Modeling.
#' 
#' @seealso \code{\link{ESM_Modeling}}, \code{\link{ESM_Projection}} and  \code{\link{ESM_Ensemble.Projection}}
#' @references 
#' Lomba, A., L. Pellissier, C.F. Randin, J. Vicente, F. Moreira, J. Honrado and A. Guisan. 2010. Overcoming the rare species 
#' modelling paradox: A novel hierarchical framework applied to an Iberian endemic plant. \emph{Biological Conservation}, \bold{143},2647-2657.
#' 
#' Breiner F.T., A. Guisan, A. Bergamini and M.P. Nobis. 2015. Overcoming limitations of modelling rare species by using ensembles of small models. \emph{Methods in Ecology and Evolution}, \bold{6},1210-1218.
#' 
#' Breiner F.T., Nobis M.P., Bergamini A., Guisan A. 2018. Optimizing ensembles of small models for predicting the distribution of species with few occurrences. \emph{Methods in Ecology and Evolution}. \doi{10.1111/2041-210X.12957}
#' 
#' Collart, F., & Guisan, A. (2023). Small to train, small to test: Dealing with low sample size in model evaluation. \emph{Ecological Informatics}. \bold{75}, 102106. \doi{10.1016/j.ecoinf.2023.102106}.
#' 
#' @export

## ESM_Ensemble.Modeling----
ESM_Ensemble.Modeling <- function(ESM.Mod,
                                 weighting.score,
                                 threshold = NULL,
                                 save.obj = TRUE){
  
  ## Check some arguments----
  SBI <- ESM.Mod$data$SBI
  if(SBI){
    if(!weighting.score %in% c("AUC", "MaxTSS", "SBI", 
                               "SomersD")) {
      stop("weighting score not supported. Choose one of the following: AUC, MaxTSS, SBI, or SomersD")
    }
  }else{
    if(!weighting.score %in% c("AUC", "MaxTSS", "Boyce", 
                               "SomersD")) {
      stop("weighting score not supported. Choose one of the following: AUC, MaxTSS, Boyce, or SomersD")
    }
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
  
  ### First generate the first level of aggregation (on algo)----
  
  models <- ESM.Mod$model.info$models
  biva.eval <- ESM.Mod$biva.evaluations
  biva.pred <- do.call(cbind,ESM.Mod$biva.predictions)
  failed.mods <- ESM.Mod$model.info$failed.mod
  cv.split.table <-ESM.Mod$cv.split.table
  run.names <- colnames(cv.split.table)
  resp <- ESM.Mod$data$resp
  pooling <- ESM.Mod$model.info$pooling
  
  for(i in 1:length(models)){
    ##Get the weights of the full model (which is = to mean across the cv)----
    w <- as.data.frame(lapply(biva.eval,FUN = function(x,model,ws){
      x[paste0("Full.",model),ws]
    },model = models[i],ws=weighting.score))
    row.names(w) = models[i]
    w[w<threshold | is.na(w)] = 0
    if(length(w)==sum(w==0)){
      warning(paste("All weights for models based on",models[i],"are null consider changing the threshold."))
    }
    
    if(i == 1){
      weights.algo <- w
    }else{
      weights.algo <- rbind(weights.algo,w)
    }

    ## Then make the ensemble----
    pred.EF.algo <- sapply(run.names,.bivaToEnsemble,model = models[i],
                           w = w,biva.pred = biva.pred) 
    pred.EF.algo <- do.call(cbind.data.frame,pred.EF.algo)
    if(i==1){
      pred.EF <- pred.EF.algo
    }else{
      pred.EF <- cbind.data.frame(pred.EF,pred.EF.algo)
      
    }
    
  }
  if(pooling){
    EF.algo.eval <- .pooling.ESM.Mod(pred.EF,
                                    resp=resp, 
                                    models=models,
                                    cv.split.table=cv.split.table[,-ncol(cv.split.table)],
                                    SBI = SBI)
  }else{
    ## Evaluate the ensembles per Run + extract w.scores----
    EF.algo.eval <- .bivaEvaluation(biva = pred.EF,
                                    resp=resp, models=models,
                                    cv.split.table=cv.split.table,
                                    SBI = SBI)
    
  }
  rownames(EF.algo.eval) = paste0(rownames(EF.algo.eval),".EF")
 

  ### Make the ensemble across the modeling techniques----
  if(length(models)>1){
    ## Take the weights for each modeling technique
    weights.EF <- EF.algo.eval[grep("Full",row.names(EF.algo.eval)),weighting.score]
    weights.EF[weights.EF<threshold | is.na(weights.EF)] = 0
    names(weights.EF) = sub("Full.","",names(weights.EF),fixed = TRUE)
    EF <- sapply(run.names,.EFToEnsemble,
                           w = weights.EF,pred.EF = pred.EF)
    colnames(EF) = paste0(colnames(EF),".EF")
    
    ##Evaluate the model for each run----
    if(pooling){
      EF.eval <- .pooling.ESM.Mod(EF,
                                 resp=resp, models="EF",
                                 cv.split.table=cv.split.table[,-ncol(cv.split.table)],
                                 SBI = SBI)
    }else{
      EF.eval <- .bivaEvaluation(biva = EF,
                                 resp=resp, models="EF",
                                 cv.split.table=cv.split.table,
                                 SBI = SBI)
    }
    
    EF.eval <- rbind(EF.algo.eval,EF.eval)
    
  }else{
    EF <- pred.EF
    EF.eval <- EF.algo.eval
    weights.EF <- weights.algo
  }
  
  ### Save the outputs----
  obj <- list(data=ESM.Mod$data,
              model.info = ESM.Mod$model.info,
              cv.split.table = cv.split.table,
              pooling = pooling,
              evaluations = EF.eval,
              weighting.score = weighting.score,
              threshold = threshold,
              EF.algo = list(pred.EF.algo = pred.EF,
                             weights.algo = weights.algo),
              EF = list(pred.EF = EF,
                        weights.EF = weights.EF))
  if(save.obj){
    save(obj,file= paste0("../ESM.EnsembleModeling.",ESM.Mod$model.info$modeling.id,".out"))
  }
  return(obj)
}

### The hidden functions ----
##.bivaToEnsemble----
## Allows to ensemble each bivariate model from a single algo
.bivaToEnsemble <- function(x,model,w,biva.pred){
  
    pred.algo <- biva.pred[,grep(paste0(x,".",model),colnames(biva.pred),fixed = T)]
    pred.EF.run <- as.data.frame(apply(pred.algo,1,stats::weighted.mean,
                                       w=w[sub(paste0(".",x,".",model), "",colnames(pred.algo),fixed = T)],na.rm=T))
    ## For the w, I make sure that it is in the same order as pred.algo
    colnames(pred.EF.run) = model
    return(pred.EF.run)
}
##.EFToEnsemble----
## Allows to ensemble each algo 
.EFToEnsemble <- function(x,pred.EF,w){
  pred.EF.algo <- pred.EF[,grep(paste0(x,"."),colnames(pred.EF),fixed = T)]
  EF <- apply(pred.EF.algo,1,stats::weighted.mean,w=w)
  return(EF)
}
