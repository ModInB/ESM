
#' @name ESM_Projection
#' @author Flavien Collart \email{flaviencollart@@hotmail.com} 
#' @title Ensemble of Small Models: Projections of Bivariate Models
#' @description Project bivariate models 
#' @param ESM.Mod The object returned by \code{\link{ESM_Modeling}}.
#' @param new.env a \code{data.frame}, \code{matrix} or a \code{SpatRaster} of the predictors where the models 
#' should be projected. Note that the colnames or the names of the predictors should be exactly the same as the one used in \code{\link{ESM_Modeling}}.
#' @param name.env a \code{character} which will be used to generate a folder to save the individual projections.
#' @param rounded \code{logical}. Should the prediction be rounded to become an integer? \emph{Default:} TRUE. 
#' @param pred.multiplier \code{numeric}. The factor to multiply the predictions before rounding them (must be a power of 10). 
#' Only needed when rounded = FALSE. \emph{Default:} 1000. 
#' @param parallel \code{logical}. Allows or not parallel job using the functions makeCluster.
#' @param n.cores \code{integer}. Number of CPU cores used to make the models.
#' @param datatype a \code{character} When new.env is a  \code{SpatRaster}, which datatype should be used to save the files? \emph{Default:} 'INT2U' if rounded = TRUE 
#' and 'FLT4S' if rounded = FALSE.
#' @param save.obj \code{logical}. Allows or not to save the final output into a '.out file'.
#' @param verbose \code{logical}. Allows or not message.
#' 
#' @return  a \code{list} containing: 
#' \itemize{
#' \item{projection.path}: a \code{character} the path to each projection.
#' \item{model.info}: a \code{list} of the models used, their options, the combination of bivariate models (which.biva),
#' the modeling ID, the path to the folder where are the stored the models (biva.path), and the failed models.
#' \item{name.env}: a \code{character} the value given to the argument name.env.
#' \item{proj.type}: a \code{character} which is either "data.frame" or "SpatRaster" depending on the class of new.env
#' }
#' @details  
#' \describe{
#' Each bivariate models are geographically projected into the new.env area. 
#' Note that if new.env is a SpatRaster, the projected maps are by default multiplied by 1000 and 
#' rounded to reduce storage space. They are also store in compress tif file (with these parameters: "COMPRESS=DEFLATE","PREDICTOR=2").
#' If the argument 'rounded' is set to FALSE and datatype is left default, the datatype will be changed to 'FLT4S' to keep decimals.
#' For the use of this function, please refer to the manual of ESM_Modeling.
#' }
#' @seealso \code{\link{ESM_Modeling}}, \code{\link{ESM_Ensemble.Modeling}}, \code{\link{ESM_Ensemble.Projection}}
#' @references Lomba, A., L. Pellissier, C.F. Randin, J. Vicente, F. Moreira, J. Honrado and A. Guisan. 2010. Overcoming the rare species 
#' modelling paradox: A novel hierarchical framework applied to an Iberian endemic plant. \emph{Biological Conservation}, \bold{143},2647-2657.
#' 
#' Breiner F.T., A. Guisan, A. Bergamini and M.P. Nobis. 2015. Overcoming limitations of modelling rare species by using ensembles of small models. \emph{Methods in Ecology and Evolution}, \bold{6},1210-1218.
#' 
#' Breiner F.T., Nobis M.P., Bergamini A., Guisan A. 2018. Optimizing ensembles of small models for predicting the distribution of species with few occurrences. \emph{Methods in Ecology and Evolution}. \doi{10.1111/2041-210X.12957}
#' 
#' @importFrom stats predict.glm predict
#' @importFrom utils read.table write.table
#' @export

ESM_Projection <- function(ESM.Mod,
                           new.env,
                           name.env,
                           rounded = TRUE,
                           pred.multiplier = 1000,
                           parallel = FALSE,
                           n.cores = 1,
                           datatype = "INT2U",
                           save.obj = T,
                           verbose = TRUE){
  iniwd <- getwd()
  on.exit(setwd(iniwd))
  setwd(ESM.Mod$model.info$biva.path)
  
  models <- ESM.Mod$model.info$models
  which.biva <- ESM.Mod$model.info$which.biva
  env.var <- ESM.Mod$data$env.var
  biva.pred = ESM.Mod$biva.predictions
  
  
  #### check new.env----
  if(is.matrix(new.env)){
    new.env  <- as.data.frame(new.env)
  }
  if(!(inherits(new.env,"data.frame")) & !(inherits(new.env,"SpatRaster"))){
    stop("new.env must be a dataframe, a matrix or a SpatRaster")
  }
  if(is.null(name.env)){
    stop("name.env cannot be null")
  }
  
  ## check other arguments ----
  if(!is.logical(rounded)){
    stop("rounded must be a logical")
  }
  
  if(!is.numeric(pred.multiplier) | pred.multiplier < 1){
    stop("pred.multiplier must be a positive numeric.")
  }
  if(round(log10(pred.multiplier)) != log10(pred.multiplier)){
    stop("pred.multiplier must be a power of 10.")
    
  }
  if(!rounded){
    pred.multiplier = 1
    if(datatype == 'INT2U'){
      datatype = 'FLT4S'
    }
  }
  
  dir.create(paste0("../",name.env))
  ## Compute the realized combinations----
  combinations <- utils::combn(colnames(env.var), 2)
  if (is.null(which.biva)) {
    which.biva <- 1:ncol(combinations)
  }
  combinations <- combinations[,which.biva]
  used.env <- unique(c(combinations[1,],combinations[2,]))

  if(is.data.frame(new.env)){
    if(sum(colnames(new.env) %in% used.env) != length(used.env) ){
      stop("new.env need to have the same variable names as used in ESM_Modeling")
    }
    proj.type = "data.frame"
    if(parallel){
      cl <- parallel::makeCluster(n.cores)
      proj <- parallel::parApply(cl,combinations, 2, .IndividualProj,
                                 new.env = new.env,models = models,
                                 name.env = name.env, biva.pred = biva.pred,
                                 datatype = datatype,
                                 rounded = rounded,
                                 pred.multiplier = pred.multiplier,
                                 verbose = verbose)
      parallel::stopCluster(cl)
      
    }else{
      proj <- apply(combinations, 2, .IndividualProj,
                    new.env = new.env,models = models,
                    name.env = name.env, biva.pred = biva.pred,
                    datatype = datatype,
                    rounded = rounded,
                    pred.multiplier = pred.multiplier,
                    verbose = verbose)
    }
   
  }else{
    if(sum(names(new.env) %in% used.env) != length(used.env) ){
      stop("new.env need to have the same variable names as used in ESM_Modeling")
    }
    proj.type = "SpatRaster"
    if(parallel){
      new.env <- terra::wrap(new.env) ## To allow a correct parallelisation
      cl <- parallel::makeCluster(n.cores)
      proj <- parallel::parApply(cl,combinations, 2, .IndividualProj,
                                 new.env = new.env,models = models,
                                 name.env = name.env,parallel = parallel, 
                                 biva.pred = biva.pred,
                                 datatype = datatype,
                                 rounded = rounded,
                                 pred.multiplier = pred.multiplier,
                                 verbose = verbose)
      parallel::stopCluster(cl)
      
    }else{
      proj <- apply(combinations, 2, .IndividualProj,
                    new.env = new.env,models = models,
                    name.env = name.env,parallel = parallel, 
                    biva.pred = biva.pred,
                    datatype = datatype,
                    rounded = rounded,
                    pred.multiplier = pred.multiplier,
                    verbose = verbose)
    }
  }
  
  obj = list(projection.path=as.list(proj),
             model.info = ESM.Mod$model.info,
             name.env = name.env,
             proj.type=proj.type,
             datatype = datatype,
             rounded = rounded,
             pred.multiplier = pred.multiplier)
  if(save.obj){
    save(obj,file=paste0("../ESM.Projection.",name.env,".",ESM.Mod$model.info$modeling.id,".out"))
  }
  return(obj)
}

## The hidden functions ----
## .IndividualProj----
## Function projecting a model onto a new environment
#' @import nnet
#' @import rpart
#' @import maxnet
.IndividualProj <- function(x,
                            new.env,
                            models,
                            name.env,
                            parallel, 
                            biva.pred,
                            datatype,
                            rounded,
                            pred.multiplier,
                            verbose){
  
  if(!is.data.frame(new.env) & verbose){
    message(paste("\n Projections of bivariate model:",x[1],x[2]))}
  done <- c()
  for(j in 1:length(models)){
    ToSkip <- anyNA(biva.pred[[paste(x[1],x[2], sep=".")]][,paste0("Full.",models[j])])
    if(ToSkip){
      if(verbose){
        message(paste("\nThe Full model",x[1],x[2],models[j]),"failed and thus won't be projected")
      }
      next
    }
    mod <- paste("ESM_Full",x[1],x[2],models[j],"model.out",sep="_")
    if(file.exists(mod)){
      mod <- get(load(mod))
      if(is.data.frame(new.env)){
        if(models[j]=="GLM"){
          pred <- pred.multiplier * predict.glm(mod, newdata = new.env, type = "response")
          if(rounded){
            pred <- round(pred)
          }
        }else if(models[j]=="GBM"){
          pred <- spsUtil::quiet(pred.multiplier * gbm::predict.gbm(mod,newdata = new.env, type = "response"))
          
          if(rounded){
            pred <- round(pred)
          }
        }else if(models[j]=="MAXNET"){
          pred <- spsUtil::quiet(pred.multiplier * predict(mod, newdata = new.env, type = "cloglog",clamp=FALSE)) ##invisible not working
          if(rounded){
            pred <- round(pred)
          }
        }else if(models[j]=="CTA"){
          pred <- spsUtil::quiet(pred.multiplier * as.data.frame(predict(mod, newdata = new.env, type = "prob")[,2]))
          if(rounded){
            pred <- round(pred)
          }
        }else if(models[j]=="GAM"){
          pred <- spsUtil::quiet(pred.multiplier * mgcv::predict.gam(mod,newdata = new.env, type = "response"))
          if(rounded){
            pred <- round(pred)
          }
        }else{
          pred <- spsUtil::quiet(pred.multiplier * predict(mod, newdata = new.env, type = "raw"))
          if(rounded){
            pred <- round(pred)
          }
        }
        done <- c(done,paste0(name.env,"/ESM_",x[1],"_",x[2],"_",models[j],".txt"))
        write.table(pred,paste0("../",name.env,"/ESM_",x[1],"_",x[2],"_",models[j],".txt"),sep="\t")
        
      }else{
        if(verbose){
          message(paste0("\n\t",models[j]))
        }
        
        
        if(parallel){
          new.env <- terra::unwrap(new.env)
        }
        if(models[j]=="GLM"){
          pred <- pred.multiplier * terra::predict(new.env, mod, type = "response",na.rm=T)
          if(rounded){
            pred <- round(pred)
          }
        }else if(models[j]=="GBM"){
          pred <- spsUtil::quiet(pred.multiplier * terra::predict(new.env,mod, fun = gbm::predict.gbm, type = "response",na.rm=TRUE)) ##need to test again
          if(rounded){
            pred <- round(pred)
          }
        }else if(models[j]=="MAXNET"){
          pred <- spsUtil::quiet(pred.multiplier * terra::predict(new.env,mod, type = "cloglog",clamp=FALSE,na.rm=TRUE)) ##invisible not working
          if(rounded){
            pred <- round(pred)
          }
        }else if(models[j]=="CTA"){
          pred <- spsUtil::quiet(pred.multiplier *  terra::predict(new.env,mod, type = "prob",na.rm = TRUE)[[2]]) ##invisible not working
          if(rounded){
            pred <- round(pred)
          }
        }else if(models[j]=="GAM"){
          pred <- spsUtil::quiet(pred.multiplier * terra::predict(new.env, mod, fun = mgcv::predict.gam, type = "response", na.rm = TRUE))
          if(rounded){
            pred <- round(pred)
          }
        }else{
          pred <- spsUtil::quiet(pred.multiplier * terra::predict(new.env,mod, type = "raw",na.rm=TRUE)) ##invisible not working
          if(rounded){
            pred <- round(pred)
          }
        }
        done <- c(done,paste0(name.env,"/ESM_",x[1],"_",x[2],"_",models[j],".tif"))
        terra::writeRaster(pred,paste0("../",name.env,"/ESM_",x[1],"_",x[2],"_",models[j],".tif"),
                           gdal=c("COMPRESS=DEFLATE","PREDICTOR=2"),datatype=datatype,overwrite=TRUE)
      }
      
    }else{
      warning(paste("\n",paste("ESM_Full",x[1],x[2],models[j],"model.out",sep="_"),"is not present and thus won't be projected"))
    }
    
    
  }
  return(done)
}
