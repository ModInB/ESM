
#' @name ESM_Projection
#' @author Flavien Collart \email{flaviencollart@hotmail.com} 
#' @title Ensemble of Small Models: Projections of Bivariate Models
#' @description Project bivariate models 
#' @param ESM.Mod The object returned by \code{\link{ESM_Modeling}}.
#' @param new.env a \code{data.frame}, \code{matrix} or a \code{SpatRaster} of the predictors where the models 
#' should be projected. Note that the colnames or the names of the predictors should be exactly the same as the one used in \code{\link{ESM_Modeling}}.
#' @param name.env a \code{character} which will be used to generate a folder to save the individual projections.
#' @param parallel \code{logical}. Allows or not parallel job using the functions makeCluster.
#' @param n.cores \code{integer}. Number of CPU cores used to make the models.
#' @param save.obj \code{logical}. Allows or not the final output.
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
#' Each bivariate models are geographically projected into the new.env area. Note that the projected maps are multiplied 
#' by 1000 and rounded to reduce storage space. 
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
                           parallel = FALSE,
                           n.cores = 1,
                           save.obj = T){
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
                                 name.env = name.env, biva.pred = biva.pred)
      parallel::stopCluster(cl)
      
    }else{
      proj <- apply(combinations, 2, .IndividualProj,
                    new.env = new.env,models = models,
                    name.env = name.env, biva.pred = biva.pred)
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
                                 biva.pred = biva.pred)
      parallel::stopCluster(cl)
      
    }else{
      proj <- apply(combinations, 2, .IndividualProj,
                    new.env = new.env,models = models,
                    name.env = name.env,parallel = parallel, 
                    biva.pred = biva.pred)
    }
  }
  
  obj = list(projection.path=as.list(proj),
             model.info = ESM.Mod$model.info,
             name.env = name.env,
             proj.type=proj.type)
  if(save.obj){
    save(obj,file=paste0("../ESM.Projection.",ESM.Mod$model.info$modeling.id,".out"))
  }
  return(obj)
}

## The hidden functions ----
## .IndividualProj----
## Function projecting a model onto a new environment
#' @import nnet
#' @import rpart
#' @import maxnet
.IndividualProj <- function(x,new.env,models,name.env,parallel, biva.pred){
  if(!is.data.frame(new.env)){
    cat(paste("\n Projections of bivariate model:",x[1],x[2]))}
  done <- c()
  for(j in 1:length(models)){
    ToSkip <- anyNA(biva.pred[[paste(x[1],x[2], sep=".")]][,paste0("Full.",models[j])])
    if(ToSkip){
      cat(paste("\nThe Full model",x[1],x[2],models[j]),"failed and thus won't be projected")
      next
    }
    mod <- paste("ESM_Full",x[1],x[2],models[j],"model.out",sep="_")
    if(file.exists(mod)){
      mod <- get(load(mod))
      if(is.data.frame(new.env)){
        if(models[j]=="GLM"){
          pred <- round(1000 * predict.glm(mod, newdata = new.env, type = "response"))
          
        }else if(models[j]=="GBM"){
          pred <- spsUtil::quiet(round(1000 * gbm::predict.gbm(mod,newdata = new.env, type = "response")))
        }else if(models[j]=="MAXNET"){
          pred <- spsUtil::quiet(round(1000 * predict(mod, newdata = new.env, type = "cloglog",clamp=FALSE))) ##invisible not working
        }else if(models[j]=="CTA"){
          pred <- spsUtil::quiet(round(1000 * as.data.frame(predict(mod, newdata = new.env, type = "prob")[,2]))) ##invisible not working
        }else{
          pred <- spsUtil::quiet(round(1000 * predict(mod, newdata = new.env, type = "raw")))
        }
        done <- c(done,paste0(name.env,"/ESM_",x[1],"_",x[2],"_",models[j],".txt"))
        write.table(pred,paste0("../",name.env,"/ESM_",x[1],"_",x[2],"_",models[j],".txt"),sep="\t")
        
      }else{
        cat(paste0("\n\t",models[j]))
        
        if(parallel){
          new.env <- terra::unwrap(new.env)
        }
        if(models[j]=="GLM"){
          pred <- round(1000 * terra::predict(new.env, mod, type = "response",na.rm=T))
        }else if(models[j]=="GBM"){
          pred <- spsUtil::quiet(round(1000 * terra::predict(new.env,mod, fun = gbm::predict.gbm, type = "response",na.rm=T))) ##need to test again
        }else if(models[j]=="MAXNET"){
          pred <- spsUtil::quiet(round(1000 * terra::predict(new.env,mod, type = "cloglog",clamp=FALSE,na.rm=T))) ##invisible not working
        }else if(models[j]=="CTA"){
          pred <- spsUtil::quiet(round(1000 *  terra::predict(new.env,mod, type = "prob",na.rm = T)[[2]])) ##invisible not working
        }else{
          pred <- spsUtil::quiet(round(1000 * terra::predict(new.env,mod, type = "raw",na.rm=T))) ##invisible not working
        }
        done <- c(done,paste0(name.env,"/ESM_",x[1],"_",x[2],"_",models[j],".tif"))
        terra::writeRaster(pred,paste0("../",name.env,"/ESM_",x[1],"_",x[2],"_",models[j],".tif"),
                           gdal=c("COMPRESS=DEFLATE","PREDICTOR=2"),datatype="INT2U",overwrite=T)
      }
      
    }else{
      cat(paste("\n",paste("ESM_Full",x[1],x[2],models[j],"model.out",sep="_"),"is not present and thus won't be projected"))
    }
    
    
  }
  return(done)
}
