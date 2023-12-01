#################################################################################################################################################
## ESM_Projection: 
##  Description:
##    Projects each bivariate model on new.env
##
##  Arguments:
## @ESM.Mod: The outputs resulted from ESM_Modeling()
## @new.env: a data.frame, a matrix or a SpatRaster containg the environment where the models should be projected. 
##           Note that the colnames or the names of the predictors should be exactly the same as the one used in ESM_Modeling()
## @name.env: a character which will be used to generate a folder to save the individual projections.
## @paralell: logical. Allows or not parallel job using the functions makeCluster
## @n.cores: numeric. Number of cores used to make the models.
##
## Values: 
##        a list containing: 
##                          projection.path: a vector of the path to each projection.
##                          model.info: contains the models used, their options, the combination of bivariate models (which.biva), the
##                                      modeling ID, the path to the folder where are the stored the models (biva.path), and the failed models
##                          name.env: the value given to the argument name.env
##                          proj.type: a character equal to data.frame or SpatRaster depending on the class of new.env
##
##        Maps or data.frame are stored in the folder named by name.env. Values are multiplied by 1000 and rounded. Note that the maps 
##        are in .tif and compressed with the option COMPRESS=DEFLATE and PREDICTOR=2
##  Authors:
##          Flavien Collart based on the previous code written by Frank Breiner with the contributions of 
##          Olivier Broennimann and Flavien Collart
#' @export
#################################################################################################################################################

ESM_Projection <- function(ESM.Mod,
                           new.env,
                           name.env,
                           parallel = FALSE,
                           n.cores = 1){
  iniwd <- getwd()
  on.exit(setwd(iniwd))
  setwd(ESM.Mod$model.info$biva.path)
  
  models <- ESM.Mod$model.info$models
  which.biva <- ESM.Mod$model.info$which.biva
  env.var <- ESM.Mod$data$env.var
  
  #### check new.env
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
  ## Compute the realized combinations
  combinations <- combn(colnames(env.var), 2)
  if (is.null(which.biva)) {
    which.biva <- 1:ncol(combinations)
  }
  combinations <- combinations[,which.biva]
  used.env <- unique(c(combinations[1,],combinations[2,]))
  if("MAXNET" %in% models){
    require(maxnet)
  }
  if(is.data.frame(new.env)){
    if(sum(colnames(new.env) %in% used.env) != length(used.env) ){
      stop("new.env need to have the same variable names as used in ESM_Modeling")
    }
    proj.type = "data.frame"
    if(parallel){
      cl <- parallel::makeCluster(n.cores)
      proj <- parallel::parApply(cl,combinations, 2, .IndividualProj,
                                 new.env = new.env,models = models,
                                 name.env = name.env)
      parallel::stopCluster(cl)
      
    }else{
      proj <- apply(combinations, 2, .IndividualProj,
                    new.env = new.env,models = models,
                    name.env = name.env)
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
                                 name.env = name.env,parallel = parallel)
      parallel::stopCluster(cl)
      
    }else{
      proj <- apply(combinations, 2, .IndividualProj,
                    new.env = new.env,models = models,
                    name.env = name.env,parallel = parallel)
    }
  }
  
  obj = list(projection.path=proj,
             model.info = ESM.Mod$model.info,
             name.env = name.env,
             proj.type=proj.type)
  if(save.obj){
    save(obj,file=paste0("../ESM.Projection.",ESM.Mod$model.info$modeling.id,".out"))
  }
  return(obj)
}

## Function projecting a model onto a new environment.
.IndividualProj <- function(x,new.env,models,name.env,parallel){
  cat(paste("\n Projections of bivariate model:",x[1],x[2]))
  done <- c()
  for(j in 1:length(models)){
    cat(paste0("\n\t",models[j]))
    mod <- paste("ESM_Full",x[1],x[2],models[j],"model.out",sep="_")
    if(file.exists(mod)){
      mod <- get(load(mod))
      if(is.data.frame(new.env)){
        if(models[j]=="GLM"){
          pred <- round(1000 * predict.glm(mod, newdata = new.env, type = "response"))
          
        }else if(models[j]=="GBM"){
          pred <- invisible(round(1000 *gbm::predict.gbm(mod,newdata = new.env, type = "response")))
        }else{
          pred <- invisible(round(1000 * predict(mod, newdata= new.env, type = "cloglog",clamp=FALSE))) ##invisible not working
        }
        done <- c(done,paste0(name.env,"/ESM_",x[1],"_",x[2],"_",models[j],".txt"))
        write.table(pred,paste0("../",name.env,"/ESM_",x[1],"_",x[2],"_",models[j],".txt"),sep="\t")
        
      }else{
        if(parallel){
          new.env <- terra::unwrap(new.env)
        }
        if(models[j]=="GLM"){
          pred <- round(1000 * terra::predict(new.env, mod, type = "response"))
        }else if(models[j]=="GBM"){
          pred <- invisible(round(1000 * terra::predict(new.env,mod, fun = gbm::predict.gbm, type = "response"))) ##invisible not working
          pred <- terra::mask(pred,subset(new.env,1)) ## was putting values where both variables had na
        }else{
          pred <- invisible(round(1000 * terra::predict(new.env,mod, type = "cloglog",clamp=FALSE,na.rm=T))) ##invisible not working
        }
        done <- c(done,paste0(name.env,"/ESM_",x[1],"_",x[2],"_",models[j],".tif"))
        terra::writeRaster(pred,paste0("../",name.env,"/ESM_",x[1],"_",x[2],"_",models[j],".tif"),
                           gdal=c("COMPRESS=DEFLATE","PREDICTOR=2"),overwrite=T)
      }
      
    }else{
      cat(paste("\n",paste("ESM_Full",x[1],x[2],models[j],"model.out",sep="_"),"is not present and thus won't be projected"))
    }
    
    
  }
  return(done)
}
