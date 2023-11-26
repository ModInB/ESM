ESM_Projection <- function(ESM.Mod,
                           new.env,
                           name.env){
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

  if(is.data.frame(new.env)){
    if(sum(colnames(new.env) %in% used.env) != length(used.env) ){
      stop("new.env need to have the same variable names as used in ESM_Modeling")
    }
    
    proj <- apply(combinations, 2, .IndividualProj,
                  new.env = new.env,models = models,
                  name.env = name.env)
  }else{
    if(sum(names(new.env) %in% used.env) != length(used.env) ){
      stop("new.env need to have the same variable names as used in ESM_Modeling")
    }
    
    proj <- apply(combinations, 2, .IndividualProj,
                  new.env = new.env,models = models,
                  name.env = name.env)
  }
  
  obj = list(projection.path=proj,
             model.info = ESM.Mod$model.info)
  if(save.obj){
    save(obj,file=paste0("../ESM.Projection.",ESM.Mod$model.info$modeling.id,".out"))
  }
  return(obj)
}
.IndividualProj <- function(x,new.env,models,name.env){
  cat(paste("\n Projections of bivariate model:",x[1],x[2]))
  done <- c()
  for(j in 1:length(models)){
    cat(paste0("\n\t",models[j]))
    mod <- paste("ESM_Full",x[1],x[2],models[j],"model.out",sep="_")
    if(file.exists(mod)){
      mod <- get(load(mod))
      if(is.data.frame(new.env)){
        if(models[j]=="GLM"){
          pred <- round(1000 * predict.glm(mod, new.env = new.env, type = "response"))
          
        }else if(models[j]=="GBM"){
          pred <- invisible(round(1000 *gbm::predict.gbm(mod,new.env = new.env, type = "response")))
        }
        done <- c(done,paste0(name.env,"/ESM_",x[1],"_",x[2],"_",models[j],".txt"))
        write.table(pred,paste0("../",name.env,"/ESM_",x[1],"_",x[2],"_",models[j],".txt"),sep="\t")
        
      }else{
        
        if(models[j]=="GLM"){
          pred <- round(1000 * terra::predict(new.env, mod, type = "response"))
          
        }else if(models[j]=="GBM"){
          pred <- invisible(round(1000 *terra::predict(new.env,mod, type = "response"))) ##invisible not working
          pred <- terra::mask(pred,subset(new.env,1)) ## was putting values where both variables had na
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