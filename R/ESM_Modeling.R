## Modeling

ESM_Modeling <- function(resp,
                          xy,
                          env,
                          sp.name = "test",
                          models,
                          models.options = NULL,
                          prevalence = 0.5,
                          cv.method = "split-sampling", #can be split, block, custom
                          cv.rep = 10,
                          cv.ratio = 0.7,
                          cv.split.table = NULL,
                          which.biva = NULL,
                          modeling.id = as.character(format(Sys.time(), "%s")),
                          pathToSaveObject = getwd()){
  
  
  if(length(resp) != nrow(xy)){
    stop("resp and xy must have the same length")
  }
  if(ncol(xy)!=2){
    stop("xy should be a two-column matrix or data.frame")
  }
  if(is.data.frame(env)){
    if(length(resp) != nrow(env) | nrow(env) != nrow(xy)){
      stop("resp, xy and env must have the same length")
    }else{
      env.var <- env
    }
    
  }else if(inherits(env,"SpatRaster")){
    xy <- as.matrix(xy)
    env.var <- terra::extract(env,xy)
  }else{
    stop("env should be either a SpatRaster or a data.frame")
  }
  if(is.null(models.options)){
    models.options = new("ESM.models.options")
  }else{
    if(!inherits(models.options,"ESM.models.options")){
      stop("models.options should be an object from ESM_Models.options")
    }
  }
  if(anyNA(resp)){
    
    warning("NAs were present in resp and were converted to 0")
    resp[is.na(resp)] = 0
    
    if(is.null(prevalence)){
      cat("\nAs NAs were present in resp, we assume that you use pseudo-absences and we thus set prevalence to 0.5")
      prevalence = 0.5
    }
  }

  ##Need to generate the diff if 
  
  
  ### Start the process
  
  # Create a folder to store objects
  iniwd <- getwd()
  on.exit(setwd(iniwd))
  dir.create(paste0(pathToSaveObject,"/ESM.BIOMOD.output_", sp.name,"/",modeling.id),recursive = T)
  newwd <- paste0("./ESM.BIOMOD.output_", sp.name,"/",modeling.id)
  setwd(newwd)
  
  # Generate all the possible combination variables
  combinations <- combn(colnames(env.var), 2)
  
  if (is.null(which.biva)) {
    which.biva <- 1:ncol(combinations)
  }
  

  if(is.null(cv.split.table)){
    cv.split.table <- .CreatingDataSplitTableV2(resp = resp,
                                                cv.rep = cv.rep,
                                                cv.method = cv.method,
                                                cv.ratio = cv.ratio)
  }
  
  combinations <- combinations[,which.biva]
  
  cat("\n################### Start Modelling ###################")
  
  biva.mods <- apply(combinations, 2, .bivaModeling,resp = resp,
                     env.var = env.var, models = models,
                     models.options = models.options,
                     cv.split.table = cv.split.table,prevalence = prevalence,
                     simplify = FALSE)
  names(biva.mods) = paste(combinations[1,],combinations[2,])
  cat("\n##################### Done #####################")
  
  cat("\n############### Start evaluations ###############")
  
  ## Evaluation
  biva.eval <- lapply(biva.mods,.bivaEvaluation,
                      resp=resp,
                      cv.split.table=cv.split.table)
  obj <- list(resp = resp,
              xy = xy,
              env.var = env.var,
              sp.name= sp.name,
              cv.split.table = cv.split.table,
              biva.predictions = biva.mods,
              biva.evaluations = biva.eval,
              modeling.id = modeling.id,
              biva.path = newwd
              )
  save(obj,file=paste0("../ESM.Modeling_",modeling.id,".out"))
  cat("\n##################### Done #####################")
  
  return(obj)
}


#Function to generate the argument DataSplitTable in the function ecospat.ESM.Modeling                                                                        
.CreatingDataSplitTableV2 <- function(resp,
                                      cv.rep,
                                      cv.method,
                                      cv.ratio){
  
  calib.Lines <- matrix(FALSE, nrow = length(resp), ncol = cv.rep)
  
  if(cv.method == "split-sampling"){

    pres <- which(resp==1)
    abs <- which(resp==0)
    
    for(i in 1:cv.rep){
      calib.Lines[sample(pres,size = round(length(pres)*cv.ratio)),i] = TRUE
      calib.Lines[sample(abs,size = round(length(abs)*cv.ratio)),i] = TRUE
    }
    
  }else{ ##Create blocks
    TODO 
  }
  calib.Lines <- cbind(calib.Lines,TRUE)
  colnames(calib.Lines) = c(paste0("_RUN",1:cv.rep),"_Full")
  return(calib.Lines)
}

.bivaModeling <- function(x,
                          resp,
                          env.var,
                          models,
                          models.options,
                          cv.split.table,
                          prevalence){
  
  cat(c("\nCombinations", as.character(x)))
  envi <- env.var[,as.character(x)]
  cv.split.table <- rbind.data.frame(colnames(cv.split.table),cv.split.table)
  d <- apply(cv.split.table,2,.doModeling,resp = resp,
             env.var = envi, models = models,
             models.options = models.options,
             prevalence = prevalence)
  return(d)
}

.doModeling <- function(x,
                        resp,
                        env.var,
                        models,
                        models.options,
                        prevalence,
                        combi){
    nameRun <- x[1]
    x <- as.logical(x[-1])
    data <- cbind.data.frame(resp = resp[x],
                             env.var[x,])
    if(!is.null(prevalence)){
      ratio.Pres.Abs <- table(data$resp)/nrow(data)
      ratio.Pres.Abs <- prevalence/ratio.Pres.Abs
      ratio.Pres.Abs <- round(ratio.Pres.Abs/min(ratio.Pres.Abs)) #so that the min weight is 1
      weights <- data$resp
      weights[data$resp==1] = ratio.Pres.Abs["1"]
      weights[data$resp==0] = ratio.Pres.Abs["0"]
    }else{
      weights <- rep(1,nrow(data))
    } 
    for(j in 1: length(models)){
      
      
      if(models[j] == "GLM"){
        if(models.options@GLM$test == "none"){
          if(is.null(models.options@GLM$myFormula)){
            formula <- as.formula(paste0("resp~", 
                                         paste0(colnames(env.var),"+I(",colnames(env.var),"^2)",collapse = "+")))
            mod <- glm(formula = formula,
                       family = models.options@GLM$family,
                       weights = weights,
                       data = data)
            
            save(mod,file=paste("ESM",nameRun,
                                colnames(env.var)[1],
                                colnames(env.var)[2],
                                models[j],"model.out",
                                sep="_"))
          } 
          
        }else if(models.options@GLM$test == "AIC"){
          step()
        }
        pred <- predict(mod,newdata = env.var,type = "response")
      }
      
      
    }
    return(pred)
}
