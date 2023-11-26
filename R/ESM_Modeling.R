## Modeling
## save.obj if FALSE, it only save the full models
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
                          pathToSaveObject = getwd(),
                          save.obj = TRUE){
  
  
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
    models.options = ESM_Models.options()
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
  combinations <- combinations[,which.biva]

  if(is.null(cv.split.table)){
    cv.split.table <- ESM.CreatingDataSplitTable(resp = resp, 
                                                 cv.rep = cv.rep,
                                                 cv.method = cv.method,
                                                 cv.ratio = cv.ratio)
  }
  
  
  cat("\n################### Start Modelling ###################")
  ## Model each bivariate model
  biva.mods <- apply(combinations, 2, .bivaModeling,resp = resp,
                     env.var = env.var, models = models,
                     models.options = models.options,
                     cv.split.table = cv.split.table,
                     prevalence = prevalence, save.obj = save.obj,
                     simplify = FALSE) 
  if(length(models)==1){
    for(i in 1:length(biva.mods)){
      colnames(biva.mods[[i]]) = paste0(colnames(cv.split.table),".",
                                        colnames(biva.mods[[i]]))
    }
    
  }
  names(biva.mods) = paste0(combinations[1,],".",combinations[2,])
  cat("\n##################### Done #####################")
  
  cat("\n############### Start evaluations ###############")
  
  ## Evaluation
  biva.eval <- lapply(biva.mods,.bivaEvaluation,
                      resp=resp, models=models,
                      cv.split.table=cv.split.table)
  
  
  obj <- list(data = list(resp = resp,
                          xy = xy,
                          env.var = env.var,
                          sp.name= sp.name),
              model.info = list(models = models,
                                models.options = models.options,
                                which.biva = which.biva,
                                modeling.id = modeling.id,
                                biva.path = newwd),
              cv.split.table = cv.split.table,
              biva.predictions = biva.mods,
              biva.evaluations = biva.eval
              )
  
  if(save.obj){
    save(obj,file=paste0("../ESM.Modeling.",modeling.id,".out"))
  }
  cat("\n##################### Done #####################")
  
  return(obj)
}


#Function to generate the argument DataSplitTable in the function ecospat.ESM.Modeling                                                                        
ESM.CreatingDataSplitTable <- function(resp,
                                       cv.method,
                                       cv.rep = NULL,
                                       cv.ratio = NULL,
                                       cv.n.blocks = NULL){
  
  calib.Lines <- matrix(FALSE, nrow = length(resp), ncol = cv.rep)
  
  if(cv.method == "split-sampling"){

    pres <- which(resp==1)
    abs <- which(resp==0)
    
    for(i in 1:cv.rep){
      calib.Lines[sample(pres,size = round(length(pres)*cv.ratio)),i] = TRUE
      calib.Lines[sample(abs,size = round(length(abs)*cv.ratio)),i] = TRUE
    }
    
    calib.Lines <- cbind(calib.Lines,TRUE)
    colnames(calib.Lines) = c(paste0("RUN",1:cv.rep),"Full")
    
  }else{ ##Create blocks
    stop("Not Yet available") 
  }
  
  return(calib.Lines)
}

.bivaModeling <- function(x,
                          resp,
                          env.var,
                          models,
                          models.options,
                          cv.split.table,
                          prevalence,
                          save.obj = TRUE){
  
  cat(c("\n\nCombinations", as.character(x)))
  envi <- env.var[,as.character(x)]
  cv.split.table <- rbind.data.frame(colnames(cv.split.table),cv.split.table)
  ## Model each run
  d <- apply(cv.split.table,2,.doModeling,resp = resp,
             env.var = envi, models = models,
             models.options = models.options,
             prevalence = prevalence, save.obj = save.obj)
  return(do.call(cbind,d))
}

.doModeling <- function(x,
                        resp,
                        env.var,
                        models,
                        models.options,
                        prevalence,
                        save.obj = TRUE){
    
  nameRun <- x[1]
  x <- as.logical(x[-1])
  data <- cbind.data.frame(resp = resp[x],
                           env.var[x,])
  if(!is.null(prevalence)){
    ratio.Pres.Abs <- table(data$resp)/nrow(data)
    ratio.Pres.Abs <- prevalence/ratio.Pres.Abs
    ratio.Pres.Abs <- round(ratio.Pres.Abs/min(ratio.Pres.Abs)) #so that the min weight is 1
    w <- data$resp
    w[data$resp==1] = ratio.Pres.Abs["1"]
    w[data$resp==0] = ratio.Pres.Abs["0"]
  }else{
    w <- rep(1,nrow(data))
  } 
  data$w = w
    for(j in 1: length(models)){
      
      
      if(models[j] == "GLM"){
        cat(paste("\nGLM", nameRun))
        
        if(is.null(models.options$GLM$myFormula)){
          if(models.options$GLM$test == "none"){
            
            formula <- .makeGLMFormula(env.var,
                                       models.options$GLM)
            mod <- glm(formula = formula,
                       family = models.options$GLM$family,
                       weights = w,
                       data = data)
            if(save.obj){
              save(mod,file=paste("ESM",nameRun,
                                  colnames(env.var)[1],
                                  colnames(env.var)[2],
                                  models[j],"model.out",
                                  sep="_"))
            }else if(nameRun == "Full"){
              save(mod,file=paste("ESM",nameRun,
                                  colnames(env.var)[1],
                                  colnames(env.var)[2],
                                  models[j],"model.out",
                                  sep="_"))
            }
            
          }else if(models.options$GLM$test == "AIC"){
            formula <- .makeGLMFormula(env.var,
                                       models.options$GLM)
            mod.full<- glm(formula = formula,
                           family = models.options$GLM$family,
                           weights = w,
                           data = data)
            mod <- step(mod.full,
                        scope = "resp~1",
                        direction = "both",
                        trace=F)
            cat(paste0("\n\tBest Formula:",deparse(mod$formula)))
            if(save.obj){
              save(mod,file=paste("ESM",nameRun,
                                  colnames(env.var)[1],
                                  colnames(env.var)[2],
                                  models[j],"model.out",
                                  sep="_"))
            }else if(nameRun == "Full"){
              save(mod,file=paste("ESM",nameRun,
                                  colnames(env.var)[1],
                                  colnames(env.var)[2],
                                  models[j],"model.out",
                                  sep="_"))
            }
          } 
          
        }else{
          mod <- glm(formula = models.options$GLM$myFormula,
                     family = models.options$GLM$family,
                     weights = w,
                     data = data)
          if(save.obj){
            save(mod,file=paste("ESM",nameRun,
                                colnames(env.var)[1],
                                colnames(env.var)[2],
                                models[j],"model.out",
                                sep="_"))
          }else if(nameRun == "Full"){
            save(mod,file=paste("ESM",nameRun,
                                colnames(env.var)[1],
                                colnames(env.var)[2],
                                models[j],"model.out",
                                sep="_"))
          }
        }
        pred <- as.data.frame(predict(mod,newdata = env.var,type = "response"))
        colnames(pred) = "GLM" 
      }
      
      if(models[j]=="GBM"){
        cat(paste("\nGBM", nameRun,"\n"))
        
        formula <- .makeGLMFormula(env.var,
                                   model.option=list(type="linear"))
        mod <- gbm::gbm(formula = formula,
                        distribution = models.options$GBM$distribution,
                        data = data,
                        weights = w,
                        n.trees = models.options$GBM$n.trees,
                        interaction.depth = models.options$GBM$interaction.depth,
                        n.minobsinnode = models.options$GBM$n.minobsinnode,
                        shrinkage = models.options$GBM$shrinkage,
                        bag.fraction = models.options$GBM$bag.fraction,
                        train.fraction = models.options$GBM$train.fraction,
                        cv.folds = models.options$GBM$cv.folds,
                        keep.data = TRUE,
                        verbose = models.options$GBM$verbose,
                        n.cores = models.options$GBM$n.cores)
        pred <- as.data.frame(gbm::predict.gbm(mod,newdata = env.var,type = "response"))
        colnames(pred) = "GBM" 
        
        if(save.obj){
          save(mod,file=paste("ESM",nameRun,
                              colnames(env.var)[1],
                              colnames(env.var)[2],
                              models[j],"model.out",
                              sep="_"))
        }else if(nameRun == "Full"){
          save(mod,file=paste("ESM",nameRun,
                              colnames(env.var)[1],
                              colnames(env.var)[2],
                              models[j],"model.out",
                              sep="_"))
        }
      }
      
      if(j == 1){
        predFin <- pred
        
      }else{
        predFin <- cbind.data.frame(predFin,pred)
      }
    }
    return(predFin)
}

.makeGLMFormula <- function(env.var = env.var,model.option){
  
  formula <- paste0("resp~", paste0(colnames(env.var),collapse = "+"))
  
  if(model.option$type == "quadratic" | model.option$type == "polynomial"){
    IsNum <- apply(env.var,2,is.numeric)
    Topaste <- paste0("I(",colnames(env.var)[IsNum],"^2)",collapse="+")
    formula <- paste(formula,Topaste,sep = "+")
    if(model.option$type == "polynomial"){
      Topaste <- paste0("I(",colnames(env.var)[IsNum],"^3)",collapse="+")
      formula <- paste(formula,Topaste,sep = "+")
    }
  }
  
  return(as.formula(formula))
  
}
