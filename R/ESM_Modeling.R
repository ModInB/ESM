## Modeling
## save.obj if FALSE, it only save the full models
ESM_Modeling <- function( resp,
                          xy,
                          env,
                          sp.name,
                          models,
                          models.options = NULL,
                          prevalence = 0.5,
                          cv.method = "split-sampling", #can be split, block, custom
                          cv.rep = 10,
                          cv.ratio = 0.7,
                          cv.n.blocks = NULL,
                          cv.split.table = NULL,
                          which.biva = NULL,
                          parallel = FALSE,
                          n.cores = 1,
                          modeling.id = as.character(format(Sys.time(), "%s")),
                          pathToSaveObject = getwd(),
                          save.obj = TRUE){
  
  ## Check resp, XY, sp.name and prevalence
  if(length(resp) != nrow(xy)){
    stop("resp and xy must have the same length")
  }
  if(anyNA(resp)){
    
    warning("NAs were present in resp and were converted to 0")
    resp[is.na(resp)] = 0
    
    if(is.null(prevalence)){
      cat("\nAs NAs were present in resp, we assume that you use pseudo-absences and we thus set prevalence to 0.5")
      prevalence = 0.5
    }
  }
  if(ncol(xy)!=2){
    stop("xy should be a two-column matrix or data.frame")
  }
  if(is.null(sp.name) | !(is.character(sp.name))){
    stop("sp.name should be a character object")
  }
  if(!is.null(prevalence) & (prevalence >= 1 | prevalence <= 0)){
    stop("prevalence must be inside ]0;1[ or null")
  }
  
  #########################################
  ## Check model names 
  
  if(any(!(models  %in% c("GLM","GBM","MAXNET")))){
    stop("models should be = to GLM, GBM, and/or MAXNET")
  }
  
  #######################################
  ## Check model options
  if(is.null(models.options)){
    models.options = ESM_Models.Options()
  }else{
    if(!is.list(models.options) | deparse(names(models.options))!= deparse(c("GLM","GBM"))){
     stop("models.options should null or formatted via ESM_Models.Options()") 
    }
  }
  
  ###############################################
  ## Check env and extract values if SpatRaster
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
  

  
  ###############################################
  ## Check split.tables and generate one if it is null
  if(sum(cv.method %in% c("split-sampling","block","custom")) != 1){
    stop("cv.method shoud be either split-sampling, block or custom ")
  }
  if(cv.method == "split-sampling"){
    if(cv.rep<1){
      stop("When cv.method=split-sampling, cv.rep must be at least 1")
    }
    if(cv.ratio > 1 | cv.ratio < 0){
      stop("When cv.method=split-sampling, cv.ratio must be comprised between 0 and 1")
      
    }
  }else if(cv.method == "block"){
    if(is.null(cv.n.blocks) | cv.n.blocks<2){
      stop("When cv.method=block, cv.n.blocks should not be null and be at least greater than 1")
    }
  }else{
    if(is.null(cv.split.table)){
      stop("When cv.method = custom, cv.split.table cannot be null")
    }
    if(sum(apply(cv.split.table, 2, is.logical)) != ncol(cv.split.table)){
      stop("All columns from cv.split.table should be logical.")
    }
    if(length(grep("RUN", colnames(cv.split.table))) == 0 | length(grep("Full", colnames(cv.split.table))) == 0){
      stop("When cv.method = custom, the colnames of cv.split.table should be RUNX, where X is a number from 1 to the needed number of replicates and
           the last column should be called Full (/!\ case sensitive).")
    }else if(sum(cv.split.table[,"Full"]) != nrow(cv.split.table)){
      stop("The column Full in cv.split.table should be entirely filled with TRUE")
    }
  }
  
  if(is.null(cv.split.table)){
    cv.split.table <- ESM.CreatingDataSplitTable(resp = resp, 
                                                 cv.rep = cv.rep,
                                                 cv.method = cv.method,
                                                 cv.ratio = cv.ratio,
                                                 cv.n.blocks = cv.n.blocks)
  }
  
  ####################################
  ## Remove NAs in env.var
  is.thereNAs <- is.na(apply(env.var,1,sum))
  if(sum(is.thereNAs)>0){
    n.presAbs <- table(resp)
    n.presAbs <- c(n.presAbs["1"],n.presAbs["0"]) ##Make sure that the order is the same
    resp <- resp[!(is.thereNAs)]
    n.presAbs.new <-  table(resp)
    n.presAbs.new <- c(n.presAbs.new["1"],n.presAbs.new["0"])
    change.n <- n.presAbs - n.presAbs.new 
    warning(paste("\nNAs were found in env and were thus removed.\n",change.n["0"],"absences",change.n["1"], "presences were removed"))
    xy <- xy[!(is.thereNAs),]
    env.var = na.omit(env.var)
  }
  # Check if presences and absences are present even after removing possible NAs
  if(sum(resp)==0 | sum(resp==0) == 0){
    stop("presences and absences/pseudo-absences should be present in resp")
  }
  
  ################################
  ### Start the process
  
  # Create a folder to store objects
  iniwd <- getwd()
  on.exit(setwd(iniwd))
  dir.create(paste0(pathToSaveObject,"/ESM.output_", sp.name,"/",modeling.id),recursive = T)
  newwd <- paste0("./ESM.output_", sp.name,"/",modeling.id)
  setwd(newwd)
  
  # Generate all the possible combination variables
  combinations <- combn(colnames(env.var), 2)
  
  if (is.null(which.biva)) {
    which.biva <- 1:ncol(combinations)
  }else if(sum(!(which.biva %in% (1:ncol(combinations))))>0){ ## Error check
    stop(paste("which.biva should be an integer vector with values inside", deparse(as.character(1:ncol(combinations)))))
  }
  
  combinations <- combinations[,which.biva]

  
 
  
  cat("\n################### Start Modelling ###################")

  if(parallel){
    cl <- parallel::makeCluster(n.cores)
    biva.mods <-  parallel::parApply(cl, combinations, 2, 
                                     .bivaModeling,resp = resp,
                       env.var = env.var, models = models,
                       models.options = models.options,
                       cv.split.table = cv.split.table,
                       prevalence = prevalence, save.obj = save.obj)
    parallel::stopCluster(cl)
  }else{
    biva.mods <- apply(combinations, 2, .bivaModeling,resp = resp,
                       env.var = env.var, models = models,
                       models.options = models.options,
                       cv.split.table = cv.split.table,
                       prevalence = prevalence, save.obj = save.obj,
                       simplify = FALSE)
  }
  
  if(length(models)==1){
    for(i in 1:length(biva.mods)){
      colnames(biva.mods[[i]]) = paste0(colnames(cv.split.table),".",
                                        colnames(biva.mods[[i]]))
    }
    
  }
  
  names(biva.mods) = paste0(combinations[1,],".",combinations[2,])
  
  cat("\n##################### Done #####################")
  
  failed.mod <- do.call(cbind, biva.mods)
  mod.pred.NAs <- apply(failed.mod, 2, anyNA)
  failed.mod <- colnames(failed.mod)[mod.pred.NAs]
  
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
                                failed.mod = failed.mod,
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
  
 
  
  pres <- which(resp==1)
  abs <- which(resp==0)
  
  if(cv.method == "split-sampling"){
    calib.Lines <- matrix(FALSE, nrow = length(resp), ncol = cv.rep)
    for(i in 1:cv.rep){
      calib.Lines[sample(pres,size = round(length(pres)*cv.ratio)),i] = TRUE
      calib.Lines[sample(abs,size = round(length(abs)*cv.ratio)),i] = TRUE
    }
    
    calib.Lines <- cbind(calib.Lines,TRUE)
    colnames(calib.Lines) = c(paste0("RUN",1:cv.rep),"Full")
    
  }else{ 
    calib.Lines <- matrix(FALSE, nrow = length(resp), ncol = cv.n.blocks)
    pres.Random <- sample(pres,size = length(pres))
    abs.Random <- sample(abs,size = length(abs))
    size.blockPres = round(length(pres)/cv.n.blocks)
    if(size.blockPres<2){
      stop("Less than 2 occurrences are present in the dataset for the evaluation. Consider reducing the number of blocks")
    }
    size.blockAbs = round(length(abs)/cv.n.blocks)
    if(size.blockAbs<2){
      stop("Less than 2 absences are present in the dataset for the evaluation. Consider reducing the number of blocks")
    }
    for(i in 0:(cv.n.blocks-1)){
      if(i < (cv.n.blocks-1)){
        calib.Lines[pres.Random[(i*size.blockPres+1):(size.blockPres*(i+1))],i] = TRUE
        calib.Lines[abs.Random[(i*size.blockAbs+1):(size.blockAbs*(i+1))],i] = TRUE
      }else{
        calib.Lines[pres.Random[(i*size.blockPres+1):length(pres.Random)],i] = TRUE
        calib.Lines[abs.Random[(i*size.blockAbs+1):length(abs.Random)],i] = TRUE
      }
      
    }
    calib.Lines <- cbind(calib.Lines,TRUE)
    colnames(calib.Lines) = c(paste0("RUN",1:cv.n.blocks),"Full")
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
  
  .doModeling <- function(x,
                          resp,
                          env.var,
                          models,
                          models.options,
                          prevalence,
                          save.obj = TRUE){
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
    data.maxnet <- data
    data$w = w
    for(j in 1: length(models)){
      err <- FALSE
      
      if(models[j] == "GLM"){
        cat(paste("\nGLM", nameRun))
        
        if(is.null(models.options$GLM$myFormula)){
          if(models.options$GLM$test == "none"){
            
            formula <- .makeGLMFormula(env.var,
                                       models.options$GLM)
            tryCatch(expr={mod <- glm(formula = formula,
                                      family = models.options$GLM$family,
                                      weights = w,
                                      data = data)}, error=function(e){
                                        cat(paste("\n model",models[j],nameRun,"failed"))
                                        err <<-TRUE
                                      })
            if(err){
              pred <- as.data.frame(rep(NA,nrow(env.var)))
              colnames(pred) = "GLM" 
            }else{
              pred <- as.data.frame(predict(mod,newdata = env.var,type = "response"))
              colnames(pred) = "GLM" 
            }
            
            if(save.obj & !(err)){
              save(mod,file=paste("ESM",nameRun,
                                  colnames(env.var)[1],
                                  colnames(env.var)[2],
                                  models[j],"model.out",
                                  sep="_"))
            }else if(nameRun == "Full" & !(err)){
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
            tryCatch(expr={
              mod <- step(mod.full,
                          scope = "resp~1",
                          direction = "both",
                          trace=F)
              cat(paste0("\n\tBest Formula:",deparse(mod$formula)))}, error=function(e){
                cat(paste("\n model",models[j],nameRun,"failed"))
                err <<-TRUE
              })
            
            if(err){
              pred <- as.data.frame(rep(NA,nrow(env.var)))
              colnames(pred) = "GLM" 
            }else{
              pred <- as.data.frame(predict(mod,newdata = env.var,type = "response"))
              colnames(pred) = "GLM" 
            }
            
            if(save.obj & !(err)){
              save(mod,file=paste("ESM",nameRun,
                                  colnames(env.var)[1],
                                  colnames(env.var)[2],
                                  models[j],"model.out",
                                  sep="_"))
            }else if(nameRun == "Full" & !(err)){
              save(mod,file=paste("ESM",nameRun,
                                  colnames(env.var)[1],
                                  colnames(env.var)[2],
                                  models[j],"model.out",
                                  sep="_"))
            }
          } 
          
        }else{
          tryCatch(expr={mod <- glm(formula = models.options$GLM$myFormula,
                                    family = models.options$GLM$family,
                                    weights = w,
                                    data = data)}, error=function(e){
                                      cat(paste("\n model",models[j],nameRun,"failed"))
                                      err <<-TRUE
                                    })
          
          if(err){
            pred <- as.data.frame(rep(NA,nrow(env.var)))
            colnames(pred) = "GLM" 
          }else{
            pred <- as.data.frame(predict(mod,newdata = env.var,type = "response"))
            colnames(pred) = "GLM" 
          }
          
          if(save.obj & !(err)){
            save(mod,file=paste("ESM",nameRun,
                                colnames(env.var)[1],
                                colnames(env.var)[2],
                                models[j],"model.out",
                                sep="_"))
          }else if(nameRun == "Full"& !(err)){
            save(mod,file=paste("ESM",nameRun,
                                colnames(env.var)[1],
                                colnames(env.var)[2],
                                models[j],"model.out",
                                sep="_"))
          }
        }
        
      }
      
      if(models[j]=="GBM"){
        cat(paste("\nGBM", nameRun,"\n"))
        formula <- .makeGLMFormula(env.var,
                                   model.option=list(type="linear"))
        tryCatch(expr = {mod <- gbm::gbm(formula = formula,
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
                                         n.cores = models.options$GBM$n.cores)}, error=function(e){
                                           cat(paste("\n model",models[j],nameRun,"failed"))
                                           err <<-TRUE
                                         })
        if(err){
          pred <- as.data.frame(rep(NA,nrow(env.var)))
          colnames(pred) = "GBM" 
        }else{
          pred <- as.data.frame(gbm::predict.gbm(mod,newdata = env.var,type = "response"))
          colnames(pred) = "GBM" 
        }
        
        if(save.obj & !(err)){
          save(mod,file=paste("ESM",nameRun,
                              colnames(env.var)[1],
                              colnames(env.var)[2],
                              models[j],"model.out",
                              sep="_"))
        }else if(nameRun == "Full"& !(err)){
          save(mod,file=paste("ESM",nameRun,
                              colnames(env.var)[1],
                              colnames(env.var)[2],
                              models[j],"model.out",
                              sep="_"))
        }
      }
      
      if(models[j] == "MAXNET"){
        cat(paste("\nMAXNET", nameRun,"\n"))
        tryCatch(expr={mod <- maxnet::maxnet(p = data.maxnet$resp,data = data.maxnet[,-1])}, error=function(e){
          cat(paste("\n model",models[j],nameRun,"failed"))
          err <<-TRUE
        })
        
        if(err){
          pred <- as.data.frame(rep(NA,nrow(env.var)))
          colnames(pred) = "MAXNET" 
        }else{
          pred <- as.data.frame(predict(mod,newdata = env.var,type = "cloglog",clamp=F))
          colnames(pred) = "MAXNET" 
        }
        
        if(save.obj & !(err)){
          save(mod,file=paste("ESM",nameRun,
                              colnames(env.var)[1],
                              colnames(env.var)[2],
                              models[j],"model.out",
                              sep="_"))
        }else if(nameRun == "Full"& !(err)){
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
  data.maxnet <- data
  data$w = w
    for(j in 1: length(models)){
      err <- FALSE
      
      if(models[j] == "GLM"){
        cat(paste("\nGLM", nameRun))
        
        if(is.null(models.options$GLM$myFormula)){
          if(models.options$GLM$test == "none"){
            
            formula <- .makeGLMFormula(env.var,
                                       models.options$GLM)
            tryCatch(expr={mod <- glm(formula = formula,
                       family = models.options$GLM$family,
                       weights = w,
                       data = data)}, error=function(e){
                         cat(paste("\n model",models[j],nameRun,"failed"))
                         err <<-TRUE
                       })
            if(err){
              pred <- as.data.frame(rep(NA,nrow(env.var)))
              colnames(pred) = "GLM" 
            }else{
              pred <- as.data.frame(predict(mod,newdata = env.var,type = "response"))
              colnames(pred) = "GLM" 
            }
            
            if(save.obj & !(err)){
              save(mod,file=paste("ESM",nameRun,
                                  colnames(env.var)[1],
                                  colnames(env.var)[2],
                                  models[j],"model.out",
                                  sep="_"))
            }else if(nameRun == "Full" & !(err)){
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
            tryCatch(expr={
              mod <- step(mod.full,
                        scope = "resp~1",
                        direction = "both",
                        trace=F)
            cat(paste0("\n\tBest Formula:",deparse(mod$formula)))}, error=function(e){
                          cat(paste("\n model",models[j],nameRun,"failed"))
                          err <<-TRUE
                        })
            
            if(err){
              pred <- as.data.frame(rep(NA,nrow(env.var)))
              colnames(pred) = "GLM" 
            }else{
              pred <- as.data.frame(predict(mod,newdata = env.var,type = "response"))
              colnames(pred) = "GLM" 
            }
            
            if(save.obj & !(err)){
              save(mod,file=paste("ESM",nameRun,
                                  colnames(env.var)[1],
                                  colnames(env.var)[2],
                                  models[j],"model.out",
                                  sep="_"))
            }else if(nameRun == "Full" & !(err)){
              save(mod,file=paste("ESM",nameRun,
                                  colnames(env.var)[1],
                                  colnames(env.var)[2],
                                  models[j],"model.out",
                                  sep="_"))
            }
          } 
          
        }else{
         tryCatch(expr={mod <- glm(formula = models.options$GLM$myFormula,
                     family = models.options$GLM$family,
                     weights = w,
                     data = data)}, error=function(e){
                       cat(paste("\n model",models[j],nameRun,"failed"))
                       err <<-TRUE
                     })
          
          if(err){
            pred <- as.data.frame(rep(NA,nrow(env.var)))
            colnames(pred) = "GLM" 
          }else{
            pred <- as.data.frame(predict(mod,newdata = env.var,type = "response"))
            colnames(pred) = "GLM" 
          }
          
          if(save.obj & !(err)){
            save(mod,file=paste("ESM",nameRun,
                                colnames(env.var)[1],
                                colnames(env.var)[2],
                                models[j],"model.out",
                                sep="_"))
          }else if(nameRun == "Full"& !(err)){
            save(mod,file=paste("ESM",nameRun,
                                colnames(env.var)[1],
                                colnames(env.var)[2],
                                models[j],"model.out",
                                sep="_"))
          }
        }
      
      }
      
      if(models[j]=="GBM"){
        cat(paste("\nGBM", nameRun,"\n"))
        formula <- .makeGLMFormula(env.var,
                                   model.option=list(type="linear"))
        tryCatch(expr = {mod <- gbm::gbm(formula = formula,
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
                        n.cores = models.options$GBM$n.cores)}, error=function(e){
                          cat(paste("\n model",models[j],nameRun,"failed"))
                          err <<-TRUE
                        })
        if(err){
          pred <- as.data.frame(rep(NA,nrow(env.var)))
          colnames(pred) = "GBM" 
        }else{
          pred <- as.data.frame(gbm::predict.gbm(mod,newdata = env.var,type = "response"))
          colnames(pred) = "GBM" 
        }
        
        if(save.obj & !(err)){
          save(mod,file=paste("ESM",nameRun,
                              colnames(env.var)[1],
                              colnames(env.var)[2],
                              models[j],"model.out",
                              sep="_"))
        }else if(nameRun == "Full"& !(err)){
          save(mod,file=paste("ESM",nameRun,
                              colnames(env.var)[1],
                              colnames(env.var)[2],
                              models[j],"model.out",
                              sep="_"))
        }
      }
      if(models[j] == "MAXNET"){
        cat(paste("\nMAXNET", nameRun,"\n"))
        tryCatch(expr={mod <- maxnet::maxnet(p = data.maxnet$resp,data = data.maxnet[,-1])}, error=function(e){
          cat(paste("\n model",models[j],nameRun,"failed"))
          err <<-TRUE
        })
        
        if(err){
          pred <- as.data.frame(rep(NA,nrow(env.var)))
          colnames(pred) = "MAXNET" 
        }else{
          pred <- as.data.frame(predict(mod,newdata = env.var,type = "cloglog",clamp=F))
          colnames(pred) = "MAXNET" 
        }
        
        if(save.obj & !(err)){
          save(mod,file=paste("ESM",nameRun,
                              colnames(env.var)[1],
                              colnames(env.var)[2],
                              models[j],"model.out",
                              sep="_"))
        }else if(nameRun == "Full"& !(err)){
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
