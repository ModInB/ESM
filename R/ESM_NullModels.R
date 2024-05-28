#' @export
ESM_Null.Models <- function(ESM.Mod,
                            ESM.ensembleMod,
                            n.rep = 99,
                            pooling = FALSE,
                            hist.plot = FALSE,
                            pathToSaveObject = getwd(),
                            parallel = FALSE,
                            n.cores = 1,
                            save.obj = FALSE){
 
  ## check info
  cv.method <- ESM.Mod$cv.method
  if(cv.method == "custom"){
    stop("Null models are not implemented when custom cross-validations are used to model the species")
  }else if(cv.method == "split-sampling"){
    cv.ratio <- round(sum(ESM.Mod$cv.split.table[,1])/length(ESM.Mod$cv.split.table[,1]),1)
    cat(paste("\ncv.ratio has been estimated to be:", cv.ratio))
  }
  models <- ESM.Mod$model.info$models
  if(parallel){
    cl <- parallel::makeCluster(n.cores)
    eval.nullModel <- parallel::parSapply(cl, 1:n.rep, 
                                          .NullModelling,
                                          ESM.Mod = ESM.Mod,
                                          ESM.ensembleMod = ESM.ensembleMod,
                                          pooling = pooling,
                                          pathToSaveObject = pathToSaveObject, 
                                          save.obj = save.obj
                                          )
    parallel::stopCluster(cl)
  }else{
    eval.nullModel <- sapply(1:n.rep, .NullModelling,
                             ESM.Mod = ESM.Mod,
                             ESM.ensembleMod = ESM.ensembleMod,
                             pooling = pooling,
                             pathToSaveObject = pathToSaveObject, 
                             save.obj = save.obj 
                             )
  }
  
  if(pooling){
    ESM.pool <- ESM_Pooling.Evaluation(ESM.Mod = ESM.Mod,
                                       ESM.ensembleMod = ESM.ensembleMod)
    
    
    if(length(models)==1){
      eval.nullModel <-cbind(ESM.pool$ESM.evaluations[models,],eval.nullModel)
    }else{
      eval.nullModel <- cbind(ESM.pool$ESM.evaluations["EF",], eval.nullModel)
    }
  }else{
    if(length(models)==1){
      eval.nullModel <-cbind(ESM.ensembleMod$evaluations[paste0("Full.",models,".EF"),],eval.nullModel)
    }else{
      eval.nullModel <- cbind(ESM.ensembleMod$evaluations["Full.EF",], eval.nullModel)
    }
  }
  
  colnames(eval.nullModel) = c("Observed",paste0("NullModel",1:n.rep))
  pval <- apply(eval.nullModel >= eval.nullModel[,1],1,sum, na.rm =T)/(apply(!is.na(eval.nullModel),1,sum))
  
  if(!save.obj){
    unlink(paste0(pathToSaveObject,"/ESM.output_",ESM.Mod$data$sp.name,".Null"),recursive = T, force = T)
  }
  
  ##Plotting the results
  if(hist.plot){
    old.par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old.par))
    graphics::par(mfrow = c(2, 2))
    for(i in 1:4){
      hist.data <- graphics::hist(eval.nullModel[i,-1],xlim = c(min(eval.nullModel[i,],na.rm=T),
                                                               max(eval.nullModel[i,],na.rm=T)+0.05),
                                  main = rownames(eval.nullModel)[i],
                                  xlab = paste(rownames(eval.nullModel)[i]), "Null")
      graphics::segments(x0=eval.nullModel[i,1], x1= eval.nullModel[i,1], y0=0,y1=max(hist.data$counts),col = "red")
      graphics::points(x = eval.nullModel[i,1], y = max(hist.data$counts), col = "red",pch = 19)
    }
  }
  
  
  obj <- list(pval = pval,
              evaluations = eval.nullModel)
  return(obj)
}

.NullModelling <- function(x,
                           ESM.Mod,
                           ESM.ensembleMod,
                           pooling,
                           pathToSaveObject,
                           save.obj){
  cat(paste("\nProcessing Null Model:", x))
  # Modeling Info
  models <- ESM.Mod$model.info$models
  models.options <- ESM.Mod$model.info$models.options
  which.biva <- ESM.Mod$model.info$which.biva
  prevalence <- ESM.Mod$model.info$prevalence
  cv.method <- ESM.Mod$cv.method
  if(cv.method == "split-sampling"){
    cv.rep <- ncol(ESM.Mod$cv.split.table)-1
    cv.ratio <- round(sum(ESM.Mod$cv.split.table[,1])/length(ESM.Mod$cv.split.table[,1]),1)
    cv.n.blocks = NULL
  }else{
    cv.n.blocks <- ncol(ESM.Mod$cv.split.table)-1
    cv.rep = NULL
    cv.ratio = NULL
  }
  
  env <- ESM.Mod$data$env.var
  resp <- ESM.Mod$data$resp
  resp.random <- sample(resp, size = length(resp))
  sp.name <- ESM.Mod$data$sp.name
  xy <- ESM.Mod$data$xy
  ## Ensemble information
  weighting.score <- ESM.ensembleMod$weighting.score
  threshold  <- ESM.ensembleMod$threshold
  
  ESM.Mod.Null <- spsUtil::quiet(ESM_Modeling(resp = resp.random,
                                              xy = xy,
                                              env = env,
                                              sp.name = paste0(sp.name,".Null"),
                                              models = models,
                                              models.options = models.options,
                                              prevalence = prevalence,
                                              cv.method = cv.method,
                                              cv.rep = cv.rep,
                                              cv.n.blocks = cv.n.blocks,
                                              which.biva = which.biva,
                                              parallel = FALSE,
                                              modeling.id = as.character(paste0("NullModel",x)),
                                              pathToSaveObject = pathToSaveObject,
                                              save.models = F,
                                              save.obj = F)
                               )
  ESM.ensembleMod.Null  <-  spsUtil::quiet(ESM_Ensemble.Modeling(ESM.Mod = ESM.Mod.Null,
                                                           weighting.score = weighting.score,
                                                           threshold = threshold,
                                                           save.obj = F
                                                           )
                                     )
  if(!save.obj){
    unlink(paste0(pathToSaveObject,"/ESM.output_",sp.name,".Null/NullModel",x),recursive = T, force = T)
  }
  
  if(pooling){
    ESM.pool.NULL <- ESM_Pooling.Evaluation(ESM.Mod = ESM.Mod.Null,
                                       ESM.ensembleMod = ESM.ensembleMod.Null)
    if(length(models)==1){
      return(ESM.pool.NULL$ESM.evaluations[models,])
    }else{
      return(ESM.pool.NULL$ESM.evaluations["EF",])
    }
    
  }else{
    if(length(models)==1){
      return(ESM.ensembleMod.Null$evaluations[paste0("Full.",models,".EF"),])
    }else{
      return(ESM.ensembleMod.Null$evaluations["Full.EF",])
    }
  }
}
