ESM_Ensemble.Projection <- function(ESM.proj,
                                    ESM.ensembleMod,
                                    chosen.models = "all",
                                    save.obj = TRUE){
  models <- ESM.proj$model.info$models
  if(chosen.models != "all"){
    if(min(chosen.models %in% models)==0){
      stop(paste("chosen.models need to be a subset of",deparse(models)))
    }else{
      models <- chosen.models
    }
  }
  
  iniwd <- getwd()
  on.exit(setwd(iniwd))
  setwd(paste0("ESM.output_",ESM.ensembleMod$data$sp.name))
  
  weights.algo <- ESM.ensembleMod$EF.algo$weights.algo
  weights.EF <- ESM.ensembleMod$EF$weights.EF
  projections <- ESM.proj$projection.path
  name.env <- ESM.proj$name.env
  
  if(ESM.proj$proj.type == "data.frame"){
    mods <- do.call(cbind,sapply(projections, function(x){read.table(x,sep = "\t")}))
    names.proj <- sub(paste0(name.env,"/ESM_"),"",as.character(projections),fixed = TRUE)
    names.proj <- sub(".txt","",names.proj,fixed = TRUE)
    names.proj <- gsub("_",".",names.proj,fixed = TRUE)
    colnames(mods) = names.proj
    
    for(j in 1: length(models)){
      to.select <- grep(paste0(".",models[j]),names.proj,fixed = TRUE)
      algo.proj <- mods[,to.select]
      w <- weights.algo[models[j],sub(paste0(".",models[j]),"",names.proj[to.select],fixed = TRUE)]
      algo.proj <- algo.proj[,w>0]
      w <- w[w>0]
      EF.algo <- apply(algo.proj, 1, weighted.mean,w=w)
      if(j==1){
        EF.toMerge <- as.data.frame(EF.algo)
      }else{
        EF.toMerge <- cbind.data.frame(EF.toMerge,EF.algo)
      }
      
    }
    colnames(EF.toMerge) = models
    
    if(length(models)>1){
      EF <- apply(EF.toMerge,1,weighted.mean,w=weights.EF[paste0(models,".EF")])
      EF <- cbind.data.frame(EF.toMerge,EF)
    }else{
      EF <- EF.toMerge
    }
    if(save.obj){
      write.table(EF,paste0("ESM_Ensemble_",name.env,".txt"),sep="\t")
    }
  }else{
    mods <- terra::rast(as.character(projections))
    names.proj <- sub(paste0(name.env,"/ESM_"),"",as.character(projections),fixed = TRUE)
    names.proj <- sub(".tif","",names.proj,fixed = TRUE)
    names.proj <- gsub("_",".",names.proj,fixed = TRUE)

    for(j in 1: length(models)){
      to.select <- grep(paste0(".",models[j]),names.proj,fixed = TRUE)
      algo.proj <- subset(mods,to.select)
      w <- weights.algo[models[j],sub(paste0(".",models[j]),"",names.proj[to.select],fixed = TRUE)]
      algo.proj <- subset(algo.proj,w>0)
      w <- w[w>0]
      EF.algo <- terra::weighted.mean(algo.proj,w=as.numeric(w))
      names(EF.algo) = models[j]
      if(j==1){
        EF.toMerge <- EF.algo 
      }else{
        EF.toMerge <- c(EF.toMerge,EF.algo)
      }
    }
    if(length(models)>1){
      EF <- terra::weighted.mean(EF.toMerge,w=weights.EF[paste0(models,".EF")])
      names(EF) = "EF"
      EF <- c(EF.toMerge,EF)
    }else{
      EF <- EF.toMerge
    }
    EF <- terra::round(EF)
    if(save.obj){
      writeRaster(EF,paste0("ESM_Ensemble_",name.env,".tif"),gdal=c("COMPRESS=DEFLATE","PREDICTOR=2"))
      
    }
  }
  return(EF)
  
}