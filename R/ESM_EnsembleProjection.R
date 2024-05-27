#' @name ESM_Ensemble.Projection
#' @author Flavien Collart \email{flaviencollart@hotmail.com} based on the previous code written by Frank Breiner 
#' and Mirko Di Febbraro with the contributions of Olivier Broennimann and Flavien Collart
#' @title Ensemble of Small Models: Projections of ESMs
#' @description Generate the ESMs into a new space
#' @param ESM.proj The object returned by \code{\link{ESM_Projection}}.
#' @param ESM.ensembleMod The object returned by \code{\link{ESM_Ensemble.Modeling}}.
#' @param chosen.models a \code{character} allowing to only make the ensemble maps on a selection of modeling techniques. \emph{Default}: "all".
#' @param save.obj   \code{logical}. If \code{TRUE}, the map or the data.frame will be saved. Maps or data.frame wiil be stored in the 
#' folder named by sp.name. Values are multiplied by 1000 and rounded. Note that the maps are in .tif and compressed with the option 
#' COMPRESS=DEFLATE and PREDICTOR=2.
#' @return a \code{SpatRaster} or \code{data.frame}. EF = the ensemble across the modeling techniques. 
#' @references Lomba, A., L. Pellissier, C.F. Randin, J. Vicente, F. Moreira, J. Honrado and A. Guisan. 2010. Overcoming the rare species 
#' modelling paradox: A novel hierarchical framework applied to an Iberian endemic plant. \emph{Biological Conservation}, \bold{143},2647-2657.
#' 
#' Breiner F.T., A. Guisan, A. Bergamini and M.P. Nobis. 2015. Overcoming limitations of modelling rare species by using ensembles of small models. \emph{Methods in Ecology and Evolution}, \bold{6},1210-1218.
#' 
#' Breiner F.T., Nobis M.P., Bergamini A., Guisan A. 2018. Optimizing ensembles of small models for predicting the distribution of species with few occurrences. \emph{Methods in Ecology and Evolution}. \doi{10.1111/2041-210X.12957}
#' 
#' @export
#################################################################################################################################################


ESM_Ensemble.Projection <- function(ESM.proj,
                                    ESM.ensembleMod,
                                    chosen.models = "all",
                                    save.obj = TRUE){
  ## models check
  models <- ESM.proj$model.info$models
  if(chosen.models != "all"){
    if(min(chosen.models %in% models)==0){
      stop(paste("chosen.models need to be a subset of",deparse(models)))
    }else{
      models <- chosen.models
    }
  }
  
  ##Set the pathway
  iniwd <- getwd()
  on.exit(setwd(iniwd))
  setwd(paste0("ESM.output_",ESM.ensembleMod$data$sp.name))
  
  ##Take the needed objects
  weights.algo <- ESM.ensembleMod$EF.algo$weights.algo
  weights.EF <- ESM.ensembleMod$EF$weights.EF
  projections <- do.call(c,ESM.proj$projection.path)
  name.env <- ESM.proj$name.env
  
  ## Make the projections----
  if(ESM.proj$proj.type == "data.frame"){
    mods <- do.call(cbind,sapply(projections, function(x){read.table(x,sep = "\t")}))
    names.proj <- sub(paste0(name.env,"/ESM_"),"",as.character(projections),fixed = TRUE)
    names.proj <- sub(".txt","",names.proj,fixed = TRUE)
    names.proj <- gsub("_",".",names.proj,fixed = TRUE)
    colnames(mods) = names.proj
    colnames(weights.algo) = gsub("_",".",colnames(weights.algo),fixed = TRUE)
    for(j in 1: length(models)){
      
      ## Projections and weights selection
      to.select <- grep(paste0(".",models[j]),names.proj,fixed = TRUE)
      algo.proj <- mods[,to.select]
      w <- weights.algo[models[j],sub(paste0(".",models[j]),"",names.proj[to.select],fixed = TRUE)] ##Here make sure of that the order is the same between weights and projections
      w[is.na(w)] = 0
      algo.proj <- algo.proj[,w>0]
      w <- w[w>0]
      ## Make the ensemble for each algo
      EF.algo <- apply(algo.proj, 1, weighted.mean,w=w)
      if(j==1){
        EF.toMerge <- as.data.frame(EF.algo)
      }else{
        EF.toMerge <- cbind.data.frame(EF.toMerge,EF.algo)
      }
      
    }
    colnames(EF.toMerge) = models
    ## Make the ensemble between algo
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
    ## if it is a SpatRaster----
    mods <- terra::rast(as.character(projections)) ## Load the maps
    names.proj <- sub(paste0(name.env,"/ESM_"),"",as.character(projections),fixed = TRUE)
    names.proj <- sub(".tif","",names.proj,fixed = TRUE)
    names.proj <- gsub("_",".",names.proj,fixed = TRUE) ## To make the selections
    colnames(weights.algo) = gsub("_",".",colnames(weights.algo),fixed = TRUE)
    for(j in 1: length(models)){
      
      ## Projections and weights selection
      to.select <- grep(paste0(".",models[j]),names.proj,fixed = TRUE)
      algo.proj <- terra::subset(mods,to.select)
      w <- weights.algo[models[j],sub(paste0(".",models[j]),"",names.proj[to.select],fixed = TRUE)]##Here make sure of that the order is the same between weights and projections
      w[is.na(w)] = 0
      algo.proj <- terra::subset(algo.proj,w>0) ## Remove the projections where the weights are equal or lower than 0
      w <- w[w>0] ## Remove the unwanted weights
      
      ## Make the ensemble for each algo----
      EF.algo <- terra::weighted.mean(algo.proj,w=as.numeric(w))
      names(EF.algo) = models[j]
      if(j==1){
        EF.toMerge <- EF.algo 
      }else{
        EF.toMerge <- c(EF.toMerge,EF.algo)
      }
    }
    ## Make the ensemble across the algo----
    if(length(models)>1){
      EF <- terra::weighted.mean(EF.toMerge,w=weights.EF[paste0(models,".EF")])
      names(EF) = "EF"
      EF <- c(EF.toMerge,EF)
    }else{
      EF <- EF.toMerge
    }
    EF <- round(EF)
    if(save.obj){
      terra::writeRaster(EF,paste0("ESM_Ensemble_",name.env,".tif"),gdal=c("COMPRESS=DEFLATE","PREDICTOR=2"), datatype = "INT2U",overwrite=TRUE)
      
    }
  }
  return(EF)
  
}