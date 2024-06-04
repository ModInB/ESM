#' @name ESM_Null.Models
#' @author Flavien Collart \email{flaviencollart@hotmail.com}
#' @title Ensemble of Small Models: Evaluation using Null Models
#' @description Performed null models to test the significativity of the ESM evaluatuion as recommended in Collart & Guisan (2023) and performed in 
#' van Proosdij et al (2016). 
#' It consists of running a series of null models, randomly shuffling the presences and absences, thus generating ‘fake’ occurrences.
#' These models are performed folowing the same methodology as for the ESM (same model parameters, same number of cross-validations and 
#' same threshold for the ensemble). The pooling evaluation can also be combined with these null models with the option pooling = TRUE 
#' (as recommended in Collart & Guisan, 2023). In addition to the pvalue, the observed evaluation metrics are then reajusted using the 
#' value of the evaluation metric at a certain quantile following the methology developped in Verdon et al (2024). 
#' 
#' @param  ESM.Mod The object returned by \code{ESM_Modeling}.
#' @param ESM.ensembleMod The object returned by \code{\link{ESM_Ensemble.Modeling}}.
#' @param n.rep \code{integer}. The number of null models. \emph{Default: 99}.
#' @param pooling \code{logical}. Should the pooling evaluation be performed? \emph{Default: FALSE}.
#' @param hist.plot \code{logical}. Should an histogram representing the null distribution be performed? \emph{Default: FALSE}.
#' @param quant \code{numeric}. The quantile used to recalibrated the evaluation metrics (between 0 and 1). \emph{Default: 0.95}.
#' @param parallel \code{logical}. Allows or not parallel job using the function parallel::makeCluster. \emph{Default: FALSE} 
#' but It is highly recommended to perform the null models in parallel as it takes quite a lot of time.
#' @param n.cores \code{integer}. Number of cores used to make the models.
#' @param pathToSaveObject The path to save the objects (temporally if save.obj = FALSE). \emph{Default: working directory}
#' @param save.obj \code{logical}. Should all the null models be kept? \emph{Default: FALSE}.
#' 
#' @return \itemize{
#' a \code{list} containing: 
#' \item{pval}: \code{numeric}. The pvalue computed for each metric.
#' \item{adj.evaluation}: \code{numeric}. Adjusted evaluation metrics for the ESM based on a certain quantile from the null distribution.
#' \item{evaluations}: \code{matrix}. Observed and null values for each metrics.
#' }
#' 
#' @references 
#' Collart, F., & Guisan, A. 2023. Small to train, small to test: Dealing with low sample size in model evaluation. 
#' \emph{Ecological Informatics}. \bold{75}, 102106. \doi{10.1016/j.ecoinf.2023.102106}.
#' 
#' van Proosdij, A.S.J., Sosef, M.S.M., Wieringa, J.J. and Raes, N. 2016. Minimum required number of specimen records to 
#' develop accurate species distribution models. \emph{Ecography}. \bold{39}, 542-552. \doi{10.1111/ecog.01509}.
#' 
#' Verdon, V., Malard, L., Collart, F., Adde, A., Yashiro, E., Lara Pandi, E., Mod, H., Singer, D., Niculita-Hirzel, H., 
#' Guex, N. and Guisan, A. 2024. Can we accurately predict the distribution of soil microorganism presence and relative abundance? 
#' \emph{Ecography}. e07086. \doi{10.1111/ecog.07086}.
#' 
#' @export
ESM_Null.Models <- function(ESM.Mod,
                            ESM.ensembleMod,
                            n.rep = 99,
                            quant = 95,
                            pooling = FALSE,
                            hist.plot = FALSE,
                            parallel = FALSE,
                            n.cores = 1,
                            pathToSaveObject = getwd(),
                            save.obj = FALSE){
 
  ## check info
  cv.method <- ESM.Mod$cv.method
  if(cv.method == "custom"){
    stop("Null models are not implemented when custom cross-validations are used to model the species")
  }else if(cv.method == "split-sampling"){
    cv.ratio <- round(sum(ESM.Mod$cv.split.table[,1])/length(ESM.Mod$cv.split.table[,1]),1)
    cat(paste("\ncv.ratio has been estimated to be:", cv.ratio))
  }
  
  if((n.rep%%1) !=0){
    stop("nrep should be an integer")
  }
  if(!is.logical(pooling)){
    stop("pooling should be a logical")
  }
  if(!is.logical(hist.plot)){
    stop("hist.plot should be a logical")
  }
  if(!is.logical(parallel)){
    stop("parallel should be a logical")
  }
  if(!is.logical(save.obj)){
    stop("save.obj should be a logical")
  }
  if(parallel & (n.cores%%1) !=0){
    stop("ncores should be an integer")
  }
  
  ##Create a check for quantile
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
  pval <- apply(eval.nullModel >= eval.nullModel[,1],1, sum, na.rm =T)/(apply(!is.na(eval.nullModel),1,sum))
  
  val.quantile <- apply(eval.nullModel[,-1],1,stats::quantile, probs = quant, na.rm = T) 
  adj.evaluation <- (eval.nullModel[,1] - val.quantile)/ (1-val.quantile)
  names(adj.evaluation) = paste0("Adjusted_",names(adj.evaluation))
  
  pval <- pval[-2]
  adj.evaluation <- adj.evaluation[-2]
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
              adj.evaluation = adj.evaluation,
              evaluations = eval.nullModel[-2,])
  return(obj)
}

### Hidden function To perform individiual Null Models ----

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
