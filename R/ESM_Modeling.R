
#' @name ESM_Modeling
#' @author Flavien Collart \email{flaviencollart@hotmail.com} based on the previous code written by Frank Breiner 
#' and Mirko Di Febbraro with the contributions of Olivier Broennimann and Flavien Collart
#' @title Ensemble of Small Models: Calibration of Bivariate Models
#' @description Model species distribution based on the method Ensemble of Small Models (ESM) Evaluate also each bivariate models.
#' 
#' @param resp \code{numeric} of 0-1. 0 the species si absent and 1 when present.
#' @param xy \code{matrix} or \code{data.frame} containing the X and Y coordinate of the species.
#' @param env \code{matrix}, \code{data.frame} or \code{SpatRaster} of the species predictors.
#' @param sp.name \code{character}. Name of the species (To generate of ESM folder with this name).
#' @param models  \code{character} of the wanted algorithm methods. Can be c("ANN","CTA","GLM","GBM","MAXNET) or a subset of these 5 techniques.
#' @param models.options \code{NULL} or the output from \code{\link{ESM_Models.Options}}
#' @param prevalence \code{NULL} or a \code{numeric} comprised between 0-1. Prevalence value is used to build 
#' 'weighted response weights'. The default is 0.5 (weighting presences equally to the absences). 
#' If \code{NULL} each observation (presence or absence) has the same weight (independent of the number of presences and absences). 
#' Note that it is not applicable for MAXNET.
#' @param cv.method \code{character}. Either "split-sampling", "block" or "custom". "split-sampling" corresponds 
#' to a repeated split-sampling cross-validations where a percentage of presences and absences, randomly selected 
#' for each run, are used to train the model and the remaining to test it. "block" corresponds to a k-fold cross-validations 
#' but where the presences (and absences) are equally split into the k blocks. If "custom", cv.split.table should be provided.
#' @param cv.rep   \code{numeric}. Number of replicates used for the split-sampling. Only applicable when cv.method="split-sampling".
#' @param cv.ratio  \code{numeric} betweem 0 and 1.Ratio of the dataset used to trained the model. Only applicable when cv.method="split-sampling".
#' @param cv.n.blocks \code{numeric}. Number of wanted blocks (k-fold cross-validation). Only applicable when cv.method = "block.
#' @param cv.split.table a \code{matrix} or a \code{data.frame} filled with TRUE/FALSE to specify which part of data must be used for models calibration (TRUE) 
#' and for models validation (FALSE). Each column corresponds to a 'RUN' and should be named "RUNX" where X correspond to the number of the run. 
#' The last column should be filled with only TRUE and named "Full" to make a full model used for the future projection. Only applicable when cv.method="custom".
#' @param which.biva \code{numeric}. which bivariate combinations should be used for modeling. \emph{Default}: \code{NULL}, 
#' meaning that all the combinations will be made.
#' @param parallel \code{logical}. Allows or not parallel job using the function parallel::makeCluster.
#' @param n.cores \code{numeric}. Number of cores used to make the models.
#' @param modeling.id  \code{character}. the ID (=name) of modeling procedure. A random number by default.
#' @param pathToSaveObject a \code{character} of a full path to store the objects. \emph{Default}: Takes the value from getwd().
#' @param save.models \code{logical}. Allows or not to save all the bivariate models. If \code{FALSE}, only the full models will be 
#' saved to make the projections possible.
#' @param save.obj \code{logical}. Allows or not to save the final output.
#' @details  
#' \describe{
#' The basic idea of ensemble of small models (ESMs) is to model a species distribution based on small, simple models, 
#' for example all possible bivariate models (i.e. models that contain only two predictors at a time out of a larger set of predictors), 
#' and then combine all possible bivariate models into an ensemble (Lomba et al. 2010; Breiner et al. 2015).
#' 
#' The ESM set of functions could be used to build ESMs using simple bivariate models which are averagedusing weights based on model performances. 
#' They provide full functionality of the approach described in Breiner et al. (2015).
#' 
#' The argument which.biva allows to split model runs, e.g. if which.biva is 1:3, only the three first bivariate variable combinations will be modeled. 
#' This allows to run different biva splits on different computers. However, it is better not to use this option if all models are run on a single computer.
#' }
#' @return 
#' \itemize{
#' a \code{list} containing: 
#' \item{data}: a \code{list} with the object resp, xy, env.var and sp.name. env.var is = to the data supplied in the argument env. 
#' If env, was a SpatRaster, it corresponds to the extracted values of these rasters.
#' \item{model.info}: a \code{list} of the models used (models), their options (model.options), the combination of bivariate models (which.biva), 
#' the failed models (failed.mod), the modeling ID (modeling.id), the prevalence argument and, the path to the folder where are the stored the models (biva.path).
#' \item{cv.split.table}: a \code{matrix} used to train and test models. See explanation of the argument cv.split.table.
#' \item{cv.method }: a \code{character} corresponding to the used cross-validation method.
#' \item{biva.predictions}: a \code{list} of the predictions of all the runs for each bivariate models.
#' \item{biva.evaluations}: a \code{list} of the evaluation of each bivariate model runs. The evaluation of the full model correspond 
#' to the mean of all the runs. Note that if one of the run has a Boyce = NA, we will consider this has a 0 when averaging.
#' \item{biva.calibration}: a \code{list} of the calibration power of each bivariate model runs including the full model.
#' }
#' 
#' @references Lomba, A., L. Pellissier, C.F. Randin, J. Vicente, F. Moreira, J. Honrado and A. Guisan. 2010. Overcoming the rare species 
#' modelling paradox: A novel hierarchical framework applied to an Iberian endemic plant. \emph{Biological Conservation}, \bold{143},2647-2657.
#' 
#' Breiner F.T., A. Guisan, A. Bergamini and M.P. Nobis. 2015. Overcoming limitations of modelling rare species by using ensembles of small models. \emph{Methods in Ecology and Evolution}, \bold{6},1210-1218.
#' 
#' Breiner F.T., Nobis M.P., Bergamini A., Guisan A. 2018. Optimizing ensembles of small models for predicting the distribution of species with few occurrences. \emph{Methods in Ecology and Evolution}. \doi{10.1111/2041-210X.12957}
#' 
#' @seealso \code{\link{ESM_Projection}}, \code{\link{ESM_Ensemble.Modeling}},   \code{\link{ESM_Ensemble.Projection}}, 
#' \code{\link{ESM_Pooling.Evaluation}}
#' 
#' @examples \donttest{library(ecospat)
#' #Loading test data
#' data(ecospat.testNiche.inv)
#' inv <- ecospat.testNiche.inv
#' #species occurrences
#' xy <- inv[,1:2]
#' resp <- inv[,11]
#' env <- inv[,3:5]
#' ### Calibration of simple bivariate models
#' my.ESM <- ESM_Modeling(resp = resp,
#'                        xy=xy,
#'                        env=env,
#'                        sp.name = "test",
#'                        models = c("GLM"),
#'                        models.options = NULL,
#'                        prevalence = 0.5,
#'                        cv.method = "split-sampling",
#'                        cv.rep = 2,
#'                        cv.ratio = 0.7,
#'                        parallel = FALSE)
#'                        
#' # Performances of each bivariate model
#' my.ESM$biva.evaluations
#' 
#' ### Ensemble models using a weighted mean based on maxTSS
#' my.ESM_EF <- ESM_Ensemble.Modeling(my.ESM,
#'                                    weighting.score=c("MaxTSS"),
#'                                    threshold=0)
#'                                    
#' ## Performances of the ensemble across the replicates
#' ## The full model evaluation corresponds to the mean value across the replicates
#' my.ESM_EF$evaluations
#' 
#' ### Evaluation of the ensemble models based on the pooling procedure 
#' ### as recommended in Collart & Guisan (2023)
#' eval <- ESM_Pooling.Evaluation(ESM.Mod = my.ESM,
#'                                ESM.ensembleMod = my.ESM_EF,
#'                                EachSmallModels = FALSE)
#'                                
#' ## Performances of the ensemble
#' eval$ESM.evaluations
#'
#' ### Predictions of each bivariate model into a new space
#' proj <- ESM_Projection(ESM.Mod = my.ESM,
#'                        new.env = env,
#'                        name.env = "current",
#'                        parallel = FALSE)
#'
#'
#' ### Ensemble predictions
#' Ens.proj <- ESM_Ensemble.Projection(ESM.proj = proj,
#'                                    ESM.ensembleMod = my.ESM_EF,
#'                                    save.obj = TRUE)
#'
#' ### thresholds to produce binary maps
#' my.ESM_thresholds <- ESM_Threshold(my.ESM_EF)
#'
#'
#' ### Binary Projection based on max TSS of calibrated ESMs into new space                                                
#' my.ESM_EFproj_current_binary <- (Ens.proj > 
#'                                 (my.ESM_thresholds$TSS.th*1000))*1
#'
#' ### get the variable contributions of ESMs
#' ESM_Variable.Contributions(my.ESM,my.ESM_EF) 
#'
#' ### get the response plots of ESMs
#' my.ESM_responsePlot<- ESM_Response.Plot(my.ESM,
#'                                         my.ESM_EF,
#'                                         fixed.var.metric = 'mean')
#'                                         
#' #To avoid a note: DO NOT RUN 
#' unlink("ESM.output_test", recursive = TRUE)
#' }
#' @export
#### ESM_Modeling----
ESM_Modeling <- function(resp,
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
                          save.models = TRUE,
                          save.obj = TRUE){
  
  ## Check resp, XY, sp.name and prevalence----
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
  
  if(!all(sort(unique(resp)) == c(0,1))){
    stop("resp should only contain 0 and 1")
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
  
  ## Check model names ----
  
  if(any(!(models  %in% c("GLM","GBM","MAXNET","ANN", "CTA")))){
    stop("models should be = to ANN, CTA, GLM, GBM, and/or MAXNET")
  }
  
  ## Check model options----
  if(is.null(models.options)){
    models.options = ESM_Models.Options()
  }else{
    if(!is.list(models.options) | deparse(names(models.options))!= deparse(c("ANN", "CTA","GLM","GBM"))){
     stop("models.options should null or formatted via ESM_Models.Options()") 
    }
  }
  
  ## Check env and extract values if SpatRaster----
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
  

  
  ## Check split.tables and generate one if it is null----
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
    cv.split.table <- .ESM.CreatingDataSplitTable(resp = resp, 
                                                 cv.rep = cv.rep,
                                                 cv.method = cv.method,
                                                 cv.ratio = cv.ratio,
                                                 cv.n.blocks = cv.n.blocks)
  }
  
  ## Remove NAs in env.var----
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
    env.var = stats::na.omit(env.var)
  }
  # Check if presences and absences are present even after removing possible NAs----
  if(sum(resp)==0 | sum(resp==0) == 0){
    stop("presences and absences/pseudo-absences should be present in resp")
  }
  
  ### Start the modeling ----
  
  # Create a folder to store objects
  iniwd <- getwd()
  on.exit(setwd(iniwd))
  dir.create(paste0(pathToSaveObject,"/ESM.output_", sp.name,"/",modeling.id),recursive = T)
  newwd <- paste0("./ESM.output_", sp.name,"/",modeling.id)
  setwd(newwd)
  
  # Generate all the possible combination variables
  combinations <- utils::combn(colnames(env.var), 2)
  
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
                       prevalence = prevalence, save.obj = save.models)
    parallel::stopCluster(cl)
  }else{
    biva.mods <- apply(combinations, 2, .bivaModeling,resp = resp,
                       env.var = env.var, models = models,
                       models.options = models.options,
                       cv.split.table = cv.split.table,
                       prevalence = prevalence, save.obj = save.models,
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

  failed.mods <- lapply(biva.mods,.checkFailedMods)
  biva.mods.filt <- lapply(1:length(biva.mods), .PutNAsFailed, 
                           biva.mods,failed.mods)
  names(biva.mods.filt) = names(biva.mods)
  lapply(1:length(biva.mods), .PrintFailedMods, 
         biva.mods,failed.mods)
  cat("\n############### Start evaluations ###############")
  
  
  ## Evaluation----
  biva.eval <- lapply(biva.mods.filt,.bivaEvaluation,
                      resp=resp, models=models,
                      cv.split.table=cv.split.table,
                      validation = TRUE) #If the full Model failed Next
  
  biva.calib <- lapply(biva.mods.filt,.bivaEvaluation,
                       resp=resp, models=models,
                       cv.split.table=!(cv.split.table),
                       validation = FALSE)

  ## Return outputs ----
  obj <- list(data = list(resp = resp,
                          xy = xy,
                          env.var = env.var,
                          sp.name= sp.name),
              model.info = list(models = models,
                                models.options = models.options,
                                which.biva = which.biva,
                                failed.mod = failed.mods,
                                modeling.id = modeling.id,
                                prevalence = prevalence,
                                biva.path = newwd),
              cv.split.table = cv.split.table,
              cv.method = cv.method,
              biva.predictions = biva.mods.filt,
              biva.calibration = biva.calib,
              biva.evaluations = biva.eval
              )
  
  if(save.obj){
    save(obj,file=paste0("../ESM.Modeling.",modeling.id,".out"))
  }
  cat("\n##################### Done #####################")
  
  return(obj)
}

#### The Hidden Functions ----
## Functions used inside ESM_Modeling
## .ESM.CreatingDataSplitTable----
# Function to generate the argument DataSplitTable                                                                     
.ESM.CreatingDataSplitTable <- function(resp,
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
        calib.Lines[pres.Random[(i*size.blockPres+1):(size.blockPres*(i+1))],(i+1)] = TRUE
        calib.Lines[abs.Random[(i*size.blockAbs+1):(size.blockAbs*(i+1))],(i+1)] = TRUE
      }else{
        calib.Lines[pres.Random[(i*size.blockPres+1):length(pres.Random)],(i+1)] = TRUE
        calib.Lines[abs.Random[(i*size.blockAbs+1):length(abs.Random)],(i+1)] = TRUE
      }
      
    }
    calib.Lines <- cbind(calib.Lines,TRUE)
    colnames(calib.Lines) = c(paste0("RUN",1:cv.n.blocks),"Full")
  }
  
  return(calib.Lines) 
}


## .bivaModeling---- 
#allow to run the functions .doModeling to model each run of a bivariate models
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

## .doModeling ----
# generates the bivariate models

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
  data.ann <- data
  data.ann$resp <- as.factor(data.ann$resp)
    for(j in 1: length(models)){
      err <- FALSE
      if(models[j] == "ANN"){
        cat(paste("\nANN", nameRun,"\n"))
        formula <- .makeGLMFormula(env.var,
                                   model.option=list(type="linear"))
        
        tryCatch(expr={mod <- nnet::nnet(formula,data = data.ann,
                                         weights = w,
                                         size = models.options$ANN$size,
                                         decay = models.options$ANN$decay,
                                         rang = models.options$ANN$rang,
                                         maxit = models.options$ANN$maxit, 
                                         trace= FALSE)}, 
                 error=function(e){
                   cat(paste("\n model",models[j],nameRun,"failed"))
                   err <<-TRUE
                 })
        
        if(err){
          pred <- as.data.frame(rep(NA,nrow(env.var)))
          colnames(pred) = "ANN" 
        }else{
          pred <- predict(mod,newdata = env.var,type="raw")
          colnames(pred) = "ANN" 
          
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
      if(models[j] == "CTA"){
        
        cat(paste("\nCTA", nameRun,"\n"))
        formula <- .makeGLMFormula(env.var,
                                   model.option=list(type="linear"))
        
        tryCatch(expr={mod <- rpart::rpart(formula = formula,
                                           data = data,
                                           weights = w,
                                           na.action = models.options$CTA$na.action,
                                           method = models.options$CTA$method,
                                           model = models.options$CTA$model,
                                           x = models.options$CTA$x,
                                           y = models.options$CTA$y,
                                           control = models.options$CTA$control
                                           )}, 
                 error=function(e){
                   cat(paste("\n model",models[j],nameRun,"failed"))
                   err <<-TRUE
                 })
        
        if(err){
          pred <- as.data.frame(rep(NA,nrow(env.var)))
          colnames(pred) = "CTA" 
        }else{
          pred <- as.data.frame(predict(mod,newdata = env.var,type="prob")[,2])
          colnames(pred) = "CTA" 
          
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
      if(models[j] == "GLM"){
        cat(paste("\nGLM", nameRun))
        
        if(is.null(models.options$GLM$myFormula)){
          if(models.options$GLM$test == "none"){
            
            formula <- .makeGLMFormula(env.var,
                                       models.options$GLM)
            tryCatch(expr={mod <- stats::glm(formula = formula,
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
            mod.full<- stats::glm(formula = formula,
                           family = models.options$GLM$family,
                           weights = w,
                           data = data)
            tryCatch(expr={
              mod <- stats::step(mod.full,
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
         tryCatch(expr={mod <- stats::glm(formula = models.options$GLM$myFormula,
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
## .makeGLMFormula ---- 
# allows to generate the formula for GLM to allow linear, quadratic or polynomial terms (and for GBM)
.makeGLMFormula <- function(env.var = env.var,model.option){
  
  formula <- paste0("resp~", paste0(colnames(env.var),collapse = "+"))
  
  if(model.option$type == "quadratic" | model.option$type == "polynomial"){
    IsNum <- sapply(env.var,is.numeric)
    Topaste <- paste0("I(",colnames(env.var)[IsNum],"^2)",collapse="+")
    formula <- paste(formula,Topaste,sep = "+")
    if(model.option$type == "polynomial"){
      Topaste <- paste0("I(",colnames(env.var)[IsNum],"^3)",collapse="+")
      formula <- paste(formula,Topaste,sep = "+")
    }
  }
  
  return(stats::as.formula(formula))
  
}
## .checkFailedMods----
# Check Failed Mods
.checkFailedMods <- function(biva.mod){
  IsNa <- apply(biva.mod, 2, anyNA)
  IsFlat <- apply(biva.mod, 2, stats::sd) == 0
  
  Failed <- IsNa | IsFlat
  return(Failed)
}
##.PutNAsFailed ----
# Transform Failed models into NAs
.PutNAsFailed <- function(biva,biva.mods,failed.mods){
  biva.mods[[biva]][failed.mods[[biva]]] = NA
  return(biva.mods[[biva]])
}
## .PrintFailedMods----
# Print Failed Mods
.PrintFailedMods <- function(biva,biva.mods,failed.mods){
  if(sum(failed.mods[[biva]])>0){
    cat(paste("\nFailed Models for combination",names(biva.mods)[biva],":",colnames(biva.mods[[biva]])[failed.mods[[biva]]] ))
  }
}
