#############################################################
## Melting pot of function and hidden function to evaluate
## models and/or check model outputs
###########################################################

#############
## ESM.Pooling.Evaluation
## Perform the pooling evaluation
#' @name ESM_Pooling.Evaluation
#' @title Ensemble of Small Models: Evaluation via the Pooling procedure
#' @author Flavien Collart \email{flaviencollart@hotmail.com}
#' @description 
#' This function evaluates the Ensemble of Small Models by pooling the different runs of the cross validation as in Collart et al (2021) 
#' and recommended for rare species by Collart & Guisan (2023).
#' @param ESM.Mod The object returned by \code{\link{ESM_Modeling}}.
#' @param ESM.ensembleMod The object returned by \code{\link{ESM_Ensemble.Modeling}}.
#' @param EachSmallModels \code{logical}. Should the individual bivariate models be evaluated by the pooling procedure?
#' 
#' @details 
#' Because a minimum sample size is needed to evaluate models (see Collart & Guisan, 2023), 
#' this function uses the approach from Collart et al.(2021), which consists to pool the suitability values of 
#' the hold-out data (evaluation dataset) across replicates. As the same data point (presence or absence or background point) 
#' is presumably sampled in several replicates, the suitability values for each data point is consequently averaged across 
#' replicates where they were sampled. This procedure generates a series of independent suitability values with a size approximately 
#' equal (as some data points may not have been sampled by chance in any of the \emph{n} replicates) to that of the number of data point.
#' 
#' @return 
#' a \code{list} containing:
#' \itemize{
#' \item{ESM.evaluations}{a \code{matrix} with the evaluation scores for the ESMs based on the different modelling algorithms and 
#' based on the consensus across the modelling algorithms (called here "ensemble").}
#' \item{ESM.fit}{a \code{matrix} of predicted values resulting from the pooling procedure and used to compute the evaluation scores. 
#' The column \emph{resp} is where the species occurs or not.}
#' \item{ESM.evaluations.bivariate.models}{a \code{list} containing a matrix of evaluation scores for each bivariate models 
#' (generated only if EachSmallModels = T).}
#' \item{ESM.fit.bivariate.models}{a \code{list} containing a matrix of of predicted values resulting from the pooling procedure 
#' for each bivariate models (generated only if EachSmallModels = T).}
#' }
#' 
#' @references 
#' Collart, F., & Guisan, A. (2023). Small to train, small to test: Dealing with low sample size in model evaluation. 
#' \emph{Ecological Informatics}. \bold{75}, 102106. \doi{10.1016/j.ecoinf.2023.102106}.
#' 
#' Collart, F., Hedenas, L., Broennimann, O., Guisan, A. and Vanderpoorten, A. 2021. Intraspecific differentiation: 
#' Implications for niche and distribution modelling. \emph{Journal of Biogeography}. \bold{48}, 415-426. \doi{10.1111/jbi.14009}.
#' 
#' @seealso \code{\link{ESM_Ensemble.Modeling}}, \code{\link{ESM_Modeling}}
#' @export

ESM_Pooling.Evaluation <- function (ESM.Mod, 
                                    ESM.ensembleMod, 
                                    EachSmallModels = FALSE) 
{
  
  if (!is.logical(EachSmallModels)) {
    stop("EachSmallModels should be logical")
  }
  
  resp <- ESM.Mod$data$resp
  modelling.techniques <- ESM.Mod$model.info$models
  nMod <- length(ESM.Mod$biva.predictions)
  nReplicate <- ncol(ESM.Mod$cv.split.table)-1
  fit <- ESM.ensembleMod$EF.algo$pred.EF.algo
  calib <- ESM.Mod$cv.split.table[, 1:nReplicate]
  
  PredFin <- NULL
  evalFin <- NULL
  for(d in 1:length(modelling.techniques)) {
    fitMod <- fit[, grep(modelling.techniques[d], colnames(fit))]
    fitMod <- fitMod[, -c(grep("Full", colnames(fitMod)))]
    fitMod <- cbind(resp,fitMod)
    Pred <- .ecospat.pooling(calib = calib, models.prediction = fitMod)
    if (d == 1) {
      PredFin <- cbind(PredFin, Pred)
    }else {
      PredFin <- cbind(PredFin, Pred[, -1])
    }
    colnames(PredFin)[ncol(PredFin)] = paste0("Fit_", modelling.techniques[d])
    evalInter <- .evaluationScores(Pred = Pred[,-1],resp=Pred[,1])
    evalFin <- rbind(evalFin, evalInter)
    rownames(evalFin)[nrow(evalFin)] = modelling.techniques[d]
  }
  if (length(modelling.techniques) > 1) {
    weights <- ESM.ensembleMod$EF$weights.EF
    PredEns <- cbind.data.frame(resp = PredFin[, 1], apply(PredFin[, 
                                                                   -1], 1, stats::weighted.mean, w = weights))
    PredFin <- cbind(PredFin, PredEns[, -1])
    colnames(PredFin)[ncol(PredFin)] = "Fit_ensemble"
    evalInter <- .evaluationScores(Pred = PredEns[,-1],
                                   resp = PredEns[, 1])
    evalFin <- rbind(evalFin, evalInter)
    rownames(evalFin)[nrow(evalFin)] = "ensemble"
  }
  output <- list(ESM.evaluations = evalFin, ESM.fit = PredFin)
  if (EachSmallModels) {
    evalBivaFin <- list()
    PredBivaFin <- list()
    for (i in 1:nMod) {
      evalBiva <- NULL
      PredBiva <- NULL
      IndivMod <- ESM.Mod$biva.predictions
      for (d in 1:length(modelling.techniques)) {
        models.prediction <- IndivMod[[i]]
        models.prediction <- models.prediction[, grep(modelling.techniques[d], colnames(models.prediction))]
        models.prediction <- models.prediction[, -c(grep("Full", colnames(models.prediction)))]
        models.prediction <- cbind.data.frame(resp = resp, 
                                              models.prediction)
        Pred <- .ecospat.pooling(calib = calib, models.prediction = models.prediction)
        PredBiva <- cbind(PredBiva, Pred)
        Pred <- stats::na.omit(Pred)
        colnames(PredBiva)[ncol(PredBiva)] = paste0("Fit_", 
                                                    modelling.techniques[d])
        evalInter <- .evaluationScores(Pred = Pred[,-1], 
                                       resp=Pred[,1])
        evalBiva <- rbind(evalBiva, evalInter)
        rownames(evalBiva)[nrow(evalBiva)] = modelling.techniques[d]
      }
      evalBivaFin[[i]] = evalBiva
      PredBivaFin[[i]] = PredBiva
    }
    names(evalBivaFin) = names(PredBivaFin) = names(ESM.Mod$biva.predictions)
    output$ESM.evaluations.bivariate.models = evalBivaFin
    output$ESM.fit.bivariate.models = PredBivaFin
  }
  return(output)
}

###############

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

#############
#' @name ESM_Threshold
#' @title Thresholds for Ensemble of Small Models
#' @author Flavien Collart \email{flaviencollart@hotmail.com} based on the scripts of Frank Breiner.
#' @description 
#' This function evaluates the full model which is used for projections and provides thresholds to produce binary maps.
#' @param ESM.ensembleMod The object returned by \code{\link{ESM_Ensemble.Modeling}}.
#' @details 
#' This function provides evaluation scores of the full model (no split sampling) and thresholds which can be used to convert suitability 
#' maps into binary maps. Various thresholds are provided: TSS (where sensitivity and specificity are maximised), MPA 1.0 (where all presences 
#' are prdicted positive), MPA 0.95 (where 95\% of all presences are predicted positive), MPA 0.90 (where 90\% of all presences are predicted positive), 
#' Boyce.th.min (the lowest suitability value where the predicted/expected ratio is >1) and Boyce.th.max (the highest suitability value where the 
#' predicted/expected ratio is =1). 
#' @return 
#' A \code{data.frame} with evaluation scores of the fit and thresholds.
#' @references 
#' Hirzel, Alexandre H., et al. Evaluating the ability of habitat suitability models to predict species presences. 
#' \emph{Ecological modelling}, \bold{199.2} (2006): 142-152.
#' 
#' Engler, Robin, Antoine Guisan, and Luca Rechsteiner. An improved approach for predicting the distribution of rare and endangered species 
#' from occurrence and pseudo-absence data. \emph{Journal of applied ecology}, \bold{41.2} (2004): 263-274.
#' 
#' Fielding, Alan H., and John F. Bell. A review of methods for the assessment of prediction errors in conservation presence/absence models. 
#' \emph{Environmental conservation}, \bold{24.1} (1997): 38-49.
#' 
#' 
#' @export
###############

ESM_Threshold <- function (ESM.ensembleMod) 
{
  models = ESM.ensembleMod$model.info$models 
  resp <-  ESM.ensembleMod$data$resp
  
  if(length(models)>1){
    Full.models <- cbind(ESM.ensembleMod$EF.algo$pred.EF.algo,ESM.ensembleMod$EF$pred.EF)
  }else{
    Full.models <- ESM.ensembleMod$EF$pred.EF
  }
  
  Full.models <- as.data.frame(Full.models[,grep("Full",colnames(Full.models))])
  if(length(models)==1){
    colnames(Full.models) = models
  }
  
  EVAL <- NULL
  for (i in 1:ncol(Full.models)) {
    
    DATA <- cbind(PlotID= 1:nrow(Full.models),
                  resp.var = resp, 
                  Full.models[, i])
    AUC <- PresenceAbsence::auc(DATA, st.dev = FALSE,which.model = 1, na.rm = TRUE)
    
    TSS <- ecospat::ecospat.max.tss(Pred = Full.models[, i],Sp.occ = resp)
    TSS.th <- TSS$max.threshold
    TSS <- TSS$max.TSS
    meva <- ecospat::ecospat.meva.table(Pred = Full.models[, i],
                               Sp.occ = resp,
                               th = TSS.th)
    EVAL1 <- t(as.data.frame(meva$EVALUATION_METRICS$Value[2:9]))
    colnames(EVAL1) = meva$EVALUATION_METRICS$Metric[2:9]
    
    SomersD <- AUC * 2 - 1
    boyce <- ecospat::ecospat.boyce(fit = Full.models[, i],
                           obs = Full.models[resp==1, i], PEplot = FALSE)
    Boyce <- boyce$cor
    MPA1.0 <- ecospat::ecospat.mpa(Full.models[resp==1, i], perc = 1)
    MPA0.95 <- ecospat::ecospat.mpa(Full.models[resp==1, i], perc = 0.95)
    MPA0.90 <- ecospat::ecospat.mpa(Full.models[resp==1, i], perc = 0.90)
    
    pos.F <- which(boyce$F.ratio > 1)
    neg.F <- which(boyce$F.ratio <= 1)
    if (max(neg.F) < min(pos.F)) {
      Boyce.th.max <- Boyce.th.min <- mean(boyce$HS[c(max(neg.F), 
                                                            min(pos.F))])
    }else {
      Boyce.th.max <- mean(boyce$HS[c(max(neg.F), 
                                      max(neg.F) + 1)])
      Boyce.th.min <- mean(boyce$HS[c(min(pos.F), 
                                      min(pos.F) - 1)])
    }
    
    EVAL1 <- cbind.data.frame(model = colnames(Full.models)[i], Boyce.th.min, Boyce.th.max,
                              MPA1.0, MPA0.95, MPA0.90, TSS.th, AUC, Boyce, SomersD, TSS, EVAL1)
    rownames(EVAL1) = NULL
    
    EVAL <- rbind(EVAL, EVAL1)
  }
  return(EVAL)
}


#############
#' @name ESM_Variable.Contributions
#' @title Variable contribution in ESMs
#' @author Olivier Broennimann \email{Olivier.Broennimann@unil.ch} with contributions of Heidi Mod and Daniel Scherrer 
#' and the adaptations of Flavien Collart \email{flaviencollart@hotmail.com}
#' @description 
#' This function evaluates the full model which is used for projections and provides thresholds to produce binary maps.
#' 
#' @param ESM.Mod The object returned by \code{\link{ESM_Modeling}}.
#' @param ESM.ensembleMod The object returned by \code{\link{ESM_Ensemble.Modeling}}.
#' 
#' @details Calculates the ratio between sum of weights of bivariate models where a focal variable was used and sum of weights of bivariate 
#' models where the focal variable was not used. The ratio is corrected for the number of models with or without the focal variable. 
#' This ratio gives an indication on the proportional contribution of the variable in the final ensemble model. A value of higher than 1 
#' indicate that the focal variable has a higher contribution than average.
#' In the case of multiple methods (e.g., GLM,...), the contributions are counted per method. For ensemble model, the contributions 
#' are then weighted means (based on the weighting score as chosen in ecospat.ESM.EnsembleModeling of single methods.
#' 
#' @return 
#' a \code{dataframe} with contribution values (i.e. proportional contribution) by variable and model
#' @seealso \code{\link{ESM_Modeling}}
#' @export
###############

ESM_Variable.Contributions <- function (ESM.Mod, 
                                        ESM.ensembleMod){
  
  var <- colnames(ESM.Mod$data$env.var)
  models <- ESM.Mod$model.info$models
  contrib <- data.frame(matrix(nrow = length(var), ncol = length(models), 
                               dimnames = list(var, c(models))))
  weights <- ESM.ensembleMod$EF.algo$weights.algo
  cb1 <- utils::combn(var, 2)[1, ]
  cb2 <- utils::combn(var, 2)[2, ]
  for (m in models) {
    for (v in var) {
      pos_models <- rownames(weights) == m
      pos_cb <- c(which(cb1 == v), which(cb2 == v))
      
      contrib[which(var == v), which(names(contrib) == m)] <- 
        sum(weights[m,pos_cb])/(sum(weights[m,setdiff(1:length(cb1), pos_cb)])) * 
        length(setdiff(1:length(cb1), pos_cb))/length(pos_cb)
    }
  }
  if (length(models) > 1) {
    EF <- matrixStats::rowWeightedMeans(x = data.matrix(contrib[, 
                                                                models]), w =ESM.ensembleMod$EF$weights.EF , na.rm = TRUE)
    contrib <- cbind(contrib, EF)
  }
  return(contrib)
}

#############
#' @name ESM_Response.Plot
#' @title Produce response plots for ESMs
#' @author Flavien Collart \email{flaviencollart@hotmail.com} from the code of Frank Breiner.
#' @description 
#' This function creates response plots (evaluation strips) for Ensembles of Small Models (ESMs).
#' 
#' @param ESM.Mod The object returned by \code{\link{ESM_Modeling}}.
#' @param ESM.ensembleMod The object returned by \code{\link{ESM_Ensemble.Modeling}}.
#' @param fixed.var.metric Either 'median' (\emph{Default}), 'mean', 'min' or 'max' specifying the statistic used 
#' to fix as constant the remaining variables when the predicted response is estimated for one of the variables.
#' 
#' @details 
#' This function plots the response curves of a model for each variable, while keeping the remianing variables constant. 
#' This is an adaptation of the Evaluation Strip method proposed by Elith et al.(2005).
#' @return 
#' A plot of the response curves is produced (red line Ensemble, other colour lines are for single algorithms) and a \code{list} with the output is provided.
#' 
#' @references 
#' Elith, J., Ferrier, S., Huettmann, FALSE. & Leathwick, J. R. 2005 The evaluation strip: A new and robust method for plotting predicted 
#' responses from species distribution models. Ecological Modelling 186, 280-289.
#' 
#' @seealso \code{\link{ESM_Modeling}}
#' @export
###############

ESM_Response.Plot <- function (ESM.Mod, 
                               ESM.ensembleMod, 
                               fixed.var.metric = "median"){
  
  models <- ESM.Mod$model.info$models
  weights <- ESM.ensembleMod$EF.algo$weights.algo
  weights.EF <- ESM.ensembleMod$EF$weights.EF
  resp <- ESM.Mod$data$resp
  data <- ESM.Mod$data$env.var
  min.data <- apply(data, 2, min)
  max.data <- apply(data, 2, max)
  data.fixed.all <- as.data.frame(t(apply(data, 2, fixed.var.metric)))
  data.fixed.all <- data.fixed.all[rep(1, each = 1000), ]
  var.names <- colnames(data)
  proj.fixed.list <- list()
  for (i in 1:ncol(data)) {
    data.fixed <- data.fixed.all
    data.fixed[, i] <- seq(min.data[i], max.data[i], (max.data[i] - 
                                                        min.data[i])/999)
    proj.fixed <- ESM_Projection(ESM.Mod,
                                 new.env = data.fixed,
                                 name.env = paste0("data.fixed",i))
    proj.fixed.list[[i]] <- ESM_Ensemble.Projection(ESM.proj = proj.fixed,
                                                    ESM.ensembleMod = ESM.ensembleMod,
                                                    save.obj = FALSE)
    proj.fixed.list[[i]] <- cbind(data.fixed[, i], proj.fixed.list[[i]])
  }
  names(proj.fixed.list) <- colnames(data)
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old.par))
  xs <- floor(sqrt(length(var.names)))
  ys <- ceiling(length(var.names)/xs)
  graphics::par(mfrow = c(xs, ys))
  ColModels <- c("#222E50", "#007991", "#439A86", "#BCD8C1", 
                 "#E9D985", "#FF6978", "#6D435A", "#352D39", "#6E8894", 
                 "#FA7921", "#FE9920")
  for (i in 1:ncol(data)) {
    if (length(models) == 1) {
      plot(proj.fixed.list[[i]][, 2] ~ proj.fixed.list[[i]][, 
                                                            1], xlab = "", main = names(proj.fixed.list)[i], 
           ylab = "predicted value", ylim = c(min(sapply(proj.fixed.list, 
                                                         function(x) {
                                                           min(x[, -1])
                                                         })), max(sapply(proj.fixed.list, function(x) {
                                                           max(x[, -1])
                                                         }))), type = "n", las = TRUE)
      graphics::points(proj.fixed.list[[i]][, 2] ~ proj.fixed.list[[i]][, 
                                                              1], xlab = names(proj.fixed.list)[i], col = "red", 
             lwd = 2, type = "l")
    }
    else {
      plot(proj.fixed.list[[i]]$EF ~ proj.fixed.list[[i]][, 
                                                          1], xlab = "", main = names(proj.fixed.list)[i], 
           ylab = "predicted value", ylim = c(min(sapply(proj.fixed.list, 
                                                         function(x) {
                                                           min(x[, -1])
                                                         })), max(sapply(proj.fixed.list, function(x) {
                                                           max(x[, -1])
                                                         }))), type = "l", lwd=2, col = "red")
      graphics::legend("topleft", legend = c("ensemble", models), 
             fill = c("red", ColModels[1:length(models)]), 
             box.lty = 0)
      for (mod.i in models) {
        graphics::points(proj.fixed.list[[i]][, mod.i] ~ proj.fixed.list[[i]][, 
                                                                    1], col = ColModels[which(models == mod.i)], 
               lwd = 2, type = "l")
        graphics::rug(data[, i], col = "black")
      }
    }
  }
  
  unlink(paste0("ESM.output_", ESM.Mod$data$sp.name,"/data.fixed",1:ncol(data)), recursive = TRUE)
  
  return(proj.fixed.list = proj.fixed.list)
}



########################################################################
## The dark side of the moon: the hidden functions

## Compute diverse metrics
.evaluationScores <- function (Pred, resp){

  pred.esmPres <- Pred[resp == 1]
  pred.esmAbs <- Pred[resp == 0]
  auc.test <- PresenceAbsence::auc(DATA=cbind(ID=1:length(Pred),
                                              resp=resp,
                                              Pred), 
                                   st.dev = FALSE,which.model = 1, 
                                   na.rm = TRUE)
  boyce.test <- ecospat::ecospat.boyce(c(pred.esmPres, pred.esmAbs), 
                              pred.esmPres, PEplot = F)$cor
  tss.test <- ecospat::ecospat.max.tss(Pred = Pred, Sp.occ = resp)[[2]]
  return(cbind(AUC = auc.test, 
               SomersD = (2 * auc.test - 
                            1),
               Boyce = boyce.test, 
               MaxTSS = tss.test))
  
}
## Allow to evaluate each runs of a bivariate model 
.bivaEvaluation <- function(biva,
                            resp,
                            models,
                            cv.split.table, 
                            validation = TRUE){
  evalBiva <- NULL
  nameBiva <- colnames(biva)
  for(i in 1: length(models)){
    eval <- NULL
    
    for(j in 1:(ncol(cv.split.table)-1)){
      
        ToDo <- grep(paste0(colnames(cv.split.table)[j],
                            ".",models[i]), 
                     nameBiva, value = TRUE)
        
        if(anyNA(biva[!(cv.split.table[,j]),ToDo])){
          if(j ==1){
            eval <- NA
          }else{
            scores <- NA
            eval <- rbind(eval,scores)
          }
          
        }else{
          scores <- .evaluationScores(Pred = biva[!(cv.split.table[,j]),ToDo],
                                      resp = resp[!(cv.split.table[,j])])
          eval <- rbind(eval,scores)
          }
        
        rownames(eval)[nrow(eval)] = ToDo
    }
    if(validation){
      if(anyNA(biva[,paste0("Full.",models[i])])){
        full = NA
      }else{
        eval2 <- eval
        eval2[is.na(eval2[,"Boyce"]),"Boyce"] <- 0
        full = apply(eval2,2,mean,na.rm=T)
      }
      
      
    }else{
      ToDo <- grep(paste0(colnames(cv.split.table)[ncol(cv.split.table)],
                          ".",models[i]), 
                   nameBiva, value = TRUE)
      
      if(anyNA(biva[!(cv.split.table[,ncol(cv.split.table)]),ToDo])){
        full = NA
      }else{
        full =  .evaluationScores(Pred = biva[!(cv.split.table[,ncol(cv.split.table)]),ToDo],
                                  resp = resp[!(cv.split.table[,ncol(cv.split.table)])])
      }
      
    }
    
    eval <- rbind(eval,full)
    ToDo <- grep(paste0(colnames(cv.split.table)[j+1],
                        ".",models[i]), 
                 nameBiva, value = TRUE)
    rownames(eval)[nrow(eval)] = ToDo
    evalBiva <- rbind(evalBiva,eval)
  }
  return(evalBiva)
}
## Perform the pooling 
.ecospat.pooling<-function (calib, models.prediction) 
{
  Pred <- NULL
  for (k in 1:nrow(calib)) {
    if (sum(!calib[k, ]) != 0) {
      valStock <- cbind(models.prediction[k, 1], mean(as.numeric(models.prediction[k, 
                                                                                   (which(!calib[k, ]) + 1)]), na.rm = T))
      colnames(valStock) = c("resp", "meanESM")
      Pred <- rbind(Pred, valStock)
    }
  }
  return(Pred)
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

