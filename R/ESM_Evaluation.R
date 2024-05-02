#############################################################
## Melting pot of function and hidden function to evaluate
## models and/or check model outputs
###########################################################

## The "user-friendly" functions that were slightly modified for ecospat package to works with the new version
## check in the help section of the ecospat package

#############
## ESM.Pooling.Evaluation
## Perform the pooling evaluation
#' @export
###############


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
                                                                   -1], 1, weighted.mean, w = weights))
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
        Pred <- na.omit(Pred)
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

#############
## ESM_Threshold
## Evaluate the fit of the ensemble models and compute thresholds
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
    
    TSS <- ecospat.max.tss(Pred = Full.models[, i],Sp.occ = resp)
    TSS.th <- TSS$max.threshold
    TSS <- TSS$max.TSS
    meva <- ecospat.meva.table(Pred = Full.models[, i],
                               Sp.occ = resp,
                               th = TSS.th)
    EVAL1 <- t(as.data.frame(meva$EVALUATION_METRICS$Value[2:9]))
    colnames(EVAL1) = meva$EVALUATION_METRICS$Metric[2:9]
    
    SomersD <- AUC * 2 - 1
    boyce <- ecospat.boyce(fit = Full.models[, i],
                           obs = Full.models[resp==1, i], PEplot = FALSE)
    Boyce <- boyce$cor
    MPA1.0 <- ecospat.mpa(Full.models[resp==1, i], perc = 1)
    MPA0.95 <- ecospat.mpa(Full.models[resp==1, i], perc = 0.95)
    MPA0.90 <- ecospat.mpa(Full.models[resp==1, i], perc = 0.90)
    
    pos.F <- which(boyce$F.ratio > 1)
    neg.F <- which(boyce$F.ratio <= 1)
    if (max(neg.F) < min(pos.F)) {
      Boyce.th.max <- EVAL1$Boyce.th.min <- mean(boyce$HS[c(max(neg.F), 
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
## ESM.Variable.Contributions
## Compute the contributions of each predictor
#' @export
###############

ESM_Variable.Contributions <- function (ESM.Mod, 
                                        ESM.ensembleMod){
  
  var <- colnames(ESM.Mod$data$env.var)
  models <- ESM.Mod$model.info$models
  contrib <- data.frame(matrix(nrow = length(var), ncol = length(models), 
                               dimnames = list(var, c(models))))
  weights <- ESM.ensembleMod$EF.algo$weights.algo
  cb1 <- combn(var, 2)[1, ]
  cb2 <- combn(var, 2)[2, ]
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
## ESM_Response.Plot
## Generates the species response curve for each predictor
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
      points(proj.fixed.list[[i]][, 2] ~ proj.fixed.list[[i]][, 
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
      legend("topleft", legend = c("ensemble", models), 
             fill = c("red", ColModels[1:length(models)]), 
             box.lty = 0)
      for (mod.i in models) {
        points(proj.fixed.list[[i]][, mod.i] ~ proj.fixed.list[[i]][, 
                                                                    1], col = ColModels[which(models == mod.i)], 
               lwd = 2, type = "l")
        rug(data[, i], col = "black")
      }
    }
  }
  
  unlink(paste0("ESM.output_", ESM.Mod$data$sp.name,"/data.fixed",1:ncol(data)), recursive = TRUE)
  
  return(proj.fixed.list = proj.fixed.list)
}

########################################################################
## The dark side of the moon: the hidden functions

## Compute diverse metrics
.evaluationScores <- function (Pred, resp) 
{

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
          next()
        }
        eval <- rbind(eval,
                          .evaluationScores(Pred = biva[!(cv.split.table[,j]),ToDo],
                                            resp = resp[!(cv.split.table[,j])]))
        rownames(eval)[nrow(eval)] = ToDo
    }
    if(nrow(eval)<=1){
      next
    }
    if(validation){
      eval2 <- eval
      eval2[is.na(eval2[,"Boyce"]),"Boyce"] <- 0
      full = apply(eval2,2,mean,na.rm=T)
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


