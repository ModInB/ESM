.evaluationScores <- function (Pred, resp) 
{

  pred.esmPres <- Pred[resp == 1]
  pred.esmAbs <- Pred[resp == 0]
  auc.test <- dismo::evaluate(p = pred.esmPres, a = pred.esmAbs)@auc
  boyce.test <- ecospat::ecospat.boyce(c(pred.esmPres, pred.esmAbs), 
                              pred.esmPres, PEplot = F)$cor
  tss.test <- ecospat::ecospat.max.tss(Pred = Pred, Sp.occ = resp)[[2]]
  return(cbind(AUC = auc.test, 
               SomersD = (2 * auc.test - 
                            1),
               Boyce = boyce.test, 
               MaxTSS = tss.test))
  
}
.bivaEvaluation <- function(biva,
                            resp,
                            models,
                            cv.split.table){
  evalBiva <- NULL
  nameBiva <- colnames(biva)
  for(i in 1: length(models)){
    eval <- NULL
    for(j in 1:(ncol(cv.split.table)-1)){
      
        ToDo <- grep(paste0(colnames(cv.split.table)[j],
                            ".",models[i]), 
                     nameBiva, value = TRUE)
        eval <- rbind(eval,
                          .evaluationScores(Pred = biva[!(cv.split.table[,j]),ToDo],
                                            resp = resp[!(cv.split.table[,j])]))
        rownames(eval)[nrow(eval)] = ToDo
    }
    full = apply(eval,2,mean)
    eval <- rbind(eval,full)
    ToDo <- grep(paste0(colnames(cv.split.table)[j+1],
                        ".",models[i]), 
                 nameBiva, value = TRUE)
    rownames(eval)[nrow(eval)] = ToDo
    evalBiva <- rbind(evalBiva,eval)
    
  }
  return(evalBiva)
}
  
  
