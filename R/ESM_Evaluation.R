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
.bivaEvaluation <- function(biva,resp,cv.split.table){
  evalBiva <- NULL
  for(j in 1:(ncol(cv.split.table)-1)){
    evalBiva <- rbind(evalBiva,
                      .evaluationScores(Pred = biva[!(cv.split.table[,j]),j],
                      resp = resp[!(cv.split.table[,j])]))
  }
  evalBiva <- t(evalBiva)
  evalBiva <- cbind(evalBiva,apply(evalBiva,1,mean))
  colnames(evalBiva) = colnames(cv.split.table)
  return(evalBiva)
}