library(ecospat)
# Loading test data
data(ecospat.testNiche.inv)
inv <- ecospat.testNiche.inv

# species occurrences
xy <- inv[,1:2]
sp_occ <- inv[11]

# env data
current <- inv[3:8]

### Formating the data with the BIOMOD_FormatingData() function from the package biomod2
sp <- 1

### Calibration of simple bivariate models
my.ESM <- ESM_Modeling(resp = sp_occ[,1],
                       xy=xy,
                       env=current,
                       sp.name = "test",
                       models = "GLM",
                       models.options = NULL,
                       prevalence = 0.5,
                       cv.method = "split-sampling", #can be split, block, custom
                       cv.rep = 2,
                       cv.ratio = 0.7,
                       cv.split.table = NULL,
                       weighting.score = "AUC",
                       which.biva = NULL,
                       modeling.id = as.character(format(Sys.time(), "%s")),
                       pathToSaveObject = getwd())

### Ensemble models
my.ESM_EF <- ecospat.ESM.EnsembleModeling(my.ESM,weighting.score=c("SomersD"),threshold=0)

### thresholds to produce binary maps
my.ESM_thresholds <- ecospat.ESM.threshold(my.ESM_EF)

### Evaluation of bivariate and ensemble models based on standard cross-validation
my.ESM_EF$ESM.evaluations
my.ESM_thresholds

### Evaluation of the ensemble models based on the pooling procedure 
my.ESM_evaluations <- ecospat.ESM.EnsembleEvaluation(ESM.modeling.output= my.ESM,
                                                     ESM.EnsembleModeling.output = my.ESM_EF,
                                                     metrics= c("AUC","MaxTSS"),
                                                     EachSmallModels = FALSE)
my.ESM_evaluations$ESM.evaluations

### Projection of simple bivariate models into new space 
my.ESM_proj_current<-ecospat.ESM.Projection(ESM.modeling.output=my.ESM,
                                            new.env=current)
### Projection of calibrated ESMs into new space 
my.ESM_EFproj_current <- ecospat.ESM.EnsembleProjection(ESM.prediction.output=my.ESM_proj_current,
                                                        ESM.EnsembleModeling.output=my.ESM_EF)
### Binary Projection based on max TSS of calibrated ESMs into new space                                                
my.ESM_EFproj_current_binary <- (my.ESM_EFproj_current > (my.ESM_thresholds$TSS.th*1000))*1

## get the variable contributions of ESMs
ecospat.ESM.VarContrib(my.ESM,my.ESM_EF)                                                      

## get the response plots of ESMs
my.ESM_responsePlot<-ecospat.ESM.responsePlot(my.ESM_EF,my.ESM,fixed.var.metric = 'mean')