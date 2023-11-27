library(ecospat)
library(terra)

source("R/ESM_Modeling.R")
source("R/ESM_Models.options.R")
source("R/ESM_EnsembleModeling.R")
source("R/ESM_Projection.R")
source("R/ESM_Evaluation.R")
source("R/ESM_EnsembleProjection.R")


test <- get(load("C:/Users/fcollart/Downloads/Code_and_data/BiomodData_Collembola_Entom_schoe.RData"))
resp <-test@data.species
env <- test@data.env.var
xy <- test@coord
sp.name <- "TestError"
models = c("GLM","GBM")
models.options = ESM_Models.Options(GLM=list(test="AIC",
                                             type="quadratic"))

### Calibration of simple bivariate models
my.ESM <- ESM_Modeling(resp = resp,
                       xy=xy,
                       env=env,
                       sp.name = sp.name,
                       models = models,
                       models.options = models.options,
                       prevalence = 0.5,
                       cv.method = "split-sampling", #can be split, block, custom
                       cv.rep = 2,
                       cv.ratio = 0.7,
                       cv.split.table = NULL,
                       which.biva = 34:38,
                       modeling.id = as.character(format(Sys.time(), "%s")),
                       pathToSaveObject = getwd(),
                       save.obj = TRUE)

my.ESM_EF <- ESM_Ensemble.Modeling(my.ESM,
                                   weighting.score=c("MaxTSS"),
                                   threshold=0,
                                   save.obj = TRUE)
my.ESM_EF$evaluations
# Loading test data
data(ecospat.testNiche.inv)
inv <- ecospat.testNiche.inv

# species occurrences
xy <- inv[,1:2]
resp <- inv[,11]

# env data
env <- inv[,3:5]
env <- terra::rast(system.file("extdata","ecospat.testEnv.tif",package="ecospat"))
xy <- ecospat.testData[,2:3]
resp <- ecospat.testData$Veronica_alpina

### Formating the data with the BIOMOD_FormatingData() function from the package biomod2
sp.name = "V.a"
models = c("GLM","GBM")
models.options = ESM_Models.Options(GLM=list(test="none",
                                             type="quadratic"))
prevalence = 0.5
cv.method = "split-sampling"
cv.rep = 2
cv.ratio = 0.7
cv.split.table = NULL
which.biva = NULL
modeling.id = as.character(format(Sys.time(), "%s"))
pathToSaveObject = getwd()
save.obj=T
### Calibration of simple bivariate models
my.ESM <- ESM_Modeling(resp = resp,
                       xy=xy,
                       env=env,
                       sp.name = sp.name,
                       models = models,
                       models.options = models.options,
                       prevalence = 0.5,
                       cv.method = "split-sampling", #can be split, block, custom
                       cv.rep = 10,
                       cv.ratio = 0.7,
                       cv.split.table = NULL,
                       which.biva = NULL,
                       parallel = TRUE,
                       n.cores = 5,
                       modeling.id = as.character(format(Sys.time(), "%s")),
                       pathToSaveObject = getwd(),
                       save.obj = TRUE)
my.ESM$biva.evaluations

### Ensemble models
my.ESM_EF <- ESM_Ensemble.Modeling(my.ESM,
                                  weighting.score=c("MaxTSS"),
                                  threshold=0,
                                  save.obj = TRUE)
my.ESM_EF$evaluations

### Predictions

proj <- ESM_Projection(ESM.Mod = my.ESM,
                       new.env = env,
                       name.env = "current",
                       parallel = T,
                       n.cores = 5)



Ens.proj <- ESM_Ensemble.Projection(ESM.proj = proj,
                                    ESM.ensembleMod = my.ESM_EF,
                                    save.obj = TRUE) #if TRUE the maps or the data.frame will be saved

####################################################################################



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