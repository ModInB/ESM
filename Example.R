library(ecospat)
library(terra)
library(ESM)
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
sp.name = "bc"
models = c("ANN","GLM")
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
                       parallel = F,
                       n.cores = 5,
                       modeling.id = as.character(format(Sys.time(), "%s")),
                       pathToSaveObject = getwd(),
                       save.models = T,
                       save.obj = T)
my.ESM$biva.evaluations
my.ESM$biva.calibration

### Ensemble models
my.ESM_EF <- ESM_Ensemble.Modeling(my.ESM,
                                  weighting.score=c("SomersD"),
                                  threshold=0,
                                  save.obj = F)
my.ESM_EF$evaluations

### Evaluation of the ensemble models based on the pooling procedure 
eval <- ESM_Pooling.Evaluation(ESM.Mod = my.ESM,
                               ESM.ensembleMod = my.ESM_EF,
                               EachSmallModels = T)

eval$ESM.evaluations
eval$ESM.evaluations.bivariate.models

### Predictions

proj <- ESM_Projection(ESM.Mod = my.ESM,
                       new.env = env,
                       name.env = "current",
                       parallel = F,
                       n.cores = 5,
                       save.obj = F)



Ens.proj <- ESM_Ensemble.Projection(ESM.proj = proj,
                                    ESM.ensembleMod = my.ESM_EF,
                                    save.obj = T) #if TRUE the maps or the data.frame will be saved

### thresholds to produce binary maps
my.ESM_thresholds <- ESM_Threshold(my.ESM_EF)

## get the variable contributions of ESMs
ESM_Variable.Contributions(my.ESM,my.ESM_EF) 

## get the response plots of ESMs
my.ESM_responsePlot<- ESM_Response.Plot(my.ESM,
                                        my.ESM_EF,
                                        fixed.var.metric = 'mean')
