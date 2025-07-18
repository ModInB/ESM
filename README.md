[![Github Version](https://img.shields.io/badge/dev%20version-0.7.2-53AA93.svg)](https://github.com/ModInB/ESM)
[![Last Commit](https://img.shields.io/github/last-commit/ModInB/ESM.svg)](https://github.com/ModInB/ESM/commits/main)
[![R-CMD-check](https://github.com/ModInB/ESM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ModInB/ESM/actions/workflows/R-CMD-check.yaml)

<img src="inst/logo/ESM_icon_V2.png" align="right" height = 180/>
<div align="center">
<b>------------------------------------------------------------<br/>
<b>ESM - Ensemble of Small Models<br/>
<b>------------------------------------------------------------<br/>

</b>
</div>


### <i class="fas fa-tools"></i> Installation



- **Development version** [![v](https://img.shields.io/badge/dev%20version-0.7.2-53AA93.svg)](https://github.com/ModInB/ESM)
```R
library(devtools)
devtools::install_github("ModInB/ESM", dependencies = TRUE)
```

### Information

Functions to perform and evaluate Ensemble of small models. These functions are made to model and predict rare species distributions. Please note that it is an alpha version and thus the code is not yet stable. Please write in the issues if you find any or if you think about possible enhancements. Thanks for using the package

### List of the functions in the package

| Category      	| Function      	| Description                                  	|
|:-----------------:|:-------------------:|:-----------------:|
| Data preparation	| Bp_Sampling	| Samples background point using 4 different methods. The 4 methods are fully random or stratified in the geographic space (“rand.geo”, “strat.geo”) and in the environmental space (“rand.env”, “strat.env”)|
|		| ESM_Models.Options	| Generates a list of model parameters|
|		| get_Chesla.Clim	| Download climatic grids from CHELSA v 2.1|
|		| get_Topography	| Download topographic grids from Amatulli et al (2017)|
| Modeling	| ESM_Modeling	| Calibrates and evaluates bivariate models|
|		| ESM_Projection	| Projects each bivariate model into a geographical space |
|		| ESM_Ensemble.Modeling	| Generates and evaluates an ensemble model (called “ESM”) |
|		| ESM_Ensemble.Projection	| Generates the ESM in the geographical space |
| Evaluation	| ESM_Pooling.Evaluation	| Evaluates each bivariate model and the ESM based on the pooling method |
|		| ESM_Null.Models	| Tests the significance of the evaluation metrics of the ESM based on null models and adjust the evaluation metrics |
|		| ESM_Variable.Contributions	| Computes the contribution of each variable in the ESM |
|		| ESM_Response.Plot	| Generates species response curve for each variable |
|   | Max_MCC | Compute the maximum value of Matthew’s Correlation Coefficient  |
|		| Smooth_CBI| Computes the Smooth continuous Boyce Index (SBI) |
| Post Modeling	| ESM_Binarize	| Binarizes probability values |
|		| ESM_Generate.ODMAP	| Generates and fills ODMAP table  |
|		| ESM_Threshold	| Computes diverse threshold to binarize ESMs |
|		| ESM_Range.Shift	| Computes range changes between two projections  |
| Test Data	| ESM_Env	| Environmental SpatRaster to perform ESM with the ESM package |
| 	| ESM_Species.Env	| Species and environmental data to perform ESM with the ESM package |
| 	| ESM_Splachnum.Data	| Species occurrence data to perform ESM |
| 	| ESM_Splachnum.Env	| Environmental SpatRaster to perform ESM in Scandinavia |
| Miscellaneous tools	| Load_ESM_Modelling	| load the ESM_Modelling object|
|   | Load_ESM_Ensemble.Modelling	| load the ESM_Ensemble.Modelling object|
|   | Load_ESM_Projection	| load the ESM_Projection object|



### Citation

To cite this package, please use the following reference:

<code> <i> Collart, F., Hotermans, A., Theunissen, K., Broennimann, O., Guisan, A. 2024. ESM: Ensemble of Small Models_. R package version 0.7.2.</code> </i>

### References

  - Breiner F.T., A. Guisan, A. Bergamini and M.P. Nobis. 2015. Overcoming limitations of modelling rare species by using ensembles of small models. Methods in Ecology and Evolution, 6,1210-1218.
  - Breiner F.T., Nobis M.P., Bergamini A., Guisan A. 2018. Optimizing ensembles of small models for predicting the distribution of species with few occurrences. Methods in Ecology and Evolution. doi:10.1111/2041-210X.12957.
  - Collart, F., & Guisan, A. (2023). Small to train, small to test: Dealing with low sample size in model evaluation. Ecological Informatics. 75, 102106. doi:10.1016/j.ecoinf.2023.102106.
  - Lomba, A., L. Pellissier, C.F. Randin, J. Vicente, F. Moreira, J. Honrado and A. Guisan. 2010. Overcoming the rare species modelling paradox: A novel hierarchical framework applied to an Iberian endemic plant. Biological Conservation, 143,2647-2657.
