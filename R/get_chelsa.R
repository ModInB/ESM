#' @name get_chelsa
#' @author Flavien Collart \email{flaviencollart@hotmail.com}
#' @title Download chelsa data
#' @description
#' Download and crop climatic variables at different time period from CHELSA v.2.1.
#' 
#' @param var.names \code{character}. Vector containing the variable names to download. See details for more information.
#' @param time \code{character}. The time period(s) wanted. Must be: '1981-2010','2011-2040','2041-2070', and/or'2071-2100'.
#' \emph{Default: '1981-2010'}
#' @param gcm \code{character}. A vector containing the Global Circulation Models for which you want to download the datasets. 
#' Only needed when future time periods are provided. Must be: 'gfdl-esm4','ipsl-cm6a-lr','mpi-esm1-2-hr','mri-esm2-0' and/or,'ukesm1-0-ll'.
#' \emph{Default: NULL}.
#' @param ssp \code{character}. A vector containing the Shared Socio-economic Pathways for which you want to download the datasets. 
#' Only needed when future time periods are provided. Must be: '126', '370' and/or '585'.\emph{Default: NULL}.
#' @param path \code{character}.The path to store the climatic grids.
#' @param extent The extent to which you want to crop your climatic datasets. extent must be either a 'SpatRaster', 'SpatVector', or'SpatExtent'.
#' @param mask \code{logical}. Do you want to mask your maps? \emph{Note that when TRUE, extent must be a 'SpatRaster' or a 'SpatVector'.}
#' @param compress \code{logical}. When extent is provided, do you want to compress the climatic grids? if TRUE, the compression
#' will be: "COMPRESS=DEFLATE", "PREDICTOR=2" and "ZLEVEL=6". \emph{Default: TRUE}.
#' @param timeout \code{integer}. The maximum downloading time allowed for each map in seconds. \emph{Default: 300L}.
#' #' @details  
#' \describe{}
#' @return 
#' \code{character}. A vector containing the path to the downloaded files.
#' @references
#' Karger, D.N., Conrad, O., Böhner, J., Kawohl, T., Kreft, H., Soria-Auza, R.W., Zimmermann, N.E., Linder, P., Kessler, M. (2017). Climatologies at high resolution for the Earth land surface areas. \emph{Scientific Data}. \bold{4}, 170122. \doi{10.1038/sdata.2017.122}.
#' 
#' Karger D.N., Conrad, O., Böhner, J., Kawohl, T., Kreft, H., Soria-Auza, R.W., Zimmermann, N.E, Linder, H.P., Kessler, M. (2018): Data from: Climatologies at high resolution for the earth’s land surface areas. \emph{EnviDat}. \doi{10.16904/envidat.228.v2.1}.
#' 
#' @examples  \donttest{
#' # get_chelsa_clim("bio1")
#' }
#' @export 

get_chelsa_clim <- function(var.names,
                            time = "1981-2010",
                            gcm = NULL,
                            ssp = NULL,
                            path = getwd(),
                            extent = NULL,
                            mask = FALSE, 
                            compress = TRUE,
                            timeout = 300L){
  
  #### Vector of variable and scenarios available in CHELSA ####
  
  var.full <- c(paste0("bio",1:19),
                     "ai",
                     "fcf","fgd","gdd0","gdd10",
                     "gdd5","gddlgd0","gddlgd5",
                     "gddlgd10","gdgfgd0","gdgfgd10",
                     "gdgfgd5","gsl","gsp","gst",
                     paste0("kg",0:5),
                     "lgd","ngd0","ngd10","ngd5",
                     "npp","scd","swb","swe")
  var.pres <- c("pet_perman_max","pet_perman_mean",
                "pet_perman_min","pet_perman_range",
                "sfcWind_max","sfcWind_mean",
                "sfcWind_min","sfcWind_range",
                "vpd_max","vpd_mean","vpd_min",
                "vpd_range",
                "hurs_max","hurs_mean","hurs_min",
                "hurs_range","clt_max","clt_mean",
                "clt_min", "clt_range",
                "cmi_max","cmi_mean",
                "cmi_min", "cmi_range",
                "rsds_1981-2010_mean",
                "rsds_1981-2010_min",
                "rsds_1981-2010_max",
                "rsds_1981-2010_range")
  
  ssp.chelsa <- c("126","370","585")
  gcm.chelsa <- c("gfdl-esm4","ipsl-cm6a-lr",
                  "mpi-esm1-2-hr","mri-esm2-0",
                  "ukesm1-0-ll")
  FilePath <- c() #File path to return at the end
  #### Check time and scenarios ####
  
  if(!all(time %in% c("1981-2010","2011-2040","2041-2070","2071-2100") )){
    stop("all elements in time must be '1981-2010','2011-2040','2041-2070','2071-2100'")
  }
  if(any(time %in% c("2011-2040","2041-2070","2071-2100"))){
    if(is.null(ssp) | is.null(gcm)){
      stop("when time is a future period, ssp and gcm must be provided")
    }else{
      ssp <- as.character(ssp)
      if(any(!(ssp %in% ssp.chelsa))){
        stop("ssp must be '126', '370' and/or '585'.")
      }
      gcm <- tolower(gcm)
      if(any(!(gcm %in% gcm.chelsa))){
        stop("gcm must be 'gfdl-esm4','ipsl-cm6a-lr','mpi-esm1-2-hr','mri-esm2-0' and/or,'ukesm1-0-ll'.")
      }
    }
  }
  
  #### Check variable names ####
  
  if(length(var.names)==1){
    if(var.names=="all"){
      var.down = c(var.full, var.pres)
      if(length(time)>1 | time != "1981-2010"){
        stop("When var.names = 'all', time must be only '1981-2010'. Please try var.names='all_future' instead.")
      }
    }else if(var.names == "bioclim"){
      var.down = paste0("bio",1:19)
    }else if(var.names == "all_fut"){
      var.down = var.full
    }else{
      if(length(time)>1){
        if(any(var.names %in% var.pres) | !all(var.names %in% var.full)){
          stop("When time contains a future period, var.names must be a group of variables available for future time period.")
        }
      }else{
        if(time == "1981-2010"){
          if(!all(var.names %in% c(var.full,var.pres))){
            stop("all elements in var.names must be present in extended bioclimate or must be either all, all_future ,or bioclim.")
          }
        }else{
          if(any(var.names %in% var.pres)){
            stop("When time contains a future period, var.names must be a group of variables available for future time period.")
          }
        }

      }
      
      var.down = var.names 
    }
  }else{
    if(length(time)>1){
      if(any(var.names %in% var.pres)){
        stop("When time contains a future period, var.names should be a group of variables available for future time period.")
      }
    }else{
      if(time == "1981-2010"){
        if(!all(var.names %in% c(var.full,var.pres))){
          stop("all elements in var.names must be present in ex bioclimate or must be either all, all_future ,or bioclim.")
        }
      }else{
        if(any(var.names %in% var.pres)){
          stop("When time contains a future period, var.names should be a group of variables available for future time period.")
        }
      }
      
    }
    
    var.down = var.names 
  }
  
  #### Check map transformation objects ####
  
  if(!is.null(extent)){
    if(!(inherits(extent,c("SpatRaster", "SpatVector","SpatExtent") ))){
      stop("extent must be either a 'SpatRaster', 'SpatVector', or'SpatExtent'.")
      if(mask & inherits(extent,c("SpatExtent") )){
        stop("when mask=TRUE, extent must be a 'SpatRaster' or a 'SpatVector'.")
      }
    }
    if(!is.logical(compress)){
      stop("compress must be a logical")
    }
  }
  
  #### Perform downloading and map cropping ####
  if(!is.integer(timeout)){
    stop("timeout must be an integer.")
  }
  options(timeout = max(timeout, getOption("timeout")))
  for(j in 1:length(time)){
    for(i in 1:length(var.down)){
      if(time[j] =="1981-2010"){
        download.file(url=paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_",
                                 var.down[i],"_1981-2010_V.2.1.tif"),
                      destfile=paste0(path,"/",var.down[i],"_CHELSA_1981-2010_V.2.1.tif"),
                      mode="wb")
        if(!is.null(extent)){
          map <- terra::rast(paste0(path,"/",var.down[i],"_CHELSA_1981-2010_V.2.1.tif"))
          map <- terra::crop(map,extent,mask=mask)
          if(compress){
            terra::writeRaster(map, 
                               paste0(path,"/",var.down[i],"_CHELSA_1981-2010_V.2.1.tif"), 
                               overwrite = TRUE,
                               gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6"))
          }else{
            terra::writeRaster(map, 
                               paste0(path,"/",var.down[i],"_CHELSA_1981-2010_V.2.1.tif"), 
                               overwrite = TRUE)
            }
        }
        FilePath <- c(FilePath, paste0(path,"/",var.down[i],"_CHELSA_1981-2010_V.2.1.tif")) 
        }else{
          for(h in 1:length(gcm)){
            for(g in 1:length(ssp)){
              download.file(url=paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/",
                                     time[j],"/",toupper(gcm[h]),"/ssp",ssp[g],
                                     "/CHELSA_",var.down[i],"_",time[j],"_",gcm[h],"_ssp",ssp[g],"_V.2.1.tif"),
                          destfile=paste0(path,"/CHELSA_",var.down[i],"_",time[j],"_",gcm[h],"_ssp",ssp[g],"_V.2.1.tif"),
                          mode="wb")
              if(!is.null(extent)){
                map <- terra::rast(paste0(path,"/CHELSA_",var.down[i],"_",time[j],"_",gcm[h],"_ssp",ssp[g],"_V.2.1.tif"))
                map <- terra::crop(map,extent,mask=mask)
                if(compress){
                  terra::writeRaster(map, 
                                     paste0(path,"/CHELSA_",var.down[i],"_",time[j],"_",gcm[h],"_ssp",ssp[g],"_V.2.1.tif"), 
                                     overwrite = TRUE,
                                     gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6"))
                }else{
                  terra::writeRaster(map, 
                                     paste0(path,"/CHELSA_",var.down[i],"_",time[j],"_",gcm[h],"_ssp",ssp[g],"_V.2.1.tif"), 
                                     overwrite = TRUE)
                }
                gc()
                
              }
              FilePath <- c(FilePath, paste0(path,"/CHELSA_",var.down[i],"_",time[j],"_",gcm[h],"_ssp",ssp[g],"_V.2.1.tif")) 
              
              }
        }
      }

      
        gc()
    }
  }
  return(FilePath)
}

#' @name get_topography
#' @author Flavien Collart \email{flaviencollart@hotmail.com}
#' @title Download topographic data
#' @description
#' Download and crop topographic variables from Amatulli et al (2017) using GMTED elevation data
#' 
#' @param var.names \code{character}. Vector containing the variable name(s) to download. Can be either 'all' to download all variables or a subset of the variables  
#' avaible in Amatulli et al (2017). See details for more information.
#' @param res \code{character}. Resolution of the data. One element from: '1KM','5KM','10KM','50KM','100KM'. \emph{Default: '1KM'}.
#' @param aggr \code{character}. Vector containing the aggregating factor(s). Can be either 'all' to download all the aggregating factors or a subset 
#' of c('md','mn','mi','ma','sd'). \emph{Default: 'md'}. See details for more information.
#' @param path \code{character}.The path to store the climatic grids.
#' @param extent The extent to which you want to crop your climatic datasets. extent must be either a 'SpatRaster', 'SpatVector', or'SpatExtent'.
#' @param mask \code{logical}. Do you want to mask your maps? \emph{Note that when TRUE, extent must be a 'SpatRaster' or a 'SpatVector'. }
#' @param compress \code{logical}. When extent is provided, do you want to compress the climatic grids. if TRUE, the compression
#' will be "COMPRESS=DEFLATE", "PREDICTOR=2" and "ZLEVEL=6". \emph{Default: TRUE}.
#' @param timeout \code{integer}. The maximum downloading time allowed for each map in seconds. \emph{Default: 300L}.
#' #' @details  
#' \describe{}
#' @return 
#' "done" when finished
#' @references
#' Amatulli, G., Domisch, S., Tuanmu, M.-N., Parmentier, B., Ranipeta, A., Malczyk, J., and Jetz, W. (2018) A suite of global, cross-scale topographic variables for environmental and biodiversity modeling. \emph{Scientific Data}. \bold{5}, 180040. \doi{10.1038/sdata.2018.40}. 
#' 
#' @examples  \donttest{
#' # get_topography("elevation")
#' }
#' @export 
#### Get Topography

get_topography <- function(var.names,
                           res = "1KM",
                           aggr = "md",
                           path = getwd(),
                           extent = NULL,
                           mask = FALSE, 
                           compress = TRUE,
                           timeout = 300L
                           ){
  aggr.full <- c("md","mn","mi","ma","sd")
  var.full <- c("elevation","slope","aspectcosine","aspectsine","eastness","northness",
                "roughness","tpi","tri","vrm","dx","dxx","dy","dyy","pcurv","tcurv")
  res.full <- c('1KM','5KM','10KM','50KM','100KM')
  ## Check variable names ####
  
  if(length(var.names)==1){
    if(var.names=="all"){
      var.down = var.full
    }else{
      if(any(var.names %in% var.full)){
        stop("var.names must be 'all' or a subset of : c('elevation','slope','aspectcosine','aspectsine','eastness','northness',
                'roughness','tpi','tri','vrm','dx','dxx','dy','dyy','pcurv','tcurv').")
      }
      var.down = var.names
    }
  }else{
    if(any(var.names %in% var.full)){
      stop("var.names must be 'all' or a subset of : c('elevation','slope','aspectcosine','aspectsine','eastness','northness',
                'roughness','tpi','tri','vrm','dx','dxx','dy','dyy','pcurv','tcurv').")
    }
    var.down = var.names 
  }
  
  ## Check aggregation names ####
  
  if(length(aggr)==1){
    if(aggr=="all"){
      aggr.down = aggr.full
    }else{
      if(any(aggr %in% aggr.full)){
        stop("aggr must be 'all' or a subset of : c('md','mn','mi','ma','sd').")
      }
      aggr.down = aggr
    }
  }else{
    if(any(aggr %in% aggr.full)){
      stop("aggr must be 'all' or a subset of : c('md','mn','mi','ma','sd').")
    }
    aggr.down = aggr 
  }
  
  ## Check resolution ####
  
  if(length(res)>1){
    stop("res must be composed of only one element")
  }else{
    if(any(aggr %in% aggr.full))
  }
  
  ## Check map transformation objects ####
  
  if(!is.null(extent)){
    if(!(inherits(extent,c("SpatRaster", "SpatVector","SpatExtent") ))){
      stop("extent must be either a 'SpatRaster', 'SpatVector', or'SpatExtent'.")
      if(mask & inherits(extent,c("SpatExtent") )){
        stop("when mask=TRUE, extent must be a 'SpatRaster' or a 'SpatVector'.")
      }
    }
    if(!is.logical(compress)){
      stop("compress must be a logical")
    }
  }
  
  ## Perform downloading and map cropping ####
  if(!is.integer(timeout)){
    stop("timeout must be an integer.")
  }
  options(timeout = max(timeout, getOption("timeout")))
  
  for(j in 1:length(aggr)){
    for(i in 1:length(var.names)){
      download.file(paste0("https://data.earthenv.org/topography/",
                           var.names[i],"_",toupper(res),aggr[j],"_GMTED",aggr[j],".tif"),
                    destfile=paste0(path,"/",var.names[i],"_",toupper(res),aggr[j],"_GMTED",aggr[j],".tif"),
                    mode="wb")
      if(!is.null(extent)){
        map <- terra::rast(paste0(path,"/",var.names[i],"_",toupper(res),aggr[j],"_GMTED",aggr[j],".tif"))
        map <- terra::crop(map,extent,mask=mask)
        if(compress){
          terra::writeRaster(map, 
                             paste0(path,"/",var.names[i],"_",toupper(res),aggr[j],"_GMTED",aggr[j],".tif"), 
                             overwrite = TRUE,
                             gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6"))
        }else{
          terra::writeRaster(map, 
                             paste0(path,"/",var.names[i],"_",toupper(res),aggr[j],"_GMTED",aggr[j],".tif"), 
                             overwrite = TRUE)
        }
      }
    }
  }
  return("done")
}


