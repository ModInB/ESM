#' @name get_Chelsa.Clim
#' @author Flavien Collart \email{flaviencollart@hotmail.com} with contributions of Adèle Hotermans
#' @title Download CHELSA data
#' @description
#' Download and crop climatic variables at different time periods from CHELSA v.2.1.
#' 
#' @param var.names \code{character}. Vector containing either 'all' to download all the variables (works only at present time), 'all_fut'
#' for all variables that are also available for future time periods, 'bioclim' to download the 19 bioclimatic variables,
#' or the variable names to download. See details for more information.
#' @param time \code{character}. The time period(s) wanted. Must be: '1981-2010','2011-2040','2041-2070',and/or '2071-2100'.
#' \emph{Default: '1981-2010'}
#' @param gcm \code{character}. A vector containing the Global Circulation Models for which you want to download the datasets. 
#' Only needed when future time periods are provided. Must be: 'gfdl-esm4','ipsl-cm6a-lr','mpi-esm1-2-hr','mri-esm2-0' and/or,'ukesm1-0-ll'.
#' \emph{Default: NULL}.
#' @param ssp \code{character}. A vector containing the Shared Socio-economic Pathways for which you want to download the datasets. 
#' Only needed when future time periods are provided. Must be: '126', '370' and/or '585'.\emph{Default: NULL}.
#' @param path \code{character}.The path to store the climatic grids.
#' @param extent The extent to which you want to crop your climatic datasets. Extent must be either a 'SpatRaster', 'SpatVector', or'SpatExtent'.
#' @param mask \code{logical}. Do you want to mask your maps? \emph{Note that when TRUE, extent must be a 'SpatRaster' or a 'SpatVector'.}
#' @param compress \code{logical}. When extent is provided, do you want to compress the climatic grids? If TRUE, the compression
#' will be: "COMPRESS=DEFLATE", "PREDICTOR=2" and "ZLEVEL=6". \emph{Default: TRUE}.
#' @param timeout \code{integer}. The maximum downloading time allowed for each map in seconds. \emph{Default: 300L}.
#' @details  
#' \describe{
#' get_chesla downloads and crops climatic variables at different time periods from CHELSA v.2.1.
#' 
#' The following list contains the different variables associated with their time period : 
#' \itemize{
#' 
#' \item{Present}: pet_perman_max, pet_perman_mean, pet_perman_min, pet_perman_range,
#' sfcWind_max, sfcWind_mean, sfcWind_min","sfcWind_range, vpd_max, vpd_mean, vpd_min, vpd_range,
#' hurs_max, hurs_mean, hurs_min, hurs_range, clt_max, clt_mean, clt_min, clt_range, cmi_max, cmi_mean,
#' cmi_min, cmi_range, rsds_1981-2010_mean, rsds_1981-2010_min, rsds_1981-2010_max, rsds_1981-2010_range.
#' 
#' \item{Present and Future}: the 19 bioclim variables (e.g. "bio1"), ai, fcf, fgd, gdd0, gdd10,
#' gddd5, gddlgd0, gddlgd5, gddlgd10, gdgfgd0, gdgfgd10, gdgfgd5, gsl, gsp, gst,
#' kg from 0 to 5 (e.g. "kg0"), lgd, ngd0, ngd10, ngd5, npp, scd, swd, swe.
#' }
#' The description of the different variables is available at the point 7 (p.11) of the following link :
#' \href{https://chelsa-climate.org/wp-admin/download-page/CHELSA_tech_specification_V2.pdf}{CHELSA_tech_specification_V2.pdf}
#' 
#' The parameter "var.names" allows you to select the variables taken into account :
#' \itemize{
#' 
#' \item{all}: All available variables are selected. \emph{Note: Only works for present time}
#' \item{all_fut}: All the variables available in present AND future time are selected.
#' \item{bioclim}: the 19 bioclimatic variables are selected.
#' }
#' The variables available in the future are to be combined with a future scenario by choosing a global circular model (gcm) and a shared-socioeconomic pathway (ssp) :
#' \itemize{
#' 
#' \item{ssp126}: ssp1 (Sustainability - Taking the green road) combined with a radiative forcing of 2.6 W/m².
#' \item{ssp370}: ssp3 (Regional rivalry - A rocky road) combined with a radiative forcing of 7 W/m².
#' \item{ssp585}: ssp5 (Fossil-fueled development - Taking the highway) with a radiative forcing of 8.5 W/m². 
#' }
#' The complete description of the ssp is available in O'Neill et al. 2017. These ssp are combined with
#' Global circulation Models (GCM) to assess future climate.
#' 
#' If the files are to heavy and/or your connection is to slow, you may want to modified the timeout parameter to a greater period of time.
#' Note that you always have to end your number with "L" to make sure your time is an integer (e.g. 300L, 500L).
#' 
#' To crop a map to a different size than the original one, it is possible to add an extent and a mask parameter. If the extent of 
#' the map is large, you may want to compress it, reducing storage space without significantly increase the reading time in R.
#' }
#' 
#' @return 
#' \code{character}. A vector containing the path to the downloaded files.
#' @references
#' Karger, D.N., Conrad, O., Böhner, J., Kawohl, T., Kreft, H., Soria-Auza, R.W., Zimmermann, N.E., Linder, P., Kessler, M. (2017). Climatologies at high resolution for the Earth land surface areas. \emph{Scientific Data}. \bold{4}, 170122. \doi{10.1038/sdata.2017.122}.
#' 
#' Karger D.N., Conrad, O., Böhner, J., Kawohl, T., Kreft, H., Soria-Auza, R.W., Zimmermann, N.E, Linder, H.P., Kessler, M. (2018): Data from: Climatologies at high resolution for the earth’s land surface areas. \emph{EnviDat}. \doi{10.16904/envidat.228.v2.1}.
#' 
#' B. C. O’Neill, E. Kriegler, K. L. Ebi, E. Kemp-Benedict, K. Riahi, D. S. Rothman, B. J. van Ruijven, D. P. van Vuuren, J. Birkmann, K. Kok, M. Levy, W. Solecki. (2017). The roads ahead: Narratives for shared socioeconomic pathways describing world futures in the 21st century.
#' \emph{Global Environmental Change}. \bold{42} 169-180. \doi{10.1016/j.gloenvcha.2015.01.004}
#'
#' @examples  \donttest{
#' # get_Chelsa.Clim("bio1")
#' }
#' @importFrom utils download.file
#' @export 

get_Chelsa.Clim <- function(var.names,
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
        stop("When var.names = 'all', time must be only '1981-2010'. Please try var.names='all_fut' instead.")
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
            stop("All elements in var.names must be present in extended bioclimate or must be either all, all_fut ,or bioclim.")
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
      if(any(var.names %in% var.pres) | !all(var.names %in% var.full)){
        stop("When time contains a future period, var.names should be a group of variables available for future time period.")
      }
    }else{
      if(time == "1981-2010"){
        if(!all(var.names %in% c(var.full,var.pres))){
          stop("All elements in var.names must be present in extended bioclimate or must be either all, all_fut ,or bioclim.")
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
        stop("When mask=TRUE, extent must be a 'SpatRaster' or a 'SpatVector'.")
      }
    }
    if(!is.logical(compress)){
      stop("compress must be a logical.")
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
        if(file.exists(paste0(path,"/CHELSA_",var.down[i],"_1981-2010_V.2.1.tif"))){
          warning(paste("File 'CHELSA_",var.down[i],"_1981-2010_V.2.1.tif' is already present in the location."))
          next
        }
        download.file(url=paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_",
                                 var.down[i],"_1981-2010_V.2.1.tif"),
                      destfile=paste0(path,"/CHELSA_",var.down[i],"_1981-2010_V.2.1.tif"),
                      mode="wb")
        if(!is.null(extent)){
          map <- terra::rast(paste0(path,"/CHELSA_",var.down[i],"_1981-2010_V.2.1.tif"))
          map <- terra::crop(map,extent,mask=mask)
          if(compress){
            terra::writeRaster(map, 
                               paste0(path,"/CHELSA_",var.down[i],"_1981-2010_V.2.1.tif"), 
                               overwrite = TRUE,
                               gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6"))
          }else{
            terra::writeRaster(map, 
                               paste0(path,"/CHELSA_",var.down[i],"_1981-2010_V.2.1.tif"), 
                               overwrite = TRUE)
            }
        }
        FilePath <- c(FilePath, paste0(path,"/CHELSA_",var.down[i],"_1981-2010_V.2.1.tif")) 
        }else{
          for(h in 1:length(gcm)){
            for(g in 1:length(ssp)){
              if(file.exists(paste0(path,"/CHELSA_",var.down[i],"_",time[j],"_",gcm[h],"_ssp",ssp[g],"_V.2.1.tif"))){
                warning(paste("File 'CHELSA_",var.down[i],"_",time[j],"_",gcm[h],"_ssp",ssp[g],"_V.2.1.tif' is already present in the location."))
                next
              }
              download.file(url=paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/",
                                     time[j],"/",toupper(gcm[h]),"/ssp",ssp[g],
                                     "/bio/CHELSA_",var.down[i],"_",time[j],"_",gcm[h],"_ssp",ssp[g],"_V.2.1.tif"),
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

#' @name get_Topography
#' @author Flavien Collart \email{flaviencollart@hotmail.com} with contributions of Adèle Hotermans
#' @title Download topographic data
#' @description
#' Download and crop topographic variables from Amatulli et al (2018) using GMTED elevation data.
#' 
#' @param var.names \code{character}. Vector containing the variable name(s) to download. Can be either 'all' to download all variables or a subset of the variables  
#' available in Amatulli et al (2018). See details for more information.
#' @param res \code{character}. Resolution of the data. One element from: '1KM','5KM','10KM','50KM','100KM'. \emph{Default: '1KM'}.
#' @param aggr \code{character}. Vector containing the aggregating factor(s). Can be either 'all' to download all the aggregating factors or a subset 
#' of c('md','mn','mi','ma','sd'). \emph{Default: 'md'}. See details for more information.
#' @param path \code{character}.The path to store the topographic grids.
#' @param extent The extent to which you want to crop your topographic datasets. Extent must be either a 'SpatRaster', 'SpatVector', or'SpatExtent'.
#' @param mask \code{logical}. Do you want to mask your maps? \emph{Note that when TRUE, extent must be a 'SpatRaster' or a 'SpatVector'. }
#' @param compress \code{logical}. When extent is provided, do you want to compress the topographic grids? If TRUE, the compression
#' will be "COMPRESS=DEFLATE", "PREDICTOR=2" and "ZLEVEL=6". \emph{Default: TRUE}.
#' @param timeout \code{integer}. The maximum downloading time allowed for each map in seconds. \emph{Default: 300L}.
#' @details  
#' \describe{
#' get_topography downloads and crops topographic variables from Amatulli et al (2018) using GMTED elevation data.
#' 
#' The description of each following variables comes from Amatulli et al (2018) : 
#' \itemize{ 
#' \item{elevation}: the elevation across the terrain (expressed in meters).
#' \item{slope}: the rate of change of elevation in the direction of the water flow line (expressed in degrees).
#' \item{aspectcosine}: cosine of the aspect (angular direction that a slope faces).
#' \item{aspectsine}: sine of the aspect (angular direction that a slope faces).
#' \item{eastness}: sine of the slope multiplied by the sine of the aspect.
#' \item{northness}: sine of the slope multiplied by the cosine of the aspect.
#' \item{roughness}: the largest inter-cell absolute difference of a focal cell and its 8 surrounding cells.
#' \item{tpi}: the topographic position index is the difference between the elevation of a focal cell and the mean of its 8 surrounding cells.
#' \item{tri}: the terrain ruggedness index  is a mean of the absolute differences in elevation between a focal cell and its 8 surrounding cells.
#' \item{vrm}: the vector ruggedness measure quantifies terrain ruggedness by measuring the variation by means of sine and cosine of the slope in the three-dimensional orientation of grid cells, within a moving window.
#' \item{dx}: the first order partial derivative (E-W slope)  the slope in an East-West direction.
#' \item{dxx}: the second order partial derivative (E-W slope) is the derivative of a slope in a East-West direction.
#' \item{dy}: the first order partial derivative (N-S slope) is the slope in a North-South direction.
#' \item{dyy}: the second order partial derivative (N-S slope) is the derivative of the slope in a North-South direction.
#' \item{pcurv}: the profile curvature measures the rate of change of a slope along a flow line, and affects the acceleration of water flow along a surface.
#' \item{tcurv}: the tangential curvature measures the rate of change perpendicular to the slope gradient and is related to the convergence and divergence of flow across a surface.
#' }
#' 
#' The aggregating factors correspond to the statistics computed from the variables.
#' Amatulli et al. used the 250m GMTED as a data source with a resolution of 250m.
#' GMTED is available in different forms derived from the different high resolution DEMs and 
#' corresponding to different statistics 
#' namely minimum (mi), maximum (ma), mean (mn), median (md) and standard deviation (sd).
#' For the elevation variable, the mi, ma, mn, md and sd elevations were obtained from the
#' source layers 250m GMTEDmi, 250 m GMTEDma, 250m GMTEDmn, 250m GMTEDmd and 250m GMTEDsd respectively
#' (note that the downloaded files finished by "md" but still correspond to the different source layers).
#' The rest of the variables were obtained from the 250m GMTEDmd and the others aggregating factors layers were calculated after.
#'
#' If the files are to heavy and/or your connection is to slow, you may want to modified the timeout parameter to a greater period of time.
#' Note that you always have to end your number with "L" to make sure your time is an integer (e.g. 300L, 500L).
#' 
#' To crop a map to a different size than the original one, it is possible to add an extent and a mask parameter. If the extent of 
#' the map is large, you may want to compress it, reducing storage space without significantly increase the reading time in R.
#' 
#' }
#' 
#' @return 
#' \code{character}. A vector containing the path to the downloaded files.
#' 
#' @references
#' Amatulli, G., Domisch, S., Tuanmu, M.-N., Parmentier, B., Ranipeta, A., Malczyk, J., and Jetz, W. (2018) A suite of global, cross-scale topographic variables for environmental and biodiversity modeling. \emph{Scientific Data}. \bold{5}, 180040. \doi{10.1038/sdata.2018.40}. 
#' 
#' @examples  \donttest{
#' # get_Topography("elevation")
#' }
#' @export 

get_Topography <- function(var.names,
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
  FilePath <- c() #File path to return at the end
  
  ## Check variable names ####
  
  if(length(var.names)==1){
    if(var.names=="all"){
      var.down = var.full
    }else{
      if(any(!(var.names %in% var.full))){
        stop("var.names must be 'all' or a subset of : c('elevation','slope','aspectcosine','aspectsine','eastness','northness',
                'roughness','tpi','tri','vrm','dx','dxx','dy','dyy','pcurv','tcurv').")
      }
      var.down = var.names
    }
  }else{
    if(any(!(var.names %in% var.full))){
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
      if(any(!(aggr %in% aggr.full))){
        stop("aggr must be 'all' or a subset of : c('md','mn','mi','ma','sd').")
      }
      aggr.down = aggr
    }
  }else{
    if(any(!(aggr %in% aggr.full))){
      stop("aggr must be 'all' or a subset of : c('md','mn','mi','ma','sd').")
    }
    aggr.down = aggr 
  }
  
  ## Check resolution ####
  
  if(length(res)>1){
    stop("res must be composed of only one element.")
  }else{
    if(!all(toupper(res) %in% res.full)){
      stop("res must be one element of c('1KM','5KM','10KM','50KM','100KM').")
    }
  }
  
  ## Check map transformation objects ####
  
  if(!is.null(extent)){
    if(!(inherits(extent,c("SpatRaster", "SpatVector","SpatExtent") ))){
      stop("extent must be either a 'SpatRaster', 'SpatVector', or'SpatExtent'.")
      if(mask & inherits(extent,c("SpatExtent") )){
        stop("When mask=TRUE, extent must be a 'SpatRaster' or a 'SpatVector'.")
      }
    }
    if(!is.logical(compress)){
      stop("compress must be a logical.")
    }
  }
  
  ## Perform downloading and map cropping ####
  if(!is.integer(timeout)){
    stop("timeout must be an integer.")
  }
  options(timeout = max(timeout, getOption("timeout")))
  
  for(j in 1:length(aggr.down)){
    for(i in 1:length(var.down)){
      if(file.exists(paste0(path,"/",var.down[i],"_",toupper(res),aggr[j],"_GMTEDmd.tif"))){
        warning(paste("File '",var.down[i],"_",toupper(res),aggr[j],"_GMTEDmd.tif' is already present in the location."))
        next
      }
      ## elevation is the only variable that is not written in the same way as the others
      if(var.down[i] == "elevation"){
        download.file(paste0("https://data.earthenv.org/topography/",
                             var.down[i],"_",toupper(res),aggr[j],"_GMTED",aggr[j],".tif"),
                      destfile=paste0(path,"/",var.down[i],"_",toupper(res),aggr[j],"_GMTEDmd.tif"),
                      mode="wb")
      }else{
        download.file(paste0("https://data.earthenv.org/topography/",
                             var.down[i],"_",toupper(res),aggr[j],"_GMTEDmd.tif"),
                      destfile=paste0(path,"/",var.down[i],"_",toupper(res),aggr[j],"_GMTEDmd.tif"),
                      mode="wb")
      }
      
      if(!is.null(extent)){
        map <- terra::rast(paste0(path,"/",var.down[i],"_",toupper(res),aggr[j],"_GMTEDmd.tif"))
        map <- terra::crop(map,extent,mask=mask)
        if(compress){
          terra::writeRaster(map, 
                             paste0(path,"/",var.down[i],"_",toupper(res),aggr[j],"_GMTEDmd.tif"), 
                             overwrite = TRUE,
                             gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6"))
        }else{
          terra::writeRaster(map, 
                             paste0(path,"/",var.down[i],"_",toupper(res),aggr[j],"_GMTEDmd.tif"), 
                             overwrite = TRUE)
        }
      }
      FilePath <- c(FilePath,paste0(path,"/",var.down[i],"_",toupper(res),aggr[j],"_GMTEDmd.tif")) 
    }
  }
  return(FilePath)
}


