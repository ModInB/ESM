### Dataset documentation

#' Species and environmental data to perform ESM with the ESM package
#' 
#' Contains sampling of 3 bryophytes in 413 plots in the Western Swiss Alps with
#' environmental values from Collart et al (2024)
#'
#' @format A data.frame with 413 plots and 9 variables:
#' \describe{
#' \item{x}{x coordinate in the CRS EPSG:2056}
#' \item{y}{y coordinate in the CRS EPSG:2056}
#' \item{Campylophyllum_halleri}{Presence/Absence data of \emph{Campylophyllum halleri}}
#' \item{Dicranella_subulata}{Presence/Absence data of \emph{Dicranella subulata}}
#' \item{Solenostoma_gracillimum}{Presence/Absence data of \emph{Solenostoma gracillimum}}
#' \item{ch_edaphic_eivdescombes_pixel_r}{Ecological indicator values reflecting 
#' Soil pH.}
#' \item{sradY}{Annual solar radiation. Sum of monthly solar radiations}
#' \item{bio3_tiso_8110_LV95}{Isothermalithy. It is equal to: Mean Diurnal 
#' Range/Temperature Annual Range * 100. The values have been averaged between 
#' 1981 and 2010}
#' \item{bio15_ps_8110_LV95}{Precipitation Seasonality. Computed using the 
#' coefficient of variation among monthly precipitation and averaged between
#' 1981-2010.}
#' }
#' 
#' @author Flavien Collart
#' @references
#' Collart, F., Kiebacher, T., Quetsch, M., Broennimann, O., Guisan, A. & Vanderpoorten, A. 2024. To what extent can we predict 
#' variation of bryophyte and tracheophyte community composition at fine spatial scale along an elevation gradient?
#' \emph{STOTEN}, 926:171741. \doi{10.1016/j.scitotenv.2024.171741}.
#' 
#' Külling, N., Adde, A., Fopp, F. et al. 2024. SWECO25: a cross-thematic raster 
#' database for ecological research in Switzerland. \emph{Sci Data}, 11:21. 
#' \doi{10.1038/s41597-023-02899-1}.
#' 
#' @examples 
#' data(ESM_species.env)
"ESM_Species.Env"


#' Environmental SpatRaster to perform ESM with the ESM package
#' 
#' Contains stack of 4 predictors at 25m resolution in the Western Swiss Alps.
#'
#' @format A SpatRaster with 4 layers:
#' \describe{
#' \item{ch_edaphic_eivdescombes_pixel_r}{Ecological indicator values reflecting 
#' Soil pH.}
#' \item{sradY}{Annual solar radiation. Sum of monthly solar radiations}
#' \item{bio3_tiso_8110_LV95}{Isothermalithy. It is equal to: Mean Diurnal 
#' Range/Temperature Annual Range * 100. The values have been averaged between 
#' 1981 and 2010}
#' \item{bio15_ps_8110_LV95}{Precipitation Seasonality. Computed using the 
#' coefficient of variation among monthly precipitation and averaged between
#' 1981-2010.}
#' }
#' @author Flavien Collart
#' @references
#' Collart, F., Kiebacher, T., Quetsch, M., Broennimann, O., Guisan, A. & Vanderpoorten, A. 2024. To what extent can we predict 
#' variation of bryophyte and tracheophyte community composition at fine spatial scale along an elevation gradient?
#' \emph{STOTEN}, 926:171741. \doi{10.1016/j.scitotenv.2024.171741}.
#' 
#' Külling, N., Adde, A., Fopp, F. et al. 2024. SWECO25: a cross-thematic raster 
#' database for ecological research in Switzerland. \emph{Sci Data}, 11:21. 
#' \doi{10.1038/s41597-023-02899-1}.
#' 
#' @examples 
#' library(terra)
#' data(ESM_Env)
#' ESM_Env <- terra::unwrap(ESM_Env)
#' 
"ESM_Env"


#' Species occurrence data to perform ESM
#' 
#' Contains a filtered dataset obtained from GBIF.org for \emph{Splachnum melanocaulon},
#' a rare endangered species in Scandinavia.
#'
#' @format A data.frame with 36 lines and 50 columns in the GBIF format:
#' \describe{
#' \item{decimalLongitude}{Longitude in the CRS EPSG:4326}
#' \item{decimalLatitude}{Latitude in the CRS EPSG:4326}
#' }
#' 
#' @author Flavien Collart
#' @details
#' The data has been downloaded from GBIF.org and afterwards filtered using the default
#' parameter of the function "clean_coordinates" from the "CoordinateCleaner" R package.
#' The temporal data is between 1981 and 2024 and coordinate uncertainty is less than 500 m.
#' 
#' @references
#' GBIF.org (16 December 2024) GBIF Occurrence Download. \doi{10.15468/dl.a5fhfv}.
#' 
#' @examples 
#' data(ESM_Splachnum.Data)
"ESM_Splachnum.Data"

#' Environmental SpatRaster to perform ESM in Scandinavia
#' 
#' Contains stack of 4 predictors at ~20km resolution in Scandinavia at present time and 3 in 2071-2100 .
#'
#' @format A SpatRaster with 5 layers:
#' \describe{
#' \item{bio1}{Annual Mean Temperature averaged between 1981-2010}
#' \item{bio3}{Isothermalithy. It is equal to: Mean Diurnal 
#' Range/Temperature Annual Range * 100. The values have been averaged between 
#' 1981 and 2010}
#' \item{bio5}{Max Temperature of Warmest Month averaged between 1981-2010}
#' \item{northness}{Northness.close to 1 corresponds to a northern exposition on a 
#' vertical slope, while a value close to -1 corresponds to a very steep southern slope}
#' \item{bio1_UKESM585}{Annual Mean Temperature averaged between 2071-2100 under the GCM 
#' UKESM 1-0-ll combined with the ssp 5-8.5}
#' \item{bio3_UKESM585}{Isothermalithy. It is equal to: Mean Diurnal 
#' Range/Temperature Annual Range * 100. The values have been averaged between 
#' 2071 and 2100 under the GCM UKESM 1-0-ll combined with the ssp 5-8.5}
#' \item{bio5_UKESM585}{Max Temperature of Warmest Month averaged between 2071-2100
#' under the GCM UKESM 1-0-ll combined with the ssp 5-8.5}
#' }
#' @details
#' Data were downloaded using get_Chelsa.Clim and get_Topography and aggregated to ~20km.
#' 
#' @author Flavien Collart
#' 
#' @references
#' Amatulli, G., Domisch, S., Tuanmu, M.-N., Parmentier, B., Ranipeta, A., Malczyk, J., and Jetz, W. (2018) A suite of global, cross-scale topographic variables for environmental and biodiversity modeling. \emph{Scientific Data}. \bold{5}, 180040. \doi{10.1038/sdata.2018.40}. 
#' 
#' Karger, D.N., Conrad, O., Böhner, J., Kawohl, T., Kreft, H., Soria-Auza, R.W., Zimmermann, N.E., Linder, P., Kessler, M. (2017). Climatologies at high resolution for the Earth land surface areas. \emph{Scientific Data}. \bold{4}, 170122. \doi{10.1038/sdata.2017.122}.
#' 
#' @examples 
#' library(terra)
#' data(ESM_Splachnum.Env)
#' ESM_Splachnum.Env <- terra::unwrap(ESM_Splachnum.Env)
#' 
"ESM_Splachnum.Env"
