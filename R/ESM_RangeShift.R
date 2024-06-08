#' @name ESM_Range.Shift
#' @author Flavien Collart \email{flaviencollart@hotmail.com}
#' @title Ensemble of Small Models: Range Shift between projections
#' @description This function compare ESM or SDM binary projections between two time periods by comparing the change 
#' in suitable/Unsuitable pixels. It also returns percentage of loss and gained
#' 
#' @param proj.curr \code{SpatRaster}. One or more ESM/SDM binary projections at current time. 
#' \emph{Note that if the number of layer is greater than one for proj.curr, the number of layer 
#' in proj.fut must be the same as comparisons will be between between map i of proj.curr and map i of proj.fut.}
#' @param proj.fut \code{SpatRaster}. One or more ESM/SDM binary projections at future time.
#' 
#' @details
#' The comparisons betwween SDM/ESM binary predictions are made following this formula: 
#' \bold{RangeShift = proj.curr + 2*proj.fut}. 
#' Therefore, when the value is equal to 0. It correspond to a pixel considered as "absent"/"unsuitable" in both 
#' time period. 
#' When it is 1, pixels considered as "presence" in proj.curr  but not in proj.fut \bold{(pixel lost)}.
#' When it is 2, pixels considered as "presence" in proj.fut  but not in proj.curr \bold{(pixel gained)}.
#' When it is 3, pixels are considered as "present" in both time periods \bold{(pixel Stable)}.
#'   
#' 
#' 
#' @return \code{list} containing:
#' \itemize{
#' \item{RangeShift}: \code{SpatRaster}. The range comparisons between current and future (\emph{see Details for the values}).
#' \item{RangeShift.table}: \code{data.frame}. Summary of the comparisons between maps. containing 9 columns:
#' \itemize{
#' \item{layer}. The layer names of proj.fut used to make the comparisons.
#' \item{Pixel.Absence.Stable}. The number of pixels that are considered as absent for both time period.
#' \item{Pixel.Presence.Current.Only}. The number of pixels considered as "presence" in proj.curr but not in proj.fut.
#' \item{Pixel.Presence.Future.Only}. The number of pixels considered as "presence" in proj.fut but not in proj.curr.
#' \item{Pixel.Presence.Stable}. The number of pixels considered as "presence" both in proj.curr and proj.fut.
#' \item{Current.Range.size}. The number of pixels considered as "presence" in proj.curr.
#' \item{Future.Range.size}. The number of pixels considered as "presence" in proj.fut assuming that species is
#' not limited by their dispersal capacities.
#' \item{Perc.Loss}. The percentage of pixel considered as "presence" in proj.curr and becoming "absent"
#' in proj.fut compared to the current range size. 
#' \item{Perc.Gain}. The percentage of pixel considered as "absent" in proj.curr and becoming 
#' "present" in proj.fut compared to the current range size. 
#' }
#' }
#' 
#' @seealso [ESM_Modeling], [ESM_Ensemble.Modeling],  [ESM_Projection], [ESM_Ensemble.Projection], [ESM_Threshold]
#' @export

ESM_Range.Shift <- function(proj.curr,
                            proj.fut){
  
  if(!(inherits(proj.curr,"SpatRaster") & inherits(proj.fut,"SpatRaster"))){
    stop("proj.curr and proj.fut must be SpatRaster a accepted")
  }

  if(terra::nlyr(proj.curr) > 1 & terra::nlyr(proj.curr) != terra::nlyr(proj.fut)){
    stop("When proj.curr has more than one layer the number of layer in proj.fut must match terra::nlyr(proj.curr)")
  }
  shift <- proj.curr + 2 * proj.fut
  results <- terra::freq(shift, wide=T)
  n.column <- ncol(results)
  while(ncol(results) < 5){
    
    results <- cbind(results,0)
  } 
  colnames(results)[(n.column+1):5] = setdiff(c("layer",0:3), colnames(results))
  results <- results[,c(1,order(colnames(results)[2:5])+1)]
  colnames(results) <-c("layer","Pixel.Absence.Stable","Pixel.Presence.Current.only",
                        "Pixel.Presence.Future.Only","Pixel.Presence.Stable")
  results$layer <- names(proj.fut)
  results$Current.Range.size <- results$Pixel.Presence.Current.only+ results$Pixel.Presence.Stable
  results$Future.Range.size <- results$Pixel.Presence.Future.Only+ results$Pixel.Presence.Stable
  results$Perc.Loss <- round(100*results$Pixel.Presence.Current.only/
                               (results$Pixel.Presence.Current.only+results$Pixel.Presence.Stable),
                             2)
  
  results$Perc.Gain <- round(100*results$Pixel.Presence.Future.Only/
                               (results$Pixel.Presence.Current.only+results$Pixel.Presence.Stable),
                             2)  
  
  return(list(RangeShift = shift,
              RangeShift.table = results))
}