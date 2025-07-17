#####################
# Load ESM objects

#' @name Load_ESM
#' @author Flavien Collart \email{flaviencollart@hotmail.com} 
#' @title Functions to load the ESM objects
#' @description These functions allow the user to easily load their ESM objects
#' @param sp.name \code{character}. The species name given in sp.name argument from ESM_Modelling()
#' @return a usable ESM object 
#' @details
#' The function allow to load the most recent ESM object. Note that the ESM folder 
#' should be in your working directory
#' @export

#' @rdname Load_ESM
#' @export
Load_ESM_Modelling <- function(sp.name){
  folder <- paste0("ESM.output_", sp.name)
  if(!dir.exists(folder)){
    stop("The ESM folder does not exist:", folder)
  }
  list.out.files <- list.files(path = folder,
                               pattern = "ESM.Modeling",
                               full.names = TRUE)
  
  if(length(list.out.files)==0){
    stop("No .out files are present in the folder. Please consider next time to turn save.obj to TRUE in order to save the object.")
  }else if(length(list.out.files)>1){
    MoreRecentFile <- which.max(file.mtime(list.out.files))
    list.out.files <- list.out.files[MoreRecentFile]
    cat(paste0("\nMore than 1 .out file is present... Taking the more recent one which is : ",list.out.files))
  }
  
  return(get(load(list.out.files)))
  
}

#' @rdname Load_ESM
#' @export
Load_ESM_Ensemble.Modelling <- function(sp.name){
  folder <- paste0("ESM.output_", sp.name)
  if(!dir.exists(folder)){
    stop("The ESM folder does not exist:", folder)
  }
  list.out.files <- list.files(path = folder,
                               pattern = "ESM.EsnsembleModeling",
                               full.names = TRUE)
  
  if(length(list.out.files)==0){
    stop("No .out files are present in the folder. Please consider next time to turn save.obj to TRUE in order to save the object.")
  }else if(length(list.out.files)>1){
    MoreRecentFile <- which.max(file.mtime(list.out.files))
    list.out.files <- list.out.files[MoreRecentFile]
    cat(paste0("\nMore than 1 .out file is present... Taking the more recent one which is : ",list.out.files))
  }
  
  return(get(load(list.out.files)))
  
}
#' @rdname Load_ESM
#' @export
Load_ESM_Projection <- function(sp.name){
  folder <- paste0("ESM.output_", sp.name)
  if(!dir.exists(folder)){
    stop("The ESM folder does not exist:", folder)
  }
  list.out.files <- list.files(path = folder,
                               pattern = "ESM.Projection",
                               full.names = TRUE)
  
  if(length(list.out.files)==0){
    stop("No .out files are present in the folder. Please consider next time to turn save.obj to TRUE in order to save the object.")
  }else if(length(list.out.files)>1){
    MoreRecentFile <- which.max(file.mtime(list.out.files))
    list.out.files <- list.out.files[MoreRecentFile]
    cat(paste0("\nMore than 1 .out file is present... Taking the more recent one which is : ",list.out.files))
  }
  
  return(get(load(list.out.files)))
  
}
