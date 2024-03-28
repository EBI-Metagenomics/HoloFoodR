#' Get metabolomic data from MetaboLights database
#'
#' @details
#' The HoloFood database primarily comprises targeted metabolomic data,
#' omitting non-targeted metabolomic information. Nonetheless, it features URLs
#' linking to studies within the MetaboLights database. This functionality
#' enables users to access non-targeted metabolomic data. The function returns
#' a structured list encompassing data frames for study metadata,
#' assay metadata, and assay.
#'
#' @param url \code{character vector} specifying the URL address of study in
#' MetaboLights database.
#'
#' @param ... optional arguments:
#' \itemize{
#'   
#'   \item \strong{cache.dir} \code{Character scalar} specifying directory
#'   where downloaded file is stored. (Default: \code{tempdir()})
#'   
#' }
#'
#' @return \code{list}
#'
#' @examples
#' 
#' # This example is not run, because the server fails to respond sometimes.
#' 
#' url <- "https://www.ebi.ac.uk/metabolights/ws/studies/MTBLS4381"
#' #res <- getMetaboLights(url)
#' #names(res)
#' #head(res[["feat_meta"]])
#' 
#' @seealso
#' \code{\link[HoloFoodR:getResult]{getResult}}
#' \code{\link[HoloFoodR:getData]{getData}}
#'
#' @name getMetaboLights
NULL

#' 
#' @rdname getMetaboLights
#' @export
#' @importFrom dplyr full_join
getMetaboLights <- function(url, ...){
    # Check urls
    temp <- .check_input(url, list("character vector"))
    #
    # Get unique urls
    url <- unique(url)
    # Loop through those unique url addresses
    res <- lapply(url, function(x) .get_metabolomic_data(x, ...))
    # Get assay, assay metadata and study metadata separately
    assay <- lapply(res, function(x) x[["assay"]])
    assay_meta <- lapply(res, function(x) x[["assay_meta"]])
    study_meta <- lapply(res, function(x) x[["study_meta"]])
    # ...and combine results from different urls
    assay <- .full_join_list(assay)
    assay_meta <- .full_join_list(assay_meta)
    study_meta <- .full_join_list(study_meta)
    # Drop duplicates
    assay <- unique(assay)
    assay_meta <- unique(assay_meta)
    study_meta <- unique(study_meta)
    # Return a list
    res <- list(assay = assay, assay_meta = assay_meta, study_meta = study_meta)
    return(res)
}

################################ HELP FUNCTIONS ################################

# This function retrieves metabolomic data from MetaboLights database for single
# URL address
#' @importFrom dplyr left_join
.get_metabolomic_data <- function(url, ...){
    # In MetaboLight, the study IDs are different than in HoloFood. Get
    # Info about the study that corresponds to this particular HoloFood study.
    study_info <- .get_study_info(url, ...)
    # Get study metadata
    study_id <- study_info[["identifier"]]
    file_name <- study_info[["filename"]]
    study_metadata <- .get_metabolights_file(study_id, file_name, ...)
    # Get assays. A Study might have multiple assays
    assays_info <- study_info[["assays"]]
    assays <- lapply(assays_info, function(assay_info){
        # Get file name associated to assay
        file_name <- assay_info[["filename"]]
        # Get metadata on assay
        assay_metadata <- .get_metabolights_file(study_id, file_name, ...)
        # Get metabolomics data, the abundance table
        file_names <- unique(assay_metadata[["Metabolite Assignment File"]])
        assay <- lapply(file_names, function(file_name){
            temp <- .get_metabolights_file(study_id, file_name, ...)
            return(temp)
        })
        # Bind tables together
        assay <- .full_join_list(assay)
        # Return a list that have metadata and abundance table
        temp <- list(assay = assay, metadata = assay_metadata)
        return(temp)
    })
    # Combine assay metadata and abundance tables
    assay <- lapply(assays, function(x) x[["assay"]])
    assay_metadata <- lapply(assays, function(x) x[["metadata"]])
    # There are columns named ununique. Use R base rbind, because it does not
    # check the names. This might file if number of columns do not match...
    assay <- .full_join_list(assay)
    assay_metadata <- .full_join_list(assay_metadata)
    # Ensure that ID columns is character
    assay[["feat_ID"]] <- as.character( assay[["feat_ID"]] )
    # Make column names unique. For some reason files include
    # non-unique column names that have unique information.
    colnames(assay) <- make.unique( colnames(assay) )
    colnames(assay_metadata) <- make.unique( colnames(assay_metadata) )
    colnames(study_metadata) <- make.unique(colnames(study_metadata))
    # Create a list of data
    res <- list(
        assay = assay, assay_meta = assay_metadata, study_meta = study_metadata)
    return(res)
}

# This function fetches info about a study
.get_study_info <- function(url, ...){
    # From the metabolights database, find associated study. Which study
    # represents this HoloFood study?
    res <- .perform_single_query(path = "metabolight", full.url = url, ...)
    # Get only relevant info
    study_info <- res[["isaInvestigation"]][["studies"]]
    return(study_info)
}

# This is a common function for downloading a file from MetaboLights database
#' @importFrom utils download.file read.delim
.get_metabolights_file <- function(
        study.id, file.name, cache.dir = tempdir(),
        metabolights.base.url = "http://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/", ...){
    # Check metabolights.base.url
    temp <- .check_input(metabolights.base.url, list("character scalar"))
    # Check study.id
    temp <- .check_input(study.id, list("character scalar"))
    # Check file.name
    temp <- .check_input(file.name, list("character scalar"))
    # Check cache.dir
    temp <- .check_input(cache.dir, list("character scalar"))
    #
    # Create url
    url <- paste0( metabolights.base.url, "/", study.id, "/", file.name)
    # Create a directory path and create the dir if it is not existing
    cache_dir <- file.path(cache.dir, "HoloFoodR_cache")
    if( !dir.exists(cache_dir) ){
        dir.create(cache_dir)
    }
    # Create a file path
    file_path <- file.path(cache_dir, file.name)
    # Check if file is already loaded. If not, download from internet.
    if( !file.exists(file_path) ){
        download.file(url, file_path, quiet = TRUE)
    }
    # Read the local file
    df <- read.delim(file_path, check.name = FALSE)
    return(df)
}
