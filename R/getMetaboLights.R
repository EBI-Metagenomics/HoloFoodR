#' Get metabolomic data from MetaboLights database
#'
#' @details
#' Non-targeted metabolomics data is in MetaboLigths...
#'
#' @param url \code{character vector} specifying the URL address of study in
#' MetaboLights database.
#'
#' @param ... optional arguments:
#' \itemize{
#'   
#'   \item \strong{cache.dir} \code{Character scalar} specifying cache
#'   directory. (Default: \code{tempdir()})
#'   
#' }
#'
#' @return \code{list}
#'
#' @examples
#'
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
    # Get assay, feature metadata and sample metadata separately
    assay <- lapply(res, function(x) x[["assay"]])
    feat_meta <- lapply(res, function(x) x[["feat_meta"]])
    sample_meta <- lapply(res, function(x) x[["sample_meta"]])
    # ...and combine results from different urls
    assay <- .full_join_list(assay)
    feat_meta <- do.call(rbind, feat_meta)
    sample_meta <- do.call(rbind, sample_meta)
    # Drop duplicates
    assay <- unique(assay)
    feat_meta <- feat_meta[ !duplicated(feat_meta[["feat_ID"]]), ]
    sample_meta <- sample_meta[ !duplicated(sample_meta[["Sample Name"]]), ]
    # Return a list
    res <- list(assay = assay, feat_meta = feat_meta, sample_meta = sample_meta)
    return(res)
}

################################ HELP FUNCTIONS ################################

# This function retrieves metabolomic data from MetaboLights database for single
# URL address
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
        assay <- bind_rows(assay)
        # Return a list that have metadata and abundance table
        temp <- list(assay = assay, metadata = assay_metadata)
        return(temp)
    })
    # Combine assay metadata and abundance tables
    assay <- lapply(assays, function(x) x[["assay"]])
    assay_metadata <- lapply(assays, function(x) x[["metadata"]])
    # There are columns named ununique. Use R base rbind, because it does not
    # check the names. This might file if number of columns do not match...
    assay <- do.call(rbind, assay)
    assay_metadata <- do.call(rbind, assay_metadata)
    # Split assay to abundance table and feature metadata
    assay_cols <- colnames(assay) %in% assay_metadata[["Sample Name"]]
    feature_metadata <- assay[ , !assay_cols, drop = FALSE]
    assay <- assay[ , assay_cols, drop = FALSE]
    assay[["feat_ID"]] <- feature_metadata[["feat_ID"]]
    # Combine assay and study metadata
    assay_metadata <- assay_metadata[
        match(study_metadata[["Sample Name"]], assay_metadata[["Sample Name"]]),
    ]
    metadata <- cbind(study_metadata, assay_metadata)
    # Create a list of data
    res <- list(
        assay = assay, feat_meta = feature_metadata, sample_meta = metadata)
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
