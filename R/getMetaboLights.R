#' Get metabolomic data from MetaboLights database
#'
#' @details
#' The HoloFood database primarily comprises targeted metabolomic data,
#' omitting non-targeted metabolomic information. Nonetheless, it features URLs
#' linking to studies within the MetaboLights database. This functionality
#' enables users to access non-targeted metabolomic data. The
#' \code{getMetaboLights} function returns
#' a structured list encompassing processed data in \code{data.frame} format
#' for study metadata, assay metadata, and assay.
#' 
#' The metadata includes the file names of spectra data. Those files can be
#' loaded with \code{getMetaboLightsFile}. Alternatively, once you've identified
#' the study and files to fetch, you can refer to this
#' [vignette](https://rformassspectrometry.github.io/MsIO/articles/MsIO.html#loading-data-from-metabolights)
#' for instructions on loading the data directly into an \code{MsExperiment}
#' object, specifically designed for metabolomics spectra data.
#'
#' @param study.id \code{character vector} specifying the study identifier of
#' data that is going to be fetched from the MetaboLights database.
#' 
#' @param file \code{character vector} specifying the files that are being
#' fetched.
#'
#' @param ... optional arguments:
#' \itemize{
#'   
#'   \item \strong{cache.dir} \code{Character scalar} specifying directory
#'   where downloaded file is stored. (Default: \code{tempdir()})
#'   
#'   \item \strong{timeout} \code{Integer scalar} specifying timeout
#'   in seconds for loading a file. (Default: \code{5*60})
#'   
#' }
#'
#' @return \code{list}
#'
#' @examples
#' 
#' # This example is not run, because the server fails to respond sometimes.
#' if( FALSE ){
#'     res <- getMetaboLights("MTBLS4381")
#'     file_paths <- getMetaLightsFile(
#'         study.id = "MTBLS4381",
#'         file = res[["assay_meta"]][["Raw Spectral Data File"]]
#'         )
#' }
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
getMetaboLights <- function(study.id, ...){
    # Check study.id
    temp <- .check_input(study.id, list("character vector"))
    #
    # Get unique urls
    study.id <- unique(study.id)
    # Loop through those unique url addresses
    res <- lapply(study.id, function(x) .get_metabolomic_data(x, ...))
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

#' @rdname getMetaboLights
#' @export
getMetaboLightsFile <- function(study.id, file, ...){
    # Check study.id
    temp <- .check_input(study.id, list("character vector"))
    # Check files
    temp <- .check_input(file, list("character vector"))
    # Check that their dimensions are correct
    if( !(length(study.id) == 1 || length(study.id) == length(file)) ){
        stop("The length of 'study.id' must be 1 or equal to length of 'file'.",
            call. = FALSE)
    }
    #
    # Create a df that stores teh study.id and file
    fetch_df <- data.frame(study_id = study.id, file = file)
    # Get unique and put each instance to columns
    fetch_df <- unique(fetch_df) |> t() |> as.data.frame()
    # Loop through files and load them
    res <- lapply(fetch_df, function(col){
        .get_metabolights_file(col[[1]], col[[2]], return.table = FALSE, ...)
    })
    res <- res |> unlist() |> unname()
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
        # Get metadata on assays
        file_names <- unique(assay_info[["filename"]])
        assay_metadata <- lapply(file_names, function(file_name){
            .get_metabolights_file(study_id, file_name, ...)
        })
        # Bind tables together
        assay_metadata <- .full_join_list(assay_metadata)
        # Get metabolomics data, the abundance table
        file_names <- unique(assay_metadata[["Metabolite Assignment File"]])
        assay <- lapply(file_names, function(file_name){
            .get_metabolights_file(study_id, file_name, ...)
        })
        # Bind tables together
        assay <- .full_join_list(assay)
        # If there are feat_ID column, ensure that it is character as it holds
        # values of identifiers. It seems that some IDs include only numeric
        # values which is why they are uncorrectly interpreted as numeric
        # values.
        if( "feat_ID" %in% colnames(assay) ){
            assay[["feat_ID"]] <- as.character( assay[["feat_ID"]] )
        }
        # Return a list that have metadata and abundance table
        temp <- list(assay = assay, metadata = assay_metadata)
        return(temp)
    })
    # Combine assay metadata and abundance tables
    assay <- lapply(assays, function(x) x[["assay"]])
    assay_metadata <- lapply(assays, function(x) x[["metadata"]])
    # Merge all data from different assays
    assay <- .full_join_list(assay)
    assay_metadata <- .full_join_list(assay_metadata)
    # Create a list of data
    res <- list(
        assay = assay, assay_meta = assay_metadata, study_meta = study_metadata)
    return(res)
}

# This function fetches info about a study
#' @importFrom httr2 url_parse
.get_study_info <- function(
        url,
        study.search.url = "https://www.ebi.ac.uk/metabolights/ws/studies",
        ...){
    # Check if study.id is already a url address. If it is not, create url
    # from study.id and base.url
    parsed_url <- url_parse(url)
    if( !(!is.null(parsed_url$scheme) && !is.null(parsed_url$hostname)) ){
        url <- paste0(study.search.url, "/", url)
    }
    # From the metabolights database, find associated study. Which study
    # represents this HoloFood study?
    res <- .perform_single_query(path = "metabolight", full.url = url, ...)
    # Get only relevant info
    study_info <- res[["isaInvestigation"]][["studies"]]
    return(study_info)
}

# This is a common function for downloading a file from MetaboLights database
#' @importFrom utils download.file read.delim
#' @importFrom httr2 url_parse
.get_metabolights_file <- function(
        study.id, file.name, cache.dir = tempdir(), unique.cols = TRUE,
        timeout = 5*60, return.table = TRUE,
        metabolights.base.url = "http://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public",
        ...){
    # Check metabolights.base.url
    temp <- .check_input(metabolights.base.url, list("character scalar"))
    # Check study.id
    temp <- .check_input(study.id, list("character scalar"))
    # Check file.name
    temp <- .check_input(file.name, list("character scalar"))
    # Check cache.dir
    temp <- .check_input(cache.dir, list("character scalar"))
    # Check unique.cols
    temp <- .check_input(unique.cols, list("logical scalar"))
    # Check timeout
    temp <- .check_input(unique.cols, list("logical scalar"))
    # Check return.table
    temp <- .check_input(return.table, list("logical scalar"))
    #
    # If the study.id is url, get the study.id from the url
    parsed_url <- url_parse(study.id)
    if( !is.null(parsed_url$scheme) && !is.null(parsed_url$hostname) ){
        temp <- strsplit(study.id, "/")[[1]]
        study.id <- temp[length(temp)]
    }
    # Create url
    url <- paste0( metabolights.base.url, "/", study.id, "/", file.name)
    # Create a directory path
    cache_dir <- file.path(cache.dir, "HoloFoodR_cache")
    # Create a file path
    file_path <- file.path(cache_dir, file.name)
    # Create the dir if it is not existing
    cache_dir <- dirname(file_path)
    if( !dir.exists(cache_dir) ){
        dir.create(cache_dir, recursive = TRUE)
    }
    # Check if file is already loaded. If not, download from internet.
    if( !file.exists(file_path) ){
        # Set timeout as user-desired time
        def_opt <- getOption("timeout")
        options(timeout = timeout)
        # Load the data
        download.file(url, file_path, quiet = FALSE, timeout = timeout)
        # Set the timeout back to default
        options(timeout = def_opt)
    }
    # By default, the loaded table is returned. However, for spectra files, we
    # do not want to return them.
    if( return.table ){
        # Read the local file
        df <- read.delim(file_path, check.name = FALSE)
        # Make column names unique if specified
        if( anyDuplicated(colnames(df)) && unique.cols ){
            colnames(df) <- make.unique(colnames(df))
        }
    } else{
        df <- file_path
    }
    return(df)
}
