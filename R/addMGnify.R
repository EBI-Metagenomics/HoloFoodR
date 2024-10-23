#' Add results from MGnifyR to HoloFoodR results
#'
#' @details
#' Metagenomic data is found in MGnify rather than HoloFoodR, and the two
#' databases use different sample identifiers. However, MGnify's sample
#' metadata includes references to the identifiers used in the HoloFood
#' database, making it straightforward to convert sample IDs for alignment
#' with HoloFood data. Despite this, HoloFood contains additional metadata
#' not available in MGnify. Moreover, integrating data into a
#' \code{MultiAssayExperiment} while maintaining accurate sample and system
#' matches can be challenging.
#' 
#' This function is designed to simplify these
#' tasks, enabling seamless integration of MGnify data with HoloFood data after
#' retrieval from the database. You need only to input the returned data from
#' \code{MGnifyR::getResult()} and \code{HoloFoodR::getResult()} functions.
#' 
#' @param x \code{SummarizedExperiment}. Results from
#' \code{MGnifyR::getResult()}.
#' 
#' @param y \code{MultiAssayExperiment} or \code{SummarizedExperiment} or
#' \code{data.frame}-like table. Results from \code{HoloFoodR::getResult()} or
#' sample metadata from it.
#' 
#' @param exp.name1 \code{Character scalar}. Specifies the name of experiment
#' that will be added to \code{y}. (Default: \code{"metagenomic"})
#' 
#' @param exp.name2 \code{Character scalar}. Specifies the name of experiment
#' from HoloFoodR results. This experiment is used to match IDs with MGnify
#' data. (Default: \code{"metagenomic_amplicon"})
#' 
#' @param id.col1 \code{Character scalar}. Specifies the name of column from
#' \code{colData(x)} that includes HoloFood identifiers.
#' (Default: \code{"sample_biosample"})
#' 
#' @param id.col2 \code{Character scalar}. Specifies the name of column
#' from \code{colData(y[[exp.name2]])} that includes HoloFood identifiers.
#' (Default: \code{"accession"})
#' 
#' @param replace \code{Logical scalar}. Whether to replace the template
#' experiment. (Default: \code{TRUE})
#'
#' @param ... optional arguments not used currently.
#'
#' @return \code{MultiAssayExperiment}
#'
#' @examples
#' 
#' \dontrun{
#' # Get data from HoloFood database
#' mae <- HoloFoodR::getResult(
#'     salmon_sample_ids,
#'     use.cache = TRUE
#' )
#' 
#' # Get data from MGnify database
#' mg <- MgnifyClient(
#'     useCache = TRUE,
#'     cacheDir = ".MGnifyR_cache"
#' )
#' tse <- MGnifyR::getResult(
#'     mg,
#'     accession = mgnify_analyses_ids,
#'     get.func = FALSE
#' )
#' 
#' # Add MGnify data to HoloFood data
#' mae <- addMGnify(tse, mae)
#' }
#'
#' @name addMGnify
NULL

#' @rdname addMGnify
#' @export
setGeneric(
    "addMGnify", signature = c("x", "y"), function(x, y, ...)
    standardGeneric("addMGnify"))

#' @rdname addMGnify
#' @export
setMethod(
    "addMGnify",
    signature = c(x = "SummarizedExperiment", y = "MultiAssayExperiment"),
    function(
        x, y, exp.name1 = "metagenomic", exp.name2 = "metagenomic_amplicon",
        replace = TRUE, ...){
        ############################## Input check #############################
        .check_input(exp.name1, list("character scalar"))
        .check_input(
            exp.name2, list("character scalar"),
            supported_values = names(y))
        .check_input(replace, list("logical scalar"))
        ############################ Input check end ###########################
        # Merge results and get SE with combined metadata
        x <- addMGnify(x, y[[exp.name2]], ...)
        # Add the combined data to MAE
        # And add to MAE in place of old experiment
        x <- .add_mgnify_to_mae(x, y, exp.name1, exp.name2, replace)
        return(x)
    }
)

#' @rdname addMGnify
#' @export
setMethod(
    "addMGnify",
    signature = c(x = "SummarizedExperiment", y = "SummarizedExperiment"),
    function(x, y, ...){
        # Get metadata of HoloFood data
        y <- colData(y)
        # Combine data
        x <- addMGnify(x, y, ...)
        return(x)
    }
)

#' @rdname addMGnify
#' @export
setMethod(
    "addMGnify",
    signature = c(x = "SummarizedExperiment", y = "ANY"),
    function(
        x, y, id.col1 = "sample_biosample", id.col2 = "accession", ...){
        ############################## Input check #############################
        .check_input(y, list("data.frame", "DFrame"))
        if( is.null(rownames(y)) ){
            stop("'rownames(y)' must be populated.", call. = FALSE)
        }
        .check_input(
            id.col1, list("character scalar"),
            supported_values = colnames(colData(x)))
        .check_input(
            id.col2, list("character scalar"),
            supported_values = colnames(y))
        if( !( is.character(colData(x)[[id.col1]]) ||
                is.factor(colData(x)[[id.col1]]) ) ){
            stop("'id.col1' must specify IDs from 'x'.", call. = FALSE)
        }
        if( !( is.character(y[[id.col2]]) || is.factor(y[[id.col2]]) ) ){
            stop("'id.col2' must specify IDs from 'y'.", call. = FALSE)
        }
        ############################ Input check end ###########################
        # Ensure that sample matadata is data.frame
        y <- data.frame(y, check.names = FALSE)
        # Add metadata to TreeSE
        x <- .add_holofood_metadata(x, y, id.col1, id.col2)
        return(x)
    }
)

################################ HELP FUNCTIONS ################################

# This function adds TreeSE to MAE
.add_mgnify_to_mae <- function(tse, mae, exp.name1, exp.name2, replace){
    # Get sample map
    sample_map <- sampleMap(mae)
    # If user wants to replace the "template" experiment, remove it
    if( replace ){
        mae <- mae[, , !names(mae) %in% c(exp.name2) ]
    }
    # Get samples that are matching with those one that are being added
    add_sample_map <- sample_map[
        match(colnames(tse), sample_map[["colname"]]), ]
    # Rename the experiment
    add_sample_map[["assay"]] <- exp.name1
    # Add to sample map'
    sample_map <- rbind(sample_map, add_sample_map)
    # Create experiment list from TreeSE being added
    tse <- ExperimentList(temp = tse)
    names(tse) <- exp.name1
    # Combine experiment lists and create MAE
    exp_list <- c(experiments(mae), tse)
    mae <- MultiAssayExperiment(
        exp_list, colData = colData(mae), sampleMap = sample_map)
    return(mae)
}

# This function adds sample metadata to TreeSE.
.add_holofood_metadata <- function(tse, metadata, id.col1, id.col2){
    # Add old sample names to column metadata
    colData(tse)[["old_mgnify_colnames"]] <- colnames(tse)
    # Rename columns based on HoloFood ID
    colnames(tse) <- colData(tse)[[id.col1]]
    # Get sample metada6a from MGnify data
    cd <- data.frame(colData(tse), check.names = FALSE)
    # Combine sample metadata
    cd <- merge(
        cd, metadata,
        by.x = id.col1,
        by.y = id.col2,
        all.x = TRUE
    )
    # Now order colData to ensure that order is correct
    cd <- cd[match(colnames(tse), cd[[id.col1]]), ]
    # Add HoloFood identifiers as sample names
    rownames(cd) <- cd[[id.col1]]
    # Add it back to TreeSE
    colData(tse) <- DataFrame(cd)
    return(tse)
}
