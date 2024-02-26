#' Get data from HoloFood database
#'
#' @details
#' With \code{getData}, you can fetch data from the database. Compared to
#' \code{getResult}, this function is more flexible since it can fetch any kind
#' of data from the database. However, this function returns the data
#' without further wrangling as \code{list} or \code{data.frame} which are not
#' optimized format for fetching data on samples.
#'
#' Search results can be filtered; for example, animals can be filtered based on
#' available samples. See [Api browser](https://www.holofooddata.org/api/docs)
#' for information on filters. You can find help on customizing queries from
#' [here](https://emg-docs.readthedocs.io/en/latest/api.html#customising-queries).
#'
#' @param type \code{NULL} or \code{character scalar} specifying the type of
#' data to query. Must be one of the following options:
#' \code{analysis-summaries}, \code{animals}, \code{genome-catalogues},
#' \code{samples}, \code{sample_metadata_markers} or \code{viral-catalogues}.
#' When genome or viral catalogues is fetched by their accession ID, the type
#' can also be \code{genomes} or \code{fragments}. (Default: \code{NULL})
#'
#' @param accession.type \code{NULL} or \code{character scalar} specifying the
#' type of accession IDs. Must be one of the following options: \code{animals},
#' \code{genome-catalogues}, \code{samples} or \code{viral-catalogues}.
#' (Default: \code{NULL})
#'
#' @param accession \code{NULL} or \code{character vector} specifying the
#' accession IDs of type \code{accession.type}. (Default: \code{NULL})
#'
#' @param flatten \code{Logical scalar} specifying whether to flatten the
#' resulting \code{data.frame}. This means that columns with multiple values
#' are separated to multiple columns. (Default: \code{FALSE})
#'
#' @param ... optional arguments:
#' \itemize{
#'
#'   \item{\code{max.hits}: }{\code{NULL} or \code{integer scalar} specifying the
#'   maximum number of results to fetch. When NULL, all results are fetched.
#'   (Default: \code{NULL})}
#'
#'   \item{\code{use.cache}: }{\code{Logical scalar} specifying whether to
#'   use.cache (Default: \code{FALSE})}
#'
#'   \item{\code{cache.dir}: }{\code{Character scalar} specifying cache directory.
#'   (Default: \code{tempdir()})}
#'
#'   \item{\code{clear.cache}: }{\code{Logical scalar} specifying whether to
#'   use.cache (Default: \code{FALSE})}
#'
#' }
#'
#' @return \code{list} or \code{data.frame}
#'
#' @examples
#'
#' # Find genome catalogues
#' catalogues <- getData(type = "genome-catalogues")
#' head(catalogues)
#'
#' # Find genomes based on certain genome catalogue iD
#' res <- getData(
#'     type = "genomes", accession.type = "genome-catalogues",
#'     accession = catalogues[1, "id"], max.hits = 100)
#' head(res)
#'
#' @seealso
#' \code{\link[HoloFoodR:getResult]{getResult}}
#'
#' @name getData
NULL

#' @rdname getData
#' @export
getData <- function(
        type = NULL, accession.type = NULL, accession = NULL, flatten = FALSE,
        ...){
    ############################### INPUT CHECK ################################
    # Check type
    supported_types <- c(
        "animals", "samples", "sample_metadata_markers", "analysis-summaries",
        "genome-catalogues", "viral-catalogues")
    if( !is.null(accession) ){
        supported_types <- c(supported_types, "genomes", "fragments")
    }
    temp <- .check_input(
        type, list(NULL, "character scalar"),
        supported_values = supported_types)
    # Check accession.type
    supported_accession_types <- c(
        "animals", "samples", "genome-catalogues", "viral-catalogues")
    temp <- .check_input(
        accession.type, list(NULL, "character scalar"),
        supported_values = supported_accession_types)
    # Check accession
    temp <- .check_input(accession, list(NULL, "character vector"))
    # Check flatten
    temp <- .check_input(flatten, list("logical vector"))
    # Every query parameter cannot be NULL
    if( is.null(type) && is.null(accession) && is.null(accession.type) ){
        stop(
            "'type', 'accession' and 'accession.type' cannot be NULL.",
            call. = FALSE)
    }
    #
    # If both accession.type and accession was not specified
    if( !((is.null(accession.type) && is.null(accession)) ||
            (!is.null(accession.type) && !is.null(accession))) ){
        stop(
            "Both 'accession.type' and 'accession' must be specified or be ",
            "NULL.", call. = FALSE)
    }
    # type can only be genomes and fragments if accession is specified -->
    # accession specifies genome or viral catalog
    if( !is.null(accession) &&
            !(type %in% c("genomes", "fragments") || is.null(type)) ){
        stop(
            "'type' must be 'genomes' or 'fragments' when accession is ",
            "specified.", call. = FALSE)
    }
    ############################### INPUT CHECK ################################
    # Create paths. If there are multiple accessions specified, this step
    # results to multiple paths
    if( !is.null(accession) ){
        path <- paste0(accession.type, "/", accession)
        if( !is.null(type) ){
            path <- paste0(path, "/", type)
        }
    } else{
        path <- type
    }
    # Loop through paths i.e. accession IDs (or if accession is not specified
    # retrieve data only once) and get data
    res <- lapply(path, function(x){ .retrieve_from_db(x, ...) })

    # If there were multiple accessions, there are multiple results -->
    # combine data from multiple results/accessions. The tables are combined so
    # that each table type is preserved separately meaning that the result is
    # not necessarily a single data.frame but a list of data.frames, each
    # including multiple results.
    if( !is.null(accession) && accession.type %in% c("animals", "samples") ){
        names(res) <- accession
        res <- .merge_data(res, ...)
        # flatten if user has specified.
        if( flatten ){
            # The result can be a list of data.frames. In order to flatten the
            # data it must be first combined into single data.frame.
            if( !is.data.frame(res) ){
                res <- .join_datatypes(res)
            }
            # Now we can flatten the data --> collapse columns that are lists
            # into multiple columns
            res <- .flatten_df(res)
        }
    } else{
        # Otherwise get the single result.
        res <- res[[1]]
    }
    return(res)
}

################################ HELP FUNCTIONS ################################

# There might be multiple data.frames in results, each representing unique
# datatype. This function combines these datatypes into single data.frame
.join_datatypes <- function(res, ...){
    # Check whether data types are empty
    not_empty <- lengths(res) > 0
    # If there are dfs that have info
    if( any(not_empty) ){
        # Get non-empty data.frames / data types
        res <- res[ not_empty ]
        datatypes <- names(res)

        # If the data type has column including analysis summaries, add
        # info to column names. That is because analysis summaries have
        # column names that are also in main data (such as title).
        analysis_tab <- "analysis_summaries"
        if( analysis_tab %in% datatypes ){
            temp <- res[[analysis_tab]]
            temp <- .add_datatype_to_colnames(temp, analysis_tab)
            res[[analysis_tab]] <- temp
        }

        # Get first data type
        tab <- res[[1]]
        tab_name <- names(res)[[1]]
        # If there were more than 1 data type
        if(length(datatypes) > 1){
            # Loop through types
            for( type in datatypes[2:length(datatypes)] ){
                temp <- res[[ type ]]
                # If the column names are already present in the data, add
                # data type to column names that are being added.
                if( sum(colnames(temp) %in% colnames(tab)) > 1 ){
                    temp <- .add_datatype_to_colnames(temp, type)
                }
                # Add data based on accession ID
                tab <- dplyr::full_join(tab, temp, by = "query_accession")
            }
        }
    } else{
        warning(
            "Cannot join data.frames since they are all empty.", call. = FALSE)
        tab <- res
    }
    return(tab)
}

# Add data type to column names. Exclude column that has accession IDs (do not
# add data type to it).
.add_datatype_to_colnames <- function(df, tab_name, ex = "query_accession"){
    #
    temp <- .check_input(ex, list("character scalar"))
    #
    nams <- colnames(df)
    nams[ !nams %in% "query_accession" ] <- paste0(
        tab_name, ".", nams[ !nams %in% "query_accession" ])
    colnames(df) <- nams
    return(df)
}

# This function loops through results and combine the tables so that each
# data type still has unique table but the table has now all accessions.
#' @importFrom dplyr bind_rows
.merge_data <- function(res, ...){
    # Remove empty elements
    res <- res[ lengths(res) > 0 ]
    # If there are still results left
    if( length(res) > 0 ){
        # Get datatypes
        datatypes <- unique(unlist(lapply(res, names)))
        # Loop through datatypes
        res <- lapply(datatypes, function(type){
            # Loop through each accession
            temp <- lapply(names(res), function(accession){
                # Get the data
                r <- res[[ accession ]][[ type ]]
                if( !is.null(r) ){
                    # If it is not data.frame, create data.frame from it
                    if( !is.data.frame(r) ){
                        r <- as.data.frame(r)
                        colnames(r) <- type
                    }
                    # Add query ID to df
                    r[["query_accession"]] <- accession
                }
                return(r)
            })
            # Combine accessions into single data.frame
            temp <- bind_rows(temp)
            temp <- as.data.frame(temp)
            return(temp)
        })
        names(res) <- datatypes
    }
    return(res)
}
