#' Search HoloFood database for animals, genome catalogues, samples,
#' or viral catalogues
#'
#' @details
#' \code{doQuery} is a flexible query function which can be utilized to search
#' available animals, genome catalogues, samples, or viral catalogues. Search
#' results can be filtered; for example, animals can be filtered based on
#' available samples. See [Api browser](https://www.holofooddata.org/api/docs)
#' for information on filters. You can find help on customizing queries from
#' [here](https://emg-docs.readthedocs.io/en/latest/api.html#customising-queries).
#'
#' @param type \code{Character scalar} specifying the type of data to
#' query. Must be one of the following options: \code{"animals"},
#' \code{"genome-catalogues"}, \code{"samples"} or \code{"viral-catalogues"}.
#'
#' @param flatten \code{Logical scalar} specifying whether to flatten the
#' resulting \code{data.frame}. This means that columns with multiple values
#' are separated to multiple columns. (Default: \code{TRUE})
#'
#' @param ... optional arguments:
#' \itemize{
#' 
#'   \item \strong{max.hits} \code{NULL} or \code{integer scalar} specifying the
#'   maximum number of results to fetch. When NULL, all results are fetched.
#'   (Default: \code{NULL})
#'   
#'   \item \strong{spread.sample.types} \code{Logical scalar} specifying whether
#'   to create spread sample types column of animals data. In animals data,
#'   sample types column might have multiple values that might be hard to
#'   explore. This argument specifies whether to create presence/absence table
#'   from sample types. (Default: \code{TRUE})
#'   
#'   \item \strong{use.cache} \code{Logical scalar} specifying whether to
#'   use cache. (Default: \code{FALSE})
#'   
#'   \item \strong{cache.dir} \code{Character scalar} specifying cache
#'   directory. (Default: \code{tempdir()})
#'   
#'   \item \strong{clear.cache} \code{Logical scalar} specifying whether to
#'   remove and clear cache (Default: \code{FALSE})
#'   
#' }
#'
#' @return \code{data.frame}
#'
#' @examples
#'
#' # Find animals results. The maximum amount of results is 100. Use filter
#' # so that only chicken is searched.
#' res <- doQuery("animals", max.hits = 100, system = "chicken")
#' head(res)
#'
#' @name doQuery
NULL

#' @rdname doQuery
#' @export
doQuery <- function(
        type, flatten = TRUE, ...){
    ############################### INPUT CHECK ################################
    # Check type
    supported_types <- c(
        "animals", "genome-catalogues", "samples", "viral-catalogues")
    temp <- .check_input(
        type, list(NULL, "character scalar"),
        supported_values = supported_types)
    # Check flatten
    temp <- .check_input(flatten, list("logical vector"))
    ############################### INPUT CHECK ################################
    # Get data
    res <- .retrieve_from_db(type, ...)

    # Create presence/absence table from the sample types (only for animals).
    col_name <- "sample_types"
    if( col_name %in% colnames(res) ){
        res <- .create_pa_table(res, col_name, ...)
    }
    # flatten if user has specified. (Each cell has only one value; it might
    # be that some columns are lists at this point)
    if( flatten ){
        res <- .flatten_df(res)
    }
    return(res)
}

################################ HELP FUNCTIONS ################################

# This function gets "sample.types" column from animals table. The column is a
# list that contains all the sample types that each animal has. Convert this
# list to table where columns are unique sample types, rows animals and cells
# are TRUE or FALSE stating whether sample type is present for certain animal.
.create_pa_table <- function(res, col_name, spread.sample.types = TRUE, ...){
    # Check spread.sample.types
    temp <- .check_input(spread.sample.types, list("logical scalar"))
    #
    # If user wants to spread sample types column
    if( spread.sample.types ){
        # Get column and remove it from the original table
        col <- res[[col_name]]
        res[[col_name]] <- NULL
        # Get unique sample types
        uniq_types <- sort(unique(unlist(col)))
        # Loop through animals
        col <- lapply(col, function(x){
            # Is certain sample types present?
            temp <- lapply(uniq_types, function(type) type %in% unlist(x) )
            temp <- unlist(temp)
            return(temp)
        })
        # Create a df from list of vectors containing boolean values
        col <- do.call(rbind, col)
        col <- as.data.frame(col)
        colnames(col) <- uniq_types
        # Add to original table
        res <- cbind(res, col)
    }
    return(res)
}
