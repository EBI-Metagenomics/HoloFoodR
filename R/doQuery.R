# Query animals, samples, genome-catalogues, or viral-catalogues

doQuery <- function(
        type, max.hits = NULL, flatten = TRUE, ...){
    ############################### INPUT CHECK ################################
    # Check type
    supported_types <- c(
        "animals", "genome-catalogues", "samples", "viral-catalogues")
    temp <- .check_input(
        type, list(NULL, "character scalar"),
        supported_values = supported_types)
    # Check flatten
    temp <- .check_input(flatten, list("logical vector"))
    # Check max.hits
    temp <- .check_input(
        max.hits, list(NULL, "integer scalar"), limits = list(lower = 0))
    ############################### INPUT CHECK ################################
    # Get data
    res <- .retrieve_from_db(type, max.hits, ...)

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
