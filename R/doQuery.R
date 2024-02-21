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

    # Creae prececnce abscence table from the sample types (only for animals)
    # PParameter for user to specify
    # hidden
    col_name <- "sample_types"
    if( col_name %in% colnames(res) ){
        res <- .create_pa_table(res, col_name, ...)
    }

    # flatten if user has specified. DO this first and rhen sampletypes, because
    # sample types is good data to test.
    if( flatten ){
        res <- .flatten_df(res)
    }

    return(res)
}

.create_pa_table <- function(res, col_name, spread.sample.types = TRUE, ...){
    # Check spread.sample.types
    temp <- .check_input(spread.sample.types, list("logical scalar"))
    #
    if( spread.sample.types ){
        col <- res[[col_name]]
        res[[col_name]] <- NULL

        uniq_types <- sort(unique(unlist(col)))

        col <- lapply(col, function(x){
            temp <- lapply(uniq_types, function(type) type %in% unlist(x) )
            temp <- unlist(temp)
            return(temp)
        })

        col <- do.call(rbind, col)
        col <- as.data.frame(res)
        colnames(col) <- uniq_types

        res <- cbind(res, col)
    }
    return(res)
}
