# Query animals, samples, genome-catalogues, or viral-catalogues

doQuery <- function(
        type, max.hits = NULL, flatten = TRUE, ...){
    #
    # Check type
    supported_types <- c(
        "animals", "genome-catalogues", "samples", "viral-catalogues")
    temp <- .check_input(
        type, list("NULL", "character scalar"),
        supported_values = supported_types)
    # Check flatten
    temp <- .check_input(flatten, list("logical vector"))
    # Check max.hits
    temp <- .check_input(
        max.hits, list("NULL", "integer scalar"), limits = list(lower = 0))
    #
    # Get data
    res <- .retrieve_from_db(type, max.hits, ...)

    # Creae prececnce abscence table from the sample types (only for animals)
    # PParameter for user to specify
    # hidden
    col_name <- "sample_types"
    if( col_name %in% colnames(res) ){
        col <- res[[col_name]]
        res[[col_name]] <- NULL
        col <- .create_pa_table(col)
        res <- cbind(res, col)
    }

    # flatten if user has specified. DO this first and rhen sampletypes, because
    # sample types is good data to test.
    if( flatten ){
        is_list <- lapply(res, is.list)
        is_list <- unlist(is_list)
        is_list <- names(is_list)[ is_list ]
        for(col_name in is_list ){
            col <- res[[col_name]]
            res[[col_name]] <- NULL
            col <- .spread_column(col)
            colnames(col) <- paste0( col_name, seq_len(ncol(col)) )
            res <- cbind(res, col)
        }
    }

    return(res)
}

.create_pa_table <- function(col){
    uniq_types <- sort(unique(unlist(col)))

    res <- lapply(col, function(x){
        temp <- lapply(uniq_types, function(type) type %in% unlist(x) )
        temp <- unlist(temp)
        return(temp)
    })

    res <- do.call(rbind, res)
    res <- as.data.frame(res)
    colnames(res) <- uniq_types
    return(res)
}

.spread_column <- function(col){
    # Spread column by unlisting values
    col <- lapply(col, function(x) as.data.frame(t(data.frame(unlist((x))))) )
    # Create a data.frame from it
    col <- bind_rows(col)
    col <- as.data.frame(col)
    rownames(col) <- NULL
    return(col)
}
