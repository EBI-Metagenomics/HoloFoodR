
getData <- function(
        type = NULL, accession.type = NULL, accession = NULL, max.hits = NULL,
        flatten = FALSE, ...){
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
    # Check max.hits
    temp <- .check_input(
        max.hits, list("NULL", "integer scalar"), limits = list(lower = 0))
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
    # Create paths
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
    res <- lapply(path, function(x){ .retrieve_from_db(x, max.hits, ...) })

    # If there were multiple results, accession was used --> combine data
    if( !is.null(accession) ){
        names(res) <- accession
        res <- .merge_data(res, ...)
    } else{
        res <- res[[1]]
    }

    # flatten if user has specified. DO this first and rhen sampletypes, because
    # sample types is good data to test.
    if( flatten ){
        if( !is.data.frame(res) ){
            res <- .join_datatypes(res)
        }

        if( is.data.frame(res) ){
            res <- .flatten_df(res)
        }
    }

    return(res)
}

############################ HELP FUNCTIONS ##################

.join_datatypes <- function(res, ...){
    not_empty <- lapply(res, function(x) nrow(x) > 0)
    not_empty <- unlist(not_empty)
    tab <- res
    if( any(not_empty) ){
        res <- res[ not_empty ]

        tab_names <- names(res)

        tab <- res[[1]]
        tab_name <- names(res)[[1]]
        if( tab_name %in% c("analysis_summaries") ){
            tab <- .add_datatype_to_colnames(tab, tab_name)
        }
        if(length(tab_names) > 1){
            for(tab_name in tab_names[2:length(tab_names)]){
                temp <- res[[ tab_name ]]
                #
                if( sum(colnames(temp) %in% colnames(tab)) > 1 ){
                    temp <- .add_datatype_to_colnames(temp, tab_name)
                }
                tab <- dplyr::full_join(tab, temp, by = "query_accession")
            }
        }
    } else{
        warning("Cannot join data.frames since they are all empty.", call. = FALSE)
    }
    return(tab)
}

.add_datatype_to_colnames <- function(df, tab_name, ex = "query_accession"){
    nams <- colnames(df)
    nams[ !nams %in% "query_accession" ] <- paste0(tab_name, ".", nams[ !nams %in% "query_accession" ])
    colnames(df) <- nams
    return(df)
}

.merge_data <- function(res, ...){
    # Remove empty elements
    res <- res[ lengths(res) > 0 ]
    # If there are still results left
    if( length(res) > 0 ){
        # Get datatypes from the first element
        datatypes <- names(res[[1]])
        # Loop through datatypes
        res <- lapply(datatypes, function(type){
            # Get the data from each element
            temp <- lapply(names(res), function(accession){
                # Get the data
                r <- res[[ accession ]][[ type ]]
                if( !is.null(r) ){
                    # If it is not data.frame, create
                    if( !is.data.frame(r) ){
                        r <- as.data.frame(r)
                        colnames(r) <- type
                    }
                    # Add query ID to df
                    r[["query_accession"]] <- accession
                }
                return(r)
            })
            # Combine
            temp <- bind_rows(temp)
            temp <- as.data.frame(temp)
            return(temp)
        })
        names(res) <- datatypes
    }
    return(res)
}
