################################### TESTING ###################################
# Methods for testing

.is_a_bool <- function(x){
    is.logical(x) && length(x) == 1L && !is.na(x)
}

.is_non_empty_character <- function(x){
    is.character(x) && all(nzchar(x))
}

.is_non_empty_string <- function(x){
    is.character(x) && length(x) == 1L
}

.is_a_numeric <- function(x){
    is.numeric(x) %% length(x) == 1L
}

.is_integer <- function(x){
    all(is.numeric(x) && x%%1==0)
}

.is_an_integer <- function(x){
    .is_integer(x) && length(x) == 1L
}

.get_name_in_parent <- function(x) {
    .safe_deparse(do.call(substitute, list(substitute(x), parent.frame())))
}

.safe_deparse <- function (expr, ...) {
    paste0(deparse(expr, width.cutoff = 500L, ...), collapse = "")
}

# This function unifies input testing. The message will always be in same format
# also it makes the code simpler in main function since testing is done here.
.check_input <- function(
        variable, supported_class, supported_values = NULL, limits = NULL,
        variable_name = .get_name_in_parent(variable)){
    # Convert supported classes to character
    classes_char <- lapply(supported_class, function(class){
        if( is.null(class) ){
            class <- "NULL"
        }
        return(class)
    })
    classes_char <- unlist(classes_char)
    # Based on length of the classes, the msg is different
    if( length(classes_char) > 2 ){
        class_txt <- paste0(
            paste(classes_char[1:(length(classes_char)-1)], collapse = ", "),
            " or ", classes_char[length(classes_char)])
    } else if( length(classes_char) == 2 ){
        class_txt <- paste0(classes_char[[1]], " or ", classes_char[[2]])
    } else{
        class_txt <- classes_char
    }
    # Create a message
    msg <- paste0("'", variable_name, "' must be ", class_txt, "." )

    # If supported values were provided
    if( !is.null(supported_values) ){
        # Convert supported values to character
        values_char <- lapply(supported_values, function(value){
            if( is.null(value) ){
                value <- "NULL"
            }
            value <- as.character(value)
            return(value)
        })
        values_char <- unlist(values_char)
        # Collapse into text
        values_txt <- paste0("'", paste(values_char, collapse = "', '"), "'")
        msg <- paste0(
            msg, " It must be one of the following options: ", values_txt)
    }

    # If limits were provided
    if( !is.null(limits) ){
        msg <- paste0(msg, " It must be ")
        # Add thresholds to message
        if( !is.null(limits$upper) ){
            msg <- paste0(msg, limits$upper, ">x")
        } else if(!is.null(limits$upper_include)){
            msg <- paste0(msg, limits$upper, ">=x")
        }
        if( !is.null(limits$lower) ){
            msg <- paste0(msg, "x>", limits$lower)
        } else if(!is.null(limits$lower_include)){
            msg <- paste0(msg, "x>=", limits$lower_include)
        }
    }

    # List all the input types. Run the check if the variable must be that type.
    # If correct type was found, change the result to TRUE.
    input_correct <- FALSE
    if( "NULL" %in% classes_char && is.null(variable) ){
        input_correct <- TRUE
    }
    if( "logical scalar" %in% classes_char && .is_a_bool(variable) ){
        input_correct <- TRUE
    }
    if( "logical vector" %in% classes_char && is.logical(variable) ){
        input_correct <- TRUE
    }
    if( "character scalar" %in% classes_char && .is_non_empty_string(
            variable) ){
        input_correct <- TRUE
    }
    if( "character vector" %in% classes_char && .is_non_empty_character(
            variable) ){
        input_correct <- TRUE
    }
    if( "numeric scalar" %in% classes_char && .is_a_numeric(variable) ){
        input_correct <- TRUE
    }
    if( "numeric vector" %in% classes_char && is.numeric(variable) ){
        input_correct <- TRUE
    }
    if( "integer vector" %in% classes_char && .is_integer(variable) ){
        input_correct <- TRUE
    }
    if( "integer scalar" %in% classes_char && .is_an_integer(variable) ){
        input_correct <- TRUE
    }
    if( "list" %in% classes_char && is.list(variable) && !is.data.frame(
            variable) ){
        input_correct <- TRUE
    }
    if( "data.frame" %in% classes_char && is.data.frame(variable) ){
        input_correct <- TRUE
    }
    if( "matrix" %in% classes_char && is.matrix(variable) ){
        input_correct <- TRUE
    }
    # If supported values were provided
    if( !is.null(supported_values) ){
        # Test that if variable is in supported values
        values_correct <- lapply(supported_values, function(value){
            res <- FALSE
            if( is.null(value) && is.null(variable) || value == variable){
                res <- TRUE
            }
            return(res)
        })
        values_correct <- unlist(values_correct)
        # If not, then give FALSE even though class checks were correct
        if( !any(values_correct) ){
            input_correct <- FALSE
        }
    }
    # If limits were provided
    if( !is.null(limits) && !is.null(variable) ){
        if( !is.null(limits$upper) && variable >= limits$upper ){
            input_correct <- FALSE
        } else if( !is.null(
                limits$upper_include) && variable > limits$upper_include ){
            input_correct <- FALSE
        }

        if( !is.null(limits$lower) && variable <= limits$lower ){
            input_correct <- FALSE
        } else if( !is.null(
                limits$upper_include) && variable < limits$upper_include ){
            input_correct <- FALSE
        }
    }
    # Give error if variable was not correct type
    if( !input_correct ){
        stop(msg, call. = FALSE)
    }
    return(input_correct)
}


############################ OTHER COMMON FUNCTIONS ############################

# This function retrieves data from the database
.retrieve_from_db <- function(path, max.hits = NULL, ...){
    # Get data
    res <- .perform_single_query(path, ...)
    # Check if there are counts info --> how many elements there are in database
    count_element <- "count"
    if( count_element %in% names(res) ){
        data_element <- names(res)[ !names(res) %in% count_element]
        # Get number of pages
        pages <- ceiling(res[[count_element]] / nrow(res[data_element][[1]]))
        # Get the data and drop pages info
        res <- res[ data_element ]
        # Get number of items that we have
        num_of_items <- nrow(res[[1]])
        # If there are multiple pages and max hits do no exceed the number that
        # we already have
        if( pages > 1 && (is.null(max.hits) || (
                !is.null(max.hits) && num_of_items < max.hits)) ){
            # Loop through pages
            for( page in 2:pages ){
                # Fetch data
                res2 <- .perform_single_query(path, page=page, ...)
                # Get only data without count information
                res2 <- res2[ data_element ]
                # Combine data
                res <- lapply(data_element, function(elem){
                    bind_rows(res[[elem]], res2[[elem]])})
                names(res) <- data_element
                # Check if we have exceeded the maximum hits
                num_of_items <- nrow(res[[1]])
                if( (!is.null(max.hits) && num_of_items >= max.hits) ){
                    break
                }
            }
        }
        # If there was only one data element get only it
        if( length(data_element) == 1 ){
            res <- res[[ data_element ]]
        }
    }
    return(res)
}

.perform_single_query <- function(path, use.cache = TRUE, cache.dir = tempdir(), clear.cache = FALSE, base.url = "https://www.holofooddata.org/api", ...){
    # Check base.url
    temp <- .check_input(base.url, list("character scalar"))
    # Check use.cache
    temp <- .check_input(use.cache, list("logical scalar"))
    # Check cache.dir
    temp <- .check_input(cache.dir, list("character scalar"))
    # Check clear.cache
    temp <- .check_input(clear.cache, list("logical scalar"))
    #
    # Create url address
    url <- paste0(base.url, "/", path)
    # Get query options
    query_params <- list(...)
    query_params <- lapply(query_params, function(x) paste(x, collapse = ","))
    query_params <- as.list(c(format="json", query_params))

    # Get cache name
    if( use.cache || clear.cache ){
        # It is built from query info
        temp <- unlist(query_params)
        cache_path <- c(cache.dir, names(temp), temp)
        cache_path <- paste(cache_path, collapse = "_")
        cache_path <- gsub(":", "_", cache_path)
        cache_path <- paste0(cache_path, ".RDS")
    }
    # Remove the file from path if specified
    if( clear.cache ){
        if( file.exists(cache_path) ){
            message("clear.cache is TRUE: deleting ", cache_path)
            unlink(cache_path)
        }
    }

    # If user wants to use the cache, try to find the file first from cache
    if( use.cache && file.exists(cache_path) ){
        res <- readRDS(cache_path)
    } else{
        # Get data from database
        req <- httr2::request(url)
        req <- httr2::req_url_query(req, !!!query_params)
        res <- httr2::req_perform(req)
        # If data was fetched successfully
        if( res$status_code == "200" ){
            # Get the data as JSON text
            res <- httr2::resp_body_string(res)
            # Convert into list of data.frames
            res <- jsonlite::fromJSON(res, flatten = TRUE)
            # Add the result to cache if specified
            if( use.cache ){
                if( !dir.exists(cache.dir) ){
                    dir.create(cache.di)
                }
                saveRDS(res, cache_path)
            }
        } else{
            # If data was not found, give warning
            warning(
                res$status_code, " error in query '", res$url, "'",
                call. = FALSE)
            res <- NULL
        }
    }
    return(res)
}
