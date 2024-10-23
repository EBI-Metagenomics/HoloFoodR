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
    # Based on number of acceptable classes, the msg is different
    class_txt <- .create_msg_from_list(classes_char)
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
        msg <- paste0(msg, " (Numeric constrains: ")
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
        msg <- paste0(msg, ")")
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
    if( "DFrame" %in% classes_char && is(variable, "DFrame") ){
        input_correct <- TRUE
    }
    if( "matrix" %in% classes_char && is.matrix(variable) ){
        input_correct <- TRUE
    }
    # If supported values were provided
    if( !is.null(supported_values) && !is.null(variable) ){
        # Test that if variable is in supported values
        values_correct <- lapply(supported_values, function(value){
            res <- FALSE
            if( is.null(value) && is.null(variable) || value %in% variable){
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

# This function creates a string from character values provided. The string
# can be used to messages. It creates a tidy list from list of values.
.create_msg_from_list <- function(classes_char, and_or = "or", ...){
    if( length(classes_char) > 2 ){
        class_txt <- paste0(
            paste(
                classes_char[seq_len(length(classes_char)-1)], collapse = ", "),
            " ", and_or, " ", classes_char[length(classes_char)])
    } else if( length(classes_char) == 2 ){
        class_txt <- paste0(
            classes_char[[1]], " ", and_or, " ", classes_char[[2]])
    } else{
        class_txt <- classes_char
    }
    return(class_txt)
}

############################ OTHER COMMON FUNCTIONS ############################

# This function retrieves data from the database
#' @importFrom dplyr bind_rows
.retrieve_from_db <- function(path, max.hits = NULL, ...){
    # Check max.hits
    temp <- .check_input(
        max.hits, list(NULL, "integer scalar"), limits = list(lower = 0))
    #
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

# This function performs a single query to database (or gets the data from
# cache) and returns a data.frame
#' @importFrom httr2 request req_url_query req_error req_perform
#' @importFrom httr2 resp_body_string
#' @importFrom jsonlite fromJSON
.perform_single_query <- function(
        path, use.cache = FALSE, cache.dir = tempdir(), clear.cache = FALSE,
        base.url = "https://www.holofooddata.org/api", full.url = NULL, ...){
    # Check base.url
    temp <- .check_input(base.url, list("character scalar"))
    # Check use.cache
    temp <- .check_input(use.cache, list("logical scalar"))
    # Check cache.dir
    temp <- .check_input(cache.dir, list("character scalar"))
    # Check clear.cache
    temp <- .check_input(clear.cache, list("logical scalar"))
    # Check full.url
    temp <- .check_input(full.url, list("NULL", "character scalar"))
    #
    # Create url address
    if( is.null(full.url) ){
        url <- paste0(base.url, "/", path)
    } else{
        url <- full.url
    }
    # Get query options
    query_params <- list(...)
    query_params <- lapply(query_params, function(x) paste(x, collapse = ","))
    query_params <- as.list(c(format="json", query_params))

    # Get cache name
    if( use.cache || clear.cache ){
        # It is built from query info
        temp <- unlist(query_params)
        # If user has specified the subdirectory, ensure that it works in any
        # system by adding correct "/".
        cache_dir <- as.list(strsplit(cache.dir, "[/\\\\]")[[1]])
        cache_dir <- do.call(file.path, cache_dir)
        # Construct directory path
        cache_dir <- file.path(cache_dir, "HoloFoodR_cache")
        # Construct file path
        cache_path <- c(path, names(temp), temp)
        cache_path <- paste(cache_path, collapse = "_")
        cache_path <- gsub(":|/", "_", cache_path)
        cache_path <- paste0(cache_path, ".RDS")
        cache_path <- file.path(cache_dir, cache_path)
    }
    # Remove the file from path if specified
    if( clear.cache ){
        if( file.exists(cache_path) ){
            message("clear.cache is TRUE: deleting file '", cache_path, "'")
            unlink(cache_path)
        }
    }

    # If user wants to use the cache, try to find the file first from cache
    if( use.cache && file.exists(cache_path) ){
        res <- readRDS(cache_path)
    } else{
        # Get data from database
        req <- request(url)
        req <- req_url_query(req, !!!query_params)
        # Suppress error
        req <- req_error(req, is_error = \(resp) FALSE)
        res <- req_perform(req)
        # If data was fetched successfully
        if( res$status_code == "200" ){
            # Get the data as JSON text
            res <- resp_body_string(res)
            # Convert into list of data.frames
            res <- fromJSON(res, flatten = TRUE)
            # Add the result to cache if specified
            if( use.cache ){
                if( !dir.exists(cache_dir) ){
                    dir.create(cache_dir, recursive = TRUE)
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

# Some columns/datatypes from database have multiple values so they are lists.
# This function goes through those columns and spread them so that each cell has
# only one value.
.flatten_df <- function(res){
    # Get which column is a list
    is_list <- lapply(res, is.list)
    is_list <- unlist(is_list)
    is_list <- names(is_list)[ is_list ]
    # Loop through those columns
    for(col_name in is_list ){
        # Get te column and remove it from the original table
        col <- res[[col_name]]
        res[[col_name]] <- NULL
        # Create a data.frame with multiple columns from the column
        col <- .spread_column(col)
        # Add info on what column was spread
        colnames(col) <- paste0(col_name, ".", colnames(col))
        res <- cbind(res, col)
    }
    return(res)
}

# This function creates a data.frame with multiple columns from column that
# is a list.
#' @importFrom dplyr bind_rows
.spread_column <- function(col){
    # Spread column by unlisting values
    col <- lapply(col, function(x){
        # Get priginal titles
        orig_names <- names(x)
        # Flatten the column
        x <- as.data.frame(t(data.frame(unlist((x)))))
        # If there were multiple values for each row, the titles are now in
        # format column* where * is an integer denoting the number of values.
        if( !is.null(orig_names) && any(!orig_names %in% colnames(x)) ){
            # Seach first values --> they have column1
            names(orig_names) <- paste0(orig_names, 1)
            # Loop through new columns and switch the names
            new_names <- colnames(x)
            new_names <- lapply(new_names, function(name){
                # Check what original name new name represents
                ind <- lapply(
                    names(orig_names), function(og_name) og_name == name )
                ind <- unlist(ind)
                # Replace
                if( any(ind) ){
                    name <- orig_names[ind]
                }
                return(name)
            })
            new_names <- unlist(new_names)
            colnames(x) <- new_names
        } else if( is.null(orig_names) ){
            # If there were no names in the beginning, columns are just named
            # like "V*" --> replace so that they have only number.
            new_names <- as.character(seq_len(ncol(x)))
            colnames(x) <- new_names
        }
        rownames(x) <- NULL
        return(x)
    })
    # Create a data.frame from it
    col <- bind_rows(col)
    col <- as.data.frame(col)
    rownames(col) <- NULL
    return(col)
}

# This function merges lists of data.frames with full_join into single
# data.frame
#' @importFrom dplyr full_join
.full_join_list <- function(res){
    # Remove empty elements
    res <- res[ lengths(res) > 0 ]
    # If there is more than one element, merge them
    if( length(res) > 1 ){
        df <- Reduce(function(df1, df2){
            # Get common columns
            common_cols <- intersect(colnames(df1), colnames(df2))
            # Merge based on common columns
            temp <- full_join(df1, df2, by = common_cols)
            return(temp)
        }, res)
    } else if( length(res) == 1 ){
        # Otherwise if there is only one element, give the element
        df <- res[[1]]
    } else{
        # If all the data.frames were without information, give NULL
        df <- NULL
    }
    return(df)
}
