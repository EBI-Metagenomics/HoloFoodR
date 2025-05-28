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
#' @param output \code{character scalar} specifying output format. Must be
#' \code{"list"}, \code{"TreeSE"} (\code{TreeSummarizedExperiment}) or
#' \code{"SE"} (\code{SummarizedExperiment}). (Default: \code{"list"})
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
#'   \item \strong{ion.mode} \code{Character scalar} specifying metabolite
#'   assignment files to fetch. If \code{"positive"} only positive ions are
#'   fetched. Similarly \code{"negative"} means that negative ions are fetched
#'   if such data exists. By selecting \code{"both"}, one can fetch both
#'   positive and negative ions. (Default: \code{"both"})
#'
#' }
#'
#' @return \code{list}, \code{SummarizedExperiment} or
#' \code{TreeSummarizedExperiment}
#'
#' @examples
#'
#' # This example is not run, because the server fails to respond sometimes.
#' if( FALSE ){
#'     res <- getMetaboLights("MTBLS4381")
#'     file_paths <- getMetaboLightsFile(
#'         study.id = "MTBLS4381",
#'         file = res[["assay_meta"]][["Raw Spectral Data File"]]
#'         )
#'     # Get data as SummarizedExperiment
#'     se <- getMetaboLights("MTBLS3540", output = "SE")
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
getMetaboLights <- function(study.id, output = "list", ...){
    # Check study.id and output
    temp <- .check_input(study.id, list("character vector"))
    temp <- .check_input(
        output, list("character scalar"),
        supported_values = c("list", "TreeSE", "SE"))
    #
    # Remove trailing spaces from study.id
    study.id <- study.id |> trimws()
    # Retrieve data as list
    res <- .retrieve_metabolights_data(study.id, ...)
    # If user wants to convert the data into SE or TreeSE
    if( output %in% c("TreeSE", "SE") ){
        res <- .construct_metabolomic_SE(res, output, ...)
    }
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
    # Remove trailing spaces from study.id
    study.id <- study.id |> trimws()
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

# This function facilitates retrieval of files from MetaboLights
.retrieve_metabolights_data <- function(study.id, ...){
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

# This function retrieves metabolomic data from MetaboLights database for single
# URL address
#' @importFrom dplyr left_join
.get_metabolomic_data <- function(url, ion.mode = "both", ...){
    # ion.mode specifies whether to fetch positive or negative ions or both
    temp <- .check_input(
        ion.mode, list("character scalar"),
        list("both", "positive", "negative"))
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
        file_names <- file_names[ !file_names %in% c("", NA, " ") ]
        # If user specified pattern to fetch by, check which file names match.
        if( ion.mode != "both" ){
            file_names <- file_names[ grepl(
                ion.mode, file_names, ignore.case = TRUE) ]
        }
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
.get_study_info <- function(
        url,
        study.search.url = "https://www.ebi.ac.uk/metabolights/ws/studies",
        ...){
    # Check if study.id is already a url address. If it is not, create url
    # from study.id and base.url
    if( !grepl("https://www.ebi.ac.uk", url, ignore.case = TRUE) ){
        url <- paste0(study.search.url, "/", url)
    }
    # From the metabolights database, find associated study. Which study
    # represents this HoloFood study?
    res <- .perform_single_query(path = "metabolight", full.url = url, ...)
    # Check if data was found
    if( is.null(res) ){
        stop("No data was found for the following URL: '", url, "'",
            call. = FALSE)
    }
    # Get only relevant info
    study_info <- res[["isaInvestigation"]][["studies"]]
    return(study_info)
}

# This is a common function for downloading a file from MetaboLights database
#' @importFrom utils download.file read.delim URLencode
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
    # If the study.id is url, get the study.id from the url. This works also
    # even if just study ID was provided without url.
    temp <- strsplit(study.id, "/")[[1]]
    study.id <- temp[length(temp)]
    # Create url
    url <- paste0(metabolights.base.url, "/", study.id, "/", file.name)
    # Some file names have spaces. Replace them with accepted character
    url <- URLencode(url)
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
        # Load the data with try catch. If the file is not found, give warning
        # instead of error.
        tryCatch({
            download.file(url, file_path, quiet = FALSE, timeout = timeout)
        }, error = function(e) {
            warning(conditionMessage(e), call. = FALSE)
        })
        # Set the timeout back to default
        options(timeout = def_opt)
    }
    # By default, the loaded table is returned. However, for spectra files, we
    # do not want to return them.
    if( file.exists(file_path) && return.table ){
        # Get the encoding of a file
        encodings <- .detect_encoding(file_path)
        encodings <- c("UTF-8", "latin1", "windows-1252", encodings) |> unique()
        # Read the local file. Try different encodings. Sometimes UTF-8 fails.
        df <- NULL
        i <- 1L
        while( (is.null(df) || nrow(df) == 0L) && i <= length(encodings) ){
            df <- tryCatch({
                read.delim(
                    file_path, check.name = FALSE, row.names = NULL,
                    fileEncoding = encodings[[i]])
                },
                error = function(e){
                    return(NULL)
                },
                warning = function(w){
                    return(NULL)
                }
            )
            # Increment to next encoding
            i <- i + 1
        }
        # If we were able to read the table
        if( !is.null(df) ){
            # Make column names unique if specified
            if( anyDuplicated(colnames(df)) && unique.cols ){
                colnames(df) <- make.unique(colnames(df))
            }
            # Add info from which file the data comes from
            df[["metabolights_url"]] <- url
            df[["file_name"]] <- basename(url)
            # If the file is metabolite assignment file, add information whether
            # the metabolite is positive or negative ion.
            if( grepl("^m_.*maf.*\\.tsv$", file.name) ){
                df[["ion_mode"]] <- .get_ion_mode(file.name)
            }
        }
    } else if( file.exists(file_path) ){
        df <- file_path
    } else{
        df <- NULL
    }
    return(df)
}

# This function detects encoding of the file
#' @importFrom stringi stri_enc_detect
.detect_encoding <- function(file_path){
    # Read characters
    characters <- rawToChar(readBin(file_path, "raw", 10000))
    # Detect the most possible encodings
    encoding <- stri_enc_detect(characters)[[1L]]
    # Return the most probable encoding
    encoding <- encoding[["Encoding"]]
    return(encoding)
}

# Identify ion mode based on filename
.get_ion_mode <- function(filename){
    res <- "unknown"
    if( grepl("MS_positive", filename, ignore.case = TRUE) ){
        res <- "positive"
    } else if( grepl("MS_negative", filename, ignore.case = TRUE) ){
        res <- "negative"
    }
    return(res)
}

# This function constucts TreeSE object from retrieved MetaboLights data
#' @importFrom SummarizedExperiment SummarizedExperiment
.construct_metabolomic_SE <- function(res, output, assay.type = "conc", ...){
    # Check assay.type
    temp <- .check_input(assay.type, list("character scalar"))

    # Check that we have all the data. Only metabolite assignment file is
    # required but give also warning on other missing files.
    missing <- c()
    if( is.null(res[["assay"]]) ){
        missing <- c(missing, "metabolite assignment ('m_' prefix)")
    }
    if( is.null(res[["assay_meta"]]) ){
        missing <- c(missing, "assay metadata ('a_' prefix)")
    }
    if( is.null(res[["study_meta"]]) ){
        missing <- c(missing, "study metadata ('s_' prefix)")
    }
    if( length(missing) > 0L ){
        msg <- paste0(
            "The experiment is missing the following file",
            ifelse(length(missing) > 1L, "s", ""), ": ",
            paste0(missing, collapse = ", ")
        )
        stop(msg, call. = FALSE)
    }

    # Match sample names
    res <- .match_sample_names(res)
    # Get tables from result list
    assay <- res[["assay"]]
    assay_meta <- res[["assay_meta"]]
    study_meta <- res[["study_meta"]]

    # Split assay to abundance table and feature metadata
    assay_cols <- colnames(assay) %in% rownames(assay_meta)
    feat_meta <- assay[ , !assay_cols, drop = FALSE]
    assay <- assay[ , assay_cols, drop = FALSE]

    # Assign feature IDs to assay
    feature_id <- .get_feature_ids(feat_meta)
    feat_meta[[feature_id]] <- feat_meta[[feature_id]] |> as.character()
    assay[[feature_id]] <- feat_meta[[feature_id]]

    # Add feature names to rownames of rowData and assay
    feat_names <- assay[[feature_id]]
    # Sometimes feature IDs are missing. Replace them with random number
    if( any(feat_names %in% c(NA, "", " ")) ){
        warning("Some features do not have IDs. Please check the data for ",
                "errors.", call. = FALSE)
        feat_names[ feat_names %in% c(NA, "", " ") ] <- "feature"
    }
    rownames(feat_meta) <- rownames(assay) <- feat_names |> make.unique()
    assay[[feature_id]] <- NULL

    # Combine study and assay metadata to be added to colData
    sample_meta <- .merge_metadata(study_meta, assay_meta)

    # Order metadatas based on assay
    feat_meta <- feat_meta[
        match(rownames(assay), rownames(feat_meta)), , drop = FALSE]
    sample_meta <- sample_meta[
        match(colnames(assay), rownames(sample_meta)), , drop = FALSE]

    # Abundance values might be in character format, and there might be trailing
    # spaces etc. Trim the values and convert to numeric.
    assay <- lapply(assay, function(x){
        if( is.character(x) ){
            x <- x |> trimws()
        }
        x <- x |> as.numeric() |> suppressWarnings()
        return(x)
        })
    assay <- do.call(cbind, assay)
    # Convert to classes supported by SE
    assay <- as.matrix(assay)
    assays <- SimpleList(assay)
    names(assays) <- assay.type
    feat_meta <- DataFrame(feat_meta, check.names = FALSE)
    sample_meta <- DataFrame(sample_meta, check.names = FALSE)
    # Create TreeSummarizedExperiment
    FUN <- if( output == "SE") SummarizedExperiment else
        TreeSummarizedExperiment
    se <- FUN(assays = assays, rowData = feat_meta, colData = sample_meta)
    return(se)
}

# This function matches sample names between tables. If they do not match,
# the function gives error.
.match_sample_names <- function(res){
    # Get all column names from assay
    assay_names <- res[["assay"]] |> colnames()
    # Get all possible sample names from assay metadata
    cols <- c(
        "Extract Name",
        "Derived Spectral Data File",
        "Derived Spectral Data File",
        "Sample Name"
    )
    meta_names <- res[["assay_meta"]][
        , colnames(res[["assay_meta"]]) %in% cols, drop = FALSE]
    # Polis sample names
    meta_names <- lapply(meta_names, function(col){
        col <- col |>
            as.character() |>
            basename() |>
            gsub(pattern = "\\.mzML$", replacement = "")
        return(col)
    })
    meta_names <- do.call(cbind.data.frame, meta_names)
    assay_names <- assay_names |>
        basename() |>
        gsub(pattern = "\\.mzML$", replacement = "")

    # Check which column from assay metadata has the highest number of mathches.
    num_matches <- vapply(meta_names, function(col){
        sum(col %in% assay_names)
    }, numeric(1L))
    # If sample anmes do not match, give error
    if( all(num_matches == 0L) ){
        stop("It seems that the files are missing columns containing ",
            "metabolite abundance data for individual samples. Each sample ",
            "should have its own column, labeled with the sample name. ",
            "Please check for errors.", call. = FALSE)
    }
    # The column with highest number of macthes includes the correct sample
    # names.
    sample_names <- meta_names[[which.max(num_matches)]]

    # Add sample names to assay and assay metadata
    res[["assay_meta"]] <- res[["assay_meta"]][
        !is.na(sample_names), , drop = FALSE]
    sample_names <- sample_names |> na.omit()
    matching <- match(assay_names, sample_names)
    cols_found <- assay_names %in% sample_names
    sample_names <- sample_names |> make.unique()
    #
    assay_names[ cols_found ] <- sample_names[ matching ] |> na.omit()
    colnames(res[["assay"]]) <- assay_names
    rownames(res[["assay_meta"]]) <- sample_names

    # Match sample names between study metadata and assay. To study metadata
    # assign names from assay.
    cols <- c("Sample Name", "sample_name", "sample id", "sample_id")
    col1 <- colnames(res[["study_meta"]]) %in% cols
    if( all(!col1) ){
        stop("Study metadata does not include sample names.", call. = FALSE)
    }
    col1 <- which(col1)[[1L]]
    col2 <- colnames(res[["assay_meta"]]) %in% cols
    if( all(!col2) ){
        stop("Assay metadata does not include sample names.", call. = FALSE)
    }
    col2 <- which(col2)[[1L]]
    res[["study_meta"]] <- res[["study_meta"]][
        match(res[["assay_meta"]][[col2]], res[["study_meta"]][[col1]]), ]
    rownames(res[["study_meta"]]) <- rownames(res[["assay_meta"]])

    return(res)
}

# This function returns a name of column from assay metadata that includes
# feature IDs
.get_feature_ids <- function(feat_meta){
    # Check which column has feature identifier
    cols <- c(
        "feat_ID",
        "metabolite_identification",
        "metabolite identification"
    )
    feature_id <- vapply(cols, function(name){
        any(grepl(name, colnames(feat_meta), ignore.case = TRUE))
    }, logical(1L))
    if( !any(feature_id) ){
        stop("No feature ID column found.", call. = FALSE)
    }
    feature_id <- names(feature_id[feature_id])[[1L]]
    return(feature_id)
}

# This function adds assay metadata to study metadata
.merge_metadata <- function(study_meta, assay_meta){
    # Store sample names
    sample_names <- study_meta |> rownames()
    # Combine assay and study metadata to metadata on samples
    common_cols <- intersect(colnames(study_meta), colnames(assay_meta))
    # Coerce first to same format
    class1 <- vapply(study_meta[common_cols], class, character(1L))
    class2 <- vapply(assay_meta[common_cols], class, character(1L))
    mod_cols <- common_cols[class1 != class2]
    study_meta[mod_cols] <- lapply(study_meta[mod_cols], as.character)
    assay_meta[mod_cols] <- lapply(assay_meta[mod_cols], as.character)
    # And then merge
    sample_meta <- left_join(study_meta, assay_meta, by = common_cols)
    # Add sample names back
    rownames(sample_meta) <- sample_names
    return(sample_meta)
}
