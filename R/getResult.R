#' Get data on samples from HoloFood database
#'
#' @details
#' With \code{getResult}, you can fetch data on samples from the HoloFood
#' database. Compared to \code{getData}, this function is more convenient for
#' fetching the samples data because it converts the data to
#' \code{MultiAssayExperiment} where different omics are stored as
#' \code{TreeSummarizedExperiment} objects which are optimized for downstream
#' analytics. Columns of returned \code{MultiAssayExperiment} are individual
#' animals. These columns are linked with individual samples that are stored in
#' \code{TreeSummarizedExperiment} objects.
#' 
#' The HoloFood database lacks non-targeted metabolomic data but fetched from
#' MetaboLights resource. The function \code{getResult} facilitates the
#' automatic retrieval of metabolomic data and its integration with other
#' datasets from HoloFood.
#' 
#' Furthermore, while the HoloFoodR database does not include metagenomic
#' assembly data, users can access such data from the MGnify database. The
#' MGnifyR package provides a convenient interface for accessing this database.
#' By employing \code{MGnifyR::getResult()}, users can obtain data formatted as
#' a \code{MultiAssayExperiment} object, containing multiple
#' \code{TreeSummarizedExperiment} objects. Consequently, data from both
#' HoloFood and MGnify databases are inherently compatible for subsequent
#' downstream analysis.
#'
#' @param accession \code{Character vector} specifying the
#' accession IDs of type samples.
#' 
#' @param get.metabolomic \code{Logical scalar} specifying whether to retrieve
#' metabolomic data from MetaboLights database. (Default: \code{FALSE})
#'
#' @param ... optional arguments:
#' \itemize{
#'   
#'   \item \strong{use.cache} \code{Logical scalar} specifying whether to
#'   use cache. Note that when \code{get.metabolomic = TRUE} is specified, the
#'   file from the MetaboLights is stored in the local system to the location
#'   specified by \code{cache.dir} despite of the value of \code{use.cache}.
#'   (Default: \code{FALSE})
#'   
#'   \item \strong{cache.dir} \code{Character scalar} specifying cache
#'   directory. (Default: \code{tempdir()})
#'   
#'   \item \strong{clear.cache} \code{Logical scalar} specifying whether to
#'   use.cache (Default: \code{FALSE})
#'   
#'   \item \strong{assay.type} \code{Character scalar} specifying the name of
#'   assay in resulting \code{TreeSummarizedExperiment} object.
#'   (Default: \code{"counts"}) 
#'   
#' }
#'
#' @return \code{MultiAssayExperiment}
#'
#' @examples
#'
#' # Find samples on certain animal
#' samples <- doQuery("samples", animal_accession = "SAMEA112904746")
#'
#' # Get the data
#' mae <- getResult(samples[["accession"]])
#' mae
#'
#' @seealso
#' \code{\link[HoloFoodR:getData]{getData}}
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' \code{\link[MultiAssayExperiment:MultiAssayExperiment-class]{MultiAssayExperiment}}
#' \code{\link[MGnifyR:getResult]{MGnifyR:getResult}}
#'
#' @name getResult
NULL

#' @rdname getResult
#' @export
getResult <- function(accession, get.metabolomic = FALSE, ...){
    # Check accession
    temp <- .check_input(accession, list("character vector"))
    # Check get.metabolomic
    temp <- .check_input(get.metabolomic, list("logical scalar"))
    #
    # If user tries to feed accession.type or type, disable them
    args <- list(...)
    args <- args[ !names(args) %in% c("type", "accession.type", "flatten")]
    # Get all the arguments
    args[["accession.type"]] <- "samples"
    args[["accession"]] <- accession
    args[["flatten"]] <- FALSE
    # Get data on samples
    sample_data <- do.call(getData, args)
    # Check if no data found
    if( length(sample_data) == 0 ){
        warning("No data found.", call. = FALSE)
        return(NULL)
    }
    # Replace accession with query accession to harmonize
    sample_data <- .query_accession_to_accession(sample_data)
    # Create a MultiAssayExperiment from the data
    omics_list <- .construct_MAE(sample_data, ...)
    mae <- omics_list[["mae"]]
    sample_metadata <- omics_list[["sample_metadata"]]
    study_metadata <- omics_list[["study_metadata"]]
    
    # If user wants to get metabolites data and retrieved sample IDs include
    # metabolite samples. It requires loading files from MetaboLights which
    # is why there is an option for not loading the data.
    metabolomics_url <- sample_metadata[["metabolomics_url"]]
    metabolomics_url <- metabolomics_url[ !is.na(metabolomics_url) ]
    if( get.metabolomic && length(metabolomics_url) > 0 ){
        # Get metabolomic data
        se_metabolomic <- .construct_metabolomic_SE(metabolomics_url, ...)
        # Add it to MAE
        mae <- .add_metabolomic_data_to_MAE(mae, se_metabolomic, accession)
    }
    
    # If there are samples that user wanted to include but are not present in
    # the data (they do not have data in HoloFood database), give warning.
    # Moreover, add them as empty experiments so that animal metadata can still
    # be included for them also.
    not_found <- accession[ !accession %in% unlist(colnames(mae)) ]
    if( length(not_found) ){
        # Create a message
        msg <- "Data for the following samples cannot be found."
        # Get the type of samples
        types <- sample_data[["sample_type"]]
        types <- types[types[["accession"]] %in% not_found, ]
        type_names <- types[["sample_type"]]
        type_names <- unique(type_names)
        # Add sample types to message
        if( !is.null(type_names) ){
            # Add info about sample types
            msg_temp <- .create_msg_from_list(type_names, "and")
            if( length(types) == 1 ){
                msg_temp2 <- "The sample type is"
            } else{
                msg_temp2 <- "The sample types are"
            }
            msg <- paste0(msg, " ", msg_temp2, " ", msg_temp, ".")
            # If metagenomic assembly was one of the samples that user wanted,
            # give information that it can be found from MGnify database.
            if( "metagenomic_assembly" %in% type_names ){
                msg <- paste0(
                    msg, " (Note that metagenomic assemblies can be ",
                    "found from the MGnify database. See MGnifyR package.)")
            }
        }
        # Add those sample IDs that cannot be found to message
        msg <- paste0(msg, ":\n'", paste(not_found, collapse = "', '"), "'")
        warning(msg, call. = FALSE)
        
        # Create empty SE objects and add them to MAE
        mae <- .add_empty_experiments(mae, types, sample_metadata)
    }
    
    # Get animal metadata for colData of MAE
    # Get all the animals present in sample metadata
    uniq_animals <- unique(sample_metadata[["animal"]])
    # Retrieve the data
    animal_data <- getData(
        accession.type = "animals", accession = uniq_animals, ...)
    # Align animal data with sample data
    animal_data <- .align_animal_and_sample_data(animal_data, sample_metadata)
    # Replace accession with query accession to harmonize
    animal_data <- .query_accession_to_accession(animal_data)
    # Remove samples and sample_types table since we already have that info
    animal_data <- animal_data[ !names(
        animal_data) %in% c("sample_types", "samples")]
    # Construct MAE from animal data
    mae_animal_list <- .construct_MAE(animal_data, ...)
    mae_animal <- mae_animal_list[["mae"]]
    study_metadata2 <- mae_animal_list[["sample_metadata"]]
    study_metadata3 <- mae_animal_list[["study_metadata"]]
    # Combine sample and animal data
    mae <- .add_animal_data_to_MAE(mae, mae_animal)
    # Add study metadata to MAE
    metadata_list <- list(study_metadata, study_metadata2, study_metadata3)
    mae <- .add_study_metadata_to_MAE(mae, metadata_list)
    # Modify MAE so that top-level has animal IDs which points to samples in
    # individual SEs
    mae <- .MAE_cols_to_animals(mae)
    return(mae)
}

################################ HELP FUNCTIONS ################################

# If accession cannot be found, animal metadata is not included for that
# accession in MAE. Since user wanted to get the data for that also, add empty
# SEs to MAE with accessions, so that animal metadata can be included for those
# accessions in colData(mae).
.add_empty_experiments <- function(mae, types, metadata){
    # Loop through types and create SE from them
    type_names <- unique(types[["sample_type"]])
    ses <- lapply(type_names, function(type){
        # Get all the accessions for that sample type
        accessions <- types[types[["sample_type"]] == type, "accession"]
        # Get metadata on those accessions
        temp <- metadata[match(accessions, metadata[["accession"]], ), ]
        rownames(temp) <- accessions
        temp <- DataFrame(temp, check.names = FALSE)
        # Create empty SE
        se <- TreeSummarizedExperiment(colData = temp)
        return(se)
    })
    names(ses) <- type_names
    # Convert to ExperimentList
    ses <- ExperimentList(ses)
    # Add to MAE
    exp_list <- c(experiments(mae), ses)
    mae <- MultiAssayExperiment(exp_list)
    return(mae)
}

# This function subsets metabolomic SE based on accession and adds the data
# to MAE.
.add_metabolomic_data_to_MAE <- function(mae, se, accession){
    # Since MetaboLights data is fetched from the files including all samples
    # from study and subsetting was not done, there might be samples that was
    # not listed in accession. --> subset
    ind <- colnames(se) %in% accession
    se <- se[, ind]
    # Create  experiment list from  metabolomic data
    experiment_list <- ExperimentList(METABOLOMIC = se)
    # Add it to existing experiment list and create MAE
    experiment_list <- c(experiments(mae), experiment_list)
    res <- MultiAssayExperiment(experiment_list)
    return(res)
}

# This function retrieves metabolomic data and constructs SE from it
.construct_metabolomic_SE <- function(urls, assay.type = "counts", ...){
    # Check assay.type
    temp <- .check_input(assay.type, list("character scalar"))
    #
    # Get unique urls
    urls <- unique(urls)
    # Get data from MetaboLigths
    res <- getMetaboLights(urls, ...)
    assay <- res[["assay"]]
    assay_meta <- res[["assay_meta"]]
    study_meta <- res[["study_meta"]]
    
    # Split assay to abundance table and feature metadata
    assay_cols <- colnames(assay) %in% assay_meta[["Sample Name"]]
    feat_meta <- assay[ , !assay_cols, drop = FALSE]
    feat_meta[["feat_ID"]] <- as.character(feat_meta[["feat_ID"]])
    assay <- assay[ , assay_cols, drop = FALSE]
    assay[["feat_ID"]] <- feat_meta[["feat_ID"]]

    # Combine assay and study metadata to metadata on samples
    common_cols <- intersect(colnames(study_meta), colnames(assay_meta))
    sample_meta <- left_join(study_meta, assay_meta, by = common_cols)

    # Add rownames to tables
    rownames(assay) <- assay[["feat_ID"]]
    assay[["feat_ID"]] <- NULL
    rownames(feat_meta) <- feat_meta[["Sample Name"]]
    
    # Order metadatas based on assay
    feat_meta <- feat_meta[match(rownames(assay), feat_meta[["feat_ID"]]), ]
    sample_meta <- sample_meta[
        match(colnames(assay), sample_meta[["Sample Name"]]), ]
    
    # Convert to classes supported by SE
    assay <- as.matrix(assay)
    assays <- SimpleList(assay)
    names(assays) <- assay.type
    feat_meta <- DataFrame(feat_meta, check.names = FALSE)
    sample_meta <- DataFrame(sample_meta, check.names = FALSE)
    # Create TreeSummarizedExperiment
    se <- TreeSummarizedExperiment(
        assays = assays, rowData = feat_meta, colData = sample_meta)
    
    # The data still has sample names from MetaboLights and not from HoloFood
    # replace them.
    colnames(se) <- colData(se)[["Comment[BioSamples accession]"]]
    return(se)
}

# This function replaces query_accession with accession --> this is to ensure
# that we merge the data based on query_accession which always point to samples.
# Accession might point to animals also.
.query_accession_to_accession <- function(sample_data){
    sample_data <- lapply(sample_data, function(x){
        x[["accession"]] <- NULL
        colnames(x)[ colnames(x) == "query_accession" ] <- "accession"
        return(x)
    })
    return(sample_data)
}

# This function constructs MAE object from the sample data.
.construct_MAE <- function(sample_data, ...){
    # Get only structured metadata table. It includes info about measurements.
    # Other tables include info for instance about animals.
    sample_tab <- sample_data[["structured_metadata"]]
    # Split results so that each metadata marker type gets own table
    f <- sample_tab[["marker.type"]]
    sample_tab <- split(sample_tab, f)

    # Some metadata markers are sample information that goes to sample metadata.
    # Some metadata is common for certain animal. That information goes to
    # common "study" metadata --> colData(mae).
    sample_metadata_types <- c("ENA Checklist", "SAMPLE")
    study_metadata_types <- c(
        "TRIAL", "TREATMENT", "TANK",
        "TOTAL FAT CONTENT", "PEN", "FAECES DIGESTIBILITY"
    )
    sample_metadata_types <-  names(sample_tab)[
        names(sample_tab) %in% sample_metadata_types]
    study_metadata_types <-  names(sample_tab)[
        names(sample_tab) %in% study_metadata_types]
    # All other datatypes are stored as unique omics
    omic_types <- names(sample_tab)[
        !names(sample_tab) %in% c(sample_metadata_types, study_metadata_types)]
    # Divide the data to metadata and omics
    sample_metadata_types <- sample_tab[ sample_metadata_types ]
    study_metadata_types <- sample_tab[ study_metadata_types ]
    omic_types <- sample_tab[ omic_types ]

    # Create metadata table from metadata types
    sample_metadata <- .construct_metadata_from_markers(
        sample_metadata_types, ...)
    study_metadata <- .construct_metadata_from_markers(
        study_metadata_types, ...)
    # Add rest of the sample data to metadata
    sample_data <- sample_data[ !names(sample_data) %in% "structured_metadata" ]
    sample_metadata <- .add_sample_data_to_metadata(
        sample_data, sample_metadata)
    # Create SE objecst from individual omics. Add metadata to SEs,
    # Add omics to MAE.
    omics <- .construct_omics_data_from_markers(
        omic_types, sample_metadata, ...)

    # Create a result list
    res <- list(
        mae = omics,
        sample_metadata = sample_metadata,
        study_metadata = study_metadata)
    return(res)
}

# We have already sample metadata from structured metadata. However, sample
# data also includes some additional tables. This function adds those additional
# tables to metadata.
#' @importFrom dplyr full_join left_join
.add_sample_data_to_metadata <- function(sample_data, metadata){
    # Get non empty tables
    sample_data <- sample_data[ lengths(sample_data) > 0 ]
    # If there is analysis summaries table, add prefixes to its names so that
    # they do not match with the main table
    if( "analysis_summaries" %in% names(sample_data) ){
        # Get colnames
        cols <- colnames(sample_data[["analysis_summaries"]])
        # Add prefix
        cols[ !cols %in% c("accession")] <- paste0(
            "analysis_summaries.", cols[ !cols %in% c("accession")])
        # Add colnames back
        colnames(sample_data[["analysis_summaries"]]) <- cols
    }
    # Merge sample data
    sample_data <- .full_join_list(sample_data)
    # Sample data has same title column that metadata --> remove
    sample_data <- sample_data[ , !colnames(
        sample_data) %in% c("title"), drop = FALSE]
    # Add metadata to sample data --> sample data includes all accessions
    if( !is.null(metadata) ){
        sample_data <- full_join(sample_data, metadata, by = "accession")
    }
    return(sample_data)
}

# This function creates a sample metadata from specific markers of structured
# metadata.
#' @importFrom dplyr bind_rows
.construct_metadata_from_markers <- function(res, ...){
    # If there are results
    if( length(res) > 0 ){
        # Convert each sample metadata type to table where each row represents
        # single sample
        res <- lapply(res, function(x) .convert_type_to_table(x, ...))
        # Combine data
        res <- bind_rows(res)
        # Convert to numeric those columns that can be converted
        res <- .convert_cols_numeric(res, ...)
    } else{
        res <- NULL
    }
    return(res)
}

# This function gets single datatype as an input. The data is in long format.
# It converts it to suitable format for sample metadata.
#' @importFrom stats reshape
#' @importFrom dplyr full_join
.convert_type_to_table <- function(type, ...){
    # Get those columns that are in long format and wide
    long_info <- c("accession", "measurement", "marker.name")
    wide_info <- c("accession", colnames(type)[!colnames(type) %in% long_info])
    long_data <- type[, long_info, drop = FALSE]
    table <- type[, wide_info, drop = FALSE]
    # Drop duplicates
    long_data <- long_data[ !duplicated(long_data), ]
    table <- table[ !duplicated(table), ]
    # Convert wide format to long
    num_col <- "measurement"
    rownames_col <- "marker.name"
    accession_col <- "accession"
    long_data <- reshape(
        long_data, direction = "wide", idvar = rownames_col,
        timevar = accession_col)
    # Add rownames and modify colnames by deleting prefix
    rownames(long_data) <- long_data[[rownames_col]]
    long_data[[rownames_col]] <- NULL
    colnames(long_data) <- gsub(paste0(num_col, "."), "", colnames(long_data))
    # Transpose to correct orientation
    long_data <- t(long_data)
    long_data <- as.data.frame(long_data)
    long_data[["accession"]] <- rownames(long_data)
    # Add only those columns that are not shared
    cols <- colnames(long_data)[ !colnames(long_data) %in% colnames(table) ]
    cols <- c(cols, "accession")
    long_data <- long_data[ , cols, drop = FALSE]
    # Add to wide data
    table <- full_join(table, long_data, by = "accession")

    # If there is units info, add it to column since now variables are as
    # columns so units are not defined correctly
    marker_col <- "marker.name"
    unit_col <- "units"
    if( all(c(marker_col, unit_col) %in% colnames(type)) ){
        # Remove units from table
        table[[unit_col]] <- NULL
        # Get unique units, drop NA values
        units <- type[ , c(marker_col, unit_col), drop = FALSE]
        units <- units[ !duplicated(units[[marker_col]]), ]
        units <- units[ !is.na(units[[unit_col]]), ]
        # If there are units info for variables
        if( nrow(units) > 0 ){
            # Modify marker name and add it to rownames
            rownames(units) <- paste0(units[[marker_col]], ", unit")
            # Create data.frame where each column represent measure name and its
            # unit
            units[[marker_col]] <- NULL
            units <- t(units)
            rownames(units) <- NULL
            units <- as.data.frame(units)
            # Replicate row so that there are as many rows as in original table
            units <- units[rep(1, nrow(table)), , drop = FALSE]
            # Add to original table
            table <- cbind(table, units)
        }
    }
    # Remove duplicates
    table <- table[ !duplicated(table), ]
    # Convert to numeric those columns that can be converted
    table <- .convert_cols_numeric(table, ...)
    return(table)
}

# If sample metadata has multiple rows for certain accession, that means that
# it has multiple values for some variables. Combine these values so that the
# variable column is converted into list.
.collapse_df <- function(table){
    # Check which columns are the problem
    dupl_col <- lapply(colnames(table), function(col){
        # Get column and accessions. Remove duplicates and check if the number
        # of unique values do not match with number of unique accessions -->
        # duplicated column
        N_accession <- length(unique( table[["accession"]] ))
        N_values <- nrow(unique( table[, c("accession", col), drop = FALSE] ))
        temp <- N_values != N_accession
        return(temp)
    })
    dupl_col <- unlist(dupl_col)
    # Get duplicated columns
    dupl_col_names <- colnames(table)[ dupl_col ]
    dupl_col <- table[, c("accession", dupl_col_names), drop = FALSE]
    # Remove duplicated columns from table and remove duplicates
    table <- table[, !colnames(table) %in% dupl_col_names, drop = FALSE]
    table <- table[ !duplicated(table), ]
    # Loop through problematic columns
    for( name in dupl_col_names ){
        # Get column
        col <- dupl_col[, c("accession", name), drop = FALSE]
        # Split column to list
        col <- split(col[[name]], col[["accession"]])
        # Loop through accessions/rows, there might be that values are NA or
        # replicated. Take only non-replicated values (or if any values are not
        # found, NA)
        col <- lapply(col, function(values){
            values <- unique(values)
            values <- values[ !is.na(values) ]
            if(length(values) == 0){
                values <- NA
            }
            return(values)
        })
        # Convert to vector if there is only one value for each accession
        if( all(lengths(col) == 1) ){
            col <- unlist(col)
        }
        # Order data based on original table where the data will be added.
        # The order can be changed during the process.
        ind <- match(table[["accession"]], names(col))
        col <- col[ ind ]
        # Add it to back to original table
        table[[name]] <- col
    }
    return(table)
}

# This function cretaes a MAE object from list of tables.
.construct_omics_data_from_markers <- function(res, metadata, ...){
    # Convert each table to TreeSummarizedExperiment
    res <- lapply(res, function(x){
        .convert_type_to_TreeSummarizedExperiment(x, metadata, ...)
    })
    # Create MultiAssayExperiment
    res <- ExperimentList(res)
    res <- MultiAssayExperiment(res)
    return(res)
}

# This function gets a single table as an input and it converts the table
# to SE.
#' @importFrom stats reshape
#' @importFrom S4Vectors SimpleList
.convert_type_to_TreeSummarizedExperiment <- function(
        type, metadata, assay.type = "counts", ...){
    # Check assay.type
    temp <- .check_input(assay.type, list("character scalar"))
    #
    # Specify columns that goes to assay and rowData. Other columns are sample
    # specific so they go to colData
    assay_info <- c("accession", "marker.name", "measurement")
    row_info <- c(
        "marker.name", "marker.type", "marker.canonical_url", "units")
    col_info <- c(
        "accession",
        colnames(type)[!colnames(type) %in% c(assay_info, row_info)])

    # Get assay
    num_col <- "measurement"
    rownames_col <- "marker.name"
    accession_col <- "accession"
    assay <- type[, assay_info, drop = FALSE]
    # From long format to wide
    assay <- reshape(
        assay, direction = "wide", idvar = accession_col,
        timevar = rownames_col)
    # Add rownames and remove prefix from column names
    rownames(assay) <- assay[[accession_col]]
    assay[[accession_col]] <- NULL
    colnames(assay) <- gsub(paste0(num_col, "."), "", colnames(assay))

    # Try to convert columns to numeric.
    assay <- .convert_cols_numeric(assay, ...)
    # Samples are in columns and features in rows
    assay <- t(assay)
    assay <- as.matrix(assay)

    # Get rowData, remove duplicates, order it based on assay and convert to DF
    row_data <- type[ , row_info, drop = FALSE]
    row_data <- row_data[ !duplicated(row_data), ]
    rownames(row_data) <- row_data[[rownames_col]]
    row_data <- row_data[ rownames(assay), ]
    row_data <- DataFrame(row_data)
    
    # Check if duplicated accession IDs
    dupl_acc <- duplicated(metadata[["accession"]])
    if( any(dupl_acc) ){
        # Create a list from columns that have multiple values for
        # certain rows. (Otherwise it could not be added to colData since
        # only one row must point to single sample)
        metadata <- .collapse_df(metadata)
    }
    rownames(metadata) <- metadata[["accession"]]
    # Get metadata fpr certain accessions
    col_data <- metadata[colnames(assay), ]
    # Remove those columns that do not have indfo
    empty <- unlist( lapply(col_data, function(x) all(is.na(x))) )
    empty <- names(empty)[ empty ]
    # Remove also unit columns if there are
    empty <- c(empty, paste0(empty, ", unit"))
    col_data <- col_data[ , !colnames(col_data) %in% empty, drop = FALSE]
    # Convert to DF
    col_data <- DataFrame(col_data, check.names = FALSE)

    # Create TreeSummarizedExperiment
    assays <- SimpleList(assay)
    names(assays) <- assay.type
    se <- TreeSummarizedExperiment(
        assays = assays, rowData = row_data,
        colData = col_data
        )
    return(se)
}

# This function gets data.frame as an input and it converts columns to numeric
# if they can be converted.
.convert_cols_numeric <- function(df, replace.lower.th = FALSE, ...){
    # Check replace.lower.th
    temp <- .check_input(replace.lower.th, list("logical scalar"))
    #
    df <- as.data.frame(df, check.names = FALSE)
    # Skip those columns that are list
    is_list <- unlist(lapply(df, is.list))
    list_cols <- df[ , is_list, drop = FALSE]
    vec_cols <- df[ , !is_list, drop = FALSE]
    # Loop through columns. Try to convert values to numeric.
    col_names <- rownames(df)
    vec_cols <- lapply(vec_cols, function(x){
        # Some numeric values have comma instead of point
        temp <- gsub(",", ".", x)
        # Some values have percentage symbol. The unit already tells that it is
        # percentage.
        temp <- gsub("%", "", temp)
        # If user has specified, replace those values that are under detection
        # threshold to 0.
        msg <- NULL
        if( replace.lower.th ){
            # Get those values that matches to give warning message
            pattern <- "<\\s*\\d+"
            matched_values <- grep(pattern, temp, value = TRUE)
            matched_values <- unique(matched_values)
            # Replace values with 0.
            if( length(matched_values) > 0 ){
                # Replace
                temp <- gsub(pattern, 0, temp)
                # Create message text
                msg <- paste0(
                    "The following values are replaced with 0: '",
                    paste0(matched_values, collapse = "', '"), "'")
            }
            
        }
        # Try to convert to numeric
        temp_num <- suppressWarnings( as.numeric(temp) )
        # Check if we lost info. If we did, then the column is not numeric. If
        # there are as many NAs in same places as before, conversion was
        # succesful and the values are numeric.
        if( all(is.na(temp_num) == is.na(x)) ){
            x <- temp_num
            # Give warning message about replacing lower threshold
            if( !is.null(msg) ){
                warning(msg, call. = FALSE)
            }
        }
        return(x)
    })
    # Convert list to data.frame --> check names so that feature names are
    # not modified
    vec_cols <- as.data.frame(vec_cols, check.names = FALSE)
    rownames(vec_cols) <- col_names
    df <- cbind(vec_cols, list_cols)
    return(df)
}

# This function orders the animal data based on samples. In input, rows
# represent animals. In putput, rows represent samples.
.align_animal_and_sample_data <- function(animal_data, sample_metadata){
    # Get all the datatypes
    datatypes <- names(animal_data)
    # Loop through datatypes
    animal_data <- lapply(datatypes, function(type){
        # Get the table (specific data type)
        tab <- animal_data[[type]]
        # If the table is not samples table, convert accessions pointing to
        # animal. (In samples table accessions points to samples)
        if( type != "samples" ){
            tab[["animal"]] <- tab[["accession"]]
        }
        # Order the data based on sample metadata table
        tab <- tab[ match(
            sample_metadata[["animal"]], tab[["query_accession"]]),  ]
        # Add "query_accession" column that the sample that the row represents
        tab[["query_accession"]] <- sample_metadata[["accession"]]
        return(tab)
    })
    names(animal_data) <- datatypes
    return(animal_data)
}

# This function adds animal data to MAE
.add_animal_data_to_MAE <- function(mae, mae_animal){
    # IF animal MultiAssayExperiment has data, add it to main MAE
    exp_list <- c( experiments(mae), experiments(mae_animal))
    mae <- MultiAssayExperiment(exp_list)
    return(mae)
}

# This functions combines metadata from samples and animals and add it to
# colData slot of MAE
.add_study_metadata_to_MAE <- function(mae, metadata_list){
    # Combine metadata
    metadata <- .full_join_list(metadata_list)
    # There might be multiple values for each accession. This means that there
    # are multiple values for each data type. Collapse values so that values
    # becomes to list
    dupl_acc <- duplicated(metadata[["accession"]])
    if( any(dupl_acc) ){
        # Create a list from columns that have multiple values for
        # certain rows. (Otherwise it could not be added to colData since
        # only one row must point to single sample)
        metadata <- .collapse_df(metadata)
    }
    rownames(metadata) <- metadata[["accession"]]
    # Create a DataFrame and drop those rows that are not present in MAE
    metadata <- DataFrame(metadata, check.names = FALSE)
    metadata <- metadata[ rownames(metadata) %in% unlist(colnames(mae)), ]
    # Add it to colData
    colData(mae) <- metadata
    return(mae)
}

# This function modifies MAE so that each column in MAE is an animal. These
# columns are pointing to samples in individual experiments.
.MAE_cols_to_animals <- function(mae){
    # Get experiments
    exp_list <- experiments(mae)
    # Get sample mapping
    sample_map <- sampleMap(mae)
    # Get animal metadata
    col_data <- colData(mae)
    # Replace MAE-level key with animal ID
    animal_ids <- col_data[ , c("accession", "animal")]
    ind <- match(sample_map[["primary"]], animal_ids[["accession"]])
    animal_ids <- animal_ids[ind, ]
    sample_map[["primary"]] <- animal_ids[["animal"]]
    # Remove samples from animal metadata
    col_data[["accession"]] <- NULL
    # Rename rows with animal IDs
    rownames(col_data) <- col_data[["animal"]]
    # Remove duplicates so that there is only one row per animal
    col_data <- unique(col_data)
    
    # Finally, create MAE
    mae <- MultiAssayExperiment(
        experiments = exp_list,
        colData = col_data,
        sampleMap = sample_map
    )
    return(mae)
}
