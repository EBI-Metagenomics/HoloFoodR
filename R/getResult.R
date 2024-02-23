# Get data with getData. Format it so that we get SE. accession must be specified and it must be samples

getResult <- function(accession, ...){
    # If user tries to feed accession.type or type, grab it
    args <- list(...)
    args <- args[ !names(args) %in% c("type", "accession.type", "flatten")]
    # Get all the arguments
    args[["accession.type"]] <- "samples"
    args[["accession"]] <- accession
    args[["flatten"]] <- FALSE
    # Get results
    res <- do.call(getData, args)
    whole_data <- res
    res <- res[["structured_metadata"]]
    colnames(res)[ colnames(res) == "query_accession" ] <- "accession"

    # Split results so that each metadata marker type gets own table
    f <- res[["marker.type"]]
    res <- split(res, f)

    # Parse results. Get those columns that go to metadata
    metadata_types <- c(
        "ENA Checklist", "TRIAL", "SAMPLE", "TREATMENT", "TANK"
        #"TOTAL FAT CONTENT", "IODINE"
        )
    metadata_types <-  names(res)[ names(res) %in% metadata_types ]
    # All other datatypes go to unique omics
    omic_types <- names(res)[ !names(res) %in% metadata_types ]
    # Divide the data to metadata and omics
    metadata_types <- res[ metadata_types ]
    omic_types <- res[ omic_types ]

    # Instead of samples, match with animal ID. Sample IDs do not match between omics
    # there might be multiple animal IDs. Matcj between animal ID and timepoint.
    metadata <- .construct_metadata(metadata_types)
    omics <- .construct_omics_data(omic_types, metadata)

    # Get animal metadata
    animal_metadata <- getData(accession.type = "animals", accession = unique(whole_data$animal$animal))
    animal_metadata <- animal_metadata$structured_metadata
    colnames(animal_metadata)[ colnames(animal_metadata) == "query_accession" ] <- "accession"
    f <- animal_metadata[["marker.type"]]
    animal_metadata <- split(animal_metadata, f)

    animal_metadata <- .construct_metadata(animal_metadata)
    animal_metadata <- animal_metadata[ match(whole_data$animal$animal, animal_metadata$accession), ]
    rownames(animal_metadata) <- whole_data$animal$query_accession

    animal_metadata <- .convert_cols_numeric(animal_metadata)

    # There are samples that are not presented in dataset
    not_found <- accession[ !accession %in% unlist(colnames(omics)) ]
    not_found <- metadata[ metadata$accession %in% not_found, ]


    if( !is.null(omics) ){
        res <- omics
        if(!is.null(metadata)){
            metadata <- DataFrame(animal_metadata, check.names = FALSE)
            colData(res) <- metadata
        }
    } else{
        res <- metadata
    }
    # Check if there are metagenomic samples. Give message about MGnifyR



    # Convert result to SE

    metadata_markers <- getData(type = "sample_metadata_markers")

    # Metagenomic cannot be included
    # The samples include metageomic that cannot eb included see MGnifyR package.
    return(res)
}

.construct_metadata <- function(res){
    # If there are results
    if( length(res) > 0 ){
        # Convert each sample metadata type to table where each row represents
        # single sample
        res <- lapply(res, .convert_to_table)
        # Combine data
        res <- Reduce(function(df1, df2) merge(df1, df2, all = TRUE), res)
        # Check if duplicated accession IDs
        dupl_acc <- duplicated(res[["accession"]])
        if( any(dupl_acc) ){
            # Create a list from columns that have multiple values for
            # certain rows.
            res <- .collapse_df(res)
        }
        # Add accession IDs as rownames
        rownames(res) <- res[["accession"]]
    } else{
        res <- NULL
    }
    return(res)
}

.construct_omics_data <- function(res, metadata){
    if( length(res) > 0 ){
        res <- lapply(res, .convert_to_SummarizedExperiment, metadata)
        res <- ExperimentList(res)
        res <- MultiAssayExperiment(res)
    } else{
        res <- NULL
    }
    return(res)
}

.convert_to_SummarizedExperiment <- function(
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
    assay <- .convert_cols_numeric(assay)
    # Samples are in columns and features in rows
    assay <- t(assay)
    assay <- as.matrix(assay)
    # Get rowData, remove duplicates, order it based on assay and convert to DF
    row_data <- type[ , row_info, drop = FALSE]
    row_data <- row_data[ !duplicated(row_data), ]
    rownames(row_data) <- row_data[[rownames_col]]
    row_data <- row_data[ rownames(assay), ]
    row_data <- DataFrame(row_data)
    # # Same for colData
    # col_data <- type[ , col_info, drop = FALSE]
    # col_data <- col_data[ !duplicated(col_data), ]
    # rownames(col_data) <- col_data[[accession_col]]
    # col_data <- col_data[colnames(assay), ]
    # # Add common metadata. Add those columns that are not yet present
    # common_cols <- colnames(metadata)[ !colnames(
    #     metadata) %in% colnames(col_data) ]
    # common_cols <- c("accession", common_cols)
    # metadata <- metadata[ , common_cols, drop = FALSE]
    # col_data <- merge(
    #     col_data, metadata, by = "accession", all.x = TRUE,
    #     suffixes = c("", ".y"))
    # # Remove those columns that do not have data
    # empty <- unlist( lapply(col_data, function(x) all(is.na(x))) )
    # empty <- names(empty)[ empty ]
    # # Remove also unit columns if there are
    # empty <- c(empty, paste0(empty, ", unit"))
    # col_data <- col_data[ , !colnames(col_data) %in% empty, drop = FALSE]
    # col_data <- DataFrame(col_data)

    col_data <- metadata[colnames(assay), ]
    empty <- unlist( lapply(col_data, function(x) all(is.na(x))) )
    empty <- names(empty)[ empty ]
    # Remove also unit columns if there are
    empty <- c(empty, paste0(empty, ", unit"))
    col_data <- col_data[ , !colnames(col_data) %in% empty, drop = FALSE]
    col_data <- DataFrame(col_data)

    # Rename samples so that they get animal accession and time point info. ################################################
    # Sample IDs are unique for omics; with animal IDs, we can match samples
    # between omics.

    # Create SummarizedExperiment
    assays <- SimpleList(assay)
    names(assays) <- assay.type
    se <- SummarizedExperiment(
        assays = assays, rowData = row_data,
        colData = col_data
        )
    return(se)
}

.convert_cols_numeric <- function(df){
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
        # Try to convert to numeric
        temp_num <-- suppressWarnings( as.numeric(temp) )
        # Check if we lost info. If we did, then the column is not numeric. If
        # there are as many NAs in same places as before, conversion was
        # succesful and the values are numeric.
        if( all(is.na(temp_num) == is.na(x)) ){
            x <- temp_num
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

.convert_to_table <- function(type){
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
    table <- merge(table, long_data, by = "accession", suffixes = c("",".y"))

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
            units <- units[rep(1, nrow(table)), ]
            # Add to original table
            table <- cbind(table, units)
        }
    }
    # Remove duplicates
    table <- table[ !duplicated(table), ]
    return(table)
}

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
        # Add it to back to original table
        table[[name]] <- col
    }
    return(table)
}

