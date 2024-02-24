# Get data with getData. Format it so that we get SE. accession must be specified and it must be samples

getResult <- function(accession, ...){
    # If user tries to feed accession.type or type, disable them
    args <- list(...)
    args <- args[ !names(args) %in% c("type", "accession.type", "flatten")]
    # Get all the arguments
    args[["accession.type"]] <- "samples"
    args[["accession"]] <- accession
    args[["flatten"]] <- FALSE
    # Get results about samples
    sample_data <- do.call(getData, args)

    omics_list <- .construct_MAE(sample_data)
    mae <- omics_list[["mae"]]
    sample_metadata <- omics_list[["metadata"]]

    # Get animal metadata
    uniq_animals <- unique(sample_metadata[["animal"]])
    animal_data <- getData(accession.type = "animals", accession = uniq_animals)

    mae <- .add_metadata_to_MAE(mae, sample_metadata, animal_data)

    # # There are samples that are not presented in dataset
    not_found <- accession[ !accession %in% unlist(colnames(mae)) ]
    if( length(not_found) ){
        types <- sample_data[["sample_type"]]
        types <- types[types[["query_accession"]] %in% accession, "sample_type"]
        msg <- "Data for following samples cannot be found."
        if( "metagenomic_assembly" %in% types ){
            msg <- paste0(
                msg, " (Note that metagenomic assemblies can be ",
                "found from MGnify database. See MGnifyR package.)")
        }
        msg <- paste0(msg, ":\n'", paste(not_found, collapse = "', '"), "'")
        warning(msg, call. = FALSE)
    }
    return(mae)
}

.construct_MAE <- function(sample_data){
    # Get only structured metadata table. It includes info about measurements.
    # Other tables include info for instance about animals.
    sample_tab <- sample_data[["structured_metadata"]]
    colnames(sample_tab)[ colnames(
        sample_tab) == "query_accession" ] <- "accession"

    # Split results so that each metadata marker type gets own table
    f <- sample_tab[["marker.type"]]
    sample_tab <- split(sample_tab, f)

    # Some metadata markers are sample information that goes to sample metadata
    metadata_types <- c(
        "ENA Checklist", "TRIAL", "SAMPLE", "TREATMENT", "TANK"
    )
    metadata_types <-  names(sample_tab)[ names(sample_tab) %in% metadata_types ]
    # All other datatypes are stored as unique omics
    omic_types <- names(sample_tab)[ !names(sample_tab) %in% metadata_types ]
    # Divide the data to metadata and omics
    metadata_types <- sample_tab[ metadata_types ]
    omic_types <- sample_tab[ omic_types ]

    # Create metadata table from metadata types
    metadata <- .construct_metadata_from_markers(metadata_types)
    # Add rest of the sample data to metadata
    sample_data <- sample_data[ !names(sample_data) %in% "structured_metadata" ]
    metadata <- .add_sample_data_to_metadata(sample_data, metadata)
    # Create SE objecst from individual omics. Add metadata to SEs,
    # Add omics to MAE.
    omics <- .construct_omics_data_from_markers(omic_types, metadata)

    # Create a result list
    res <- list(mae = omics, metadata = metadata)
    return(res)
}

.add_metadata_to_MAE <- function(mae, metadata, animal_data){

    animal_tab <- animal_data[["structured_metadata"]]
    colnames(animal_tab)[ colnames(animal_tab) == "query_accession" ] <- "accession"
    f <- animal_tab[["marker.type"]]
    animal_tab <- split(animal_tab, f)

    animal_tab <- .construct_metadata_from_markers(animal_tab)
    animal_tab <- animal_tab[ match(metadata[["animal"]], animal_tab[["accession"]]), ]
    rownames(animal_tab) <- metadata[["accession"]]

    # Convert accesison to animal which format is used in sample data also
    colnames(animal_tab)[ colnames(animal_tab) == "accession" ] <- "animal"

    animal_tab <- .convert_cols_numeric(animal_tab)

    metadata <- DataFrame(animal_tab, check.names = FALSE)
    metadata <- metadata[ rownames(metadata) %in% unlist(colnames(mae)), ]

    colData(mae) <- metadata
    # Check if there are metagenomic samples. Give message about MGnifyR

    # Metagenomic cannot be included
    return(mae)
}

.add_sample_data_to_metadata <- function(sample_data, metadata){
    # Get non empty tables
    sample_data <- sample_data[ lengths(sample_data) > 0 ]
    # If there is analysis summaries table, add prefixes to its names so that
    # they do not match with the main table
    if( "analysis_summaries" %in% names(sample_data) ){
        # Get colnames
        cols <- colnames(sample_data[["analysis_summaries"]])
        # Add prefix
        cols[ !cols %in% c("query_accession")] <- paste0(
            "analysis_summaries.", cols[ !cols %in% c("query_accession")])
        # Add colnames back
        colnames(sample_data[["analysis_summaries"]]) <- cols
    }
    # Merge sample data
    sample_data <- Reduce(
        function(df1, df2) merge(df1, df2, all = TRUE), sample_data)
    # Sample data has same title column that metadata --> remove
    sample_data <- sample_data[ , !colnames(sample_data) %in% c("title"), drop = FALSE]
    # Add to metadata
    metadata <- merge(metadata, sample_data, by = "accession", all = TRUE)
    # Remove duplicated accessions.
    metadata <- metadata[ !duplicated(metadata[["accession"]]), ]
    rownames(metadata) <- metadata[["accession"]]
    return(metadata)
}

.construct_metadata_from_markers <- function(res){
    # If there are results
    if( length(res) > 0 ){
        # Convert each sample metadata type to table where each row represents
        # single sample
        res <- lapply(res, .convert_type_to_table)
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

.construct_omics_data_from_markers <- function(res, metadata){
    if( length(res) > 0 ){
        res <- lapply(res, .convert_type_to_SummarizedExperiment, metadata)
        res <- ExperimentList(res)
        res <- MultiAssayExperiment(res)
    } else{
        res <- NULL
    }
    return(res)
}

.convert_type_to_SummarizedExperiment <- function(
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

.convert_type_to_table <- function(type){
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

