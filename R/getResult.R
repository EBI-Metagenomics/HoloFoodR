# Get data with getData. Format it so that we get SE. accession must be specified and it must be samples

getResult <- function(accession, ...){
    # If user tries to feed accession.type or type, disable them
    args <- list(...)
    args <- args[ !names(args) %in% c("type", "accession.type", "flatten")]
    # Get all the arguments
    args[["accession.type"]] <- "samples"
    args[["accession"]] <- accession
    args[["flatten"]] <- FALSE
    # Get data on samples
    sample_data <- do.call(getData, args)
    # Replace accession with query accession to harmonize
    sample_data <- .query_accession_to_accession(sample_data)
    # Create a MultiAssayExperiment from the data
    omics_list <- .construct_MAE(sample_data)
    mae <- omics_list[["mae"]]
    sample_metadata <- omics_list[["metadata"]]

    # Get animal metadata for colData of MAE
    # Get all the animals present in sample metadata
    uniq_animals <- unique(sample_metadata[["animal"]])
    # Retrieve the data
    animal_data <- getData(accession.type = "animals", accession = uniq_animals)
    # Align animal data with sample data
    animal_data <- .align_animal_and_sample_data(animal_data, sample_metadata)
    # Replace accession with query accession to harmonize
    animal_data <- .query_accession_to_accession(animal_data)
    # Remove samples and sample_types table since we already have that info
    animal_data <- animal_data[ !names(
        animal_data) %in% c("sample_types", "samples")]
    # Construct MAE from animal data
    mae_animal_list <- .construct_MAE(animal_data)
    mae_animal <- mae_animal_list[["mae"]]
    animal_metadata <- mae_animal_list[["metadata"]]
    # Combine sample and animal data
    mae <- .add_animal_data_to_MAE(mae, mae_animal, animal_metadata)

    # If there are samples that user wanted to include but are not present in
    # the data (they do not have data in HoloFood database), give warning
    not_found <- accession[ !accession %in% unlist(colnames(mae)) ]
    if( length(not_found) ){
        # Create a message
        msg <- "Data for the following samples cannot be found."
        # Get the type of samples
        types <- sample_data[["sample_type"]]
        types <- types[types[["query_accession"]] %in% accession, "sample_type"]
        # If metagenomic assembly was one of the samples that user wanted, give
        # information that it can be found from MGnify database.
        if( "metagenomic_assembly" %in% types ){
            msg <- paste0(
                msg, " (Note that metagenomic assemblies can be ",
                "found from the MGnify database. See MGnifyR package.)")
        }
        # Add those sample IDs that cannot be found to message
        msg <- paste0(msg, ":\n'", paste(not_found, collapse = "', '"), "'")
        warning(msg, call. = FALSE)
    }
    return(mae)
}

################################ HELP FUNCTIONS ################################

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
.construct_MAE <- function(sample_data){
    # Get only structured metadata table. It includes info about measurements.
    # Other tables include info for instance about animals.
    sample_tab <- sample_data[["structured_metadata"]]
    # colnames(sample_tab)[ colnames(
    #     sample_tab) == "query_accession" ] <- "accession"

    # Split results so that each metadata marker type gets own table
    f <- sample_tab[["marker.type"]]
    sample_tab <- split(sample_tab, f)

    # Some metadata markers are sample information that goes to sample metadata
    metadata_types <- c(
        "ENA Checklist", "TRIAL", "SAMPLE", "TREATMENT", "TANK",
        "TOTAL FAT CONTENT", "PEN", "FAECES DIGESTIBILITY"
    )
    metadata_types <-  names(sample_tab)[ names(sample_tab) %in% metadata_types ]
    # All other datatypes are stored as unique omics
    omic_types <- names(sample_tab)[ !names(sample_tab) %in% metadata_types ]
    # Divide the data to metadata and omics
    metadata_types <- sample_tab[ metadata_types ]
    omic_types <- sample_tab[ omic_types ]

    # Create metadata table from metadata types
    metadata <- .construct_metadata_from_markers(metadata_types)
    # Convert to numeric those columns that can be converted
    metadata <- .convert_cols_numeric(metadata)
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

# We have already sample metadata from structured metadata. However, sample
# data also includes some additional tables. This function adds those additional
# tables to metadata.
.add_sample_data_to_metadata <- function(sample_data, metadata){
    # Get non empty tables
    sample_data <- sample_data[ lengths(sample_data) > 0 ]
    # Similarly to main table, query_accession is used as an ID. The column name
    # is changed because accession column can also point to animal.
    # query_accession point always to sample. However, if we would point just
    # to query_accession without modifying accession, it might be that we leave
    # sample_data <- lapply(sample_data, function(x){
    #     x[["accession"]] <- NULL
    #     colnames(x)[ colnames(x) == "query_accession" ] <- "accession"
    #     return(x)
    # })

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
    sample_data <- sample_data[ , !colnames(
        sample_data) %in% c("title"), drop = FALSE]
    # Add to metadata
    metadata <- merge(metadata, sample_data, by = "accession", all = TRUE)
    # Remove duplicated accessions.
    metadata <- metadata[ !duplicated(metadata[["accession"]]), ]
    rownames(metadata) <- metadata[["accession"]]
    return(metadata)
}

# This function creates a sample metadata from specific markers of structured
# metadata.
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
            # certain rows. (Otherwise it could not be added to colData since
            # only one row must point to single sample)
            res <- .collapse_df(res)
        }
        # Add accession IDs as rownames
        rownames(res) <- res[["accession"]]
    } else{
        res <- NULL
    }
    return(res)
}

# This function gets single datatype as an input. The data is in long format.
# It converts it to suitable format for sample metadata.
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
    # Convert to numeric those columns that can be converted
    table <- .convert_cols_numeric(table)
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
        # Add it to back to original table
        table[[name]] <- col
    }
    return(table)
}

# This function cretaes a MAE object from list of tables.
.construct_omics_data_from_markers <- function(res, metadata){
    # Convert each table to SummarizedExperiment
    res <- lapply(res, .convert_type_to_SummarizedExperiment, metadata)
    # Create MultiAssayExperiment
    res <- ExperimentList(res)
    res <- MultiAssayExperiment(res)
    return(res)
}

# This function gets a single table as an input and it converts the table
# to SE.
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

    # Create SummarizedExperiment
    assays <- SimpleList(assay)
    names(assays) <- assay.type
    se <- SummarizedExperiment(
        assays = assays, rowData = row_data,
        colData = col_data
        )
    return(se)
}

# This function gets data.frame as an input and it converts columns to numeric
# if they can be converted.
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

# This function modifies animal metadata and adds it to MAE's colData slot.
.add_animal_data_to_MAE <- function(mae, mae_animal, metadata){
    # IF animal MultiAssayExperiment has data, add it to main MAE
    experiments(mae) <- c( experiments(mae), experiments(mae_animal))
    # Create a DataFrame and drop those rows that are not present in MAE
    metadata <- DataFrame(metadata, check.names = FALSE)
    metadata <- metadata[ rownames(metadata) %in% unlist(colnames(mae)), ]
    # Add it to colData
    colData(mae) <- metadata
    return(mae)
}
