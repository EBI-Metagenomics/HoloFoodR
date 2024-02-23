# Get data with getData. Format it so that we get SE. accession must be specified and it must be samples

getResult <- function(accession, ...){
    # If user tries to feed accession.type or type, grab it
    args <- list(...)
    args <- args[ !names(args) %in% c("type", "accession.type", "flatten")]
    args[["accession.type"]] <- "samples"
    args[["accession"]] <- accession
    args[["flatten"]] <- TRUE
    res <- do.call(getData, args)

    # Check if there are metagenomic samples. Give message about MGnifyR



    # Convert result to SE

    metadata_markers <- getData(type = "sample_metadata_markers")

    # Metagenomic cannot be included
    # The samples include metageomic that cannot eb included see MGnifyR package.
    return(res)
}


.convert_to_SummarizedExperiment <- function(type, assay.type = "counts", ...){
    # Check assay.type
    temp <- .check_input(assay.type, list("character scalar"))
    #
    # Specify columns that goes to assay and rowData. Other columns are sample
    # specific so they go to colData
    assay_info <- c("accession", "marker.name", "measurement")
    row_info <- c(
        "marker.name", "marker.type", "marker.canonical_url", "measurement",
        "units")
    col_info <- c(
        "accession",
        colnames(type)[!colnames(type) %in% c(assay_info, row_info)])

    # Get assay
    num_col <- "measurement"
    rownames_col <- "marker.name"
    accession_col <- "accession"
    assay <- type[, assay_info, drop = FALSE]
    # Convert values to numeric
    num_values <- gsub(",", ".", assay[[num_col]])
    num_values <- suppressWarnings( as.numeric(num_values) )
    assay[[num_col]] <- num_values
    # From long format to wide
    assay <- reshape(
        assay, direction = "wide", idvar = rownames_col,
        timevar = accession_col)
    # Add rownames and remove prefix from column names
    rownames(assay) <- assay[[rownames_col]]
    assay[[rownames_col]] <- NULL
    colnames(assay) <- gsub(paste0(num_col, "."), "", colnames(assay))
    # Convert to matrix
    assay <- as.matrix(assay)

    # Get rowData, remove duplicates, order it based on assay and convert to DF
    row_data <- type[ , row_info, drop = FALSE]
    row_data <- row_data[ !duplicated(row_data), ]
    rownames(row_data) <- row_data[[rownames_col]]
    row_data <- row_data[ rownames(assay), ]
    row_data <- DataFrame(row_data)
    # Same for colData
    col_data <- type[ , col_info, drop = FALSE]
    col_data <- col_data[ !duplicated(col_data), ]
    rownames(col_data) <- col_data[[accession_col]]
    col_data <- col_data[colnames(assay), ]
    col_data <- DataFrame(col_data)

    # Create SummarizedExperiment
    assays <- SimpleList(assay)
    names(assays) <- assay.type
    se <- SummarizedExperiment(
        assays = assays, rowData = row_data, colData = col_data)
    return(se)
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
    # Add to wide data
    long_data[["accession"]] <- rownames(long_data)
    table <- merge(table, long_data, by = "accession")

    # Remove duplicates
    table <- table[ !duplicated(table), ]
    # Check if duplicated accession IDs
    dupl_acc <- duplicated(table[["accession"]])
    if( any(dupl_acc) ){
        # Create a list from columns that have multiple values for certain rows.
        table <- .collapse_df(table)
    }


    # Check if there are analysis summaries tables, t
    # ind <- grepl("analysis_summaries", colnames(table))
    # if( any(ind) ){
    #     col_names <- colnames(table)[ind]
    #     # # Split columns based on accession
    #     # col <- split(table[, col_names, drop = FALSE], table[["accession"]])
    #     # # Remove column and the duplicates
    #     # table <- table[ , !colnames(table) %in% col_names, drop = FALSE]
    #     # table <- table[ !duplicated(table), ]
    #     # # Spread column, from long to wide
    #     # col <- .spread_column(col)
    #     # # Add back
    #     # table <- cbind(table, col)
    #
    #     ########################################################### Give warning if flattened --> cretae unflatten_col function
    #     # Loop through column names. Take every value related to certain
    #     # accession --> create a list and assign the list to column
    #     col_list <- lapply(col_names, function(x){
    #         split(table[[x]], table[["accession"]])
    #     })
    #     names(col_list) <- col_names
    #     # Remove column and the duplicates
    #     table <- table[ , !colnames(table) %in% col_names, drop = FALSE]
    #     table <- table[ !duplicated(table), ]
    #     #
    #     for( name in col_names ){
    #         table[[name]] <- col_list[[name]]
    #     }
    # }
    return(table)
}

.collapse_df <- function(table){
    # Check which columns are the problem
    dupl_col <- lapply(colnames(table), function(col){
        # Get column and accessons. Check if there are only unique values -->
        # then this column is the problematic column
        any(duplicated( table[, c("accession", col)] ))
    })
    dupl_col <- unlist(dupl_col)
    # Get duplicated columns
    dupl_col_names <- colnames(table)[ !dupl_col ]
    dupl_col <- table[, c("accession", dupl_col_names), drop = FALSE]
    # Remove duplicated columns from table and remove duplicates
    table <- table[, !colnames(table) %in% dupl_col_names, drop = FALSE]
    table <- table[ !duplicated(table), ]
    # Loop through problematic columns
    for( name in dupl_col_names ){
        # Split column to list
        col <- split(dupl_col[[name]], dupl_col[["accession"]])
        # Add it to back to origibal table
        table[[name]] <- col
    }
    return(table)
}

