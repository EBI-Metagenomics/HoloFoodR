context("getResult")
test_that("getResult", {
    # Expect errors when input is wrong
    expect_error( getResult() )
    expect_error( getResult(1) )
    expect_error( getResult(TRUE) )
    expect_error( getResult(c(TRUE, 1)) )
    expect_error( getResult(NULL) )

    # Require internet access
    skip_if_offline(host = "holofooddata.org")

    # Get some data
    samples <- c("SAMEA112950905", "SAMEA112952305", "SAMEA112949948")
    res <- getResult(samples)
    # The data must be MAE that have SEs inside
    expect_s4_class(res, "MultiAssayExperiment")
    expect_s4_class(res[[1]], "SummarizedExperiment")

    # This function fetches corresponding values from MAE (getResult) and list of
    # data.frames (getData)
    .get_values_to_test <- function(res, ref, num = 50){
        ref_temp <- ref[["structured_metadata"]]

        # Get all types that can be found from
        types <- intersect(ref_temp[["marker.type"]], names(res))
        # Loop through types
        result <- lapply(types, function(type){
            temp1 <- res[[type]]
            temp2 <- ref_temp[ ref_temp[["marker.type"]] == type, ]

            # Take samples from features
            samples <- unique( temp2[["query_accession"]] )
            features <- sample(
                unique(temp2[["marker.name"]]), num, replace = TRUE )
            test_set <- expand.grid(sample = samples, feature = features)
            # Loop through test set
            result <- apply(test_set, 1, function(row){
                sample <- row[["sample"]]
                feature <- row[["feature"]]
                # Get value from MAE
                val1 <- assay(temp1)[feature, sample]
                # Get value from reference table
                ind <- temp2[["marker.name"]] %in% feature & temp2[[
                    "query_accession"]] %in% sample
                val2 <- temp2[ ind, "measurement"]
                result <- c(mae = val1, ref = val2, slot = "assay")
                return(result)
            })
            result <- as.data.frame(result)
            result <- t(result)
            result <- as.data.frame(result)
            result[["type"]] <- type
            # Get also unit info from rowData
            val2 <- temp2[ temp2[["marker.name"]] %in% features, ]
            val1 <- rowData(temp1)[val2[["marker.name"]], "units"]
            val2 <- val2[["units"]]
            temp <- data.frame(
                mae = val1, ref = val2, slot = "rowData", type = type)
            result <- rbind(result, temp)
            # From colData, get value from ENA Checklist
            type <- "ENA Checklist"
            temp2 <- ref_temp[ ref_temp[["marker.type"]] == type, ]
            cols <- "project"
            val1 <- unique(colData(temp1)[[cols]])
            val2 <- unique(unlist(
                temp2[ temp2[["marker.name"]] %in% cols, "measurement"]))
            temp <- data.frame(
                mae = val1, ref = val2, slot = "colData", type = type)
            result <- rbind(result, temp)
            return(result)
        })
        result <- do.call(rbind, result)
        return(result)
    }

    # Get reference data with getData
    ref <- getData(accession.type = "samples", accession = samples)
    # Get corresponding values
    values <- .get_values_to_test(res, ref)
    # Some numeric values in original table have commas instead of points
    num_vals <- suppressWarnings( as.numeric(gsub(",", ".", values[["mae"]])) )
    num_vals_ref <- suppressWarnings( as.numeric(gsub(",", ".", values[["ref"]])) )
    char_vals <- values[is.na(num_vals), "mae"]
    char_vals_ref <- values[is.na(num_vals_ref), "ref"]
    # Check that values are correct
    expect_true(
        all(num_vals == num_vals_ref | (is.na(num_vals) & is.na(num_vals_ref))))
    expect_true(
        all(char_vals == char_vals_ref |
                (is.na(char_vals) & is.na(char_vals_ref))))
    
    # Check that data from MetaboLights is correct
    samples <- c("SAMEA112952704", "SAMEA112952705")
    res <- getResult(samples)
    # The data must be MAE that have SEs inside
    expect_s4_class(res, "MultiAssayExperiment")
    expect_s4_class(res[[1]], "SummarizedExperiment")
    
    ref <- getMetaboLights(
        "https://www.ebi.ac.uk/metabolights/ws/studies/MTBLS4381")
    # Get assay and test that the values are correct
    assay <- assay(res[["METABOLOMIC"]])
    assay_ref <- ref[["assay"]]
    rownames(assay_ref) <- assay_ref[["feat_ID"]]
    assay_ref[["feat_ID"]] <- NULL
    # The column names are HoloFood IDs in getResult
    col_names <- ref[["sample_meta"]][["Comment[BioSamples accession]"]][
        match(colnames(assay_ref), ref[["sample_meta"]][["Sample Name"]])]
    colnames(assay_ref) <- col_names
    assay_ref <- assay_ref[ , colnames(assay)]
    assay_ref <- as.matrix(assay_ref)
    expect_equal(assay, assay_ref)
    
    # Check sample meta
    coldata <- colData(res[["METABOLOMIC"]])
    coldata_ref <- ref[["sample_meta"]]
    # Esnure that the order is correct
    coldata_ref <- coldata_ref[
        match(
            rownames(coldata),
            coldata_ref[["Comment[BioSamples accession]"]]), ]
    rownames(coldata) <- NULL
    coldata <- data.frame(coldata, check.names = FALSE)
    expect_equal(coldata, coldata_ref)
    
    # Check feature meta
    rowdata <- rowData(res[["METABOLOMIC"]])
    rowdata_ref <- ref[["feat_meta"]]
    # Remove rownames and convert to data.frame
    rownames(rowdata) <- NULL
    rowdata <- data.frame(rowdata, check.names = FALSE)
    expect_equal(rowdata, rowdata_ref)
})
