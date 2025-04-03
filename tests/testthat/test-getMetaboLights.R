context("getMetaboLigths")
test_that("getMetaboLigths", {
    # Expect errors when input is wrong
    expect_error( getMetaboLigths() )
    expect_error( getMetaboLigths(1) )
    expect_error( getMetaboLigths(TRUE) )
    expect_error( getMetaboLigths(c(1, TRUE)) )
    expect_error( getMetaboLigths(NULL) )
    
    # Require internet access
    skip_if_offline(host = "ebi.ac.uk")
    
    # Get the data
    study_id <- "MTBLS4381"
    res <- getMetaboLights(study_id)
    
    # There should be certain named data.frames
    expect_true( all(c("assay", "assay_meta", "study_meta") %in% names(res)) )
    # Expect that all sample names are shared between sampple meta and study
    # meta and assay
    expect_true(all(
        res[["assay_meta"]][["Sample Name"]] %in%
        res[["study_meta"]][["Sample Name"]]))
    expect_true(all(
        res[["assay_meta"]][["Sample Name"]] %in% colnames(res[["assay"]])))
    
    # Get data in SE format
    res2 <- getMetaboLights(study_id, output = "SE")
    expect_s4_class(res2, "SummarizedExperiment")
    expect_equal(rownames(res2), res[[1L]][["feat_ID"]])
    
    # Get data in TreeSE format
    res3 <- getMetaboLights(study_id, output = "TreeSE")
    expect_s4_class(res3, "TreeSummarizedExperiment")
    mat <- assay(res3)
    ref <- res[[1]][, colnames(mat)] |> as.matrix()
    expect_equal(mat, ref, check.attributes = FALSE)
})
