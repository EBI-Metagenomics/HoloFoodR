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
    expect_true( all(c("assay", "feat_meta", "sample_meta") %in% names(res)) )
    # The feature names and sample names should match
    expect_equal(res[["assay"]][["feat_ID"]], res[["feat_meta"]][["feat_ID"]])
    test_names <- colnames(res[["assay"]])
    ref_names <- c(res[["sample_meta"]][["Sample Name"]], "feat_ID")
    expect_true(
        all(test_names %in% ref_names) && all(ref_names %in% test_names))
})
