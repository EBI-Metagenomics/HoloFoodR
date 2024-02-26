context("doQuery")
test_that("doQuery", {
    # Expect errors when input is wrong
    expect_error( doQuery() )
    expect_error( doQuery("test") )
    expect_error( doQuery(TRUE) )
    expect_error( doQuery(1) )
    expect_error( doQuery(c("samples", "animals")) )
    expect_error( doQuery("samples", max.hits = TRUE) )

    # Require internet access
    skip_if_offline(host = "holofooddata.org")

    # Get result and check that it is correct
    res <- doQuery("samples", max.hits = 100)
    # The result must be data.frame
    expect_s3_class(res, "data.frame")
    # It must include these columns
    column_names <- c(
        "accession", "title", "sample_type", "animal", "canonical_url",
        "metagenomics_url", "metabolomics_url")
    expect_true(all(column_names %in% colnames(res)))
    # The accession column must not include empty values
    expect_true(!any(is.na(res$accession)))
})
