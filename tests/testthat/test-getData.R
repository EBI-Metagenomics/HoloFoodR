context("getData")
test_that("getData", {
    # Expect errors when input is wrong
    expect_error( getData() )
    expect_error( getData(type = "test") )
    expect_error( getData(c("samples", "animals")) )
    expect_error( getData("samples", accession = "SAMEA10104950") )
    expect_error( getData(
        "samples", accession = "SAMEA10104950", accession.type = "samples") )
    expect_error( getData(
        accession = "SAMEA10104950", accession.type = "samples",
        flatten = "test") )

    # Require internet access
    skip_if_offline(host = "holofooddata.org")

    # Get result and check that it is correct
    res <- getData(accession = "SAMEA10104950", accession.type = "samples")
    # Result must be a list
    expect_type(res, "list")
    expect_s3_class(res, NA)
    # Result must be now data.frame when flatten is used
    res <- getData(
        accession = "SAMEA10104950", accession.type = "samples", flatten = TRUE)
    expect_s3_class(res, "data.frame")
    # Check that result is the same when daa is fetched with doQuery
    ref <- doQuery("samples", max.hits = 100)
    res <- getData("samples", max.hits = 100)
    expect_true( all(res == ref | (is.na(res) & is.na(ref))) )
})
