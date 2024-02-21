# Query animals, samples, genome-catalogues, or viral-catalogues

doQuery <- function(
        type, max.hits = NULL, flatten = TRUE, ...){
    #
    # Check type
    supported_types <- c(
        "animals", "genome-catalogues", "samples", "viral-catalogues")
    temp <- .check_input(
        type, list("NULL", "character scalar"),
        supported_values = supported_types)
    # Check flatten
    temp <- .check_input(flatten, list("logical vector"))
    # Check max.hits
    temp <- .check_input(
        max.hits, list("NULL", "integer scalar"), limits = list(lower = 0))
    #
    # Get data
    res <- .retrieve_from_db(type, max.hits, ...)

    # Creae prececnce abscence table from the sample types (only for animals)
    # PParameter for user to specify
    # hidden

    # flatten if user has specified. DO this first and rhen sampletypes, because
    # sample types is good data to test.

    return(res)
}
