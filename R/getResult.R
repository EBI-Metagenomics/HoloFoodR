# Get data with getData. Format it so that we get SE. accession must be specified and it must be samples

getResult <- function(accession, ...){
    # Of user tries to feed accession.type or type, grab it
    args <- list(...)
    args <- args[ !names(args) %in% c("type", "accession.type", "flatten")]
    args[["accession.type"]] <- "samples"
    args[["accession"]] <- accession
    args[["flatten"]] <- TRUE
    res <- do.call(getData, args)



    # Convert result to SE

    metadata_markers <- getData(type = "sample_metadata_markers")

    # Metagenomic cannot be included
    # The samples include metageomic that cannot eb included see MGnifyR package.
    return(res)
}
