library(HoloFoodR)
library(MGnifyR)
#
animals <- doQuery("animals", max.hits = 100)
animals <- animals[ animals[["histological"]], ]
animal_data <- getData(
    accession.type = "animals", accession = animals[["accession"]])
samples <- animal_data[["samples"]]
sample_ids <- unique(samples[["accession"]])
mae <- getResult(sample_ids)
#
mg <- MgnifyClient(useCache = TRUE)
metagenomic_samples <- samples[
    samples[["sample_type"]] == "metagenomic_assembly", ]
analysis_ids <- searchAnalysis(
    mg, type = "samples", metagenomic_samples[["accession"]])
path <- system.file("extdata", "analysis_ids.rds", package = "HoloFoodR")
saveRDS(analysis_ids, path)
mae_metagenomic <- MGnifyR::getResult(mg, analysis_ids)
path <- system.file("extdata", "mae_metagenomic.rds", package = "HoloFoodR")
mae_metagenomic <- saveRDS(mae_metagenomic, path)
