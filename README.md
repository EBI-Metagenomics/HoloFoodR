# HoloFoodR <img src="inst/extdata/holofood_logo.png" align="right" width="120" />

An R package developed for streamlining the integration and analysis of EBI
HoloFood data. This utility package simplifies access to the resource, enabling
direct loading of data into formats tailored for downstream analytics. With
this tool, users can efficiently search and retrieve data from the EBI
HoloFood resource.

The retrieved data is available in `MultiAssayExperiment` /
`SummarizedExperiment` formats, enabling effortless integration with
metagenomics data acquired through the `MGnifyR` package from the MGnify
database. This compatibility ensures users can seamlessly combine and
analyze datasets from both sources, leading to valuable insights into intricate
biological systems.

## Installation

### Bioc-release

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("HoloFoodR")
```

### Bioc-devel

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("HoloFoodR")
```

### GitHub

```
remotes::install_github("TuomasBorman/HoloFoodR")
```

## Basic usage
For more detailed instructions read the associated function help and vignette (`vignette("MGNifyR")`)

```
library(HoloFoodR)

# Search samples
samples <- doQuery("samples")

# Search animals
animals <- doQuery("animal")

# Fetch data on certain sample
samples <- c("ACCESSION_ID")
sample_data <- getData(accession.type = "samples", accession = samples)

# Fetch data on genome catalogues
genome_catalogues <- getData(type = "genome-catalogues")

# Fetch data on genomes in certain genome catalogue
catalogues <- c("ACCESSION_ID")
genomes <- getData(
    type = "genomes",
    accession.type = "genome-catalogues",
    accession = catalogues)

# Fetch data as MultiAssayExperiment
samples <- c("ACCESSION_ID")
mae <- getResult(accession = samples)
```

