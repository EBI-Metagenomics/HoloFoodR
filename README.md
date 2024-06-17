# HoloFoodR <img src="man/figures/holofoodr_logo.png" align="right" width="120" />

An R package developed for streamlining the integration and analysis of EMBL-EBI
HoloFood data. This utility package simplifies access to the resource, enabling
direct loading of data into formats tailored for multiomics downstream analytics.
With this tool, users can efficiently search and retrieve data from the EBI
HoloFood resource.

The retrieved data is available in `MultiAssayExperiment` /
`TreeSummarizedExperiment` format similarly to the data acquired with the
`MGnifyR` package from the MGnify database. This compatibility ensures users
can seamlessly combine and analyze datasets from both sources, leading to
valuable insights into intricate biological systems.

This research has received funding from the Horizon 2020 Programme of the
European Union within the FindingPheno project under grant agreement No 952914.

## FindingPheno <img src="man/figures/findingpheno_logo.png" align="right" width="160" />

FindingPheno, an EU-funded project, is dedicated to developing computational
tools and methodologies for the integration and analysis of multi-omics data.
Its primary objective is to deepen our understanding of the interactions
between hosts and their microbiomes.

- [FindingPheno website](https://findingpheno.eu/)

## HoloFood <img src="man/figures/holofood_logo.png" align="right" width="60" />

HoloFood, a project funded under EU's Horizon 2020 programme, employed a
holistic, "hologenomic", approach to enhance the efficiency of food production 
systems. This involved exploring the biomolecular and physiological processes 
triggered by the incorporation of feed additives and novel sustainable feeds in
farmed animals.

The HoloFood database, hosted by European Bioinformatics Institute (EMBL-EBI), 
houses data gathered during the project, encompassing multiple omics, including
metabolomics and various other biomolecular measurements. Notably, it does not
include data on metagenomic and  untargeted metabolomic analyses. However,
metagenomic data from the project can be accessed through the MGnify database,
while untargeted metabolomic data is stored in the MetaboLights database. To
explore available datasets in HoloFood, you can utilize the API browser.

- [HoloFood website](https://www.holofood.eu/)
- [HoloFood publications](https://www.holofood.eu/publications.html)
- [API browser](https://www.holofooddata.org/)
- [API documentation](https://docs.holofooddata.org/api.html)

## MGnify <img src="man/figures/mgnify_logo.jpg" align="right" width="60" />

EMBL-EBI's MGnify serves as a repository for microbiome data, offering a wide array
of analyses encompassing metabarcoding, metatranscriptomic, and metagenomic
datasets from diverse environments. This platform provides comprehensive
taxonomic assignments and functional annotations for these datasets. The data
can be accessed through MGnifyR package.

- [MGnify website](https://www.ebi.ac.uk/metagenomics)
- [MGnifyR](https://github.com/EBI-Metagenomics/MGnifyR)

## MetaboLights <img src="man/figures/metabolights_logo.jpg" align="right" width="160" />

MetaboLights, managed by EMBL-EBI, serves as a repository for metabolomic data.
It can be accessed through HoloFoodR package.

- [MetaboLights website](https://www.ebi.ac.uk/metabolights/)

## miaverse <img src="man/figures/mia_logo.png" align="right" width="60" />

miaverse is an R/Bioconductor framework specialized for microbiome downstream
data analysis, leveraging the TreeSummarizedExperiment class. It offers a
comprehensive suite of tools for microbiome bioinformatics. Additionally,
miaverse includes the tutorial book Orchestrating Microbiome Analysis (OMA),
which aims to guide users in conducting microbiome data analysis and sharing
best practices in microbiome data science.

- [miaverse website](https://microbiome.github.io/)
- [Orchestrating Microbiome Analysis (OMA)](https://microbiome.github.io/OMA/docs/devel/)

--------------------------------------------------------------------------------

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
BiocManager::install(version="devel")

BiocManager::install("HoloFoodR")
```

### GitHub

```
remotes::install_github("EBI-Metagenomics/HoloFoodR")
```

## Basic usage
For more detailed instructions read the associated function help,
[function reference page](https://EBI-Metagenomics.github.io/HoloFoodR/) and
vignette (`vignette("HoloFoodR")`)

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

# Fetch data on untargeted metabolites
metabolites <- getMetaboLights(url)

# Fetch data as MultiAssayExperiment
samples <- c("ACCESSION_ID")
mae <- getResult(accession = samples)
```
