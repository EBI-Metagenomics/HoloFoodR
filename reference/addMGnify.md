# Add results from MGnifyR to HoloFoodR results

Add results from MGnifyR to HoloFoodR results

## Usage

``` r
addMGnify(x, y, ...)

# S4 method for class 'SummarizedExperiment,MultiAssayExperiment'
addMGnify(
  x,
  y,
  exp.name1 = "metagenomic",
  exp.name2 = "metagenomic_amplicon",
  replace = TRUE,
  ...
)

# S4 method for class 'SummarizedExperiment,SummarizedExperiment'
addMGnify(x, y, ...)

# S4 method for class 'SummarizedExperiment,ANY'
addMGnify(x, y, id.col1 = "sample_biosample", id.col2 = "accession", ...)
```

## Arguments

- x:

  `SummarizedExperiment`. Results from
  [`MGnifyR::getResult()`](https://rdrr.io/pkg/MGnifyR/man/getResult.html).

- y:

  `MultiAssayExperiment` or `SummarizedExperiment` or `data.frame`-like
  table. Results from [`HoloFoodR::getResult()`](getResult.md) or sample
  metadata from it.

- ...:

  optional arguments not used currently.

- exp.name1:

  `Character scalar`. Specifies the name of experiment that will be
  added to `y`. (Default: `"metagenomic"`)

- exp.name2:

  `Character scalar`. Specifies the name of experiment from HoloFoodR
  results. This experiment is used to match IDs with MGnify data.
  (Default: `"metagenomic_amplicon"`)

- replace:

  `Logical scalar`. Whether to replace the template experiment.
  (Default: `TRUE`)

- id.col1:

  `Character scalar`. Specifies the name of column from `colData(x)`
  that includes HoloFood identifiers. (Default: `"sample_biosample"`)

- id.col2:

  `Character scalar`. Specifies the name of column from
  `colData(y[[exp.name2]])` that includes HoloFood identifiers.
  (Default: `"accession"`)

## Value

`MultiAssayExperiment`

## Details

Metagenomic data is found in MGnify rather than HoloFoodR, and the two
databases use different sample identifiers. However, MGnify's sample
metadata includes references to the identifiers used in the HoloFood
database, making it straightforward to convert sample IDs for alignment
with HoloFood data. Despite this, HoloFood contains additional metadata
not available in MGnify. Moreover, integrating data into a
`MultiAssayExperiment` while maintaining accurate sample and system
matches can be challenging.

This function is designed to simplify these tasks, enabling seamless
integration of MGnify data with HoloFood data after retrieval from the
database. You need only to input the returned data from
[`MGnifyR::getResult()`](https://rdrr.io/pkg/MGnifyR/man/getResult.html)
and [`HoloFoodR::getResult()`](getResult.md) functions.

## Examples

``` r

if (FALSE) { # \dontrun{
# Get data from HoloFood database
mae <- HoloFoodR::getResult(
    salmon_sample_ids,
    use.cache = TRUE
)

# Get data from MGnify database
mg <- MgnifyClient(
    useCache = TRUE,
    cacheDir = ".MGnifyR_cache"
)
tse <- MGnifyR::getResult(
    mg,
    accession = mgnify_analyses_ids,
    get.func = FALSE
)

# Add MGnify data to HoloFood data
mae <- addMGnify(tse, mae)
} # }
```
