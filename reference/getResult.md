# Get data on samples from HoloFood database

Get data on samples from HoloFood database

## Usage

``` r
getResult(accession, ...)
```

## Arguments

- accession:

  `Character vector` specifying the accession IDs of type samples.

- ...:

  optional arguments:

  - **use.cache** `Logical scalar` specifying whether to use cache. Note
    that when `get.metabolomic = TRUE` is specified, the file from the
    MetaboLights is stored in the local system to the location specified
    by `cache.dir` despite of the value of `use.cache`. (Default:
    `FALSE`)

  - **cache.dir** `Character scalar` specifying cache directory.
    (Default: [`tempdir()`](https://rdrr.io/r/base/tempfile.html))

  - **clear.cache** `Logical scalar` specifying whether to use.cache
    (Default: `FALSE`)

  - **assay.type** `Character scalar` specifying the name of assay in
    resulting `TreeSummarizedExperiment` object. (Default: `"counts"`)

  - **get.metabolomic** `Logical scalar` specifying whether to retrieve
    processed metabolomic data from MetaboLights database. For
    retrieving spectra data, refer to
    [`getMetaboLights`](getMetaboLights.md) documentation. (Default:
    `FALSE`)

## Value

`MultiAssayExperiment`

## Details

With `getResult`, you can fetch data on samples from the HoloFood
database. Compared to `getData`, this function is more convenient for
fetching the samples data because it converts the data to
`MultiAssayExperiment` where different omics are stored as
`TreeSummarizedExperiment` objects which are optimized for downstream
analytics. Columns of returned `MultiAssayExperiment` are individual
animals. These columns are linked with individual samples that are
stored in `TreeSummarizedExperiment` objects.

The HoloFood database lacks non-targeted metabolomic data but they can
be fetched from MetaboLights resource. Certain datasets include
processed features. Those datasets can be retrieved with the function
`getResult` which integrates metabolomic data with other datasets from
HoloFood.

Furthermore, while the HoloFoodR database does not include metagenomic
assembly data, users can access such data from the MGnify database. The
MGnifyR package provides a convenient interface for accessing this
database. By employing
[`MGnifyR::getResult()`](https://rdrr.io/pkg/MGnifyR/man/getResult.html),
users can obtain data formatted as a `MultiAssayExperiment` object,
containing multiple `TreeSummarizedExperiment` objects. Consequently,
data from both HoloFood and MGnify databases are inherently compatible
for subsequent downstream analysis.

## See also

[`getData`](getData.md)
[`TreeSummarizedExperiment`](https://rdrr.io/pkg/TreeSummarizedExperiment/man/TreeSummarizedExperiment-class.html)
[`MultiAssayExperiment`](https://github.com/waldronlab/MultiAssayExperiment/reference/MultiAssayExperiment-class.html)
[`MGnifyR:getResult`](https://rdrr.io/pkg/MGnifyR/man/getResult.html)

## Examples

``` r

# Find samples on certain animal
samples <- doQuery("samples", animal_accession = "SAMEA112904746")

# Get the data
mae <- getResult(samples[["accession"]])
#> Warning: Data for the following samples cannot be found. The sample types are metagenomic_assembly, host_genomic, transcriptomic and metatranscriptomic (Note that metagenomic assemblies can be found from the MGnify database. See MGnifyR package.):
#> 'SAMEA10130039', 'SAMEA13389496', 'SAMEA13389497', 'SAMEA13901618', 'SAMEA13901619', 'SAMEA13929785', 'SAMEA7571815'
mae
#> A MultiAssayExperiment object of 8 listed
#>  experiments with user-defined names and respective classes.
#>  Containing an ExperimentList class object of length 8:
#>  [1] BIOGENIC AMINES: TreeSummarizedExperiment with 7 rows and 2 columns
#>  [2] FATTY ACIDS: TreeSummarizedExperiment with 19 rows and 2 columns
#>  [3] HISTOLOGY: TreeSummarizedExperiment with 20 rows and 2 columns
#>  [4] INFLAMMATORY MARKERS: TreeSummarizedExperiment with 14 rows and 2 columns
#>  [5] metagenomic_assembly: TreeSummarizedExperiment with 0 rows and 2 columns
#>  [6] host_genomic: TreeSummarizedExperiment with 0 rows and 2 columns
#>  [7] transcriptomic: TreeSummarizedExperiment with 0 rows and 2 columns
#>  [8] metatranscriptomic: TreeSummarizedExperiment with 0 rows and 1 columns
#> Functionality:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DataFrame
#>  sampleMap() - the sample coordination DataFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DataFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
#>  exportClass() - save data to flat files
```
