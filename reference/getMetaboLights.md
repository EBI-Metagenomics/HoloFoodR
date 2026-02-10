# Get metabolomic data from MetaboLights database

Get metabolomic data from MetaboLights database

## Usage

``` r
getMetaboLights(study.id, output = "list", ...)

getMetaboLightsFile(study.id, file, ...)
```

## Arguments

- study.id:

  `character vector` specifying the study identifier of data that is
  going to be fetched from the MetaboLights database.

- output:

  `character scalar` specifying output format. Must be `"list"`,
  `"TreeSE"` (`TreeSummarizedExperiment`) or `"SE"`
  (`SummarizedExperiment`). (Default: `"list"`)

- ...:

  optional arguments:

  - **cache.dir** `Character scalar` specifying directory where
    downloaded file is stored. (Default:
    [`tempdir()`](https://rdrr.io/r/base/tempfile.html))

  - **timeout** `Integer scalar` specifying timeout in seconds for
    loading a file. (Default: `5*60`)

  - **ion.mode** `Character scalar` specifying metabolite assignment
    files to fetch. If `"positive"` only positive ions are fetched.
    Similarly `"negative"` means that negative ions are fetched if such
    data exists. By selecting `"both"`, one can fetch both positive and
    negative ions. (Default: `"both"`)

- file:

  `character vector` specifying the files that are being fetched.

## Value

`list`, `SummarizedExperiment` or `TreeSummarizedExperiment`

## Details

The HoloFood database primarily comprises targeted metabolomic data,
omitting non-targeted metabolomic information. Nonetheless, it features
URLs linking to studies within the MetaboLights database. This
functionality enables users to access non-targeted metabolomic data. The
`getMetaboLights` function returns a structured list encompassing
processed data in `data.frame` format for study metadata, assay
metadata, and assay.

The metadata includes the file names of spectra data. Those files can be
loaded with `getMetaboLightsFile`. Alternatively, once you've identified
the study and files to fetch, you can refer to this
\[vignette\](https://rformassspectrometry.github.io/MsIO/articles/MsIO.html#loading-data-from-metabolights)
for instructions on loading the data directly into an `MsExperiment`
object, specifically designed for metabolomics spectra data.

## See also

[`getResult`](getResult.md) [`getData`](getData.md)

## Examples

``` r

# This example is not run, because the server fails to respond sometimes.
if( FALSE ){
    res <- getMetaboLights("MTBLS4381")
    file_paths <- getMetaboLightsFile(
        study.id = "MTBLS4381",
        file = res[["assay_meta"]][["Raw Spectral Data File"]]
        )
    # Get data as SummarizedExperiment
    se <- getMetaboLights("MTBLS3540", output = "SE")
}
```
