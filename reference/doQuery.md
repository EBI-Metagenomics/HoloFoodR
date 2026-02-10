# Search HoloFood database for animals, genome catalogues, samples, or viral catalogues

Search HoloFood database for animals, genome catalogues, samples, or
viral catalogues

## Usage

``` r
doQuery(type, flatten = TRUE, ...)
```

## Arguments

- type:

  `Character scalar` specifying the type of data to query. Must be one
  of the following options: `"animals"`, `"genome-catalogues"`,
  `"samples"` or `"viral-catalogues"`.

- flatten:

  `Logical scalar` specifying whether to flatten the resulting
  `data.frame`. This means that columns with multiple values are
  separated to multiple columns. (Default: `TRUE`)

- ...:

  optional arguments:

  - **max.hits** `NULL` or `integer scalar` specifying the maximum
    number of results to fetch. When NULL, all results are fetched.
    (Default: `NULL`)

  - **spread.sample.types** `Logical scalar` specifying whether to
    create spread sample types column of animals data. In animals data,
    sample types column might have multiple values that might be hard to
    explore. This argument specifies whether to create presence/absence
    table from sample types. (Default: `TRUE`)

  - **use.cache** `Logical scalar` specifying whether to use cache.
    (Default: `FALSE`)

  - **cache.dir** `Character scalar` specifying cache directory.
    (Default: [`tempdir()`](https://rdrr.io/r/base/tempfile.html))

  - **clear.cache** `Logical scalar` specifying whether to remove and
    clear cache (Default: `FALSE`)

## Value

`data.frame`

## Details

`doQuery` is a flexible query function which can be utilized to search
available animals, genome catalogues, samples, or viral catalogues.
Search results can be filtered; for example, animals can be filtered
based on available samples. See \[Api
browser\](https://www.holofooddata.org/api/docs) for information on
filters. You can find help on customizing queries from
\[here\](https://emg-docs.readthedocs.io/en/latest/api.html#customising-queries).

## Examples

``` r

# Find animals results. The maximum amount of results is 100. Use filter
# so that only chicken is searched. (See details on customizing queries)
res <- doQuery("animals", max.hits = 100, system = "chicken")
head(res)
#>        accession  system                                   canonical_url
#> 1 SAMEA112904733 chicken https://www.ebi.ac.uk/biosamples/SAMEA112904733
#> 2 SAMEA112904734 chicken https://www.ebi.ac.uk/biosamples/SAMEA112904734
#> 3 SAMEA112904735 chicken https://www.ebi.ac.uk/biosamples/SAMEA112904735
#> 4 SAMEA112904736 chicken https://www.ebi.ac.uk/biosamples/SAMEA112904736
#> 5 SAMEA112904737 chicken https://www.ebi.ac.uk/biosamples/SAMEA112904737
#> 6 SAMEA112904738 chicken https://www.ebi.ac.uk/biosamples/SAMEA112904738
#>   histological host_genomic inflammatory_markers metabolomic
#> 1        FALSE         TRUE                FALSE       FALSE
#> 2         TRUE         TRUE                 TRUE       FALSE
#> 3         TRUE         TRUE                 TRUE       FALSE
#> 4        FALSE         TRUE                FALSE       FALSE
#> 5         TRUE         TRUE                 TRUE       FALSE
#> 6        FALSE         TRUE                FALSE       FALSE
#>   metabolomic_targeted metagenomic_assembly metatranscriptomic transcriptomic
#> 1                FALSE                 TRUE              FALSE          FALSE
#> 2                 TRUE                 TRUE               TRUE           TRUE
#> 3                 TRUE                 TRUE              FALSE           TRUE
#> 4                FALSE                 TRUE              FALSE          FALSE
#> 5                 TRUE                 TRUE               TRUE           TRUE
#> 6                 TRUE                 TRUE              FALSE          FALSE
```
