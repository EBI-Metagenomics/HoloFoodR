# Get data from HoloFood database

Get data from HoloFood database

## Usage

``` r
getData(
  type = NULL,
  accession.type = NULL,
  accession = NULL,
  flatten = FALSE,
  ...
)
```

## Arguments

- type:

  `NULL` or `character scalar` specifying the type of data to query.
  Must be one of the following options: `"analysis-summaries"`,
  `"animals"`, `"genome-catalogues"`, `"samples"`,
  `"sample_metadata_markers"` or `"viral-catalogues"`. When genome or
  viral catalogues is fetched by their accession ID, the type can also
  be `"genomes"` or `"fragments"`. (Default: `NULL`)

- accession.type:

  `NULL` or `character scalar` specifying the type of accession IDs.
  Must be one of the following options: `"animals"`,
  `"genome-catalogues"`, `"samples"` or `"viral-catalogues"`. (Default:
  `NULL`)

- accession:

  `NULL` or `character vector` specifying the accession IDs of type
  `accession.type`. (Default: `NULL`)

- flatten:

  `Logical scalar` specifying whether to flatten the resulting
  `data.frame`. This means that columns with multiple values are
  separated to multiple columns. (Default: `FALSE`)

- ...:

  optional arguments:

  - **max.hits** `NULL` or `integer scalar` specifying the maximum
    number of results to fetch. When NULL, all results are fetched.
    (Default: `NULL`)

  - **use.cache** `Logical scalar` specifying whether to use cache
    (Default: `FALSE`)

  - **cache.dir** `Character scalar` specifying cache directory.
    (Default: [`tempdir()`](https://rdrr.io/r/base/tempfile.html))

  - **clear.cache** `Logical scalar` specifying whether to remove and
    clear cache (Default: `FALSE`)

## Value

`list` or `data.frame`

## Details

With `getData`, you can fetch data from the database. Compared to
`getResult`, this function is more flexible since it can fetch any kind
of data from the database. However, this function returns the data
without further wrangling as `list` or `data.frame` which are not
optimized format for fetching data on samples.

Search results can be filtered; for example, animals can be filtered
based on available samples. See \[Api
browser\](https://www.holofooddata.org/api/docs) for information on
filters. You can find help on customizing queries from
\[here\](https://emg-docs.readthedocs.io/en/latest/api.html#customising-queries).

## See also

[`getResult`](getResult.md)

## Examples

``` r

# Find genome catalogues
catalogues <- getData(type = "genome-catalogues")
head(catalogues)
#>                 id                   title
#> 1  salmon-gut-v2-0  HoloFood Salmon Gut v2
#> 2 chicken-gut-v2-0 HoloFood Chicken Gut v2
#>                                         biome related_mag_catalogue_id  system
#> 1  root:Host-associated:Fish:Digestive system  non-model-fish-gut-v2-0  salmon
#> 2 root:Host-associated:Birds:Digestive system       chicken-gut-v1-0-1 chicken
#>   analysis_summaries
#> 1       c("HoloF....
#> 2       c("HoloF....

# Find genomes based on certain genome catalogue iD
res <- getData(
    type = "genomes", accession.type = "genome-catalogues",
    accession = catalogues[1, "id"], max.hits = 100)
# See the data.
head(res)
#>       accession cluster_representative
#> 1 MGYG000307500          MGYG000307500
#> 2 MGYG000307501          MGYG000307501
#> 3 MGYG000307502          MGYG000307501
#> 4 MGYG000307503          MGYG000299622
#> 5 MGYG000307504          MGYG000299579
#> 6 MGYG000307505          MGYG000299579
#>                                                                                                                      taxonomy
#> 1                                                    Bacteria > Firmicutes_A > Clostridia > Oscillospirales > Ruminococcaceae
#> 2 Bacteria > Proteobacteria > Gammaproteobacteria > Pseudomonadales > Pseudomonadaceae > Pseudomonas > Pseudomonas aeruginosa
#> 3 Bacteria > Proteobacteria > Gammaproteobacteria > Pseudomonadales > Pseudomonadaceae > Pseudomonas > Pseudomonas aeruginosa
#> 4          Bacteria > Proteobacteria > Gammaproteobacteria > Enterobacterales > Aeromonadaceae > Aeromonas > Aeromonas sobria
#> 5                                                       Bacteria > Firmicutes > Bacilli > Mycoplasmatales > Mycoplasmoidaceae
#> 6                                                       Bacteria > Firmicutes > Bacilli > Mycoplasmatales > Mycoplasmoidaceae
#>                                                representative_url metadata1
#> 1 https://www.ebi.ac.uk/metagenomics/api/v1/genomes/MGYG000307500       196
#> 2 https://www.ebi.ac.uk/metagenomics/api/v1/genomes/MGYG000307501       197
#> 3 https://www.ebi.ac.uk/metagenomics/api/v1/genomes/MGYG000307501       198
#> 4 https://www.ebi.ac.uk/metagenomics/api/v1/genomes/MGYG000299622       199
#> 5 https://www.ebi.ac.uk/metagenomics/api/v1/genomes/MGYG000299579       200
#> 6 https://www.ebi.ac.uk/metagenomics/api/v1/genomes/MGYG000299579       201
#>   metadata.Genome_type metadata.Length metadata.N_contigs metadata.N50
#> 1                  MAG         1319783                 19       153761
#> 2                  MAG         6676386                136        83355
#> 3                  MAG         4766025               4741         1104
#> 4                  MAG         2381842                367         8519
#> 5                  MAG          681326                 38        37281
#> 6                  MAG          627706                 26        43119
#>   metadata.GC_content metadata.Completeness metadata.Contamination
#> 1                27.3                 85.79                    0.0
#> 2               66.27                 98.37                  0.192
#> 3               64.82                 68.36                  2.053
#> 4               59.04                 58.95                   0.09
#> 5               25.14                 96.02                  0.384
#> 6               25.04                 95.25                  0.384
#>   metadata.rRNA_5S metadata.rRNA_16S metadata.rRNA_23S metadata.tRNAs
#> 1              0.0               0.0               0.0             18
#> 2              0.0               0.0               0.0             17
#> 3            93.28              36.2             99.04             15
#> 4             91.6               0.0               0.0             14
#> 5              0.0             99.74             99.59             19
#> 6              0.0             99.74             99.59             19
#>   metadata.Genome_accession metadata.Sample_accession metadata.Study_accession
#> 1               ERZ15182294            SAMEA112246717                ERP136460
#> 2               ERZ15182298            SAMEA112246686                ERP136460
#> 3               ERZ15182309            SAMEA112246683                ERP136460
#> 4               ERZ15182314            SAMEA112246713                ERP136460
#> 5               ERZ15233654            SAMEA112264417                ERP125469
#> 6               ERZ15233655            SAMEA112264472                ERP125469
#>   metadata.Country metadata.Continent
#> 1           Norway             Europe
#> 2           Norway             Europe
#> 3           Norway             Europe
#> 4           Norway             Europe
#> 5           Norway             Europe
#> 6           Norway             Europe
#>                                                                                                                                       metadata.FTP_download
#> 1 ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/non-model-fish-gut/v2.0/all_genomes/MGYG0003075/MGYG000307500/genomes1/MGYG000307500.gff.gz
#> 2 ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/non-model-fish-gut/v2.0/all_genomes/MGYG0003075/MGYG000307501/genomes1/MGYG000307501.gff.gz
#> 3 ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/non-model-fish-gut/v2.0/all_genomes/MGYG0003075/MGYG000307501/genomes1/MGYG000307502.gff.gz
#> 4 ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/non-model-fish-gut/v2.0/all_genomes/MGYG0002996/MGYG000299622/genomes1/MGYG000307503.gff.gz
#> 5 ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/non-model-fish-gut/v2.0/all_genomes/MGYG0002995/MGYG000299579/genomes1/MGYG000307504.gff.gz
#> 6 ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/non-model-fish-gut/v2.0/all_genomes/MGYG0002995/MGYG000299579/genomes1/MGYG000307505.gff.gz
#>   metadata.Host_species annotations.cazy.GH annotations.cazy.PL
#> 1           Salmo salar                   7                   0
#> 2           Salmo salar                  24                   0
#> 3           Salmo salar                  24                   0
#> 4           Salmo salar                  27                   0
#> 5           Salmo salar                   0                   0
#> 6           Salmo salar                   0                   0
#>   annotations.cazy.CE annotations.cazy.AA annotations.cazy.CB
#> 1                   0                   0                   1
#> 2                   3                   1                   4
#> 3                   3                   1                   4
#> 4                   1                   1                   6
#> 5                   0                   0                   0
#> 6                   0                   0                   0
#>   annotations.cazy.GT annotations.cazy.CL
#> 1                   7                   0
#> 2                  34                   0
#> 3                  34                   0
#> 4                  20                   0
#> 5                   2                   0
#> 6                   2                   0
# It includes for instance summary of the CAZy
# (Carbohydrate-Active enZymes) annotations as a counts per category
cazy <- res[ , grepl("annotations.cazy", colnames(res)), drop = FALSE]
head(cazy)
#>   annotations.cazy.GH annotations.cazy.PL annotations.cazy.CE
#> 1                   7                   0                   0
#> 2                  24                   0                   3
#> 3                  24                   0                   3
#> 4                  27                   0                   1
#> 5                   0                   0                   0
#> 6                   0                   0                   0
#>   annotations.cazy.AA annotations.cazy.CB annotations.cazy.GT
#> 1                   0                   1                   7
#> 2                   1                   4                  34
#> 3                   1                   4                  34
#> 4                   1                   6                  20
#> 5                   0                   0                   2
#> 6                   0                   0                   2
#>   annotations.cazy.CL
#> 1                   0
#> 2                   0
#> 3                   0
#> 4                   0
#> 5                   0
#> 6                   0
# Moreover, it includes a sample list. This sample list represents a 
# collection of samples where the MAG was identified. Thr data has also the
# completeness of MAG in a sample.
head(res[ c("metadata.Sample_accession", "metadata.Completeness")])
#>   metadata.Sample_accession metadata.Completeness
#> 1            SAMEA112246717                 85.79
#> 2            SAMEA112246686                 98.37
#> 3            SAMEA112246683                 68.36
#> 4            SAMEA112246713                 58.95
#> 5            SAMEA112264417                 96.02
#> 6            SAMEA112264472                 95.25
```
