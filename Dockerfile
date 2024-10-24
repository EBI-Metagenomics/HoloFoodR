ARG BIOC_VERSION
FROM bioconductor/bioconductor_docker:${BIOC_VERSION}

# Copy repository to Docker image
COPY --chown=rstudio:rstudio . /opt/pkg

# Install essentials
RUN apt-get install -y python3 python3-setuptools python3-dev python3-pip

# Install HoloFoodR dependencies
RUN Rscript -e 'repos <- BiocManager::repositories(); \
    remotes::install_local(path = "/opt/pkg/", repos=repos, \
    dependencies=TRUE, build_vignettes=FALSE, upgrade=TRUE); \
    sessioninfo::session_info(installed.packages()[,"Package"], \
    include_base = TRUE)'

# Istall CRAN packages for case study
RUN Rscript -e 'install.packages(c("DT", "patchwork", "reticulate", "reshape", "shadowtext", "shadowtext", \
    "scater", "ggsignif", "stringr", "ggpubr", "GGally", "ggplot2", "knitr", "latex2exp", "UpSetR"))'

# Install Bioconductor packages for case study
RUN R -e 'BiocManager::install(c("basilisk", "biomformat", \
    "ComplexHeatmap", "MGnifyR", "mia", "miaViz", "MOFA2", "tidyverse"))'

# Install latest development version of HoloFoodR
RUN R -e 'remotes::install_github("EBI-Metagenomics/HoloFoodR")'

# Install mofapy2 for case study
RUN python3 -m pip install 'https://github.com/bioFAM/mofapy2/tarball/master'

# Internal port for RStudio server is 8787
EXPOSE 8787
