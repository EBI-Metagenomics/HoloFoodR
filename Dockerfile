ARG BIOC_VERSION
FROM bioconductor/bioconductor_docker:${BIOC_VERSION}

# Copy repository to Docker image
COPY . /opt/pkg

# Install essentials
RUN apt-get install -y python3 python3-setuptools python3-dev python3-pip

# Install HoloFoodR
RUN Rscript -e 'repos <- BiocManager::repositories(); \
    remotes::install_local(path = "/opt/pkg/", repos=repos, \
    dependencies=TRUE, build_vignettes=FALSE, upgrade=TRUE); \
    sessioninfo::session_info(installed.packages()[,"Package"], \
    include_base = TRUE)'

# Istall CRAN packages
RUN Rscript -e 'install.packages(c("patchwork", "randomForest", "reticulate", "reshape", "rJava", "shadowtext", "shadowtext", \
    "scater", "ggsignif", "SuperLearner", "mcmcplots", "stringr", "ggpubr", \
    "bartMachine", "GGally", "ggplot2", "knitr"))'

# Install Bioconductor packages
RUN R -e 'BiocManager::install(c("basilisk", "biomformat", \
    "ComplexHeatmap", "MGnifyR", "mia", "miaViz", "MOFA2", "tidyverse"))'

# Install IntegratedLearner
RUN R -e 'devtools::install_github("himelmallick/IntegratedLearner")'

# Install mofapy2
RUN python3 -m pip install 'https://github.com/bioFAM/mofapy2/tarball/master'

# Internal port for RStudio server is 8787
EXPOSE 8787
