FROM bioconductor/bioconductor_docker:RELEASE_3_19

#Install essentials
RUN apt-get install -y python3 python3-setuptools python3-dev python3-pip

# Istall CRAN packages
RUN R -e 'install.packages(c("bartMachine", "GGally", "ggplot2", "knitr", "patchwork", "randomForest", "reticulate", "reshape", "rJava", "shadowtext", "ggsignif"))'

# Install Bioconductor packages
RUN R -e 'BiocManager::install(c("basilisk", "biomformat", "ComplexHeatmap", "MGnifyR", "mia", "miaViz", "MOFA2"))'

# Install GitHub packages
RUN R -e 'devtools::install_github("EBI-Metagenomics/HoloFoodR")'
RUN R -e 'devtools::install_github("himelmallick/IntegratedLearner")'

# Install Quarto
RUN TEMP_DEB="$(mktemp)" && wget -O "$TEMP_DEB" 'https://github.com/quarto-dev/quarto-cli/releases/download/v1.5.47/quarto-1.5.47-linux-amd64.deb' && sudo dpkg -i "$TEMP_DEB"
RUN rm -f "$TEMP_DEB"

# Install mofapy2
RUN python3 -m pip install 'https://github.com/bioFAM/mofapy2/tarball/master'

# Create directory for HolofoodR 
WORKDIR /root/holofoodr

# Internal port for RStudio server is 8787
EXPOSE 8787

