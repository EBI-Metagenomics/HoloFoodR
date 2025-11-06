FROM ghcr.io/bioconductor/bioconductor_docker:RELEASE_3_22-R-4.5.2

WORKDIR /project
COPY DESCRIPTION DESCRIPTION
RUN R -e 'install.packages("remotes"); \
    remotes::install_deps(pkgdir = ".", dependencies = "Suggests")'

