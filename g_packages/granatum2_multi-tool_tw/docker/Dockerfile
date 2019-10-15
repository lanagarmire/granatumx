## Emacs, make this -*- mode: sh; -*-
 
# [tkwolf 20180712] Using specific dated version
FROM debian:testing-20180426

## This handle reaches Carl and Dirk
MAINTAINER "Carl Boettiger and Dirk Eddelbuettel" rocker-maintainers@eddelbuettel.com

## Set a default user. Available via runtime flag `--user docker` 
## Add user to 'staff' group, granting them write privileges to /usr/local/lib/R/site.library
## User should also have & own a home directory (for rstudio or linked volumes to work properly). 
RUN useradd docker \
	&& mkdir /home/docker \
	&& chown docker:docker /home/docker \
	&& addgroup docker staff

RUN apt-get update \ 
	&& apt-get install -y --no-install-recommends \
		ed \
		less \
		locales \
		vim-tiny \
		wget \
		ca-certificates \
		fonts-texgyre \
		libcurl4-openssl-dev \
		libc6-dev \
		libcairo2-dev \
		libssl-dev \
		libxt-dev \
		r-cran-slam \
		libxml2-dev \
	&& rm -rf /var/lib/apt/lists/*

## Configure default locale, see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
	&& locale-gen en_US.utf8 \
	&& /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

## Use Debian unstable via pinning -- new style via APT::Default-Release
# NOTE: "sid" is always "unstable"; better to set a specific version, e.g., "stretch" [tkwolf 20180712]
#RUN echo "deb http://http.debian.net/debian sid main" > /etc/apt/sources.list.d/debian-unstable.list \
RUN echo "deb http://http.debian.net/debian stretch main" > /etc/apt/sources.list.d/debian-testing.list \
	&& echo 'APT::Default-Release "testing";' > /etc/apt/apt.conf.d/default

# Changes to Debian repository led to loss of 3.4.4 in sid; looks like stretch has 3.3.3 and 3.5.1.
# Versions 3.3.3 and 3.5.1 have problems with the current configuration - fgsea didn't install [tkwolf 20180712]
#ENV R_BASE_VERSION 3.4.4
#ENV R_BASE_VERSION 3.3.3-1
ENV R_BASE_VERSION 3.5.1-1

## Now install R and littler, and create a link for littler in /usr/local/bin
## Also set a default CRAN repo, and make sure littler knows about it too
RUN apt-get update \
	#&& apt-get install -t unstable -y --no-install-recommends \
	&& apt-get install -t testing -y --no-install-recommends \
		littler \
                r-cran-littler \
		r-base=${R_BASE_VERSION} \
		r-base-dev=${R_BASE_VERSION} \
		r-recommended=${R_BASE_VERSION} \
        && echo 'options(repos = c(CRAN = "https://cran.rstudio.com/"), download.file.method = "libcurl")' >> /etc/R/Rprofile.site \
        && echo 'source("/etc/R/Rprofile.site")' >> /etc/littler.r \
	&& ln -s /usr/share/doc/littler/examples/install.r /usr/local/bin/install.r \
	&& ln -s /usr/share/doc/littler/examples/install2.r /usr/local/bin/install2.r \
	&& ln -s /usr/share/doc/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
	&& ln -s /usr/share/doc/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
	&& install.r docopt \
	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
	&& rm -rf /var/lib/apt/lists/*

# twolf 20180221
# To build an image with this Dockerfile:
# $ docker build -t my_image .
# To run the image:
# $ docker run -it --rm --name my_container my_image
# FROM rocker/verse:3.4.2

WORKDIR /usr/local/src/myscripts/

# Installs based on Granatum v1 depenencies install and likely have redundancies,
#   but using the full set of install commands with versions indicated to be clear.
# Performing multiple RUN calls to allow caching to speedup additing future
#   dependencies at the end of the Dockerfile.
# BASH install section
# RUN ["bash","./install.sh"]
# Note: warning about apt-utils likely can be ignored, according to forum: https://github.com/phusion/baseimage-docker/issues/319

# R install section
# RUN ["Rscript","./install.R"]
RUN R -e 'install.packages("devtools")'
# devtools_1.13.5 (above)
RUN R -e 'devtools::install_version("NMF", version = "0.20.6", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("Rtsne", version = "0.13", repos = "http://cran.us.r-project.org")'
# Perhaps unlikely to need shiny, etc.; but keeping for now (below)
RUN R -e 'devtools::install_version("shiny", version = "1.0.5", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("shinyjs", version = "0.9.1", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("htmltools", version = "0.3.6", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("htmlwidgets", version = "0.9", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("plotly", version = "4.7.1", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("Rlof", version = "1.1.1", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("scales", version = "0.5.0", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("reshape2", version = "1.4.2", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("purrr", version = "0.2.4", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("readr", version = "1.1.1", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("stringr", version = "1.2.0", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("tidyr", version = "0.7.2", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("tibble", version = "1.3.4", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("forcats", version = "0.2.0", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("dplyr", version = "0.7.4", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("igraph", version = "1.1.2", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("visNetwork", version = "2.0.1", repos = "http://cran.us.r-project.org")'

# From Bioconductor version 3.6 (Release)
# Note: some of these packages may end up in archive folders, at which time the links should be changed, e.g.,
#   from: http://bioconductor.org/packages/release/bioc/src/contrib/limma_3.34.4.tar.gz
#   to:   http://bioconductor.org/packages/3.6/bioc/src/contrib/Archive/limma/limma_3.34.4.tar.gz
# Could use biocLite in the future, but here we specify versions explicitly.
RUN R -e 'install.packages("https://www.bioconductor.org/packages/3.6/bioc/src/contrib/preprocessCore_1.40.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'devtools::install_version("futile.logger", version = "1.4.3", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("snow", version = "0.4-2", repos = "http://cran.us.r-project.org")'
# RUN R -e 'install.packages("http://bioconductor.org/packages/release/bioc/src/contrib/BiocParallel_1.12.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/bioc/src/contrib/BiocParallel_1.12.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'devtools::install_version("gridExtra", version = "2.3", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("fastmatch", version = "1.1-0", repos = "http://cran.us.r-project.org")'
# Remove below:
# RUN R -e 'install.packages("http://bioconductor.org/packages/release/bioc/src/contrib/fgsea_1.4.0.tar.gz", repos = NULL, type = "source")'
#
#RUN R -e 'install.packages("https://www.bioconductor.org/packages/release/bioc/src/contrib/BiocGenerics_0.24.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/bioc/src/contrib/BiocGenerics_0.24.0.tar.gz", repos = NULL, type = "source")'
#RUN R -e 'install.packages("https://www.bioconductor.org/packages/release/bioc/src/contrib/Biobase_2.38.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/bioc/src/contrib/Biobase_2.38.0.tar.gz", repos = NULL, type = "source")'
#RUN R -e 'install.packages("http://bioconductor.org/packages/release/bioc/src/contrib/S4Vectors_0.16.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/bioc/src/contrib/S4Vectors_0.16.0.tar.gz", repos = NULL, type = "source")'
#RUN R -e 'install.packages("https://www.bioconductor.org/packages/release/bioc/src/contrib/IRanges_2.12.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/bioc/src/contrib/IRanges_2.12.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'devtools::install_version("DBI", version = "0.7", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("RSQLite", version = "2.0", repos = "http://cran.us.r-project.org")'
#RUN R -e 'install.packages("https://www.bioconductor.org/packages/release/bioc/src/contrib/AnnotationDbi_1.40.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/bioc/src/contrib/AnnotationDbi_1.40.0.tar.gz", repos = NULL, type = "source")'
#RUN R -e 'install.packages("https://bioconductor.org/packages/release/data/annotation/src/contrib/KEGG.db_3.2.3.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/data/annotation/src/contrib/KEGG.db_3.2.3.tar.gz", repos = NULL, type = "source")'
#RUN R -e 'install.packages("https://bioconductor.org/packages/release/data/annotation/src/contrib/GO.db_3.5.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/data/annotation/src/contrib/GO.db_3.5.0.tar.gz", repos = NULL, type = "source")'
#RUN R -e 'install.packages("https://bioconductor.org/packages/release/data/annotation/src/contrib/org.Hs.eg.db_3.5.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/data/annotation/src/contrib/org.Hs.eg.db_3.5.0.tar.gz", repos = NULL, type = "source")'
#RUN R -e 'install.packages("https://bioconductor.org/packages/release/data/annotation/src/contrib/org.Mm.eg.db_3.5.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/data/annotation/src/contrib/org.Mm.eg.db_3.5.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'devtools::install_version("RCurl", version = "1.95-4.8", repos = "http://cran.us.r-project.org")'
#RUN R -e 'install.packages("http://bioconductor.org/packages/release/data/annotation/src/contrib/GenomeInfoDbData_1.0.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/data/annotation/src/contrib/GenomeInfoDbData_1.0.0.tar.gz", repos = NULL, type = "source")'
#RUN R -e 'install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/GenomeInfoDb_1.14.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/bioc/src/contrib/GenomeInfoDb_1.14.0.tar.gz", repos = NULL, type = "source")'
#RUN R -e 'install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/zlibbioc_1.24.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/bioc/src/contrib/zlibbioc_1.24.0.tar.gz", repos = NULL, type = "source")'
#RUN R -e 'install.packages("http://bioconductor.org/packages/release/bioc/src/contrib/XVector_0.18.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/bioc/src/contrib/XVector_0.18.0.tar.gz", repos = NULL, type = "source")'
# Remove below:
# RUN R -e 'install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/GenomicRanges_1.30.1.tar.gz", repos = NULL, type = "source")'
#
RUN R -e 'devtools::install_version("matrixStats", version = "0.52.2", repos = "http://cran.us.r-project.org")'
#RUN R -e 'install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/DelayedArray_0.4.1.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/bioc/src/contrib/DelayedArray_0.4.1.tar.gz", repos = NULL, type = "source")'
# Remove below:
# RUN R -e 'install.packages("http://bioconductor.org/packages/release/bioc/src/contrib/SummarizedExperiment_1.8.1.tar.gz", repos = NULL, type = "source")'
#
RUN R -e 'devtools::install_version("XML", version = "3.98-1.9", repos = "http://cran.us.r-project.org")'
#RUN R -e 'install.packages("https://www.bioconductor.org/packages/release/bioc/src/contrib/annotate_1.56.1.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/bioc/src/contrib/annotate_1.56.2.tar.gz", repos = NULL, type = "source")'
#RUN R -e 'install.packages("http://bioconductor.org/packages/release/bioc/src/contrib/genefilter_1.60.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/bioc/src/contrib/genefilter_1.60.0.tar.gz", repos = NULL, type = "source")'
# Remove below:
# RUN R -e 'install.packages("https://www.bioconductor.org/packages/release/bioc/src/contrib/sva_3.26.0.tar.gz", repos = NULL, type = "source")'
#
#RUN R -e 'install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/preprocessCore_1.40.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/bioc/src/contrib/preprocessCore_1.40.0.tar.gz", repos = NULL, type = "source")'
#RUN R -e 'install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/impute_1.52.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/bioc/src/contrib/impute_1.52.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'devtools::install_version("locfit", version = "1.5-9.1", repos = "http://cran.us.r-project.org")'
# Remove below:
# RUN R -e 'install.packages("http://bioconductor.org/packages/release/bioc/src/contrib/edgeR_3.20.2.tar.gz", repos = NULL, type = "source")'
#
RUN R -e 'devtools::install_version("MetaDE", version = "1.0.5", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("VGAM", version = "1.0-4", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("DDRTree", version = "0.1.5", repos = "http://cran.us.r-project.org")'
#RUN R -e 'install.packages("https://bioconductor.org/packages/release/data/experiment/src/contrib/HSMMSingleCell_0.112.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/data/experiment/src/contrib/HSMMSingleCell_0.112.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'devtools::install_version("fastICA", version = "1.2-1", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("densityClust", version = "0.3", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("qlcMatrix", version = "0.9.5", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("pheatmap", version = "1.0.8", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("proxy", version = "0.4-19", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("viridis", version = "0.4.0", repos = "http://cran.us.r-project.org")'
#RUN R -e 'install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/graph_1.56.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/bioc/src/contrib/graph_1.56.0.tar.gz", repos = NULL, type = "source")'
#RUN R -e 'install.packages("http://bioconductor.org/packages/release/bioc/src/contrib/RBGL_1.54.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/bioc/src/contrib/RBGL_1.54.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'devtools::install_version("RUnit", version = "0.4.31", repos = "http://cran.us.r-project.org")'
#RUN R -e 'install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/biocViews_1.46.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/bioc/src/contrib/biocViews_1.46.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'devtools::install_version("RANN", version = "2.5.1", repos = "http://cran.us.r-project.org")'
RUN R -e 'install.packages("http://bioconductor.org/packages/3.6/bioc/src/contrib/Archive/limma/limma_3.34.8.tar.gz", repos = NULL, type = "source")'
# Remove below:
# RUN R -e 'install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/monocle_2.6.1.tar.gz", repos = NULL, type = "source")'
#
RUN R -e 'devtools::install_version("shinythemes", version = "1.1.1", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_github("mohuangx/SAVER", ref="v0.3.0")'
RUN R -e 'devtools::install_version("flexmix", version = "2.3-13", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("RcppArmadillo", version = "0.8.100.1.0", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("Rook", version = "1.1-1", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("rjson", version = "0.2.15", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("Cairo", version = "1.5-9", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("quantreg", version = "5.34", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("RMTstat", version = "0.3", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("extRemes", version = "2.0-8", repos = "http://cran.us.r-project.org")'
#RUN R -e 'install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/pcaMethods_1.70.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/bioc/src/contrib/pcaMethods_1.70.0.tar.gz", repos = NULL, type = "source")'
# Remove below:
# RUN R -e 'install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/scde_2.6.0.tar.gz", repos = NULL, type = "source")'
#
RUN R -e 'devtools::install_version("jsonlite", version = "1.5")'
RUN R -e 'devtools::install_version("rmarkdown", version = "1.7", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("base64enc", version = "0.1-3")'

# Updated/added installs
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/bioc/src/contrib/Archive/edgeR/edgeR_3.20.2.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/bioc/src/contrib/Archive/fgsea/fgsea_1.4.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/bioc/src/contrib/Archive/GenomicRanges/GenomicRanges_1.30.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/bioc/src/contrib/Archive/SummarizedExperiment/SummarizedExperiment_1.8.0.tar.gz", repos = NULL, type = "source")'
#RUN R -e 'install.packages("https://www.bioconductor.org/packages/release/bioc/src/contrib/sva_3.26.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/bioc/src/contrib/sva_3.26.0.tar.gz", repos = NULL, type = "source")'
#RUN R -e 'install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/scde_2.6.0.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/bioc/src/contrib/scde_2.6.0.tar.gz", repos = NULL, type = "source")'
RUN wget https://www.dropbox.com/s/pno78mmlj0exv7s/NODES_0.0.0.9010.tar.gz
RUN R -e 'install.packages("NODES_0.0.0.9010.tar.gz", repos = NULL, type = "source")'
RUN R -e 'install.packages("https://bioconductor.org/packages/3.6/bioc/src/contrib/Archive/monocle/monocle_2.6.1.tar.gz", repos = NULL, type = "source")'
# Would need addl dependencies: RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/tidyverse/tidyverse_1.1.1.tar.gz", repos = NULL, type = "source")'

# Cedric's modified scImpute for Granatum v1 on GitHub is kept here for reference,
#   but doesn't compile at the moment.
# RUN R -e 'devtools::install_github('lanagarmire/scImpute', ref='ca_modified')'

# Pseudotime (monocle)
#COPY input_test_pseudo-time.json .
#COPY run_pseudo-time.sh .
#COPY run_monocle_from_json_2018012500.R .
#COPY run_monocle_from_json_Xun_modified.R .
#
## Normalization
#COPY input_test_normalization.json .
#COPY run_normalization.sh .
#COPY normalize.R .
#COPY granatum_sdk.R .
#COPY normalize_xzhu.R .
#
## Differential expression
## COPY input_test_de.json .
## COPY input_test_de_scde_small.json .
## COPY input_test_de_3_groups.json .
#COPY input_with_genes_de.json .
#COPY run_de.sh .
#COPY de_dyn.R .
#COPY run_de_dyn.sh .
#COPY de.R .
#COPY run_de.sh .
#
## Gene filtering
#COPY input_with_genes_de_plus.json .
#COPY filter_genes.R .
#COPY run_filter_genes.sh .
#
## Enrichment
#COPY input_test_enrichment.json .
#COPY enrichment.R .
#COPY run_enrichment.sh .
#
#COPY usage.R .
##CMD Rscript ./usage.R
#COPY run_usage.sh .

# Protein-protein interaction (PPI) network
## Running here at the end to take advantage of existing caching - can move up to front later
#RUN pip3 install --upgrade pip setuptools wheel \
#RUN pip3 install setuptools wheel \
#	&& pip3 install python-igraph
RUN apt-get update \ 
	&& apt-get install -t testing -y --no-install-recommends \
		python3-igraph=0.7.1.post6-5 \
		python3-pip \
                python3-dev \
		libcairo2-dev
RUN pip3 install --upgrade pip setuptools wheel \
	pandas=='0.23.3'
	# pycairo=='1.17.1'
RUN pip3 install pycairo=='1.17.1'
RUN pip3 install scanpy=='1.2.2' # For granatum_sdk
RUN apt-get install -t testing -y --no-install-recommends \
	python3-tk # For granatum_sdk
RUN R -e 'install.packages("pROC")'

COPY . .

CMD ./run_usage.sh
