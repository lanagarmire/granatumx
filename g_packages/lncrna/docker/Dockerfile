FROM rocker/verse:3.4.3

MAINTAINER "Tianying Zhao" tianying@hawaii.edu

WORKDIR /usr/src/app

RUN R --slave -e 'devtools::install_version("gbm", version = "2.1.4", repos = "http://cran.us.r-project.org")'

# RUN R --slave -e 'devtools::install_version("NMF", version = "0.20.6", repos = "http://cran.us.r-project.org")'
# RUN R --slave -e 'install.packages("https://www.bioconductor.org/packages/3.6/bioc/src/contrib/preprocessCore_1.40.0.tar.gz", repos = NULL, type = "source")'

COPY . .
CMD [ "Rscript", "./greet.R" ]
