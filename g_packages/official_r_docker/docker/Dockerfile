FROM continuumio/anaconda3:2019.03

MAINTAINER "Xun Zhu" zhuxun2@gmail.com

WORKDIR /usr/src/app

RUN apt-get install -y vim

RUN conda install -y -c r r-base rpy2

# -------

RUN pip install git+https://gitlab.com/xz/granatum_sdk.git@1602ed8510cb84a24725fb966446d36b90e0f070
RUN pip install scanpy==1.4.1

COPY ./install_from_CRAN.R .

RUN Rscript ./install_from_CRAN.R \
            BiocManager

COPY ./install_from_bioconductor.R .

RUN apt-get install -y --no-install-recommends \
    gcc \
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
		libxml2-dev

RUN conda install -y -c r r-rcurl
RUN conda install -y -c r r-xml
# RUN Rscript ./install_from_bioconductor.R GenomeInfoDb
# RUN Rscript ./install_from_bioconductor.R GenomicRanges
# RUN Rscript ./install_from_bioconductor.R BiocParallel
# RUN conda install -y -c conda-forge gcc

RUN Rscript ./install_from_bioconductor.R SummarizedExperiment
RUN Rscript ./install_from_bioconductor.R genefilter
RUN Rscript ./install_from_bioconductor.R geneplotter
RUN Rscript ./install_from_bioconductor.R DESeq2

COPY . .

CMD [ "Rscript", "./greet.R" ]
