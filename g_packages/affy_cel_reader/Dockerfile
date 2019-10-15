FROM rocker/verse:3.4.3

MAINTAINER "Xun Zhu" zhuxun2@gmail.com

WORKDIR /usr/src/app

RUN install2.r --deps TRUE --error base64enc

COPY ./install_from_bioconductor.R .
RUN Rscript ./install_from_bioconductor.R

COPY . .

CMD [ "Rscript", "./greet.R" ]
