FROM r-base:3.4.3

WORKDIR /usr/app/scimpute_module

RUN apt-get update -y \
      && apt-get install -y python3 python3-pip

COPY ./python_requirements.txt .
RUN pip3 install --trusted-host pypi.python.org -r python_requirements.txt

COPY ./install_R_packages_from_CRAN.R .
RUN Rscript ./install_R_packages_from_CRAN.R \
       foreach \
       Matrix \
       parallel \
       stats \
       utils \
       glmnet

#COPY ./install_R_packages_from_Bioconductor.R .
#RUN Rscript ./install_R_packages_from_Bioconductor.R \
#       BiocInstaller \
#       littler \
#       MASS \
#       mgcv \
#       rpart

COPY ./update_biocLite.R .
RUN Rscript ./update_biocLite.R

COPY ./additional_requirements.txt .
RUN pip3 install --trusted-host pypi.python.org -r additional_requirements.txt

COPY . .

CMD ["python3", "greet.py"]

