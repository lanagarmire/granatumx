FROM python:3.6.4

MAINTAINER "Arr N. Eh" rneh@gmail.com

WORKDIR /usr/src/app

RUN pip install pandas==0.23.4

RUN pip install git+https://gitlab.com/xz/granatum_sdk@master

COPY . .

CMD [ "python", "./greet.py" ]
