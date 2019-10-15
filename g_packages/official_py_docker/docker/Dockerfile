FROM python:3.6.8
# FROM python:3.6.4
# FROM python:2.7.14

MAINTAINER "Xun Zhu" zhuxun2@gmail.com

WORKDIR /usr/src/app

# RUN apt-get -y update

# RUN apt-get -y install vim less

RUN pip install cmake==3.11.0

COPY ./requirements.txt .
RUN pip install --no-cache-dir -r ./requirements.txt

COPY ./granatum_clustering/ ./granatum_clustering
RUN cd ./granatum_clustering && pip install -e .

COPY ./granatum_deeplearning/ ./granatum_deeplearning
RUN cd ./granatum_deeplearning && pip install -e .

COPY ./scanpy/ ./scanpy
RUN cd ./scanpy && pip install -e .

RUN pip install git+https://github.com/DmitryUlyanov/Multicore-TSNE.git@682531fe21db7e10c1f7b0a783b7be86128273bc

# COPY ./granatum_sdk/ ./granatum_sdk
# RUN cd ./granatum_sdk && pip install -e .

COPY ./zgsea ./zgsea
RUN cd ./zgsea && pip install -e .

COPY ./jammit/ ./jammit
RUN cd ./jammit && pip install -e .

RUN pip install git+https://gitlab.com/xz/granatum_sdk.git@2172359

RUN pip install git+https://github.com/chriscainx/mnnpy.git@750bba42ef7e26ac1eafe75e772440b12160814c


COPY ./keras.json /root/.keras/keras.json
COPY ./matplotlibrc /root/.config/matplotlib/matplotlibrc
COPY . .

CMD [ "python", "./greet.py" ]
