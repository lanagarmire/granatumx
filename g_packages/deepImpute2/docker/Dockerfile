FROM python:3.6.4

MAINTAINER "Cedric Arisdakessian" cedric.arisdakessian@gmail.com

WORKDIR /usr/src/app

COPY ./requirements.txt .
RUN pip install --upgrade pip
RUN pip install --no-cache-dir -r ./requirements.txt

COPY ./deepimpute /tmp/deepimpute
COPY ./setup.py /tmp/deepimpute
RUN cd /tmp/deepimpute && pip install .

COPY ./matplotlibrc /root/.config/matplotlib/matplotlibrc

COPY . .

CMD [ "python", "./greet.py" ]
