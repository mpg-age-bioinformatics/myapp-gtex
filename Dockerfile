# Copyright (c) Bioinformatics Core Facility of the Max Planck Institute for Biology of Ageing.
# Distributed under the terms of the Modified BSD License.
ARG MYAPP_IMAGE=mpgagebioinformatics/myapp:latest
FROM $MYAPP_IMAGE

LABEL maintainer "bioinformatics@age.mpg.de"

COPY app.py /myapp/myapp/routes/index.py
COPY _app.py /myapp/myapp/routes/_app.py

RUN mkdir -p /flaski_private/aarnaseqlake/

COPY requirements.txt /requirements.txt
RUN pip3 install -r /requirements.txt

USER myapp