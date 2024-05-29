# Copyright (c) Bioinformatics Core Facility of the Max Planck Institute for Biology of Ageing.
# Distributed under the terms of the Modified BSD License.
ARG MYAPP_IMAGE=mpgagebioinformatics/flaski:latest
FROM $MYAPP_IMAGE

LABEL maintainer "bioinformatics@age.mpg.de"

COPY app.py /myapp/myapp/routes/index.py
COPY _app.py /myapp/myapp/routes/_app.py

USER myapp