# Copyright (c) Bioinformatics Core Facility of the Max Planck Institute for Biology of Ageing.
# Distributed under the terms of the Modified BSD License.
ARG MYAPP_IMAGE=mpgagebioinformatics/myapp:latest
FROM $MYAPP_IMAGE

LABEL maintainer "bioinformatics@age.mpg.de"

COPY app.py /myapp/myapp/routes/index.py

USER myapp