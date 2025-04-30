FROM python:3.10-slim

WORKDIR /root
ENV VENV /opt/venv
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8
ENV PYTHONPATH /root

RUN apt-get update && apt-get install -y build-essential

ENV VENV /opt/venv

RUN python3 -m venv ${VENV}
ENV PATH="${VENV}/bin:$PATH"

# Install from PyPI
RUN pip install --upgrade pip
RUN pip install crispr-ambiguous-mapping==0.0.200
