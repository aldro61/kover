FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

# Install system requirements
RUN apt-get update
RUN apt-get -y install build-essential cmake curl git hdf5-tools nano software-properties-common unzip wget zlib1g-dev

# Install Python 2 (needed for Kover)
RUN add-apt-repository universe && \
    apt-get -y install python2.7 python-dev && \
    update-alternatives --install /usr/bin/python python /usr/bin/python2.7 30 && \
    curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py && \
    python ./get-pip.py

# Install Kover and requirements
ADD https://github.com/aldro61/kover/zipball/kover2 kover.zip
RUN pip install --upgrade pip cython numpy scipy && \
    unzip kover.zip && mv aldro61* kover && \
    cd kover && \
    sh ./install.sh
ENV PATH="/kover/bin:${PATH}"
