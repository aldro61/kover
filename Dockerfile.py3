FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

# Install system requirements
RUN apt-get update && apt-get -y install build-essential cmake curl git hdf5-tools nano python3.8 python3-pip software-properties-common unzip wget zlib1g-dev
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.8 30 && \
    update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 30    

# Install Kover and requirements
COPY . kover
RUN pip install --upgrade pip cython numpy scipy && \
    cd kover && \
    sh ./install.sh
ENV PATH="/kover/bin:${PATH}"

