FROM ubuntu:18.04

# Install system requirements
RUN apt-get update
RUN apt-get -y install cmake git hdf5-tools python2.7 python-pip unzip wget zlib1g-dev

# Install python requirements
RUN pip install --upgrade pip
RUN pip install cython numpy scipy

# Install Kover
ADD https://github.com/aldro61/kover/zipball/kover2 kover.zip
RUN unzip kover.zip && mv aldro61* kover
RUN cd kover; sh ./install.sh
ENV PATH="/kover/bin:${PATH}"
