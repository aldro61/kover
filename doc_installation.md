---
title: Installation
keywords: installation, install, kover, machine learning, genomics
last_updated: November 5, 2015
tags: [getting-started]
summary: "This page will walk you through the installation of Kover."
---

## Prerequisites
 
Kover is written in Python and Cython. Most of the dependecies listed below should be installed automatically when running the installer:

* [H5py](http://docs.h5py.org/en/latest/build.html) (and the HDF5 library)
* [Numpy](http://docs.scipy.org/doc/numpy/user/install.html)
* [Pandas](http://pandas.pydata.org/pandas-docs/stable/install.html#installing-pandas)
* [Progressbar](https://pypi.python.org/pypi/progressbar)
* [Python 2.7.x](https://www.python.org/download/releases/2.7/)

## Linux and Mac

Download the latest version of Kover from the GitHub repository:

```
 wget https://github.com/aldro61/kover/archive/master.zip
```

or alternatively

```
git clone https://github.com/aldro61/kover.git
```

Then, in the kover directory, run the installer:

```
sh ./install.sh
```

This will build and install Kover and its dependencies. A *bin* directory containing the Kover executable will be created. You can add the *bin* directory to your $PATH environment variable:

```
export PATH=[PATH_TO_KOVER_DIRECTORY]/bin:$PATH
```

Thats it! You're done.