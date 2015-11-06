---
title: Installation
keywords: installation, install, kover, machine learning, genomics
last_updated: November 5, 2015
tags: [getting-started]
summary: "This page will walk you through the installation of Kover."
---

## Prerequisites
 
Kover is written in Python and Cython. Most of its dependencies should be installed automatically when 
running the installer. You can also install these dependencies manually or use your operating system's package manager.

What will not be installed automatically:

* [GNU C++ compiler (g++)]()
* [GNU Fortran (gfortran)]()
* [The HDF5 library](https://www.hdfgroup.org/HDF5/release/obtain5.html)
* [Python 2.7.x](https://www.python.org/download/releases/2.7/)
* [Python header files](https://docs.python.org/2/c-api/intro.html)

What should be installed automatically:

* [Cython](http://docs.cython.org/src/quickstart/install.html)
* [H5py](http://docs.h5py.org/en/latest/build.html)
* [Numpy](http://docs.scipy.org/doc/numpy/user/install.html)
* [Pandas](http://pandas.pydata.org/pandas-docs/stable/install.html#installing-pandas)
* [Progressbar](https://pypi.python.org/pypi/progressbar)

## Linux and Mac

Download the latest version of Kover from the GitHub repository:

```
 wget https://github.com/aldro61/kover/archive/master.zip
```

or alternatively

```
git clone https://github.com/aldro61/kover.git
```

Then, in the *kover* directory, run the installer:

```
sh ./install.sh
```

This will build and install Kover and its dependencies. A *bin* directory containing the Kover executable will be created.
You can add the *bin* directory to your $PATH environment variable by adding this command to your ~/.bashrc file:

```
export PATH=[PATH_TO_KOVER_DIRECTORY]/bin:$PATH
```

Now change directories and test the installation:

```
cd
kover -v
```

This should print the version of you Kover installation.

Thats it! You're done.

## Windows

Windows is currently not supported, but will be eventually.