---
title: Installation
keywords: installation, install, kover, machine learning, genomics
last_updated: May 21, 2016
tags: [getting-started]
summary: "This page will walk you through the installation of Kover."
---

## Prerequisites
 
Kover is written in Python and Cython. We provide an installer that should install most of the required Python packages 
automatically. However, some dependencies must be installed manually or using **your operating system's 
package manager** (e.g.: [apt-get](http://linux.die.net/man/8/apt-get) on Linux and [Homebrew](http://brew.sh/) on Mac).
For each dependency, a link to further installation instructions is provided.

### Need to install yourself

* [GNU C++ compiler (g++)](https://gcc.gnu.org/)
* [GNU Fortran (gfortran)](https://gcc.gnu.org/wiki/GFortran)
* [The HDF5 library](https://www.hdfgroup.org/HDF5/release/obtain5.html)
* [NumPy](http://docs.scipy.org/doc/numpy/user/install.html)*
* [Python 2.7.x](https://www.python.org/download/releases/2.7/)
* [Python development headers](https://docs.python.org/2/c-api/intro.html)
* [SciPy](https://github.com/scipy/scipy/releases)*

### Will be installed automatically

* [Cython](http://docs.cython.org/src/quickstart/install.html)
* [H5py >= 2.4.0](http://docs.h5py.org/en/latest/build.html)
* [Numpy](http://docs.scipy.org/doc/numpy/user/install.html)*
* [Pandas](http://pandas.pydata.org/pandas-docs/stable/install.html#installing-pandas)
* [Progressbar](https://pypi.python.org/pypi/progressbar)
* [SciPy](https://github.com/scipy/scipy/releases)*

*Numpy and Scipy can be installed automatically, but we recommend a manual installation for optimal performances.


## Linux and Mac

Download the latest version of Kover from the GitHub repository by clicking [here](https://github.com/aldro61/kover/archive/master.zip), or alternatively:

```
 wget https://github.com/aldro61/kover/archive/master.zip
```

or

```
git clone --recursive https://github.com/aldro61/kover.git
```

Then, in the *kover* directory, run the installer (you might have to run it as super user):

```
./install.sh
```

This will build and install Kover and its dependencies. A *bin* directory containing the Kover executable will be created.
You can add the *bin* directory to your $PATH environment variable by adding this command to your ~/.bashrc file:

```
export PATH=[PATH_TO_KOVER_DIRECTORY]/bin:$PATH
```

Now change directories and test the installation:

```
cd
kover --version
```

This should print the version of the command line interface (cli-x.y.z) and the core (core-x.y.z).

Thats it! You're done.

## Windows

Windows is currently not supported.
