try:
    from numpy import get_include as get_numpy_include
except:
    print "Numpy 1.6.2 or greater is required to install this package. Please install Numpy and try again."
    exit()

from setuptools import setup, find_packages, Extension
setup(
    name = "kover",
    version = "0.1",
    packages = find_packages(),

    install_requires = ['h5py', 'numpy', 'pandas', 'progressbar'],

    author = "Alexandre Drouin",
    author_email = "aldro61@gmail.com",
    description = "Kover: Learn computational phenotyping models from k-merized genomic data",
    license = "GPLv3",
    keywords = "genomics, machine learning, set covering machine, kmer",
    url = "http://github.com/aldro61/kover",

    # Cython Extension
    ext_modules = [Extension("kover/core/learning/set_covering_machine/popcount", ["kover/core/learning/set_covering_machine/popcount.c"], include_dirs=[get_numpy_include()], extra_compile_args=["-march=native"])]
)
