from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext as _build_ext

class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)
        # Prevent numpy from thinking it is still in its setup process:
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())

setup(
    name = "kover",
    version = "1.1.3",
    packages = find_packages(),

    cmdclass={'build_ext':build_ext},
    setup_requires = ['numpy'],
    install_requires = ['h5py>=2.4.0', 'numpy', 'pandas', 'progressbar', 'scipy'],

    author = "Alexandre Drouin",
    author_email = "aldro61@gmail.com",
    description = "Kover: Learn computational phenotyping models from k-merized genomic data",
    license = "GPLv3",
    keywords = "genomics, machine learning, set covering machine, kmer",
    url = "http://github.com/aldro61/kover",
    
    # Cython Extension
    ext_modules = [Extension("kover/learning/set_covering_machine/popcount", ["kover/learning/set_covering_machine/popcount.c"], extra_compile_args=["-march=native"])],
    data_files=[('kover/dataset/tools/contigs_tools/', ['kover/dataset/tools/contigs_tools/build/bin/dsk2kover', 
																  'kover/dataset/tools/contigs_tools/build/bin/multidsk', 
																  'kover/dataset/tools/contigs_tools/build/bin/dsk'])]
)
