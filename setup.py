from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
import numpy as np
from glob import glob

# Define Cython extensions
extensions = [
    Extension(
        "taggd.core.demultiplex_core_functions",
        ["taggd/core/demultiplex_core_functions.pyx"],
    ),
    Extension(
        "taggd.core.demultiplex_sub_functions",
        ["taggd/core/demultiplex_sub_functions.pyx"],
    ),
    Extension(
        "taggd.core.demultiplex_search_functions",
        ["taggd/core/demultiplex_search_functions.pyx"],
    ),
    Extension("taggd.core.match", ["taggd/core/match.pyx"]),
    Extension("taggd.core.match_type", ["taggd/core/match_type.pyx"]),
    Extension("taggd.core.statistics", ["taggd/core/statistics.pyx"]),
    Extension("taggd.misc.distance_metrics", ["taggd/misc/distance_metrics.pyx"]),
    Extension("taggd.misc.kmer_utils", ["taggd/misc/kmer_utils.pyx"]),
    Extension("taggd.io.fastq_utils", ["taggd/io/fastq_utils.pyx"]),
    Extension("taggd.io.barcode_utils", ["taggd/io/barcode_utils.pyx"]),
    Extension("taggd.io.record", ["taggd/io/record.pyx"]),
    Extension("taggd.io.sam_record", ["taggd/io/sam_record.pyx"]),
    Extension("taggd.io.fasta_record", ["taggd/io/fasta_record.pyx"]),
    Extension("taggd.io.fastq_record", ["taggd/io/fastq_record.pyx"]),
    Extension("taggd.io.reads_reader_writer", ["taggd/io/reads_reader_writer.pyx"]),
]

setup(
    name="taggd",
    version="0.3.7",
    author="Joel Sjostrand, Jose Fernandez",
    author_email="joel.sjostrand@scilifelab.se, jc.fernandez.navarro@gmail.com",
    license="BSD-3-Clause",
    description="Bioinformatics genetic barcode demultiplexing",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/jfnavarro/taggd",
    download_url="https://github.com/jfnavarro/taggd/0.3.7",
    scripts=glob("scripts/*.py"),
    packages=find_packages(include=["taggd", "taggd.*"]),
    setup_requires=["cython", "numpy"],
    install_requires=[
        "setuptools",
        "pysam",
        "numpy",
    ],
    python_requires=">=3.6",
    classifiers=[
        "Development Status :: 5 - Stable",
        "Programming Language :: Python :: 3.10"
        "Programming Language :: Python :: 3.11"
        "Programming Language :: Python :: 3.12"
        "License :: OSI Approved :: BSD License",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    test_suite="tests",
    ext_modules=cythonize(extensions, language_level="3"),
    include_dirs=[np.get_include()],
    keywords=["bioinformatics", "demultiplexing"],
)
