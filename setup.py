from setuptools import setup, find_packages, Extension  # type: ignore
from Cython.Build import cythonize  # type: ignore
import numpy as np

# Define Cython extensions
extensions = [
    Extension(
        "taggd.core.demultiplex_sub_functions",
        ["taggd/core/demultiplex_sub_functions.pyx"],
    ),
    Extension(
        "taggd.core.demultiplex_search_functions",
        ["taggd/core/demultiplex_search_functions.pyx"],
    ),
    Extension("taggd.misc.distance_metrics", ["taggd/misc/distance_metrics.pyx"]),
    Extension("taggd.misc.kmer_utils", ["taggd/misc/kmer_utils.pyx"]),
    Extension("taggd.io.barcode_utils", ["taggd/io/barcode_utils.pyx"]),
]

setup(
    name="taggd",
    version="0.4.0",
    author="Joel Sjostrand, Jose Fernandez",
    author_email="joel.sjostrand@scilifelab.se, jc.fernandez.navarro@gmail.com",
    license="BSD-3-Clause",
    description="Bioinformatics genetic barcode demultiplexing",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/jfnavarro/taggd",
    download_url="https://github.com/jfnavarro/taggd/0.4.0",
    packages=find_packages(include=["scripts", "scripts.*", "taggd", "taggd.*"]),
    entry_points={
        "console_scripts": [
            "taggd_demultiplex=scripts.taggd_demultiplex:main",
        ],
    },
    setup_requires=["cython", "numpy"],
    install_requires=["setuptools", "pysam", "numpy", "dnaio"],
    python_requires=">=3.10",
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
