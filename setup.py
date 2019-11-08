#!/usr/bin/env python3

from setuptools import find_packages, setup
import sys


if sys.version_info.major < 3:
    raise SystemExit("Python 3 is required!")

setup(
    name="nddia",
    version=0.01,
    description="N/D diabatization in python.",
    url="https://github.com/eljost/pysisyphus",
    maintainer="Johannes Steinmetzer",
    maintainer_email="johannes.steinmetzer@uni-jena.de",
    license="GPL 3",
    platforms=["unix"],
    packages=find_packages(),
    install_requires=[
        "numpy",
    ],
)
