#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""Setup for PyPI usage."""

import os.path as op
from glob import glob
from setuptools import setup
from setup_helper import SetupHelper

import versioneer

NAME = "goatools"
CLASSIFIERS = [
    "Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Programming Language :: Python",
    "Programming Language :: Python :: 2",
    "Programming Language :: Python :: 3",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

# Use the helper
HLPR = SetupHelper(initfile="goatools/__init__.py", readmefile="README.md")

SETUP_DIR = op.abspath(op.dirname(__file__))
REQUIREMENTS = [
    x.strip() for x in open(op.join(SETUP_DIR, "requirements.txt")).readlines()
]

setup(
    name=NAME,
    version=versioneer.get_version(),
    author=HLPR.author,
    author_email=HLPR.email,
    license=HLPR.license,
    long_description=HLPR.long_description,
    long_description_content_type="text/markdown",
    cmdclass=versioneer.get_cmdclass(),
    packages=[
        NAME,
        NAME + ".godag",
        NAME + ".gosubdag",
        NAME + ".gosubdag.plot",
        NAME + ".gosubdag.rpt",
        NAME + ".test_data",
        NAME + ".test_data.sections",
        NAME + ".test_data.cli",
        NAME + ".cli",
        NAME + ".rpt",
        NAME + ".anno",
        NAME + ".anno.init",
        NAME + ".anno.extensions",
        NAME + ".goea",
        NAME + ".grouper",
        NAME + ".parsers",
    ],
    include_package_data=True,
    package_data={"goatools.test_data.nbt_3102": ["*.*"]},
    scripts=glob("scripts/*.py"),
    classifiers=CLASSIFIERS,
    url="http://github.com/tanghaibao/goatools",
    description="Python scripts to find enrichment of GO terms",
    install_requires=REQUIREMENTS,
)
