#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""Setup for PyPI usage."""

import glob
import os.path as op
import sys

from setuptools import setup
from setuptools.command.test import test as testCommand
from setup_helper import SetupHelper

import versioneer


class PyTest(testCommand):
    """Allow testing to be run from setuptools."""

    def initialize_options(self):
        testCommand.initialize_options(self)
        self.test_args = []

    def finalize_options(self):
        testCommand.finalize_options(self)
        self.test_args += ["--cov", "goatools", "tests"]

    def run_tests(self):
        # pylint:disable=import-outside-toplevel
        import pytest

        errno = pytest.main(self.test_args)
        sys.exit(errno)


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

PACKAGES = [
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
    NAME + ".semsim",
    NAME + ".semsim.termwise",
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
    packages=PACKAGES,
    include_package_data=True,
    package_data={"goatools.test_data.nbt_3102": ["*.*"]},
    scripts=glob.glob("scripts/*.py"),
    classifiers=CLASSIFIERS,
    url="http://github.com/tanghaibao/goatools",
    description="Python scripts to find enrichment of GO terms",
    install_requires=REQUIREMENTS,
    tests_require=["pytest", "pytest-cov", "nose"],
)
