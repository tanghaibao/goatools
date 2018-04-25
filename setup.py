#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os.path as op

from glob import glob
from setuptools import setup
from setup_helper import SetupHelper


name = "goatools"
classifiers = [
    'Development Status :: 4 - Beta',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Programming Language :: Python',
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 3',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    ]

# Use the helper
h = SetupHelper(initfile="goatools/__init__.py", readmefile="README.md")

setup_dir = op.abspath(op.dirname(__file__))
requirements = ['wget'] + [x.strip() for x in
                           open(op.join(setup_dir, 'requirements.txt')).readlines()]

setup(
    name=name,
    version=h.version,
    author=h.author,
    author_email=h.email,
    license=h.license,
    long_description=h.long_description,
    packages=[
        name,
        name + ".godag",
        name + ".test_data",
        name + ".cli",
        name + ".anno",
        name + ".anno.extensions",
        name + ".parsers"],
    include_package_data=True,
    package_data={"goatools.test_data.nbt_3102": ["*.*"]},
    scripts=glob('scripts/*.py'),
    classifiers=classifiers,
    url='http://github.com/tanghaibao/goatools',
    description="Python scripts to find enrichment of GO terms",
    install_requires=requirements
)
