#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from setuptools import setup
from glob import glob
import os

setup_dir = os.path.abspath(os.path.dirname(__file__))
classifiers = [
    'Development Status :: 4 - Beta',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Programming Language :: Python',
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 3',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    ]

requirements = ['wget'] + [x.strip() for x in
                           open(os.path.join(setup_dir, 'requirements.txt')).readlines()]

exec(open(os.path.join(setup_dir, "goatools", "version.py")).read())
setup(
    name="goatools",
    version=__version__,
    author='Haibao Tang',
    author_email='tanghaibao@gmail.com',
    packages=['goatools'],
    scripts=glob('scripts/*.py'),
    license='BSD',
    classifiers=classifiers,
    url='http://github.com/tanghaibao/goatools',
    description="Python scripts to find enrichment of GO terms",
    long_description=open(os.path.join(setup_dir, "README.md")).read(),
    install_requires=requirements
    )
