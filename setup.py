#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from setuptools import setup
from glob import glob

classifiers = [
    'Development Status :: 4 - Beta',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Programming Language :: Python',
    'Programming Language :: Python :: 2',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    ]

setup(
    name="goatools",
    version='0.4.7',
    author='Haibao Tang',
    author_email='tanghaibao@gmail.com',
    packages=['goatools'],
    scripts=glob('scripts/*.py'),
    license='BSD',
    classifiers=classifiers,
    url='http://github.com/tanghaibao/goatools',
    description="Python scripts to find enrichment of GO terms",
    long_description=open("README.rst").read(),
    install_requires=['fisher', 'pygraphviz']
    )
