#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from setuptools import setup
from glob import glob

setup(
    name="goatools",
    version='0.4.7',
    author='Haibao Tang',
    author_email='tanghaibao@gmail.com',
    packages=['goatools'],
    scripts=glob('scripts/*.py'),
    license='LICENSE',
    url='http://pypi.python.org/pypi/goatools/',
    description="Python scripts to find enrichment of GO terms",
    long_description=open("README.rst").read(),
    install_requires=['fisher', 'pygraphviz']
    )
