#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from distutils.core import setup
from glob import glob

setup(
    name="goatools",
    packages=['goatools'],
    scripts=glob('scripts/*.py'),
    requires=['fisher', 'pygraphviz']
    )
