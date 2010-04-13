#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from distutils.core import setup #, Extension


setup(
      name="goatools",
    #ext_modules=[fisher_module],
      packages=['goatools'],
      scripts = ['scripts/find_enrichment.py'],
      requires=['fisher'],
      )
