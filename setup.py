#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from distutils.core import setup, Extension

fisher_module = Extension('goatools/_fisher', ['src/fisher.c','src/fisher_wrap.c'])

setup(
      name="goatools",
      ext_modules=[fisher_module],
      packages=['goatools'],
      scripts = ['goatools/go_enrichment.py'],
      )
