#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from distutils.core import setup, Extension

fisher_module = Extension('_fisher', ['fisher.c','fisher_wrap.c'])

setup(ext_modules=[fisher_module],
      py_modules=['fisher'],
      name="goatools",
      )
