#!/usr/bin/env python2

"""
setup.py file for SWIG sa module
"""

from distutils.core import setup, Extension

sa_module = Extension('_sa',
                      sources=['sa_wrap.c'],
                      depends=['../sa.c'],
                      include_dirs=['.', '..'],
                      )
setup(name='sa',
      version='0.1',
      author="Ben Kovitz and Dave Bender",
      description="""Spreading Activation genetic algorithm library""",
      ext_modules=[sa_module],
      py_modules=["sa"],
      )
