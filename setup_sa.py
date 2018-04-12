#!/usr/bin/env python2

"""
setup.py file for SWIG sa module
"""

from distutils.core import setup, Extension
import distutils.sysconfig

# Remove the "-Wstrict-prototypes" compiler option so we don't have to add
# void when there's no parameters
cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")

sa_module = Extension('_sa',
                      sources=['sa_wrap.c', 'sds.c', 'coordset.c'],
                      depends=['sa.i', 'sa.c', 'sds.c', 'coordset.c', 'coordset.h', 'Makefile'],
                      include_dirs=['.'],
                      )
setup(name='sa',
      version='0.1',
      author="Ben Kovitz and Dave Bender",
      description="""Genetic algorithm using spreading activation""",
      ext_modules=[sa_module],
      py_modules=["sa"],
      )
