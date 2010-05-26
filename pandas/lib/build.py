#/usr/bin/env python

from distutils.extension import Extension
from numpy.distutils.core import setup
import numpy
from Cython.Distutils import build_ext

pyx_ext = Extension('tseries', ['src/tseries.pyx'],
                    include_dirs=[numpy.get_include()])

setup(name='pandas.lib.tseries', description='Nothing',
      ext_modules=[pyx_ext],
      cmdclass = {
          'build_ext' : build_ext
      })
