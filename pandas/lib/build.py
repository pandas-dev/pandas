#/usr/bin/env python

from distutils.extension import Extension
from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup
import numpy

config = Configuration('pandas', parent_package=None, top_path=None)

from Cython.Distutils import build_ext


pyx_ext = Extension('tseries', ['pandas/lib/src/tseries.pyx',
                                'pandas/lib/src/wirth.c'],
                    include_dirs=[numpy.get_include(),
                                  'pandas/lib/include/'])


setup(name='pandas.lib.tseries', description='Nothing',
      ext_modules=[pyx_ext],
      cmdclass = {
          'build_ext' : build_ext
      })
