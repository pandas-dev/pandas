#/usr/bin/env python

from distutils.extension import Extension
from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup
import numpy

config = Configuration('pandas', parent_package=None, top_path=None)

from Cython.Distutils import build_ext


pyx_ext = Extension('tseries', ['src/tseries.pyx',
                                'src/wirth.c'],
                    include_dirs=[numpy.get_include(),
                                  'include/'])


dates_ext = Extension('tdates', ['src/tdates.c'],
                      include_dirs=[numpy.get_include()])

setup(name='tdates', description='Nothing', ext_modules=[dates_ext])

setup(name='pandas.lib.tseries', description='Nothing',
      ext_modules=[pyx_ext],
      cmdclass = {
          'build_ext' : build_ext
      })
