#/usr/bin/env python

from distutils.core import Extension
from numpy.distutils.misc_util import Configuration
from numpy.distutils.system_info import get_info
import numpy
import os
import sys

config = Configuration('pandas', parent_package=None, top_path=None)

cython_ext = Extension('pandas.lib.tseries', ['pandas/lib/src/tseries.c'], 
                       include_dirs=[numpy.get_include(),
                                     'pandas/lib/include/']) 

dates_ext = Extension('pandas.lib.tdates', ['pandas/lib/src/tdates.c'])

config_dict = config.todict()
try:
    config_dict.pop('packages')
except:
    pass

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(version="0.1",
          description="Panel and time series data analysis toolkit",
          author="AQR Capital Management, LLC",
          author_email='wesmckinn@gmail.com',
          url="pandas.googlecode.com",
          license="BSD License",
          classifiers=[
            'Development Status :: 2 - Pre-Alpha',
            'Environment :: Console',
            'Operating System :: OS Independent',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python',
            'Programming Language :: Cython',
            'Topic :: Scientific/Engineering',
            ],
          requires=['NumPy (>=1.2)',],
          platforms='any',
          long_description="""
          """,
          packages=["pandas", "pandas.core", "pandas.core.tests", "pandas.io",
                    "pandas.lib"],
          ext_modules=[cython_ext, dates_ext],
          **(config_dict))
