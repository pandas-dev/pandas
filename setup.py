#/usr/bin/env python

from distutils.core import Extension
from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup
import numpy

config = Configuration('pandas', parent_package=None, top_path=None)

def get_cython_ext():
    from Cython.Distutils import build_ext
    from distutils.extension import Extension

    pyx_ext = Extension('tseries', ['pandas/lib/src/tseries.pyx',
                                    'pandas/lib/src/wirth.c'],
                        include_dirs=[numpy.get_include(),
                                      'pandas/lib/include/'])


    setup(name='pandas.lib.tseries', description='Nothing',
          ext_modules=[pyx_ext],
          cmdclass = {
              'build_ext' : build_ext
          })

# get_cython_ext()
# sys.exit()

cython_ext = Extension('pandas.lib.tseries', ['pandas/lib/src/tseries.c',
                                              'pandas/lib/src/wirth.c'],
                       include_dirs=[numpy.get_include(),
                                     'pandas/lib/include/'])

dates_ext = Extension('pandas.lib.tdates', ['pandas/lib/src/tdates.c'])

config_dict = config.todict()
try:
    config_dict.pop('packages')
except:
    pass

if __name__ == '__main__':
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
          packages=["pandas",
                    "pandas.core",
                    "pandas.core.tests",
                    "pandas.io",
                    "pandas.lib",
                    "pandas.stats"],
          ext_modules=[cython_ext, dates_ext],
          **(config_dict))
