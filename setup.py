#/usr/bin/env python

from distutils.core import Extension

from numpy.distutils.misc_util import Configuration
import setuptools
from numpy.distutils.core import setup
import numpy

DESCRIPTION = "Cross-section and time series data analysis toolkit"
LONG_DESCRIPTION = """
Pandas provides data structures and statistical tools for common
time-series and cross-sectional data sets.
"""

DISTNAME = 'pandas'
LICENSE = 'BSD'
MAINTAINER = "AQR Capital Management, LLC"
MAINTAINER_EMAIL = "wesmckinn@gmail.com"
URL = "pandas.googlecode.com"
DOWNLOAD_URL = ''
CLASSIFIERS = [
    'Development Status :: 3 - Alpha',
    'Environment :: Console',
    'Operating System :: OS Independent',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python',
    'Programming Language :: Cython',
    'Topic :: Scientific/Engineering',
]

MAJOR = 0
MINOR = 1

def get_version():
    return '%d.%d' % (MAJOR, MINOR)

def get_cython_ext():
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

def configuration(parent_package='', top_path=None, package_name=DISTNAME):
    config = Configuration(None, parent_package, top_path,
                           name=DISTNAME,
                           version=get_version(),
                           maintainer =MAINTAINER,
                           maintainer_email=MAINTAINER_EMAIL,
                           description=DESCRIPTION,
                           license=LICENSE,
                           url=URL,
                           download_url=DOWNLOAD_URL,
                           long_description=LONG_DESCRIPTION)

    config.add_extension('lib.tseries',
                         sources=['pandas/lib/src/tseries.c',
                                  'pandas/lib/src/wirth.c'],
                         include_dirs=[numpy.get_include(),
                                       'pandas/lib/include/'])
    config.add_extension('lib.tdates',
                         sources=['pandas/lib/src/tdates.c'])

    config.set_options(
            ignore_setup_xxx_py=True,
            assume_default_configuration=True,
            delegate_options_to_subpackages=True,
            quiet=False,
            )

    return config

if __name__ == '__main__':
    setup(configuration=configuration,
          packages=setuptools.find_packages(),
          classifiers=CLASSIFIERS,
          requires=['numpy', 'scikits.statsmodels', 'dateutil'],
          platforms='any',
          test_suite='nose.collector',
          zip_safe=False)
