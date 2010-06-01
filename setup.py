#/usr/bin/env python

from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup

DESCRIPTION = "Cross-section and time series data analysis toolkit"
LONG_DESCRIPTION = """
pandas provides NumPy-based data structures and statistical tools for
common time series and cross-sectional data sets. It is intended to
accomplish the following:

* Simplify working with possibly labeled 1, 2, and 3 dimensional
  heterogeneous data sets commonly found in statistics, finance, and
  econometrics.

* Provide IO utilities for getting data in and out of pandas

* Implement common statistical models with a convenient interface,
  handling missing data and other common problems associated with
  messy statistical data sets

Note
----
Windows binaries built against NumPy 1.3.0
"""

DISTNAME = 'pandas'
LICENSE = 'BSD'
MAINTAINER = "AQR Capital Management, LLC"
MAINTAINER_EMAIL = "wesmckinn@gmail.com"
URL = "http://pandas.sourceforge.net"
DOWNLOAD_URL = ''
CLASSIFIERS = [
    'Development Status :: 4 - Beta',
    'Environment :: Console',
    'Operating System :: OS Independent',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python',
    'Programming Language :: Cython',
    'Topic :: Scientific/Engineering',
]

MAJOR = 0
MINOR = 2

def get_version():
    return '%s.%s' % (MAJOR, MINOR)


def configuration(parent_package='', top_path=None):
    config = Configuration(None, parent_package, top_path,
                           version=get_version())
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)

    config.add_subpackage('pandas')
    return config

if __name__ == '__main__':
    setup(name=DISTNAME,
          maintainer=MAINTAINER,
          maintainer_email=MAINTAINER_EMAIL,
          description=DESCRIPTION,
          license=LICENSE,
          url=URL,
          download_url=DOWNLOAD_URL,
          long_description=LONG_DESCRIPTION,
          classifiers=CLASSIFIERS,
          platforms='any',
          configuration=configuration)
