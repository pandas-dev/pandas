#/usr/bin/env python

import os
import re
import subprocess
import warnings

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

* Provide tools for working with dates, fixed-frequency custom date
  ranges, and other necessary time series functionality

* Provide IO utilities for getting data in and out of pandas

* Implement common statistical and time series functionality with a
  convenient interface, handling missing data and other common
  problems associated with messy statistical data sets

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
MICRO = 1
ISRELEASED = False
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)

# klau'd from SciPy's setup.py

def svn_version():
    env = os.environ.copy()
    env['LC_ALL'] = 'C'
    try:
        out = subprocess.Popen(['svn', 'info'], stdout=subprocess.PIPE,
                env=env).communicate()[0]
    except OSError:
        warnings.warn(" --- Could not run svn info --- ")
        return ""

    r = re.compile('Revision: ([0-9]+)')
    svnver = None
    for line in out.split('\n'):
        m = r.match(line)
        if m:
            svnver = m.group(1)

    if not svnver:
        raise ValueError("Error while parsing svn version ?")
    return svnver

FULLVERSION = VERSION
if not ISRELEASED:
    FULLVERSION += '.dev'
    # If in git or something, bypass the svn rev
    if os.path.exists('.svn'):
        FULLVERSION += svn_version()

def write_version_py(filename='pandas/version.py'):
    cnt = """\
# THIS FILE IS GENERATED FROM PANDAS setup.py
short_version='%(version)s'
version='%(version)s'
release=%(isrelease)s

if not release:
    version += '.dev'
    import os
    svn_version_file = os.path.join(os.path.dirname(__file__),
                                    '__svn_version__.py')
    if os.path.isfile(svn_version_file):
        import imp
        svn = imp.load_module('pandas.__svn_version__',
                              open(svn_version_file),
                              svn_version_file,
                              ('.py','U',1))
        version += svn.version
"""
    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': VERSION, 'isrelease': str(ISRELEASED)})
    finally:
        a.close()

def configuration(parent_package='', top_path=None):
    write_version_py()

    config = Configuration(None, parent_package, top_path,
                           version=FULLVERSION)
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
