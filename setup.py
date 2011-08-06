#/usr/bin/env python

"""
Parts of this file were taken from the pyzmq project
(https://github.com/zeromq/pyzmq) and hence are subject to the terms of the
Lesser GPU General Public License.
"""
# use setuptools if available
# try:
#     from setuptools import setup
#     _have_setuptools = True
# except ImportError:
#     _have_setuptools = False

from datetime import datetime
from glob import glob
import os
import sys
import shutil

import numpy as np

# from numpy.distutils.core import setup

from distutils.core import setup, Command
from distutils.extension import Extension
from distutils.command.build import build
from distutils.command.build_ext import build_ext
from distutils.command.sdist import sdist

from os.path import splitext, basename, join as pjoin

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
Windows binaries built against NumPy 1.5.1
"""

DISTNAME = 'pandas'
LICENSE = 'BSD'
AUTHOR = "AQR Capital Management, LLC"
MAINTAINER = "Wes McKinney"
MAINTAINER_EMAIL = "wesmckinn@gmail.com"
URL = "http://github.com/wesm/pandas"
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
MINOR = 4
MICRO = 0
ISRELEASED = False
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)

FULLVERSION = VERSION
if not ISRELEASED:
    FULLVERSION += '.dev'
    try:
        import subprocess
        pipe = subprocess.Popen(["git", "rev-parse", "--short", "HEAD"],
                                stdout=subprocess.PIPE).stdout
        rev = pipe.read().strip()
        FULLVERSION += "-%s" % rev
    except:
        print "WARNING: Couldn't get git revision"

def write_version_py(filename='pandas/version.py'):
    cnt = """\
from datetime import datetime

version = '%s'
"""
    a = open(filename, 'w')
    try:
        a.write(cnt % FULLVERSION)
    finally:
        a.close()

class CleanCommand(Command):
    """Custom distutils command to clean the .so and .pyc files."""

    user_options = [("all", "a", "") ]

    def initialize_options(self):
        self.all = True
        self._clean_me = []
        self._clean_trees = []
        for root, dirs, files in list(os.walk('pandas')):
            for f in files:
                if os.path.splitext(f)[-1] in ('.pyc', '.so', '.o', '.pyd'):
                    self._clean_me.append(pjoin(root, f))
            for d in dirs:
                if d == '__pycache__':
                    self._clean_trees.append(pjoin(root, d))

        for d in ('build',):
            if os.path.exists(d):
                self._clean_trees.append(d)

    def finalize_options(self):
        pass

    def run(self):
        for clean_me in self._clean_me:
            try:
                os.unlink(clean_me)
            except Exception:
                pass
        for clean_tree in self._clean_trees:
            try:
                shutil.rmtree(clean_tree)
            except Exception:
                pass

class CheckSDist(sdist):
    """Custom sdist that ensures Cython has compiled all pyx files to c."""

    _pyxfiles = ['pandas/src/tseries.pyx'
                 'pandas/src/sparse.pyx']

    def initialize_options(self):
        sdist.initialize_options(self)

        '''
        self._pyxfiles = []
        for root, dirs, files in os.walk('pandas'):
            for f in files:
                if f.endswith('.pyx'):
                    self._pyxfiles.append(pjoin(root, f))
        '''

    def run(self):
        if 'cython' in cmdclass:
            self.run_command('cython')
        else:
            for pyxfile in self._pyxfiles:
                cfile = pyxfile[:-3]+'c'
                msg = "C-source file '%s' not found."%(cfile)+\
                " Run 'setup.py cython' before sdist."
                assert os.path.isfile(cfile), msg
        sdist.run(self)

class CheckingBuildExt(build_ext):
    """Subclass build_ext to get clearer report if Cython is neccessary."""

    def check_cython_extensions(self, extensions):
        for ext in extensions:
          for src in ext.sources:
            if not os.path.exists(src):
                print """Cython-generated file '%s' not found.
                Cython is required to compile pandas from a development branch.
                Please install Cython or download a release package of pandas.
                """ % src

    def build_extensions(self):
        self.check_cython_extensions(self.extensions)
        self.check_extensions_list(self.extensions)

        for ext in self.extensions:
            self.build_extension(ext)

cmdclass = {'clean': CleanCommand,
            'build': build}

try:
    from Cython.Distutils import build_ext
    cython=True
except ImportError:
    cython=False
    suffix = '.c'
    cmdclass['build_ext'] = CheckingBuildExt
else:
    suffix = '.pyx'
    class CythonCommand(build_ext):
        """Custom distutils command subclassed from Cython.Distutils.build_ext
        to compile pyx->c, and stop there. All this does is override the
        C-compile method build_extension() with a no-op."""
        def build_extension(self, ext):
            pass

    class DummyBuildSrc(Command):
        """ numpy's build_src command interferes with Cython's build_ext.
        """
        user_options = []
        def initialize_options(self):
            self.py_modules_dict = {}
        def finalize_options(self):
            pass
        def run(self):
            pass

    cmdclass['build_src'] = DummyBuildSrc
    cmdclass['cython'] = CythonCommand
    cmdclass['build_ext'] =  build_ext
    cmdclass['sdist'] =  CheckSDist

tseries_depends = ['reindex', 'io', 'common', 'groupby'
                   'skiplist', 'isnull', 'moments', 'operators']

def srcpath(name=None, suffix='.pyx', subdir='src'):
    return pjoin('pandas', subdir, name+suffix)

tseries_ext = Extension('pandas._tseries',
                        sources=[srcpath('tseries', suffix=suffix)],
                        # depends=[srcpath(f, suffix='.pyx')
                        #          for f in tseries_depends],
                        include_dirs=[np.get_include()])
sparse_ext = Extension('pandas._sparse',
                       sources=[srcpath('sparse', suffix=suffix)],
                       include_dirs=[np.get_include()])
extensions = [tseries_ext,
              sparse_ext]

setuptools_args = {}

# if _have_setuptools:
#     setuptools_args["test_suite"] = "nose.collector"

setup(name=DISTNAME,
      version=FULLVERSION,
      maintainer=MAINTAINER,
      packages=['pandas',
                'pandas.core',
                'pandas.io',
                'pandas.stats',
                'pandas.util'],
      package_data={'pandas' : ['tests/*.py'],
                    'pandas.io' : ['tests/*.py',
                                   'tests/*.h5',
                                   'tests/*.csv',
                                   'tests/*.xls'],
                    'pandas.stats' : ['tests/*.py']},
      ext_modules=extensions,
      maintainer_email=MAINTAINER_EMAIL,
      description=DESCRIPTION,
      license=LICENSE,
      cmdclass = cmdclass,
      url=URL,
      download_url=DOWNLOAD_URL,
      long_description=LONG_DESCRIPTION,
      classifiers=CLASSIFIERS,
      platforms='any',
      **setuptools_args)
