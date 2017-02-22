#!/usr/bin/env python

"""
Parts of this file were taken from the pyzmq project
(https://github.com/zeromq/pyzmq) which have been permitted for use under the
BSD license. Parts are from lxml (https://github.com/lxml/lxml)
"""

import os
import sys
import shutil
import warnings
import re
import platform
from distutils.version import LooseVersion

def is_platform_windows():
    return sys.platform == 'win32' or sys.platform == 'cygwin'

def is_platform_linux():
    return sys.platform == 'linux2'

def is_platform_mac():
    return sys.platform == 'darwin'

# versioning
import versioneer
cmdclass = versioneer.get_cmdclass()

min_cython_ver = '0.23'
try:
    import Cython
    ver = Cython.__version__
    _CYTHON_INSTALLED = ver >= LooseVersion(min_cython_ver)
except ImportError:
    _CYTHON_INSTALLED = False

try:
    import pkg_resources
    from setuptools import setup, Command
    _have_setuptools = True
except ImportError:
    # no setuptools installed
    from distutils.core import setup, Command
    _have_setuptools = False

setuptools_kwargs = {}
min_numpy_ver = '1.7.0'
if sys.version_info[0] >= 3:

    setuptools_kwargs = {
                         'zip_safe': False,
                         'install_requires': ['python-dateutil >= 2',
                                              'pytz >= 2011k',
                                              'numpy >= %s' % min_numpy_ver],
                         'setup_requires': ['numpy >= %s' % min_numpy_ver],
                         }
    if not _have_setuptools:
        sys.exit("need setuptools/distribute for Py3k"
                 "\n$ pip install distribute")

else:
    setuptools_kwargs = {
        'install_requires': ['python-dateutil',
                            'pytz >= 2011k',
                             'numpy >= %s' % min_numpy_ver],
        'setup_requires': ['numpy >= %s' % min_numpy_ver],
        'zip_safe': False,
    }

    if not _have_setuptools:
        try:
            import numpy
            import dateutil
            setuptools_kwargs = {}
        except ImportError:
            sys.exit("install requires: 'python-dateutil < 2','numpy'."
                     "  use pip or easy_install."
                     "\n   $ pip install 'python-dateutil < 2' 'numpy'")

from distutils.extension import Extension
from distutils.command.build import build
from distutils.command.build_ext import build_ext as _build_ext

try:
    if not _CYTHON_INSTALLED:
        raise ImportError('No supported version of Cython installed.')
    try:
        from Cython.Distutils.old_build_ext import old_build_ext as _build_ext
    except ImportError:
        # Pre 0.25
        from Cython.Distutils import build_ext as _build_ext
    cython = True
except ImportError:
    cython = False


if cython:
    try:
        try:
            from Cython import Tempita as tempita
        except ImportError:
            import tempita
    except ImportError:
        raise ImportError('Building pandas requires Tempita: '
                          'pip install Tempita')


from os.path import join as pjoin


_pxipath = pjoin('pandas', 'src')
_pxi_dep_template = {
    'algos': ['algos_common_helper.pxi.in', 'algos_groupby_helper.pxi.in',
              'algos_take_helper.pxi.in', 'algos_rank_helper.pxi.in'],
    '_join': ['join_helper.pxi.in', 'joins_func_helper.pxi.in'],
    'hashtable': ['hashtable_class_helper.pxi.in',
                  'hashtable_func_helper.pxi.in'],
    'index': ['index_class_helper.pxi.in'],
    '_sparse': ['sparse_op_helper.pxi.in']
}
_pxifiles = []
_pxi_dep = {}
for module, files in _pxi_dep_template.items():
    pxi_files = [pjoin(_pxipath, x) for x in files]
    _pxifiles.extend(pxi_files)
    _pxi_dep[module] = pxi_files


class build_ext(_build_ext):
    def build_extensions(self):

        # if builing from c files, don't need to
        # generate template output
        if cython:
            for pxifile in _pxifiles:
                # build pxifiles first, template extention must be .pxi.in
                assert pxifile.endswith('.pxi.in')
                outfile = pxifile[:-3]

                if (os.path.exists(outfile) and
                    os.stat(pxifile).st_mtime < os.stat(outfile).st_mtime):
                    # if .pxi.in is not updated, no need to output .pxi
                    continue

                with open(pxifile, "r") as f:
                    tmpl = f.read()
                pyxcontent = tempita.sub(tmpl)

                with open(outfile, "w") as f:
                    f.write(pyxcontent)

        numpy_incl = pkg_resources.resource_filename('numpy', 'core/include')

        for ext in self.extensions:
            if hasattr(ext, 'include_dirs') and not numpy_incl in ext.include_dirs:
                ext.include_dirs.append(numpy_incl)
        _build_ext.build_extensions(self)


DESCRIPTION = ("Powerful data structures for data analysis, time series,"
               "and statistics")
LONG_DESCRIPTION = """
**pandas** is a Python package providing fast, flexible, and expressive data
structures designed to make working with structured (tabular, multidimensional,
potentially heterogeneous) and time series data both easy and intuitive. It
aims to be the fundamental high-level building block for doing practical,
**real world** data analysis in Python. Additionally, it has the broader goal
of becoming **the most powerful and flexible open source data analysis /
manipulation tool available in any language**. It is already well on its way
toward this goal.

pandas is well suited for many different kinds of data:

  - Tabular data with heterogeneously-typed columns, as in an SQL table or
    Excel spreadsheet
  - Ordered and unordered (not necessarily fixed-frequency) time series data.
  - Arbitrary matrix data (homogeneously typed or heterogeneous) with row and
    column labels
  - Any other form of observational / statistical data sets. The data actually
    need not be labeled at all to be placed into a pandas data structure

The two primary data structures of pandas, Series (1-dimensional) and DataFrame
(2-dimensional), handle the vast majority of typical use cases in finance,
statistics, social science, and many areas of engineering. For R users,
DataFrame provides everything that R's ``data.frame`` provides and much
more. pandas is built on top of `NumPy <http://www.numpy.org>`__ and is
intended to integrate well within a scientific computing environment with many
other 3rd party libraries.

Here are just a few of the things that pandas does well:

  - Easy handling of **missing data** (represented as NaN) in floating point as
    well as non-floating point data
  - Size mutability: columns can be **inserted and deleted** from DataFrame and
    higher dimensional objects
  - Automatic and explicit **data alignment**: objects can be explicitly
    aligned to a set of labels, or the user can simply ignore the labels and
    let `Series`, `DataFrame`, etc. automatically align the data for you in
    computations
  - Powerful, flexible **group by** functionality to perform
    split-apply-combine operations on data sets, for both aggregating and
    transforming data
  - Make it **easy to convert** ragged, differently-indexed data in other
    Python and NumPy data structures into DataFrame objects
  - Intelligent label-based **slicing**, **fancy indexing**, and **subsetting**
    of large data sets
  - Intuitive **merging** and **joining** data sets
  - Flexible **reshaping** and pivoting of data sets
  - **Hierarchical** labeling of axes (possible to have multiple labels per
    tick)
  - Robust IO tools for loading data from **flat files** (CSV and delimited),
    Excel files, databases, and saving / loading data from the ultrafast **HDF5
    format**
  - **Time series**-specific functionality: date range generation and frequency
    conversion, moving window statistics, moving window linear regressions,
    date shifting and lagging, etc.

Many of these principles are here to address the shortcomings frequently
experienced using other languages / scientific research environments. For data
scientists, working with data is typically divided into multiple stages:
munging and cleaning data, analyzing / modeling it, then organizing the results
of the analysis into a form suitable for plotting or tabular display. pandas is
the ideal tool for all of these tasks.

Note
----
Windows binaries built against NumPy 1.8.1
"""

DISTNAME = 'pandas'
LICENSE = 'BSD'
AUTHOR = "The PyData Development Team"
EMAIL = "pydata@googlegroups.com"
URL = "http://pandas.pydata.org"
DOWNLOAD_URL = ''
CLASSIFIERS = [
    'Development Status :: 5 - Production/Stable',
    'Environment :: Console',
    'Operating System :: OS Independent',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python',
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Cython',
    'Topic :: Scientific/Engineering',
]

class CleanCommand(Command):
    """Custom distutils command to clean the .so and .pyc files."""

    user_options = [("all", "a", "")]

    def initialize_options(self):
        self.all = True
        self._clean_me = []
        self._clean_trees = []

        base = pjoin('pandas','src')
        dt = pjoin(base,'datetime')
        src = base
        util = pjoin('pandas','util')
        parser = pjoin(base,'parser')
        ujson_python = pjoin(base,'ujson','python')
        ujson_lib = pjoin(base,'ujson','lib')
        self._clean_exclude = [pjoin(dt,'np_datetime.c'),
                               pjoin(dt,'np_datetime_strings.c'),
                               pjoin(src,'period_helper.c'),
                               pjoin(parser,'tokenizer.c'),
                               pjoin(parser,'io.c'),
                               pjoin(ujson_python,'ujson.c'),
                               pjoin(ujson_python,'objToJSON.c'),
                               pjoin(ujson_python,'JSONtoObj.c'),
                               pjoin(ujson_lib,'ultrajsonenc.c'),
                               pjoin(ujson_lib,'ultrajsondec.c'),
                               pjoin(util,'move.c'),
                               ]

        for root, dirs, files in os.walk('pandas'):
            for f in files:
                filepath = pjoin(root, f)
                if filepath in self._clean_exclude:
                    continue

                if os.path.splitext(f)[-1] in ('.pyc', '.so', '.o',
                                               '.pyo',
                                               '.pyd', '.c', '.orig'):
                    self._clean_me.append(filepath)
            for d in dirs:
                if d == '__pycache__':
                    self._clean_trees.append(pjoin(root, d))

        # clean the generated pxi files
        for pxifile in _pxifiles:
            pxifile = pxifile.replace(".pxi.in", ".pxi")
            self._clean_me.append(pxifile)

        for d in ('build', 'dist'):
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


# we need to inherit from the versioneer
# class as it encodes the version info
sdist_class = cmdclass['sdist']

class CheckSDist(sdist_class):
    """Custom sdist that ensures Cython has compiled all pyx files to c."""

    _pyxfiles = ['pandas/lib.pyx',
                 'pandas/hashtable.pyx',
                 'pandas/tslib.pyx',
                 'pandas/index.pyx',
                 'pandas/algos.pyx',
                 'pandas/join.pyx',
                 'pandas/window.pyx',
                 'pandas/parser.pyx',
                 'pandas/src/period.pyx',
                 'pandas/src/sparse.pyx',
                 'pandas/src/testing.pyx',
                 'pandas/src/hash.pyx',
                 'pandas/io/sas/saslib.pyx']

    def initialize_options(self):
        sdist_class.initialize_options(self)

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
                cfile = pyxfile[:-3] + 'c'
                msg = "C-source file '%s' not found." % (cfile) +\
                    " Run 'setup.py cython' before sdist."
                assert os.path.isfile(cfile), msg
        sdist_class.run(self)


class CheckingBuildExt(build_ext):
    """
    Subclass build_ext to get clearer report if Cython is necessary.

    """

    def check_cython_extensions(self, extensions):
        for ext in extensions:
            for src in ext.sources:
                if not os.path.exists(src):
                    raise Exception("""Cython-generated file '%s' not found.
                Cython is required to compile pandas from a development branch.
                Please install Cython or download a release package of pandas.
                """ % src)

    def build_extensions(self):
        self.check_cython_extensions(self.extensions)
        build_ext.build_extensions(self)


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

cmdclass.update({'clean': CleanCommand,
                 'build': build})

try:
    from wheel.bdist_wheel import bdist_wheel

    class BdistWheel(bdist_wheel):
        def get_tag(self):
            tag = bdist_wheel.get_tag(self)
            repl = 'macosx_10_6_intel.macosx_10_9_intel.macosx_10_9_x86_64'
            if tag[2] == 'macosx_10_6_intel':
                tag = (tag[0], tag[1], repl)
            return tag
    cmdclass['bdist_wheel'] = BdistWheel
except ImportError:
    pass

if cython:
    suffix = '.pyx'
    cmdclass['build_ext'] = CheckingBuildExt
    cmdclass['cython'] = CythonCommand
else:
    suffix = '.c'
    cmdclass['build_src'] = DummyBuildSrc
    cmdclass['build_ext'] = CheckingBuildExt

lib_depends = ['reduce', 'inference', 'properties']


def srcpath(name=None, suffix='.pyx', subdir='src'):
    return pjoin('pandas', subdir, name + suffix)

if suffix == '.pyx':
    lib_depends = [srcpath(f, suffix='.pyx') for f in lib_depends]
    lib_depends.append('pandas/src/util.pxd')
else:
    lib_depends = []
    plib_depends = []

common_include = ['pandas/src/klib', 'pandas/src']


def pxd(name):
    return os.path.abspath(pjoin('pandas', name + '.pxd'))

# args to ignore warnings
if is_platform_windows():
    extra_compile_args=[]
else:
    extra_compile_args=['-Wno-unused-function']

lib_depends = lib_depends + ['pandas/src/numpy_helper.h',
                             'pandas/src/parse_helper.h']


tseries_depends = ['pandas/src/datetime/np_datetime.h',
                   'pandas/src/datetime/np_datetime_strings.h',
                   'pandas/src/datetime_helper.h',
                   'pandas/src/period_helper.h',
                   'pandas/src/datetime.pxd']


# some linux distros require it
libraries = ['m'] if not is_platform_windows() else []

ext_data = dict(
    lib={'pyxfile': 'lib',
         'pxdfiles': [],
         'depends': lib_depends},
    hashtable={'pyxfile': 'hashtable',
               'pxdfiles': ['hashtable'],
               'depends': (['pandas/src/klib/khash_python.h']
                           + _pxi_dep['hashtable'])},
    tslib={'pyxfile': 'tslib',
           'depends': tseries_depends,
           'sources': ['pandas/src/datetime/np_datetime.c',
                       'pandas/src/datetime/np_datetime_strings.c',
                       'pandas/src/period_helper.c']},
    _period={'pyxfile': 'src/period',
             'depends': tseries_depends,
             'sources': ['pandas/src/datetime/np_datetime.c',
                         'pandas/src/datetime/np_datetime_strings.c',
                         'pandas/src/period_helper.c']},
    index={'pyxfile': 'index',
           'sources': ['pandas/src/datetime/np_datetime.c',
                       'pandas/src/datetime/np_datetime_strings.c'],
           'pxdfiles': ['src/util', 'hashtable'],
           'depends': _pxi_dep['index']},
    algos={'pyxfile': 'algos',
           'pxdfiles': ['src/util', 'hashtable'],
           'depends': _pxi_dep['algos']},
    _join={'pyxfile': 'src/join',
           'pxdfiles': ['src/util', 'hashtable'],
           'depends': _pxi_dep['_join']},
    _window={'pyxfile': 'window',
             'pxdfiles': ['src/skiplist', 'src/util'],
             'depends': ['pandas/src/skiplist.pyx',
                         'pandas/src/skiplist.h']},
    parser={'pyxfile': 'parser',
            'depends': ['pandas/src/parser/tokenizer.h',
                        'pandas/src/parser/io.h',
                        'pandas/src/numpy_helper.h'],
            'sources': ['pandas/src/parser/tokenizer.c',
                        'pandas/src/parser/io.c']},
    _sparse={'pyxfile': 'src/sparse',
             'depends': ([srcpath('sparse', suffix='.pyx')] +
                         _pxi_dep['_sparse'])},
    _testing={'pyxfile': 'src/testing',
              'depends': [srcpath('testing', suffix='.pyx')]},
    _hash={'pyxfile': 'src/hash',
           'depends': [srcpath('hash', suffix='.pyx')]},
)

ext_data["io.sas.saslib"] = {'pyxfile': 'io/sas/saslib'}

extensions = []

for name, data in ext_data.items():
    sources = [srcpath(data['pyxfile'], suffix=suffix, subdir='')]
    pxds = [pxd(x) for x in data.get('pxdfiles', [])]
    if suffix == '.pyx' and pxds:
        sources.extend(pxds)

    sources.extend(data.get('sources', []))

    include = data.get('include', common_include)

    obj = Extension('pandas.%s' % name,
                    sources=sources,
                    depends=data.get('depends', []),
                    include_dirs=include,
                    extra_compile_args=extra_compile_args)

    extensions.append(obj)


#----------------------------------------------------------------------
# msgpack

if sys.byteorder == 'big':
    macros = [('__BIG_ENDIAN__', '1')]
else:
    macros = [('__LITTLE_ENDIAN__', '1')]

packer_ext = Extension('pandas.msgpack._packer',
                        depends=['pandas/src/msgpack/pack.h',
                                 'pandas/src/msgpack/pack_template.h'],
                        sources = [srcpath('_packer',
                                   suffix=suffix if suffix == '.pyx' else '.cpp',
                                   subdir='msgpack')],
                        language='c++',
                        include_dirs=['pandas/src/msgpack'] + common_include,
                        define_macros=macros,
                        extra_compile_args=extra_compile_args)
unpacker_ext = Extension('pandas.msgpack._unpacker',
                        depends=['pandas/src/msgpack/unpack.h',
                                 'pandas/src/msgpack/unpack_define.h',
                                 'pandas/src/msgpack/unpack_template.h'],
                        sources = [srcpath('_unpacker',
                                   suffix=suffix if suffix == '.pyx' else '.cpp',
                                   subdir='msgpack')],
                        language='c++',
                        include_dirs=['pandas/src/msgpack'] + common_include,
                        define_macros=macros,
                        extra_compile_args=extra_compile_args)
extensions.append(packer_ext)
extensions.append(unpacker_ext)

#----------------------------------------------------------------------
# ujson

if suffix == '.pyx' and 'setuptools' in sys.modules:
    # undo dumb setuptools bug clobbering .pyx sources back to .c
    for ext in extensions:
        if ext.sources[0].endswith(('.c','.cpp')):
            root, _ = os.path.splitext(ext.sources[0])
            ext.sources[0] = root + suffix

ujson_ext = Extension('pandas.json',
                      depends=['pandas/src/ujson/lib/ultrajson.h',
                               'pandas/src/datetime_helper.h',
                               'pandas/src/numpy_helper.h'],
                      sources=['pandas/src/ujson/python/ujson.c',
                               'pandas/src/ujson/python/objToJSON.c',
                               'pandas/src/ujson/python/JSONtoObj.c',
                               'pandas/src/ujson/lib/ultrajsonenc.c',
                               'pandas/src/ujson/lib/ultrajsondec.c',
                               'pandas/src/datetime/np_datetime.c',
                               'pandas/src/datetime/np_datetime_strings.c'],
                      include_dirs=['pandas/src/ujson/python',
                                    'pandas/src/ujson/lib',
                                    'pandas/src/datetime'] + common_include,
                      extra_compile_args=['-D_GNU_SOURCE'] + extra_compile_args)


extensions.append(ujson_ext)

#----------------------------------------------------------------------
# util
# extension for pseudo-safely moving bytes into mutable buffers
_move_ext = Extension('pandas.util._move',
                      depends=[],
                      sources=['pandas/util/move.c'])
extensions.append(_move_ext)


if _have_setuptools:
    setuptools_kwargs["test_suite"] = "nose.collector"

# The build cache system does string matching below this point.
# if you change something, be careful.

setup(name=DISTNAME,
      maintainer=AUTHOR,
      version=versioneer.get_version(),
      packages=['pandas',
                'pandas.api',
                'pandas.api.types',
                'pandas.compat',
                'pandas.compat.numpy',
                'pandas.computation',
                'pandas.core',
                'pandas.indexes',
                'pandas.io',
                'pandas.io.json',
                'pandas.io.sas',
                'pandas.formats',
                'pandas.sparse',
                'pandas.stats',
                'pandas.util',
                'pandas.tests',
                'pandas.tests.api',
                'pandas.tests.computation',
                'pandas.tests.frame',
                'pandas.tests.indexes',
                'pandas.tests.indexes.datetimes',
                'pandas.tests.indexes.timedeltas',
                'pandas.tests.indexes.period',
                'pandas.tests.io',
                'pandas.tests.io.json',
                'pandas.tests.io.parser',
                'pandas.tests.io.sas',
                'pandas.tests.groupby',
                'pandas.tests.series',
                'pandas.tests.formats',
                'pandas.tests.msgpack',
                'pandas.tests.scalar',
                'pandas.tests.sparse',
                'pandas.tests.tseries',
                'pandas.tests.tools',
                'pandas.tests.types',
                'pandas.tests.plotting',
                'pandas.tools',
                'pandas.tseries',
                'pandas.types',
                'pandas.msgpack',
                'pandas.util.clipboard'
                ],
      package_data={'pandas.tests': ['data/*.csv'],
                    'pandas.tests.formats': ['data/*.csv'],
                    'pandas.tests.indexes': ['data/*.pickle'],
                    'pandas.tests.io': ['data/legacy_hdf/*.h5',
                                        'data/legacy_pickle/*/*.pickle',
                                        'data/legacy_msgpack/*/*.msgpack',
                                        'data/*.csv*',
                                        'data/*.dta',
                                        'data/*.pickle',
                                        'data/*.txt',
                                        'data/*.xls',
                                        'data/*.xlsx',
                                        'data/*.xlsm',
                                        'data/*.table',
                                        'parser/data/*.csv',
                                        'parser/data/*.gz',
                                        'parser/data/*.bz2',
                                        'parser/data/*.txt',
                                        'sas/data/*.csv',
                                        'sas/data/*.xpt',
                                        'sas/data/*.sas7bdat',
                                        'data/*.html',
                                        'data/html_encoding/*.html',
                                        'json/data/*.json'],
                    'pandas.tests.tools': ['data/*.csv'],
                    'pandas.tests.tseries': ['data/*.pickle']
                    },
      ext_modules=extensions,
      maintainer_email=EMAIL,
      description=DESCRIPTION,
      license=LICENSE,
      cmdclass=cmdclass,
      url=URL,
      download_url=DOWNLOAD_URL,
      long_description=LONG_DESCRIPTION,
      classifiers=CLASSIFIERS,
      platforms='any',
      **setuptools_kwargs)
