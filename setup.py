#!/usr/bin/env python

"""
Parts of this file were taken from the pyzmq project
(https://github.com/zeromq/pyzmq) which have been permitted for use under the
BSD license. Parts are from lxml (https://github.com/lxml/lxml)
"""

import glob
import os
import os.path as osp
import sys
import shutil
import re
import platform
from distutils.version import LooseVersion
from distutils import sysconfig

from os.path import join as pjoin

# versioning
import versioneer
cmdclass = versioneer.get_cmdclass()

min_cython_ver = '0.19.1'
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
            import numpy  # noqa
            import dateutil  # noqa
            setuptools_kwargs = {}
        except ImportError:
            sys.exit("install requires: 'python-dateutil < 2','numpy'."
                     "  use pip or easy_install."
                     "\n   $ pip install 'python-dateutil < 2' 'numpy'")


# Check if we're running 64-bit Python
is_64_bit = sys.maxsize > 2**32

# Check if this is a debug build of Python.
if hasattr(sys, 'gettotalrefcount'):
    build_type = 'Debug'
else:
    build_type = 'Release'


from distutils.extension import Extension
from distutils.command.build import build
from distutils.command.build_ext import build_ext as _build_ext

try:
    if not _CYTHON_INSTALLED:
        raise ImportError('No supported version of Cython installed.')
    from Cython.Distutils import build_ext as _build_ext  # noqa
    cython = True
except ImportError:
    cython = False


class build_ext(_build_ext):

    def build_extensions(self):
        numpy_incl = pkg_resources.resource_filename('numpy', 'core/include')

        for ext in self.extensions:
            if (hasattr(ext, 'include_dirs') and
                    numpy_incl not in ext.include_dirs):
                ext.include_dirs.append(numpy_incl)
        _build_ext.build_extensions(self)

    def run(self):
        self._run_cmake()
        _build_ext.run(self)

    # adapted from cmake_build_ext in dynd-python
    # github.com/libdynd/dynd-python

    description = "Build the C-extension for libpandas and regular pandas"
    user_options = ([('extra-cmake-args=', None,
                      'extra arguments for CMake')] +
                    _build_ext.user_options)

    def initialize_options(self):
        _build_ext.initialize_options(self)
        self.extra_cmake_args = ''

    def _run_cmake(self):
        # The directory containing this setup.py
        source = osp.dirname(osp.abspath(__file__))

        # The staging directory for the module being built
        build_temp = pjoin(os.getcwd(), self.build_temp)

        # Change to the build directory
        saved_cwd = os.getcwd()
        if not os.path.isdir(self.build_temp):
            self.mkpath(self.build_temp)
        os.chdir(self.build_temp)

        # Detect if we built elsewhere
        if os.path.isfile('CMakeCache.txt'):
            cachefile = open('CMakeCache.txt', 'r')
            cachedir = re.search('CMAKE_CACHEFILE_DIR:INTERNAL=(.*)',
                                 cachefile.read()).group(1)
            cachefile.close()
            if (cachedir != build_temp):
                return

        pyexe_option = '-DPYTHON_EXECUTABLE=%s' % sys.executable
        static_lib_option = ''
        build_tests_option = ''

        if sys.platform != 'win32':
            cmake_command = ['cmake', self.extra_cmake_args, pyexe_option,
                             build_tests_option,
                             static_lib_option, source]

            self.spawn(cmake_command)
            self.spawn(['make'])
        else:
            import shlex
            cmake_generator = 'Visual Studio 14 2015'
            if is_64_bit:
                cmake_generator += ' Win64'
            # Generate the build files
            extra_cmake_args = shlex.split(self.extra_cmake_args)
            cmake_command = (['cmake'] + extra_cmake_args +
                             [source, pyexe_option,
                              static_lib_option,
                              build_tests_option,
                             '-G', cmake_generator])
            if "-G" in self.extra_cmake_args:
                cmake_command = cmake_command[:-2]

            self.spawn(cmake_command)
            # Do the build
            self.spawn(['cmake', '--build', '.', '--config', build_type])

        if self.inplace:
            # a bit hacky
            build_lib = saved_cwd
        else:
            build_lib = pjoin(os.getcwd(), self.build_lib)

        # Move the built libpandas library to the place expected by the Python
        # build
        if sys.platform != 'win32':
            name, = glob.glob('libpandas.*')
            try:
                os.makedirs(pjoin(build_lib, 'pandas'))
            except OSError:
                pass
            shutil.move(name, pjoin(build_lib, 'pandas', name))
        else:
            shutil.move(pjoin(build_type, 'pandas.dll'),
                        pjoin(build_lib, 'pandas', 'pandas.dll'))

        # Move the built C-extension to the place expected by the Python build
        self._found_names = []
        for name in self.get_cmake_cython_names():
            built_path = self.get_ext_built(name)
            if not os.path.exists(built_path):
                raise RuntimeError('libpandas C-extension failed to build:',
                                   os.path.abspath(built_path))

            ext_path = pjoin(build_lib, self._get_cmake_ext_path(name))
            if os.path.exists(ext_path):
                os.remove(ext_path)
            self.mkpath(os.path.dirname(ext_path))
            print('Moving built libpandas C-extension', built_path,
                  'to build path', ext_path)
            shutil.move(self.get_ext_built(name), ext_path)
            self._found_names.append(name)

        os.chdir(saved_cwd)

    def _get_inplace_dir(self):
        pass

    def _get_cmake_ext_path(self, name):
        # Get the package directory from build_py
        build_py = self.get_finalized_command('build_py')
        package_dir = build_py.get_package_dir('pandas')
        # This is the name of the pandas C-extension
        suffix = sysconfig.get_config_var('EXT_SUFFIX')
        if suffix is None:
            suffix = sysconfig.get_config_var('SO')
        filename = name + suffix
        return pjoin(package_dir, filename)

    def get_ext_built(self, name):
        if sys.platform == 'win32':
            head, tail = os.path.split(name)
            suffix = sysconfig.get_config_var('SO')
            return pjoin(head, build_type, tail + suffix)
        else:
            suffix = sysconfig.get_config_var('SO')
            return name + suffix

    def get_cmake_cython_names(self):
        return ['native', osp.join('internals', 'config')]

    def get_names(self):
        return self._found_names

    def get_outputs(self):
        # Just the C extensions
        cmake_exts = [self._get_cmake_ext_path(name)
                      for name in self.get_names()]
        regular_exts = _build_ext.get_outputs(self)
        return regular_exts + cmake_exts


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
    'Programming Language :: Python :: 2.6',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
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

        base = pjoin('pandas', 'src')
        dt = pjoin(base, 'datetime')
        src = base
        parser = pjoin(base, 'parser')
        ujson_python = pjoin(base, 'ujson', 'python')
        ujson_lib = pjoin(base, 'ujson', 'lib')
        self._clean_exclude = [pjoin(dt, 'np_datetime.c'),
                               pjoin(dt, 'np_datetime_strings.c'),
                               pjoin(src, 'period_helper.c'),
                               pjoin(parser, 'tokenizer.c'),
                               pjoin(parser, 'io.c'),
                               pjoin(ujson_python, 'ujson.c'),
                               pjoin(ujson_python, 'objToJSON.c'),
                               pjoin(ujson_python, 'JSONtoObj.c'),
                               pjoin(ujson_lib, 'ultrajsonenc.c'),
                               pjoin(ujson_lib, 'ultrajsondec.c'),
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
                 'pandas/parser.pyx',
                 'pandas/src/period.pyx',
                 'pandas/src/sparse.pyx',
                 'pandas/src/testing.pyx']

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
    Also, add some platform based compiler flags.
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
        self.add_gnu_inline_flag(self.extensions)
        build_ext.build_extensions(self)

    def add_gnu_inline_flag(self, extensions):
        '''
        Add CFLAGS `-fgnu89-inline` for clang on FreeBSD 10+
        '''
        if not platform.system() == 'FreeBSD':
            return

        try:
            bsd_release = float(platform.release().split('-')[0])
        except ValueError:  # unknow freebsd version
            return

        if bsd_release < 10:  # 9 or earlier still using gcc42
            return

        for ext in extensions:
            ext.extra_compile_args += ['-fgnu89-inline']


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


lib_depends = lib_depends + ['pandas/src/numpy_helper.h',
                             'pandas/src/parse_helper.h']


tseries_depends = ['pandas/src/datetime/np_datetime.h',
                   'pandas/src/datetime/np_datetime_strings.h',
                   'pandas/src/period_helper.h']


# some linux distros require it
libraries = ['m'] if 'win32' not in sys.platform else []

ext_data = dict(
    lib={'pyxfile': 'lib',
         'pxdfiles': [],
         'depends': lib_depends},
    hashtable={'pyxfile': 'hashtable',
               'pxdfiles': ['hashtable'],
               'depends': ['pandas/src/klib/khash_python.h']},
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
                       'pandas/src/datetime/np_datetime_strings.c']},
    algos={'pyxfile': 'algos',
           'pxdfiles': ['src/skiplist'],
           'depends': [srcpath('generated', suffix='.pyx'),
                       srcpath('join', suffix='.pyx'),
                       'pandas/src/skiplist.pyx',
                       'pandas/src/skiplist.h']},
    parser={'pyxfile': 'parser',
            'depends': ['pandas/src/parser/tokenizer.h',
                        'pandas/src/parser/io.h',
                        'pandas/src/numpy_helper.h'],
            'sources': ['pandas/src/parser/tokenizer.c',
                        'pandas/src/parser/io.c']}
)

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
                    include_dirs=include)

    extensions.append(obj)


sparse_ext = Extension('pandas._sparse',
                       sources=[srcpath('sparse', suffix=suffix)],
                       include_dirs=[],
                       libraries=libraries)

extensions.extend([sparse_ext])

testing_ext = Extension('pandas._testing',
                        sources=[srcpath('testing', suffix=suffix)],
                        include_dirs=[],
                        libraries=libraries)

extensions.extend([testing_ext])

# ---------------------------------------------------------------------
# msgpack stuff here

if sys.byteorder == 'big':
    macros = [('__BIG_ENDIAN__', '1')]
else:
    macros = [('__LITTLE_ENDIAN__', '1')]

packer_ext = Extension('pandas.msgpack._packer',
                       depends=['pandas/src/msgpack/pack.h',
                                'pandas/src/msgpack/pack_template.h'],
                       sources=[
                           srcpath('_packer', suffix=suffix
                                   if suffix == '.pyx' else '.cpp',
                                   subdir='msgpack')],
                       language='c++',
                       include_dirs=['pandas/src/msgpack'] + common_include,
                       define_macros=macros)

unpacker_ext = Extension('pandas.msgpack._unpacker',
                         depends=['pandas/src/msgpack/unpack.h',
                                  'pandas/src/msgpack/unpack_define.h',
                                  'pandas/src/msgpack/unpack_template.h'],
                         sources=[
                             srcpath('_unpacker',
                                     suffix=suffix
                                     if suffix == '.pyx' else '.cpp',
                                     subdir='msgpack')],
                         language='c++',
                         include_dirs=['pandas/src/msgpack'] + common_include,
                         define_macros=macros)
extensions.append(packer_ext)
extensions.append(unpacker_ext)

if suffix == '.pyx' and 'setuptools' in sys.modules:
    # undo dumb setuptools bug clobbering .pyx sources back to .c
    for ext in extensions:
        if ext.sources[0].endswith(('.c', '.cpp')):
            root, _ = os.path.splitext(ext.sources[0])
            ext.sources[0] = root + suffix

ujson_ext = Extension('pandas.json',
                      depends=['pandas/src/ujson/lib/ultrajson.h',
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
                      extra_compile_args=['-D_GNU_SOURCE'])


extensions.append(ujson_ext)


if _have_setuptools:
    setuptools_kwargs["test_suite"] = "nose.collector"

# The build cache system does string matching below this point.
# if you change something, be careful.

setup(name=DISTNAME,
      maintainer=AUTHOR,
      version=versioneer.get_version(),
      packages=['pandas',
                'pandas.compat',
                'pandas.computation',
                'pandas.computation.tests',
                'pandas.core',
                'pandas.io',
                'pandas.rpy',
                'pandas.sandbox',
                'pandas.sparse',
                'pandas.sparse.tests',
                'pandas.stats',
                'pandas.util',
                'pandas.tests',
                'pandas.tests.test_msgpack',
                'pandas.tools',
                'pandas.tools.tests',
                'pandas.tseries',
                'pandas.tseries.tests',
                'pandas.io.tests',
                'pandas.io.tests.test_json',
                'pandas.stats.tests',
                'pandas.msgpack'
                ],
      package_data={'pandas.io': ['tests/data/legacy_hdf/*.h5',
                                  'tests/data/legacy_pickle/*/*.pickle',
                                  'tests/data/legacy_msgpack/*/*.msgpack',
                                  'tests/data/*.csv*',
                                  'tests/data/*.XPT',
                                  'tests/data/*.dta',
                                  'tests/data/*.txt',
                                  'tests/data/*.xls',
                                  'tests/data/*.xlsx',
                                  'tests/data/*.xlsm',
                                  'tests/data/*.table',
                                  'tests/data/*.html',
                                  'tests/data/html_encoding/*.html',
                                  'tests/test_json/data/*.json'],
                    'pandas.tools': ['tests/*.csv'],
                    'pandas.tests': ['data/*.pickle',
                                     'data/*.csv'],
                    'pandas.tseries.tests': ['data/*.pickle',
                                             'data/*.csv']
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
