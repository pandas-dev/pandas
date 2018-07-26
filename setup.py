#!/usr/bin/env python

"""
Parts of this file were taken from the pyzmq project
(https://github.com/zeromq/pyzmq) which have been permitted for use under the
BSD license. Parts are from lxml (https://github.com/lxml/lxml)
"""

import os
from os.path import join as pjoin

import pkg_resources
import sys
import shutil
from distutils.version import LooseVersion
from setuptools import setup, Command, find_packages

# versioning
import versioneer
cmdclass = versioneer.get_cmdclass()


def is_platform_windows():
    return sys.platform == 'win32' or sys.platform == 'cygwin'


min_numpy_ver = '1.9.0'
setuptools_kwargs = {
    'install_requires': [
        'python-dateutil >= 2.5.0',
        'pytz >= 2011k',
        'numpy >= {numpy_ver}'.format(numpy_ver=min_numpy_ver),
    ],
    'setup_requires': ['numpy >= {numpy_ver}'.format(numpy_ver=min_numpy_ver)],
    'zip_safe': False,
}


min_cython_ver = '0.28.2'
try:
    import Cython
    ver = Cython.__version__
    from Cython.Build import cythonize
    _CYTHON_INSTALLED = ver >= LooseVersion(min_cython_ver)
except ImportError:
    _CYTHON_INSTALLED = False
    cythonize = lambda x, *args, **kwargs: x  # dummy func

# The import of Extension must be after the import of Cython, otherwise
# we do not get the appropriately patched class.
# See https://cython.readthedocs.io/en/latest/src/reference/compilation.html
from distutils.extension import Extension  # noqa:E402
from distutils.command.build import build  # noqa:E402

try:
    if not _CYTHON_INSTALLED:
        raise ImportError('No supported version of Cython installed.')
    from Cython.Distutils.old_build_ext import old_build_ext as _build_ext
    cython = True
except ImportError:
    from distutils.command.build_ext import build_ext as _build_ext
    cython = False
else:
    try:
        try:
            from Cython import Tempita as tempita
        except ImportError:
            import tempita
    except ImportError:
        raise ImportError('Building pandas requires Tempita: '
                          'pip install Tempita')


_pxi_dep_template = {
    'algos': ['_libs/algos_common_helper.pxi.in',
              '_libs/algos_take_helper.pxi.in',
              '_libs/algos_rank_helper.pxi.in'],
    'groupby': ['_libs/groupby_helper.pxi.in'],
    'join': ['_libs/join_helper.pxi.in', '_libs/join_func_helper.pxi.in'],
    'reshape': ['_libs/reshape_helper.pxi.in'],
    'hashtable': ['_libs/hashtable_class_helper.pxi.in',
                  '_libs/hashtable_func_helper.pxi.in'],
    'index': ['_libs/index_class_helper.pxi.in'],
    'sparse': ['_libs/sparse_op_helper.pxi.in'],
    'interval': ['_libs/intervaltree.pxi.in']}

_pxifiles = []
_pxi_dep = {}
for module, files in _pxi_dep_template.items():
    pxi_files = [pjoin('pandas', x) for x in files]
    _pxifiles.extend(pxi_files)
    _pxi_dep[module] = pxi_files


class build_ext(_build_ext):
    @classmethod
    def render_templates(cls, pxifiles):
        for pxifile in pxifiles:
            # build pxifiles first, template extension must be .pxi.in
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

    def build_extensions(self):
        # if building from c files, don't need to
        # generate template output
        if cython:
            self.render_templates(_pxifiles)

        numpy_incl = pkg_resources.resource_filename('numpy', 'core/include')

        for ext in self.extensions:
            if (hasattr(ext, 'include_dirs') and
                    numpy_incl not in ext.include_dirs):
                ext.include_dirs.append(numpy_incl)
        _build_ext.build_extensions(self)


DESCRIPTION = ("Powerful data structures for data analysis, time series, "
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
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Cython',
    'Topic :: Scientific/Engineering']


class CleanCommand(Command):
    """Custom distutils command to clean the .so and .pyc files."""

    user_options = [("all", "a", "")]

    def initialize_options(self):
        self.all = True
        self._clean_me = []
        self._clean_trees = []

        base = pjoin('pandas', '_libs', 'src')
        dt = pjoin(base, 'datetime')
        src = base
        util = pjoin('pandas', 'util')
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
                               pjoin(util, 'move.c'),
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

    _pyxfiles = ['pandas/_libs/lib.pyx',
                 'pandas/_libs/hashtable.pyx',
                 'pandas/_libs/tslib.pyx',
                 'pandas/_libs/index.pyx',
                 'pandas/_libs/internals.pyx',
                 'pandas/_libs/algos.pyx',
                 'pandas/_libs/join.pyx',
                 'pandas/_libs/indexing.pyx',
                 'pandas/_libs/interval.pyx',
                 'pandas/_libs/hashing.pyx',
                 'pandas/_libs/missing.pyx',
                 'pandas/_libs/reduction.pyx',
                 'pandas/_libs/testing.pyx',
                 'pandas/_libs/skiplist.pyx',
                 'pandas/_libs/sparse.pyx',
                 'pandas/_libs/ops.pyx',
                 'pandas/_libs/parsers.pyx',
                 'pandas/_libs/tslibs/ccalendar.pyx',
                 'pandas/_libs/tslibs/period.pyx',
                 'pandas/_libs/tslibs/strptime.pyx',
                 'pandas/_libs/tslibs/np_datetime.pyx',
                 'pandas/_libs/tslibs/timedeltas.pyx',
                 'pandas/_libs/tslibs/timestamps.pyx',
                 'pandas/_libs/tslibs/timezones.pyx',
                 'pandas/_libs/tslibs/conversion.pyx',
                 'pandas/_libs/tslibs/fields.pyx',
                 'pandas/_libs/tslibs/offsets.pyx',
                 'pandas/_libs/tslibs/frequencies.pyx',
                 'pandas/_libs/tslibs/resolution.pyx',
                 'pandas/_libs/tslibs/parsing.pyx',
                 'pandas/_libs/writers.pyx',
                 'pandas/io/sas/sas.pyx']

    _cpp_pyxfiles = ['pandas/_libs/window.pyx',
                     'pandas/io/msgpack/_packer.pyx',
                     'pandas/io/msgpack/_unpacker.pyx']

    def initialize_options(self):
        sdist_class.initialize_options(self)

    def run(self):
        if 'cython' in cmdclass:
            self.run_command('cython')
        else:
            # If we are not running cython then
            # compile the extensions correctly
            pyx_files = [(self._pyxfiles, 'c'), (self._cpp_pyxfiles, 'cpp')]

            for pyxfiles, extension in pyx_files:
                for pyxfile in pyxfiles:
                    sourcefile = pyxfile[:-3] + extension
                    msg = ("{extension}-source file '{source}' not found.\n"
                           "Run 'setup.py cython' before sdist.".format(
                               source=sourcefile, extension=extension))
                    assert os.path.isfile(sourcefile), msg
        sdist_class.run(self)


class CheckingBuildExt(build_ext):
    """
    Subclass build_ext to get clearer report if Cython is necessary.
    """

    def check_cython_extensions(self, extensions):
        for ext in extensions:
            for src in ext.sources:
                if not os.path.exists(src):
                    print("{}: -> [{}]".format(ext.name, ext.sources))
                    raise Exception("""Cython-generated file '{src}' not found.
                Cython is required to compile pandas from a development branch.
                Please install Cython or download a release package of pandas.
                """.format(src=src))

    def build_extensions(self):
        self.check_cython_extensions(self.extensions)
        build_ext.build_extensions(self)


class CythonCommand(build_ext):
    """
    Custom distutils command subclassed from Cython.Distutils.build_ext
    to compile pyx->c, and stop there. All this does is override the
    C-compile method build_extension() with a no-op.
    """
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

# ----------------------------------------------------------------------
# Preparation of compiler arguments

if sys.byteorder == 'big':
    endian_macro = [('__BIG_ENDIAN__', '1')]
else:
    endian_macro = [('__LITTLE_ENDIAN__', '1')]


if is_platform_windows():
    extra_compile_args = []
else:
    # args to ignore warnings
    extra_compile_args = ['-Wno-unused-function']


# enable coverage by building cython files by setting the environment variable
# "PANDAS_CYTHON_COVERAGE" (with a Truthy value) or by running build_ext
# with `--with-cython-coverage`enabled
linetrace = os.environ.get('PANDAS_CYTHON_COVERAGE', False)
if '--with-cython-coverage' in sys.argv:
    linetrace = True
    sys.argv.remove('--with-cython-coverage')

# Note: if not using `cythonize`, coverage can be enabled by
# pinning `ext.cython_directives = directives` to each ext in extensions.
# github.com/cython/cython/wiki/enhancements-compilerdirectives#in-setuppy
directives = {'linetrace': False}
macros = []
if linetrace:
    # https://pypkg.com/pypi/pytest-cython/f/tests/example-project/setup.py
    directives['linetrace'] = True
    macros = [('CYTHON_TRACE', '1'), ('CYTHON_TRACE_NOGIL', '1')]


# ----------------------------------------------------------------------
# Specification of Dependencies

# TODO: Need to check to see if e.g. `linetrace` has changed and possibly
# re-compile.
def maybe_cythonize(extensions, *args, **kwargs):
    """
    Render tempita templates before calling cythonize
    """
    if len(sys.argv) > 1 and 'clean' in sys.argv:
        # Avoid running cythonize on `python setup.py clean`
        # See https://github.com/cython/cython/issues/1495
        return extensions

    numpy_incl = pkg_resources.resource_filename('numpy', 'core/include')
    # TODO: Is this really necessary here?
    for ext in extensions:
        if (hasattr(ext, 'include_dirs') and
                numpy_incl not in ext.include_dirs):
            ext.include_dirs.append(numpy_incl)

    if cython:
        build_ext.render_templates(_pxifiles)
        return cythonize(extensions, *args, **kwargs)
    else:
        return extensions


lib_depends = ['inference']


def srcpath(name=None, suffix='.pyx', subdir='src'):
    return pjoin('pandas', subdir, name + suffix)


if suffix == '.pyx':
    lib_depends = [srcpath(f, suffix='.pyx', subdir='_libs/src')
                   for f in lib_depends]
else:
    lib_depends = []

common_include = ['pandas/_libs/src/klib', 'pandas/_libs/src']


lib_depends = lib_depends + ['pandas/_libs/src/numpy_helper.h',
                             'pandas/_libs/src/parse_helper.h',
                             'pandas/_libs/src/compat_helper.h']

np_datetime_headers = ['pandas/_libs/src/datetime/np_datetime.h',
                       'pandas/_libs/src/datetime/np_datetime_strings.h']
np_datetime_sources = ['pandas/_libs/src/datetime/np_datetime.c',
                       'pandas/_libs/src/datetime/np_datetime_strings.c']

tseries_depends = np_datetime_headers


ext_data = {
    '_libs.algos': {
        'pyxfile': '_libs/algos',
        'depends': _pxi_dep['algos']},
    '_libs.groupby': {
        'pyxfile': '_libs/groupby',
        'depends': _pxi_dep['groupby']},
    '_libs.hashing': {
        'pyxfile': '_libs/hashing'},
    '_libs.hashtable': {
        'pyxfile': '_libs/hashtable',
        'depends': (['pandas/_libs/src/klib/khash_python.h'] +
                    _pxi_dep['hashtable'])},
    '_libs.index': {
        'pyxfile': '_libs/index',
        'depends': _pxi_dep['index'],
        'sources': np_datetime_sources},
    '_libs.indexing': {
        'pyxfile': '_libs/indexing'},
    '_libs.internals': {
        'pyxfile': '_libs/internals'},
    '_libs.interval': {
        'pyxfile': '_libs/interval',
        'depends': _pxi_dep['interval']},
    '_libs.join': {
        'pyxfile': '_libs/join',
        'depends': _pxi_dep['join']},
    '_libs.lib': {
        'pyxfile': '_libs/lib',
        'depends': lib_depends + tseries_depends},
    '_libs.missing': {
        'pyxfile': '_libs/missing',
        'depends': tseries_depends},
    '_libs.parsers': {
        'pyxfile': '_libs/parsers',
        'depends': ['pandas/_libs/src/parser/tokenizer.h',
                    'pandas/_libs/src/parser/io.h',
                    'pandas/_libs/src/numpy_helper.h'],
        'sources': ['pandas/_libs/src/parser/tokenizer.c',
                    'pandas/_libs/src/parser/io.c']},
    '_libs.reduction': {
        'pyxfile': '_libs/reduction'},
    '_libs.ops': {
        'pyxfile': '_libs/ops'},
    '_libs.properties': {
        'pyxfile': '_libs/properties',
        'include': []},
    '_libs.reshape': {
        'pyxfile': '_libs/reshape',
        'depends': _pxi_dep['reshape']},
    '_libs.skiplist': {
        'pyxfile': '_libs/skiplist',
        'depends': ['pandas/_libs/src/skiplist.h']},
    '_libs.sparse': {
        'pyxfile': '_libs/sparse',
        'depends': _pxi_dep['sparse']},
    '_libs.tslib': {
        'pyxfile': '_libs/tslib',
        'depends': tseries_depends,
        'sources': np_datetime_sources},
    '_libs.tslibs.ccalendar': {
        'pyxfile': '_libs/tslibs/ccalendar'},
    '_libs.tslibs.conversion': {
        'pyxfile': '_libs/tslibs/conversion',
        'depends': tseries_depends,
        'sources': np_datetime_sources},
    '_libs.tslibs.fields': {
        'pyxfile': '_libs/tslibs/fields',
        'depends': tseries_depends,
        'sources': np_datetime_sources},
    '_libs.tslibs.frequencies': {
        'pyxfile': '_libs/tslibs/frequencies'},
    '_libs.tslibs.nattype': {
        'pyxfile': '_libs/tslibs/nattype'},
    '_libs.tslibs.np_datetime': {
        'pyxfile': '_libs/tslibs/np_datetime',
        'depends': np_datetime_headers,
        'sources': np_datetime_sources},
    '_libs.tslibs.offsets': {
        'pyxfile': '_libs/tslibs/offsets',
        'depends': tseries_depends,
        'sources': np_datetime_sources},
    '_libs.tslibs.parsing': {
        'pyxfile': '_libs/tslibs/parsing'},
    '_libs.tslibs.period': {
        'pyxfile': '_libs/tslibs/period',
        'depends': tseries_depends + ['pandas/_libs/src/period_helper.h'],
        'sources': np_datetime_sources + ['pandas/_libs/src/period_helper.c']},
    '_libs.tslibs.resolution': {
        'pyxfile': '_libs/tslibs/resolution',
        'depends': tseries_depends,
        'sources': np_datetime_sources},
    '_libs.tslibs.strptime': {
        'pyxfile': '_libs/tslibs/strptime',
        'depends': tseries_depends,
        'sources': np_datetime_sources},
    '_libs.tslibs.timedeltas': {
        'pyxfile': '_libs/tslibs/timedeltas',
        'depends': np_datetime_headers,
        'sources': np_datetime_sources},
    '_libs.tslibs.timestamps': {
        'pyxfile': '_libs/tslibs/timestamps',
        'depends': tseries_depends,
        'sources': np_datetime_sources},
    '_libs.tslibs.timezones': {
        'pyxfile': '_libs/tslibs/timezones'},
    '_libs.testing': {
        'pyxfile': '_libs/testing'},
    '_libs.window': {
        'pyxfile': '_libs/window',
        'language': 'c++',
        'suffix': '.cpp'},
    '_libs.writers': {
        'pyxfile': '_libs/writers'},
    'io.sas._sas': {
        'pyxfile': 'io/sas/sas'},
    'io.msgpack._packer': {
        'macros': endian_macro + macros,
        'depends': ['pandas/_libs/src/msgpack/pack.h',
                    'pandas/_libs/src/msgpack/pack_template.h'],
        'include': ['pandas/_libs/src/msgpack'] + common_include,
        'language': 'c++',
        'suffix': '.cpp',
        'pyxfile': 'io/msgpack/_packer',
        'subdir': 'io/msgpack'},
    'io.msgpack._unpacker': {
        'depends': ['pandas/_libs/src/msgpack/unpack.h',
                    'pandas/_libs/src/msgpack/unpack_define.h',
                    'pandas/_libs/src/msgpack/unpack_template.h'],
        'macros': endian_macro + macros,
        'include': ['pandas/_libs/src/msgpack'] + common_include,
        'language': 'c++',
        'suffix': '.cpp',
        'pyxfile': 'io/msgpack/_unpacker',
        'subdir': 'io/msgpack'
    }
}

extensions = []

for name, data in ext_data.items():
    source_suffix = suffix if suffix == '.pyx' else data.get('suffix', '.c')

    sources = [srcpath(data['pyxfile'], suffix=source_suffix, subdir='')]

    sources.extend(data.get('sources', []))

    include = data.get('include', common_include)

    obj = Extension('pandas.{name}'.format(name=name),
                    sources=sources,
                    depends=data.get('depends', []),
                    include_dirs=include,
                    language=data.get('language', 'c'),
                    define_macros=data.get('macros', macros),
                    extra_compile_args=extra_compile_args)

    extensions.append(obj)

# ----------------------------------------------------------------------
# ujson

if suffix == '.pyx':
    # undo dumb setuptools bug clobbering .pyx sources back to .c
    for ext in extensions:
        if ext.sources[0].endswith(('.c', '.cpp')):
            root, _ = os.path.splitext(ext.sources[0])
            ext.sources[0] = root + suffix

ujson_ext = Extension('pandas._libs.json',
                      depends=['pandas/_libs/src/ujson/lib/ultrajson.h'],
                      sources=(['pandas/_libs/src/ujson/python/ujson.c',
                                'pandas/_libs/src/ujson/python/objToJSON.c',
                                'pandas/_libs/src/ujson/python/JSONtoObj.c',
                                'pandas/_libs/src/ujson/lib/ultrajsonenc.c',
                                'pandas/_libs/src/ujson/lib/ultrajsondec.c'] +
                               np_datetime_sources),
                      include_dirs=['pandas/_libs/src/ujson/python',
                                    'pandas/_libs/src/ujson/lib',
                                    'pandas/_libs/src/datetime'],
                      extra_compile_args=(['-D_GNU_SOURCE'] +
                                          extra_compile_args),
                      define_macros=macros)


extensions.append(ujson_ext)

# ----------------------------------------------------------------------
# util
# extension for pseudo-safely moving bytes into mutable buffers
_move_ext = Extension('pandas.util._move',
                      depends=[],
                      sources=['pandas/util/move.c'],
                      define_macros=macros)
extensions.append(_move_ext)

# The build cache system does string matching below this point.
# if you change something, be careful.

setup(name=DISTNAME,
      maintainer=AUTHOR,
      version=versioneer.get_version(),
      packages=find_packages(include=['pandas', 'pandas.*']),
      package_data={'': ['templates/*', '_libs/*.dll']},
      ext_modules=maybe_cythonize(extensions, compiler_directives=directives),
      maintainer_email=EMAIL,
      description=DESCRIPTION,
      license=LICENSE,
      cmdclass=cmdclass,
      url=URL,
      download_url=DOWNLOAD_URL,
      long_description=LONG_DESCRIPTION,
      classifiers=CLASSIFIERS,
      platforms='any',
      python_requires='>=2.7,!=3.0.*,!=3.1.*,!=3.2.*,!=3.3.*,!=3.4.*',
      **setuptools_kwargs)
