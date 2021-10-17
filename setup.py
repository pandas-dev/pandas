#!/usr/bin/env python3

"""
Parts of this file were taken from the pyzmq project
(https://github.com/zeromq/pyzmq) which have been permitted for use under the
BSD license. Parts are from lxml (https://github.com/lxml/lxml)
"""

import argparse
import multiprocessing
import os
from os.path import join as pjoin
import platform
import shutil
import sys
from sysconfig import get_config_vars

import numpy
from pkg_resources import parse_version
from setuptools import (
    Command,
    Extension,
    setup,
)
from setuptools.command.build_ext import build_ext as _build_ext

import versioneer

cmdclass = versioneer.get_cmdclass()


def is_platform_windows():
    return sys.platform == "win32" or sys.platform == "cygwin"


def is_platform_mac():
    return sys.platform == "darwin"


min_cython_ver = "0.29.24"  # note: sync with pyproject.toml

try:
    from Cython import (
        Tempita,
        __version__ as _CYTHON_VERSION,
    )
    from Cython.Build import cythonize

    _CYTHON_INSTALLED = parse_version(_CYTHON_VERSION) >= parse_version(min_cython_ver)
except ImportError:
    _CYTHON_VERSION = None
    _CYTHON_INSTALLED = False
    cythonize = lambda x, *args, **kwargs: x  # dummy func


_pxi_dep_template = {
    "algos": ["_libs/algos_common_helper.pxi.in", "_libs/algos_take_helper.pxi.in"],
    "hashtable": [
        "_libs/hashtable_class_helper.pxi.in",
        "_libs/hashtable_func_helper.pxi.in",
        "_libs/khash_for_primitive_helper.pxi.in",
    ],
    "index": ["_libs/index_class_helper.pxi.in"],
    "sparse": ["_libs/sparse_op_helper.pxi.in"],
    "interval": ["_libs/intervaltree.pxi.in"],
}

_pxifiles = []
_pxi_dep = {}
for module, files in _pxi_dep_template.items():
    pxi_files = [pjoin("pandas", x) for x in files]
    _pxifiles.extend(pxi_files)
    _pxi_dep[module] = pxi_files


class build_ext(_build_ext):
    @classmethod
    def render_templates(cls, pxifiles):
        for pxifile in pxifiles:
            # build pxifiles first, template extension must be .pxi.in
            assert pxifile.endswith(".pxi.in")
            outfile = pxifile[:-3]

            if (
                os.path.exists(outfile)
                and os.stat(pxifile).st_mtime < os.stat(outfile).st_mtime
            ):
                # if .pxi.in is not updated, no need to output .pxi
                continue

            with open(pxifile) as f:
                tmpl = f.read()
            pyxcontent = Tempita.sub(tmpl)

            with open(outfile, "w") as f:
                f.write(pyxcontent)

    def build_extensions(self):
        # if building from c files, don't need to
        # generate template output
        if _CYTHON_INSTALLED:
            self.render_templates(_pxifiles)

        super().build_extensions()


class CleanCommand(Command):
    """Custom command to clean the .so and .pyc files."""

    user_options = [("all", "a", "")]

    def initialize_options(self):
        self.all = True
        self._clean_me = []
        self._clean_trees = []

        base = pjoin("pandas", "_libs", "src")
        tsbase = pjoin("pandas", "_libs", "tslibs", "src")
        dt = pjoin(tsbase, "datetime")
        util = pjoin("pandas", "util")
        parser = pjoin(base, "parser")
        ujson_python = pjoin(base, "ujson", "python")
        ujson_lib = pjoin(base, "ujson", "lib")
        self._clean_exclude = [
            pjoin(dt, "np_datetime.c"),
            pjoin(dt, "np_datetime_strings.c"),
            pjoin(parser, "tokenizer.c"),
            pjoin(parser, "io.c"),
            pjoin(ujson_python, "ujson.c"),
            pjoin(ujson_python, "objToJSON.c"),
            pjoin(ujson_python, "JSONtoObj.c"),
            pjoin(ujson_python, "date_conversions.c"),
            pjoin(ujson_lib, "ultrajsonenc.c"),
            pjoin(ujson_lib, "ultrajsondec.c"),
            pjoin(util, "move.c"),
        ]

        for root, dirs, files in os.walk("pandas"):
            for f in files:
                filepath = pjoin(root, f)
                if filepath in self._clean_exclude:
                    continue

                if os.path.splitext(f)[-1] in (
                    ".pyc",
                    ".so",
                    ".o",
                    ".pyo",
                    ".pyd",
                    ".c",
                    ".cpp",
                    ".orig",
                ):
                    self._clean_me.append(filepath)
            for d in dirs:
                if d == "__pycache__":
                    self._clean_trees.append(pjoin(root, d))

        # clean the generated pxi files
        for pxifile in _pxifiles:
            pxifile = pxifile.replace(".pxi.in", ".pxi")
            self._clean_me.append(pxifile)

        for d in ("build", "dist"):
            if os.path.exists(d):
                self._clean_trees.append(d)

    def finalize_options(self):
        pass

    def run(self):
        for clean_me in self._clean_me:
            try:
                os.unlink(clean_me)
            except OSError:
                pass
        for clean_tree in self._clean_trees:
            try:
                shutil.rmtree(clean_tree)
            except OSError:
                pass


# we need to inherit from the versioneer
# class as it encodes the version info
sdist_class = cmdclass["sdist"]


class CheckSDist(sdist_class):
    """Custom sdist that ensures Cython has compiled all pyx files to c."""

    _pyxfiles = [
        "pandas/_libs/arrays.pyx",
        "pandas/_libs/lib.pyx",
        "pandas/_libs/hashtable.pyx",
        "pandas/_libs/tslib.pyx",
        "pandas/_libs/index.pyx",
        "pandas/_libs/internals.pyx",
        "pandas/_libs/algos.pyx",
        "pandas/_libs/join.pyx",
        "pandas/_libs/indexing.pyx",
        "pandas/_libs/interval.pyx",
        "pandas/_libs/hashing.pyx",
        "pandas/_libs/missing.pyx",
        "pandas/_libs/reduction.pyx",
        "pandas/_libs/testing.pyx",
        "pandas/_libs/sparse.pyx",
        "pandas/_libs/ops.pyx",
        "pandas/_libs/parsers.pyx",
        "pandas/_libs/tslibs/base.pyx",
        "pandas/_libs/tslibs/ccalendar.pyx",
        "pandas/_libs/tslibs/dtypes.pyx",
        "pandas/_libs/tslibs/period.pyx",
        "pandas/_libs/tslibs/strptime.pyx",
        "pandas/_libs/tslibs/np_datetime.pyx",
        "pandas/_libs/tslibs/timedeltas.pyx",
        "pandas/_libs/tslibs/timestamps.pyx",
        "pandas/_libs/tslibs/timezones.pyx",
        "pandas/_libs/tslibs/conversion.pyx",
        "pandas/_libs/tslibs/fields.pyx",
        "pandas/_libs/tslibs/offsets.pyx",
        "pandas/_libs/tslibs/parsing.pyx",
        "pandas/_libs/tslibs/tzconversion.pyx",
        "pandas/_libs/tslibs/vectorized.pyx",
        "pandas/_libs/window/indexers.pyx",
        "pandas/_libs/writers.pyx",
        "pandas/io/sas/sas.pyx",
    ]

    _cpp_pyxfiles = [
        "pandas/_libs/window/aggregations.pyx",
    ]

    def initialize_options(self):
        sdist_class.initialize_options(self)

    def run(self):
        if "cython" in cmdclass:
            self.run_command("cython")
        else:
            # If we are not running cython then
            # compile the extensions correctly
            pyx_files = [(self._pyxfiles, "c"), (self._cpp_pyxfiles, "cpp")]

            for pyxfiles, extension in pyx_files:
                for pyxfile in pyxfiles:
                    sourcefile = pyxfile[:-3] + extension
                    msg = (
                        f"{extension}-source file '{sourcefile}' not found.\n"
                        "Run 'setup.py cython' before sdist."
                    )
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
                    print(f"{ext.name}: -> [{ext.sources}]")
                    raise Exception(
                        f"""Cython-generated file '{src}' not found.
                Cython is required to compile pandas from a development branch.
                Please install Cython or download a release package of pandas.
                """
                    )

    def build_extensions(self):
        self.check_cython_extensions(self.extensions)
        build_ext.build_extensions(self)


class CythonCommand(build_ext):
    """
    Custom command subclassed from Cython.Distutils.build_ext
    to compile pyx->c, and stop there. All this does is override the
    C-compile method build_extension() with a no-op.
    """

    def build_extension(self, ext):
        pass


class DummyBuildSrc(Command):
    """numpy's build_src command interferes with Cython's build_ext."""

    user_options = []

    def initialize_options(self):
        self.py_modules_dict = {}

    def finalize_options(self):
        pass

    def run(self):
        pass


cmdclass["clean"] = CleanCommand
cmdclass["build_ext"] = CheckingBuildExt

if _CYTHON_INSTALLED:
    suffix = ".pyx"
    cmdclass["cython"] = CythonCommand
else:
    suffix = ".c"
    cmdclass["build_src"] = DummyBuildSrc

# ----------------------------------------------------------------------
# Preparation of compiler arguments

debugging_symbols_requested = "--with-debugging-symbols" in sys.argv
if debugging_symbols_requested:
    sys.argv.remove("--with-debugging-symbols")


if sys.byteorder == "big":
    endian_macro = [("__BIG_ENDIAN__", "1")]
else:
    endian_macro = [("__LITTLE_ENDIAN__", "1")]


extra_compile_args = []
extra_link_args = []
if is_platform_windows():
    if debugging_symbols_requested:
        extra_compile_args.append("/Z7")
        extra_link_args.append("/DEBUG")
else:
    # PANDAS_CI=1 is set by ci/setup_env.sh
    if os.environ.get("PANDAS_CI", "0") == "1":
        extra_compile_args.append("-Werror")
    if debugging_symbols_requested:
        extra_compile_args.append("-g")
        extra_compile_args.append("-UNDEBUG")
        extra_compile_args.append("-O0")

# Build for at least macOS 10.9 when compiling on a 10.9 system or above,
# overriding CPython distuitls behaviour which is to target the version that
# python was built for. This may be overridden by setting
# MACOSX_DEPLOYMENT_TARGET before calling setup.py
if is_platform_mac():
    if "MACOSX_DEPLOYMENT_TARGET" not in os.environ:
        current_system = platform.mac_ver()[0]
        python_target = get_config_vars().get(
            "MACOSX_DEPLOYMENT_TARGET", current_system
        )
        target_macos_version = "10.9"
        parsed_macos_version = parse_version(target_macos_version)
        if (
            parse_version(str(python_target)) < parsed_macos_version
            and parse_version(current_system) >= parsed_macos_version
        ):
            os.environ["MACOSX_DEPLOYMENT_TARGET"] = target_macos_version

    if sys.version_info[:2] == (3, 8):  # GH 33239
        extra_compile_args.append("-Wno-error=deprecated-declarations")

    # https://github.com/pandas-dev/pandas/issues/35559
    extra_compile_args.append("-Wno-error=unreachable-code")

# enable coverage by building cython files by setting the environment variable
# "PANDAS_CYTHON_COVERAGE" (with a Truthy value) or by running build_ext
# with `--with-cython-coverage`enabled
linetrace = os.environ.get("PANDAS_CYTHON_COVERAGE", False)
if "--with-cython-coverage" in sys.argv:
    linetrace = True
    sys.argv.remove("--with-cython-coverage")

# Note: if not using `cythonize`, coverage can be enabled by
# pinning `ext.cython_directives = directives` to each ext in extensions.
# github.com/cython/cython/wiki/enhancements-compilerdirectives#in-setuppy
directives = {"linetrace": False, "language_level": 3}
macros = []
if linetrace:
    # https://pypkg.com/pypi/pytest-cython/f/tests/example-project/setup.py
    directives["linetrace"] = True
    macros = [("CYTHON_TRACE", "1"), ("CYTHON_TRACE_NOGIL", "1")]

# silence build warnings about deprecated API usage
# we can't do anything about these warnings because they stem from
# cython+numpy version mismatches.
macros.append(("NPY_NO_DEPRECATED_API", "0"))


# ----------------------------------------------------------------------
# Specification of Dependencies

# TODO: Need to check to see if e.g. `linetrace` has changed and possibly
# re-compile.
def maybe_cythonize(extensions, *args, **kwargs):
    """
    Render tempita templates before calling cythonize. This is skipped for

    * clean
    * sdist
    """
    if "clean" in sys.argv or "sdist" in sys.argv:
        # See https://github.com/cython/cython/issues/1495
        return extensions

    elif not _CYTHON_INSTALLED:
        # GH#28836 raise a helfpul error message
        if _CYTHON_VERSION:
            raise RuntimeError(
                f"Cannot cythonize with old Cython version ({_CYTHON_VERSION} "
                f"installed, needs {min_cython_ver})"
            )
        raise RuntimeError("Cannot cythonize without Cython installed.")

    # reuse any parallel arguments provided for compilation to cythonize
    parser = argparse.ArgumentParser()
    parser.add_argument("--parallel", "-j", type=int, default=1)
    parsed, _ = parser.parse_known_args()

    kwargs["nthreads"] = parsed.parallel
    build_ext.render_templates(_pxifiles)
    return cythonize(extensions, *args, **kwargs)


def srcpath(name=None, suffix=".pyx", subdir="src"):
    return pjoin("pandas", subdir, name + suffix)


lib_depends = ["pandas/_libs/src/parse_helper.h"]

klib_include = ["pandas/_libs/src/klib"]

tseries_depends = [
    "pandas/_libs/tslibs/src/datetime/np_datetime.h",
    "pandas/_libs/tslibs/src/datetime/np_datetime_strings.h",
]

ext_data = {
    "_libs.algos": {
        "pyxfile": "_libs/algos",
        "include": klib_include,
        "depends": _pxi_dep["algos"],
    },
    "_libs.arrays": {"pyxfile": "_libs/arrays"},
    "_libs.groupby": {"pyxfile": "_libs/groupby"},
    "_libs.hashing": {"pyxfile": "_libs/hashing", "depends": []},
    "_libs.hashtable": {
        "pyxfile": "_libs/hashtable",
        "include": klib_include,
        "depends": (
            ["pandas/_libs/src/klib/khash_python.h", "pandas/_libs/src/klib/khash.h"]
            + _pxi_dep["hashtable"]
        ),
    },
    "_libs.index": {
        "pyxfile": "_libs/index",
        "include": klib_include,
        "depends": _pxi_dep["index"],
    },
    "_libs.indexing": {"pyxfile": "_libs/indexing"},
    "_libs.internals": {"pyxfile": "_libs/internals"},
    "_libs.interval": {
        "pyxfile": "_libs/interval",
        "include": klib_include,
        "depends": _pxi_dep["interval"],
    },
    "_libs.join": {"pyxfile": "_libs/join", "include": klib_include},
    "_libs.lib": {
        "pyxfile": "_libs/lib",
        "depends": lib_depends + tseries_depends,
        "include": klib_include,  # due to tokenizer import
        "sources": ["pandas/_libs/src/parser/tokenizer.c"],
    },
    "_libs.missing": {"pyxfile": "_libs/missing", "depends": tseries_depends},
    "_libs.parsers": {
        "pyxfile": "_libs/parsers",
        "include": klib_include + ["pandas/_libs/src"],
        "depends": [
            "pandas/_libs/src/parser/tokenizer.h",
            "pandas/_libs/src/parser/io.h",
        ],
        "sources": [
            "pandas/_libs/src/parser/tokenizer.c",
            "pandas/_libs/src/parser/io.c",
        ],
    },
    "_libs.reduction": {"pyxfile": "_libs/reduction"},
    "_libs.ops": {"pyxfile": "_libs/ops"},
    "_libs.ops_dispatch": {"pyxfile": "_libs/ops_dispatch"},
    "_libs.properties": {"pyxfile": "_libs/properties"},
    "_libs.reshape": {"pyxfile": "_libs/reshape", "depends": []},
    "_libs.sparse": {"pyxfile": "_libs/sparse", "depends": _pxi_dep["sparse"]},
    "_libs.tslib": {"pyxfile": "_libs/tslib", "depends": tseries_depends},
    "_libs.tslibs.base": {"pyxfile": "_libs/tslibs/base"},
    "_libs.tslibs.ccalendar": {"pyxfile": "_libs/tslibs/ccalendar"},
    "_libs.tslibs.dtypes": {"pyxfile": "_libs/tslibs/dtypes"},
    "_libs.tslibs.conversion": {
        "pyxfile": "_libs/tslibs/conversion",
        "depends": tseries_depends,
        "sources": ["pandas/_libs/tslibs/src/datetime/np_datetime.c"],
    },
    "_libs.tslibs.fields": {
        "pyxfile": "_libs/tslibs/fields",
        "depends": tseries_depends,
    },
    "_libs.tslibs.nattype": {"pyxfile": "_libs/tslibs/nattype"},
    "_libs.tslibs.np_datetime": {
        "pyxfile": "_libs/tslibs/np_datetime",
        "depends": tseries_depends,
        "sources": [
            "pandas/_libs/tslibs/src/datetime/np_datetime.c",
            "pandas/_libs/tslibs/src/datetime/np_datetime_strings.c",
        ],
    },
    "_libs.tslibs.offsets": {
        "pyxfile": "_libs/tslibs/offsets",
        "depends": tseries_depends,
    },
    "_libs.tslibs.parsing": {
        "pyxfile": "_libs/tslibs/parsing",
        "include": klib_include,
        "depends": ["pandas/_libs/src/parser/tokenizer.h"],
        "sources": ["pandas/_libs/src/parser/tokenizer.c"],
    },
    "_libs.tslibs.period": {
        "pyxfile": "_libs/tslibs/period",
        "depends": tseries_depends,
        "sources": ["pandas/_libs/tslibs/src/datetime/np_datetime.c"],
    },
    "_libs.tslibs.strptime": {
        "pyxfile": "_libs/tslibs/strptime",
        "depends": tseries_depends,
    },
    "_libs.tslibs.timedeltas": {
        "pyxfile": "_libs/tslibs/timedeltas",
        "depends": tseries_depends,
    },
    "_libs.tslibs.timestamps": {
        "pyxfile": "_libs/tslibs/timestamps",
        "depends": tseries_depends,
    },
    "_libs.tslibs.timezones": {"pyxfile": "_libs/tslibs/timezones"},
    "_libs.tslibs.tzconversion": {
        "pyxfile": "_libs/tslibs/tzconversion",
        "depends": tseries_depends,
    },
    "_libs.tslibs.vectorized": {"pyxfile": "_libs/tslibs/vectorized"},
    "_libs.testing": {"pyxfile": "_libs/testing"},
    "_libs.window.aggregations": {
        "pyxfile": "_libs/window/aggregations",
        "language": "c++",
        "suffix": ".cpp",
        "depends": ["pandas/_libs/src/skiplist.h"],
    },
    "_libs.window.indexers": {"pyxfile": "_libs/window/indexers"},
    "_libs.writers": {"pyxfile": "_libs/writers"},
    "io.sas._sas": {"pyxfile": "io/sas/sas"},
}

extensions = []

for name, data in ext_data.items():
    source_suffix = suffix if suffix == ".pyx" else data.get("suffix", ".c")

    sources = [srcpath(data["pyxfile"], suffix=source_suffix, subdir="")]

    sources.extend(data.get("sources", []))

    include = data.get("include", [])
    include.append(numpy.get_include())

    undef_macros = []

    if (
        sys.platform == "zos"
        and data.get("language") == "c++"
        and os.path.basename(os.environ.get("CXX", "/bin/xlc++")) in ("xlc", "xlc++")
    ):
        data.get("macros", macros).append(("__s390__", "1"))
        extra_compile_args.append("-qlanglvl=extended0x:nolibext")
        undef_macros.append("_POSIX_THREADS")

    obj = Extension(
        f"pandas.{name}",
        sources=sources,
        depends=data.get("depends", []),
        include_dirs=include,
        language=data.get("language", "c"),
        define_macros=data.get("macros", macros),
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
        undef_macros=undef_macros,
    )

    extensions.append(obj)

# ----------------------------------------------------------------------
# ujson

if suffix == ".pyx":
    # undo dumb setuptools bug clobbering .pyx sources back to .c
    for ext in extensions:
        if ext.sources[0].endswith((".c", ".cpp")):
            root, _ = os.path.splitext(ext.sources[0])
            ext.sources[0] = root + suffix

ujson_ext = Extension(
    "pandas._libs.json",
    depends=[
        "pandas/_libs/src/ujson/lib/ultrajson.h",
        "pandas/_libs/src/ujson/python/date_conversions.h",
    ],
    sources=(
        [
            "pandas/_libs/src/ujson/python/ujson.c",
            "pandas/_libs/src/ujson/python/objToJSON.c",
            "pandas/_libs/src/ujson/python/date_conversions.c",
            "pandas/_libs/src/ujson/python/JSONtoObj.c",
            "pandas/_libs/src/ujson/lib/ultrajsonenc.c",
            "pandas/_libs/src/ujson/lib/ultrajsondec.c",
        ]
        + [
            "pandas/_libs/tslibs/src/datetime/np_datetime.c",
            "pandas/_libs/tslibs/src/datetime/np_datetime_strings.c",
        ]
    ),
    include_dirs=[
        "pandas/_libs/src/ujson/python",
        "pandas/_libs/src/ujson/lib",
        "pandas/_libs/src/datetime",
        numpy.get_include(),
    ],
    extra_compile_args=(["-D_GNU_SOURCE"] + extra_compile_args),
    extra_link_args=extra_link_args,
    define_macros=macros,
)


extensions.append(ujson_ext)

# ----------------------------------------------------------------------


if __name__ == "__main__":
    # Freeze to support parallel compilation when using spawn instead of fork
    multiprocessing.freeze_support()
    setup(
        version=versioneer.get_version(),
        ext_modules=maybe_cythonize(extensions, compiler_directives=directives),
        cmdclass=cmdclass,
    )
