from setuptools import distutils as dutils
from setuptools.command import build_ext
from setuptools.extension import Extension

import os
import shutil
import sys
import tempfile

from numba.core import typing, sigutils
from numba.core.compiler_lock import global_compiler_lock
from numba.pycc.compiler import ModuleCompiler, ExportEntry
from numba.pycc.platform import Toolchain
from numba import cext


dir_util = dutils.dir_util
log = dutils.log
extension_libs = cext.get_extension_libs()


class CC(object):
    """
    An ahead-of-time compiler to create extension modules that don't
    depend on Numba.
    """

    # NOTE: using ccache can speed up repetitive builds
    # (especially for the mixin modules)

    _mixin_sources = ['modulemixin.c',]  + extension_libs

    # -flto strips all unused helper functions, which 1) makes the
    # produced output much smaller and 2) can make the linking step faster.
    # (the Windows linker seems to do this by default, judging by the results)

    _extra_cflags = {
        # Comment out due to odd behavior with GCC 4.9+ with LTO
        # 'posix': ['-flto'],
        }

    _extra_ldflags = {
        # Comment out due to odd behavior with GCC 4.9+ with LTO
        # 'posix': ['-flto'],
        }

    def __init__(self, extension_name, source_module=None):
        if '.' in extension_name:
            raise ValueError("basename should be a simple module name, not "
                             "qualified name")

        self._basename = extension_name
        self._init_function = 'pycc_init_' + extension_name
        self._exported_functions = {}
        # Resolve source module name and directory
        f = sys._getframe(1)
        if source_module is None:
            dct = f.f_globals
            source_module = dct['__name__']
        elif hasattr(source_module, '__name__'):
            dct = source_module.__dict__
            source_module = source_module.__name__
        else:
            dct = sys.modules[source_module].__dict__

        self._source_path = dct.get('__file__', '')
        self._source_module = source_module
        self._toolchain = Toolchain()
        self._verbose = False
        # By default, output in directory of caller module
        self._output_dir = os.path.dirname(self._source_path)
        self._output_file = self._toolchain.get_ext_filename(extension_name)
        self._use_nrt = True
        self._target_cpu = ''

    @property
    def name(self):
        """
        The name of the extension module to create.
        """
        return self._basename

    @property
    def output_file(self):
        """
        The specific output file (a DLL) that will be generated.
        """
        return self._output_file

    @output_file.setter
    def output_file(self, value):
        self._output_file = value

    @property
    def output_dir(self):
        """
        The directory the output file will be put in.
        """
        return self._output_dir

    @output_dir.setter
    def output_dir(self, value):
        self._output_dir = value

    @property
    def use_nrt(self):
        return self._use_nrt

    @use_nrt.setter
    def use_nrt(self, value):
        self._use_nrt = value

    @property
    def target_cpu(self):
        """
        The target CPU model for code generation.
        """
        return self._target_cpu

    @target_cpu.setter
    def target_cpu(self, value):
        self._target_cpu = value

    @property
    def verbose(self):
        """
        Whether to display detailed information when compiling.
        """
        return self._verbose

    @verbose.setter
    def verbose(self, value):
        self._verbose = value

    def export(self, exported_name, sig):
        """
        Mark a function for exporting in the extension module.
        """
        fn_args, fn_retty = sigutils.normalize_signature(sig)
        sig = typing.signature(fn_retty, *fn_args)
        if exported_name in self._exported_functions:
            raise KeyError("duplicated export symbol %s" % (exported_name))

        def decorator(func):
            entry = ExportEntry(exported_name, sig, func)
            self._exported_functions[exported_name] = entry
            return func

        return decorator

    @property
    def _export_entries(self):
        return sorted(self._exported_functions.values(),
                      key=lambda entry: entry.symbol)

    def _get_mixin_sources(self):
        here = os.path.dirname(__file__)
        mixin_sources = self._mixin_sources[:]
        if self._use_nrt:
            mixin_sources.append('../core/runtime/nrt.cpp')
        return [os.path.join(here, f) for f in mixin_sources]

    def _get_mixin_defines(self):
        # Macro definitions required by modulemixin.c
        return [
            ('PYCC_MODULE_NAME', self._basename),
            ('PYCC_USE_NRT', int(self._use_nrt)),
            ]

    def _get_extra_cflags(self):
        extra_cflags = self._extra_cflags.get(sys.platform, [])
        if not extra_cflags:
            extra_cflags = self._extra_cflags.get(os.name, [])
        return extra_cflags

    def _get_extra_ldflags(self):
        extra_ldflags = self._extra_ldflags.get(sys.platform, [])
        if not extra_ldflags:
            extra_ldflags = self._extra_ldflags.get(os.name, [])
        # helperlib uses pthread on linux. make sure we are linking to it.
        if sys.platform.startswith("linux"):
            if "-pthread" not in extra_ldflags:
                extra_ldflags.append('-pthread')
        return extra_ldflags

    def _compile_mixins(self, build_dir):
        sources = self._get_mixin_sources()
        macros = self._get_mixin_defines()
        include_dirs = self._toolchain.get_python_include_dirs()

        extra_cflags = self._get_extra_cflags()
        # XXX distutils creates a whole subtree inside build_dir,
        # e.g. /tmp/test_pycc/home/antoine/numba/numba/pycc/modulemixin.o
        objects = self._toolchain.compile_objects(sources, build_dir,
                                                  include_dirs=include_dirs,
                                                  macros=macros,
                                                  extra_cflags=extra_cflags)
        return objects

    @global_compiler_lock
    def _compile_object_files(self, build_dir):
        compiler = ModuleCompiler(self._export_entries, self._basename,
                                self._use_nrt, cpu_name=self._target_cpu)
        compiler.external_init_function = self._init_function
        temp_obj = os.path.join(build_dir,
                                os.path.splitext(self._output_file)[0] + '.o')
        log.info("generating LLVM code for '%s' into %s",
                self._basename, temp_obj)
        compiler.write_native_object(temp_obj, wrap=True)
        return [temp_obj], compiler.dll_exports

    @global_compiler_lock
    def compile(self):
        """
        Compile the extension module.
        """
        self._toolchain.verbose = self.verbose
        build_dir = tempfile.mkdtemp(prefix='pycc-build-%s-' % self._basename)

        # Compile object file
        objects, dll_exports = self._compile_object_files(build_dir)

        # Compile mixins
        objects += self._compile_mixins(build_dir)

        # Then create shared library
        extra_ldflags = self._get_extra_ldflags()
        output_dll = os.path.join(self._output_dir, self._output_file)
        libraries = self._toolchain.get_python_libraries()
        library_dirs = self._toolchain.get_python_library_dirs()
        self._toolchain.link_shared(output_dll, objects,
                                    libraries, library_dirs,
                                    export_symbols=dll_exports,
                                    extra_ldflags=extra_ldflags)

        shutil.rmtree(build_dir)

    def distutils_extension(self, **kwargs):
        """
        Create a distutils extension object that can be used in your
        setup.py.
        """
        macros = kwargs.pop('macros', []) + self._get_mixin_defines()
        depends = kwargs.pop('depends', []) + [self._source_path]
        extra_compile_args = (kwargs.pop('extra_compile_args', [])
                              + self._get_extra_cflags())
        extra_link_args = (kwargs.pop('extra_link_args', [])
                           + self._get_extra_ldflags())
        include_dirs = (kwargs.pop('include_dirs', [])
                        + self._toolchain.get_python_include_dirs())
        libraries = (kwargs.pop('libraries', [])
                     + self._toolchain.get_python_libraries())
        library_dirs = (kwargs.pop('library_dirs', [])
                        + self._toolchain.get_python_library_dirs())
        python_package_path = self._source_module[:self._source_module.rfind('.')+1]

        ext = _CCExtension(name=python_package_path + self._basename,
                           sources=self._get_mixin_sources(),
                           depends=depends,
                           define_macros=macros,
                           include_dirs=include_dirs,
                           libraries=libraries,
                           library_dirs=library_dirs,
                           extra_compile_args=extra_compile_args,
                           extra_link_args=extra_link_args,
                           **kwargs)
        ext.monkey_patch_distutils()
        ext._cc = self
        return ext


class _CCExtension(Extension):
    """
    A Numba-specific Extension subclass to LLVM-compile pure Python code
    to an extension module.
    """

    _cc = None
    _distutils_monkey_patched = False

    def _prepare_object_files(self, build_ext):
        cc = self._cc
        dir_util.mkpath(os.path.join(build_ext.build_temp, *self.name.split('.')[:-1]))
        objects, _ = cc._compile_object_files(build_ext.build_temp)
        # Add generated object files for linking
        self.extra_objects = objects

    @classmethod
    def monkey_patch_distutils(cls):
        """
        Monkey-patch distutils with our own build_ext class knowing
        about pycc-compiled extensions modules.
        """
        if cls._distutils_monkey_patched:
            return

        _orig_build_ext = build_ext.build_ext

        class _CC_build_ext(_orig_build_ext):

            def build_extension(self, ext):
                if isinstance(ext, _CCExtension):
                    ext._prepare_object_files(self)

                _orig_build_ext.build_extension(self, ext)

        build_ext.build_ext = _CC_build_ext

        cls._distutils_monkey_patched = True
