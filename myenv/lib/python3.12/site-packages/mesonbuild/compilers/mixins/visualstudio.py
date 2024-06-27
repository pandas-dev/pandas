# Copyright 2019 The meson development team
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
from __future__ import annotations

"""Abstractions to simplify compilers that implement an MSVC compatible
interface.
"""

import abc
import os
import typing as T

from ... import arglist
from ... import mesonlib
from ... import mlog

if T.TYPE_CHECKING:
    from ...environment import Environment
    from ...dependencies import Dependency
    from .clike import CLikeCompiler as Compiler
else:
    # This is a bit clever, for mypy we pretend that these mixins descend from
    # Compiler, so we get all of the methods and attributes defined for us, but
    # for runtime we make them descend from object (which all classes normally
    # do). This gives up DRYer type checking, with no runtime impact
    Compiler = object

vs32_instruction_set_args = {
    'mmx': ['/arch:SSE'], # There does not seem to be a flag just for MMX
    'sse': ['/arch:SSE'],
    'sse2': ['/arch:SSE2'],
    'sse3': ['/arch:AVX'], # VS leaped from SSE2 directly to AVX.
    'sse41': ['/arch:AVX'],
    'sse42': ['/arch:AVX'],
    'avx': ['/arch:AVX'],
    'avx2': ['/arch:AVX2'],
    'neon': None,
}  # T.Dicst[str, T.Optional[T.List[str]]]

# The 64 bit compiler defaults to /arch:avx.
vs64_instruction_set_args = {
    'mmx': ['/arch:AVX'],
    'sse': ['/arch:AVX'],
    'sse2': ['/arch:AVX'],
    'sse3': ['/arch:AVX'],
    'ssse3': ['/arch:AVX'],
    'sse41': ['/arch:AVX'],
    'sse42': ['/arch:AVX'],
    'avx': ['/arch:AVX'],
    'avx2': ['/arch:AVX2'],
    'neon': None,
}  # T.Dicst[str, T.Optional[T.List[str]]]

msvc_optimization_args = {
    'plain': [],
    '0': ['/Od'],
    'g': [], # No specific flag to optimize debugging, /Zi or /ZI will create debug information
    '1': ['/O1'],
    '2': ['/O2'],
    '3': ['/O2', '/Gw'],
    's': ['/O1', '/Gw'],
}  # type: T.Dict[str, T.List[str]]

msvc_debug_args = {
    False: [],
    True: ['/Zi']
}  # type: T.Dict[bool, T.List[str]]


class VisualStudioLikeCompiler(Compiler, metaclass=abc.ABCMeta):

    """A common interface for all compilers implementing an MSVC-style
    interface.

    A number of compilers attempt to mimic MSVC, with varying levels of
    success, such as Clang-CL and ICL (the Intel C/C++ Compiler for Windows).
    This class implements as much common logic as possible.
    """

    std_warn_args = ['/W3']
    std_opt_args = ['/O2']
    ignore_libs = arglist.UNIXY_COMPILER_INTERNAL_LIBS + ['execinfo']
    internal_libs = []  # type: T.List[str]

    crt_args = {
        'none': [],
        'md': ['/MD'],
        'mdd': ['/MDd'],
        'mt': ['/MT'],
        'mtd': ['/MTd'],
    }  # type: T.Dict[str, T.List[str]]

    # /showIncludes is needed for build dependency tracking in Ninja
    # See: https://ninja-build.org/manual.html#_deps
    # Assume UTF-8 sources by default, but self.unix_args_to_native() removes it
    # if `/source-charset` is set too.
    # It is also dropped if Visual Studio 2013 or earlier is used, since it would
    # not be supported in that case.
    always_args = ['/nologo', '/showIncludes', '/utf-8']
    warn_args = {
        '0': [],
        '1': ['/W2'],
        '2': ['/W3'],
        '3': ['/W4'],
        'everything': ['/Wall'],
    }  # type: T.Dict[str, T.List[str]]

    INVOKES_LINKER = False

    def __init__(self, target: str):
        self.base_options = {mesonlib.OptionKey(o) for o in ['b_pch', 'b_ndebug', 'b_vscrt']} # FIXME add lto, pgo and the like
        self.target = target
        self.is_64 = ('x64' in target) or ('x86_64' in target)
        # do some canonicalization of target machine
        if 'x86_64' in target:
            self.machine = 'x64'
        elif '86' in target:
            self.machine = 'x86'
        elif 'aarch64' in target:
            self.machine = 'arm64'
        elif 'arm' in target:
            self.machine = 'arm'
        else:
            self.machine = target
        if mesonlib.version_compare(self.version, '>=19.28.29910'): # VS 16.9.0 includes cl 19.28.29910
            self.base_options.add(mesonlib.OptionKey('b_sanitize'))
        assert self.linker is not None
        self.linker.machine = self.machine

    # Override CCompiler.get_always_args
    def get_always_args(self) -> T.List[str]:
        # TODO: use ImmutableListProtocol[str] here instead
        return self.always_args.copy()

    def get_pch_suffix(self) -> str:
        return 'pch'

    def get_pch_name(self, header: str) -> str:
        chopped = os.path.basename(header).split('.')[:-1]
        chopped.append(self.get_pch_suffix())
        pchname = '.'.join(chopped)
        return pchname

    def get_pch_base_name(self, header: str) -> str:
        # This needs to be implemented by inheriting classes
        raise NotImplementedError

    def get_pch_use_args(self, pch_dir: str, header: str) -> T.List[str]:
        base = self.get_pch_base_name(header)
        pchname = self.get_pch_name(header)
        return ['/FI' + base, '/Yu' + base, '/Fp' + os.path.join(pch_dir, pchname)]

    def get_preprocess_only_args(self) -> T.List[str]:
        return ['/EP']

    def get_preprocess_to_file_args(self) -> T.List[str]:
        return ['/EP', '/P']

    def get_compile_only_args(self) -> T.List[str]:
        return ['/c']

    def get_no_optimization_args(self) -> T.List[str]:
        return ['/Od', '/Oi-']

    def sanitizer_compile_args(self, value: str) -> T.List[str]:
        if value == 'none':
            return []
        if value != 'address':
            raise mesonlib.MesonException('VS only supports address sanitizer at the moment.')
        return ['/fsanitize=address']

    def get_output_args(self, target: str) -> T.List[str]:
        if self.mode == 'PREPROCESSOR':
            return ['/Fi' + target]
        if target.endswith('.exe'):
            return ['/Fe' + target]
        return ['/Fo' + target]

    def get_buildtype_args(self, buildtype: str) -> T.List[str]:
        return []

    def get_debug_args(self, is_debug: bool) -> T.List[str]:
        return msvc_debug_args[is_debug]

    def get_optimization_args(self, optimization_level: str) -> T.List[str]:
        args = msvc_optimization_args[optimization_level]
        if mesonlib.version_compare(self.version, '<18.0'):
            args = [arg for arg in args if arg != '/Gw']
        return args

    def linker_to_compiler_args(self, args: T.List[str]) -> T.List[str]:
        return ['/link'] + args

    def get_pic_args(self) -> T.List[str]:
        return [] # PIC is handled by the loader on Windows

    def gen_vs_module_defs_args(self, defsfile: str) -> T.List[str]:
        if not isinstance(defsfile, str):
            raise RuntimeError('Module definitions file should be str')
        # With MSVC, DLLs only export symbols that are explicitly exported,
        # so if a module defs file is specified, we use that to export symbols
        return ['/DEF:' + defsfile]

    def gen_pch_args(self, header: str, source: str, pchname: str) -> T.Tuple[str, T.List[str]]:
        objname = os.path.splitext(source)[0] + '.obj'
        return objname, ['/Yc' + header, '/Fp' + pchname, '/Fo' + objname]

    def openmp_flags(self) -> T.List[str]:
        return ['/openmp']

    def openmp_link_flags(self) -> T.List[str]:
        return []

    # FIXME, no idea what these should be.
    def thread_flags(self, env: 'Environment') -> T.List[str]:
        return []

    @classmethod
    def unix_args_to_native(cls, args: T.List[str]) -> T.List[str]:
        result: T.List[str] = []
        for i in args:
            # -mms-bitfields is specific to MinGW-GCC
            # -pthread is only valid for GCC
            if i in {'-mms-bitfields', '-pthread'}:
                continue
            if i.startswith('-LIBPATH:'):
                i = '/LIBPATH:' + i[9:]
            elif i.startswith('-L'):
                i = '/LIBPATH:' + i[2:]
            # Translate GNU-style -lfoo library name to the import library
            elif i.startswith('-l'):
                name = i[2:]
                if name in cls.ignore_libs:
                    # With MSVC, these are provided by the C runtime which is
                    # linked in by default
                    continue
                else:
                    i = name + '.lib'
            elif i.startswith('-isystem'):
                # just use /I for -isystem system include path s
                if i.startswith('-isystem='):
                    i = '/I' + i[9:]
                else:
                    i = '/I' + i[8:]
            elif i.startswith('-idirafter'):
                # same as -isystem, but appends the path instead
                if i.startswith('-idirafter='):
                    i = '/I' + i[11:]
                else:
                    i = '/I' + i[10:]
            # -pthread in link flags is only used on Linux
            elif i == '-pthread':
                continue
            # cl.exe does not allow specifying both, so remove /utf-8 that we
            # added automatically in the case the user overrides it manually.
            elif (i.startswith('/source-charset:')
                    or i.startswith('/execution-charset:')
                    or i == '/validate-charset-'):
                try:
                    result.remove('/utf-8')
                except ValueError:
                    pass
            result.append(i)
        return result

    @classmethod
    def native_args_to_unix(cls, args: T.List[str]) -> T.List[str]:
        result = []
        for arg in args:
            if arg.startswith(('/LIBPATH:', '-LIBPATH:')):
                result.append('-L' + arg[9:])
            elif arg.endswith(('.a', '.lib')) and not os.path.isabs(arg):
                result.append('-l' + arg)
            else:
                result.append(arg)
        return result

    def get_werror_args(self) -> T.List[str]:
        return ['/WX']

    def get_include_args(self, path: str, is_system: bool) -> T.List[str]:
        if path == '':
            path = '.'
        # msvc does not have a concept of system header dirs.
        return ['-I' + path]

    def compute_parameters_with_absolute_paths(self, parameter_list: T.List[str], build_dir: str) -> T.List[str]:
        for idx, i in enumerate(parameter_list):
            if i[:2] == '-I' or i[:2] == '/I':
                parameter_list[idx] = i[:2] + os.path.normpath(os.path.join(build_dir, i[2:]))
            elif i[:9] == '/LIBPATH:':
                parameter_list[idx] = i[:9] + os.path.normpath(os.path.join(build_dir, i[9:]))

        return parameter_list

    # Visual Studio is special. It ignores some arguments it does not
    # understand and you can't tell it to error out on those.
    # http://stackoverflow.com/questions/15259720/how-can-i-make-the-microsoft-c-compiler-treat-unknown-flags-as-errors-rather-t
    def has_arguments(self, args: T.List[str], env: 'Environment', code: str, mode: str) -> T.Tuple[bool, bool]:
        warning_text = '4044' if mode == 'link' else '9002'
        with self._build_wrapper(code, env, extra_args=args, mode=mode) as p:
            if p.returncode != 0:
                return False, p.cached
            return not (warning_text in p.stderr or warning_text in p.stdout), p.cached

    def get_compile_debugfile_args(self, rel_obj: str, pch: bool = False) -> T.List[str]:
        pdbarr = rel_obj.split('.')[:-1]
        pdbarr += ['pdb']
        args = ['/Fd' + '.'.join(pdbarr)]
        return args

    def get_instruction_set_args(self, instruction_set: str) -> T.Optional[T.List[str]]:
        if self.is_64:
            return vs64_instruction_set_args.get(instruction_set, None)
        return vs32_instruction_set_args.get(instruction_set, None)

    def _calculate_toolset_version(self, version: int) -> T.Optional[str]:
        if version < 1310:
            return '7.0'
        elif version < 1400:
            return '7.1' # (Visual Studio 2003)
        elif version < 1500:
            return '8.0' # (Visual Studio 2005)
        elif version < 1600:
            return '9.0' # (Visual Studio 2008)
        elif version < 1700:
            return '10.0' # (Visual Studio 2010)
        elif version < 1800:
            return '11.0' # (Visual Studio 2012)
        elif version < 1900:
            return '12.0' # (Visual Studio 2013)
        elif version < 1910:
            return '14.0' # (Visual Studio 2015)
        elif version < 1920:
            return '14.1' # (Visual Studio 2017)
        elif version < 1930:
            return '14.2' # (Visual Studio 2019)
        elif version < 1940:
            return '14.3' # (Visual Studio 2022)
        mlog.warning(f'Could not find toolset for version {self.version!r}')
        return None

    def get_toolset_version(self) -> T.Optional[str]:
        # See boost/config/compiler/visualc.cpp for up to date mapping
        try:
            version = int(''.join(self.version.split('.')[0:2]))
        except ValueError:
            return None
        return self._calculate_toolset_version(version)

    def get_default_include_dirs(self) -> T.List[str]:
        if 'INCLUDE' not in os.environ:
            return []
        return os.environ['INCLUDE'].split(os.pathsep)

    def get_crt_compile_args(self, crt_val: str, buildtype: str) -> T.List[str]:
        if crt_val in self.crt_args:
            return self.crt_args[crt_val]
        assert crt_val in {'from_buildtype', 'static_from_buildtype'}
        dbg = 'mdd'
        rel = 'md'
        if crt_val == 'static_from_buildtype':
            dbg = 'mtd'
            rel = 'mt'
        # Match what build type flags used to do.
        if buildtype == 'plain':
            return []
        elif buildtype == 'debug':
            return self.crt_args[dbg]
        elif buildtype == 'debugoptimized':
            return self.crt_args[rel]
        elif buildtype == 'release':
            return self.crt_args[rel]
        elif buildtype == 'minsize':
            return self.crt_args[rel]
        else:
            assert buildtype == 'custom'
            raise mesonlib.EnvironmentException('Requested C runtime based on buildtype, but buildtype is "custom".')

    def has_func_attribute(self, name: str, env: 'Environment') -> T.Tuple[bool, bool]:
        # MSVC doesn't have __attribute__ like Clang and GCC do, so just return
        # false without compiling anything
        return name in {'dllimport', 'dllexport'}, False

    def get_argument_syntax(self) -> str:
        return 'msvc'

    def symbols_have_underscore_prefix(self, env: 'Environment') -> bool:
        '''
        Check if the compiler prefixes an underscore to global C symbols.

        This overrides the Clike method, as for MSVC checking the
        underscore prefix based on the compiler define never works,
        so do not even try.
        '''
        # Try to consult a hardcoded list of cases we know
        # absolutely have an underscore prefix
        result = self._symbols_have_underscore_prefix_list(env)
        if result is not None:
            return result

        # As a last resort, try search in a compiled binary
        return self._symbols_have_underscore_prefix_searchbin(env)


class MSVCCompiler(VisualStudioLikeCompiler):

    """Specific to the Microsoft Compilers."""

    id = 'msvc'

    def __init__(self, target: str):
        super().__init__(target)

        # Visual Studio 2013 and earlier don't support the /utf-8 argument.
        # We want to remove it. We also want to make an explicit copy so we
        # don't mutate class constant state
        if mesonlib.version_compare(self.version, '<19.00') and '/utf-8' in self.always_args:
            self.always_args = [r for r in self.always_args if r != '/utf-8']

    def get_compile_debugfile_args(self, rel_obj: str, pch: bool = False) -> T.List[str]:
        args = super().get_compile_debugfile_args(rel_obj, pch)
        # When generating a PDB file with PCH, all compile commands write
        # to the same PDB file. Hence, we need to serialize the PDB
        # writes using /FS since we do parallel builds. This slows down the
        # build obviously, which is why we only do this when PCH is on.
        # This was added in Visual Studio 2013 (MSVC 18.0). Before that it was
        # always on: https://msdn.microsoft.com/en-us/library/dn502518.aspx
        if pch and mesonlib.version_compare(self.version, '>=18.0'):
            args = ['/FS'] + args
        return args

    # Override CCompiler.get_always_args
    # We want to drop '/utf-8' for Visual Studio 2013 and earlier
    def get_always_args(self) -> T.List[str]:
        return self.always_args

    def get_instruction_set_args(self, instruction_set: str) -> T.Optional[T.List[str]]:
        if self.version.split('.')[0] == '16' and instruction_set == 'avx':
            # VS documentation says that this exists and should work, but
            # it does not. The headers do not contain AVX intrinsics
            # and they cannot be called.
            return None
        return super().get_instruction_set_args(instruction_set)

    def get_pch_base_name(self, header: str) -> str:
        return os.path.basename(header)


class ClangClCompiler(VisualStudioLikeCompiler):

    """Specific to Clang-CL."""

    id = 'clang-cl'

    def __init__(self, target: str):
        super().__init__(target)

        # Assembly
        self.can_compile_suffixes.add('s')
        self.can_compile_suffixes.add('sx')

    def has_arguments(self, args: T.List[str], env: 'Environment', code: str, mode: str) -> T.Tuple[bool, bool]:
        if mode != 'link':
            args = args + ['-Werror=unknown-argument', '-Werror=unknown-warning-option']
        return super().has_arguments(args, env, code, mode)

    def get_toolset_version(self) -> T.Optional[str]:
        # XXX: what is the right thing to do here?
        return '14.1'

    def get_pch_base_name(self, header: str) -> str:
        return header

    def get_include_args(self, path: str, is_system: bool) -> T.List[str]:
        if path == '':
            path = '.'
        return ['/clang:-isystem' + path] if is_system else ['-I' + path]

    def get_dependency_compile_args(self, dep: 'Dependency') -> T.List[str]:
        if dep.get_include_type() == 'system':
            converted = []
            for i in dep.get_compile_args():
                if i.startswith('-isystem'):
                    converted += ['/clang:' + i]
                else:
                    converted += [i]
            return converted
        else:
            return dep.get_compile_args()
