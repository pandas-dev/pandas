# SPDX-License-Identifier: Apache-2.0
# Copyright 2012-2021 The Meson development team
# Copyright © 2021 Intel Corporation
from __future__ import annotations

import collections
import enum
import functools
import os
import itertools
import typing as T

from .. import build
from .. import coredata
from .. import dependencies
from .. import mesonlib
from .. import mlog
from ..compilers import SUFFIX_TO_LANG
from ..compilers.compilers import CompileCheckMode
from ..interpreterbase import (ObjectHolder, noPosargs, noKwargs,
                               FeatureNew, disablerIfNotFound,
                               InterpreterException)
from ..interpreterbase.decorators import ContainerTypeInfo, typed_kwargs, KwargInfo, typed_pos_args
from ..mesonlib import OptionKey
from .interpreterobjects import (extract_required_kwarg, extract_search_dirs)
from .type_checking import REQUIRED_KW, in_set_validator, NoneType

if T.TYPE_CHECKING:
    from ..interpreter import Interpreter
    from ..compilers import Compiler, RunResult
    from ..interpreterbase import TYPE_var, TYPE_kwargs
    from .kwargs import ExtractRequired, ExtractSearchDirs
    from .interpreter.interpreter import SourceOutputs

    from typing_extensions import TypedDict, Literal

    class GetSupportedArgumentKw(TypedDict):

        checked: Literal['warn', 'require', 'off']

    class AlignmentKw(TypedDict):

        prefix: str
        args: T.List[str]
        dependencies: T.List[dependencies.Dependency]

    class CompileKW(TypedDict):

        name: str
        no_builtin_args: bool
        include_directories: T.List[build.IncludeDirs]
        args: T.List[str]
        dependencies: T.List[dependencies.Dependency]

    class CommonKW(TypedDict):

        prefix: str
        no_builtin_args: bool
        include_directories: T.List[build.IncludeDirs]
        args: T.List[str]
        dependencies: T.List[dependencies.Dependency]

    class ComputeIntKW(CommonKW):

        guess: T.Optional[int]
        high: T.Optional[int]
        low: T.Optional[int]

    class HeaderKW(CommonKW, ExtractRequired):
        pass

    class FindLibraryKW(ExtractRequired, ExtractSearchDirs):

        disabler: bool
        has_headers: T.List[str]
        static: bool

        # This list must be all of the `HeaderKW` values with `header_`
        # prepended to the key
        header_args: T.List[str]
        header_dependencies: T.List[dependencies.Dependency]
        header_include_directories: T.List[build.IncludeDirs]
        header_no_builtin_args: bool
        header_prefix: str
        header_required: T.Union[bool, coredata.UserFeatureOption]

    class PreprocessKW(TypedDict):
        output: str
        compile_args: T.List[str]
        include_directories: T.List[build.IncludeDirs]
        dependencies: T.List[dependencies.Dependency]


class _TestMode(enum.Enum):

    """Whether we're doing a compiler or linker check."""

    COMPILER = 0
    LINKER = 1


class TryRunResultHolder(ObjectHolder['RunResult']):
    def __init__(self, res: 'RunResult', interpreter: 'Interpreter'):
        super().__init__(res, interpreter)
        self.methods.update({'returncode': self.returncode_method,
                             'compiled': self.compiled_method,
                             'stdout': self.stdout_method,
                             'stderr': self.stderr_method,
                             })

    @noPosargs
    @noKwargs
    def returncode_method(self, args: T.List['TYPE_var'], kwargs: 'TYPE_kwargs') -> int:
        return self.held_object.returncode

    @noPosargs
    @noKwargs
    def compiled_method(self, args: T.List['TYPE_var'], kwargs: 'TYPE_kwargs') -> bool:
        return self.held_object.compiled

    @noPosargs
    @noKwargs
    def stdout_method(self, args: T.List['TYPE_var'], kwargs: 'TYPE_kwargs') -> str:
        return self.held_object.stdout

    @noPosargs
    @noKwargs
    def stderr_method(self, args: T.List['TYPE_var'], kwargs: 'TYPE_kwargs') -> str:
        return self.held_object.stderr


_ARGS_KW: KwargInfo[T.List[str]] = KwargInfo(
    'args',
    ContainerTypeInfo(list, str),
    listify=True,
    default=[],
)
_DEPENDENCIES_KW: KwargInfo[T.List['dependencies.Dependency']] = KwargInfo(
    'dependencies',
    ContainerTypeInfo(list, dependencies.Dependency),
    listify=True,
    default=[],
)
_INCLUDE_DIRS_KW: KwargInfo[T.List[build.IncludeDirs]] = KwargInfo(
    'include_directories',
    ContainerTypeInfo(list, build.IncludeDirs),
    default=[],
    listify=True,
)
_PREFIX_KW: KwargInfo[str] = KwargInfo(
    'prefix',
    (str, ContainerTypeInfo(list, str)),
    default='',
    since_values={list: '1.0.0'},
    convertor=lambda x: '\n'.join(x) if isinstance(x, list) else x)

_NO_BUILTIN_ARGS_KW = KwargInfo('no_builtin_args', bool, default=False)
_NAME_KW = KwargInfo('name', str, default='')

# Many of the compiler methods take this kwarg signature exactly, this allows
# simplifying the `typed_kwargs` calls
_COMMON_KWS: T.List[KwargInfo] = [_ARGS_KW, _DEPENDENCIES_KW, _INCLUDE_DIRS_KW, _PREFIX_KW, _NO_BUILTIN_ARGS_KW]

# Common methods of compiles, links, runs, and similar
_COMPILES_KWS: T.List[KwargInfo] = [_NAME_KW, _ARGS_KW, _DEPENDENCIES_KW, _INCLUDE_DIRS_KW, _NO_BUILTIN_ARGS_KW]

_HEADER_KWS: T.List[KwargInfo] = [REQUIRED_KW.evolve(since='0.50.0', default=False), *_COMMON_KWS]

class CompilerHolder(ObjectHolder['Compiler']):
    preprocess_uid: T.Dict[str, itertools.count] = collections.defaultdict(itertools.count)

    def __init__(self, compiler: 'Compiler', interpreter: 'Interpreter'):
        super().__init__(compiler, interpreter)
        self.environment = self.env
        self.methods.update({'compiles': self.compiles_method,
                             'links': self.links_method,
                             'get_id': self.get_id_method,
                             'get_linker_id': self.get_linker_id_method,
                             'compute_int': self.compute_int_method,
                             'sizeof': self.sizeof_method,
                             'get_define': self.get_define_method,
                             'check_header': self.check_header_method,
                             'has_header': self.has_header_method,
                             'has_header_symbol': self.has_header_symbol_method,
                             'run': self.run_method,
                             'has_function': self.has_function_method,
                             'has_member': self.has_member_method,
                             'has_members': self.has_members_method,
                             'has_type': self.has_type_method,
                             'alignment': self.alignment_method,
                             'version': self.version_method,
                             'cmd_array': self.cmd_array_method,
                             'find_library': self.find_library_method,
                             'has_argument': self.has_argument_method,
                             'has_function_attribute': self.has_func_attribute_method,
                             'get_supported_function_attributes': self.get_supported_function_attributes_method,
                             'has_multi_arguments': self.has_multi_arguments_method,
                             'get_supported_arguments': self.get_supported_arguments_method,
                             'first_supported_argument': self.first_supported_argument_method,
                             'has_link_argument': self.has_link_argument_method,
                             'has_multi_link_arguments': self.has_multi_link_arguments_method,
                             'get_supported_link_arguments': self.get_supported_link_arguments_method,
                             'first_supported_link_argument': self.first_supported_link_argument_method,
                             'symbols_have_underscore_prefix': self.symbols_have_underscore_prefix_method,
                             'get_argument_syntax': self.get_argument_syntax_method,
                             'preprocess': self.preprocess_method,
                             })

    @property
    def compiler(self) -> 'Compiler':
        return self.held_object

    def _dep_msg(self, deps: T.List['dependencies.Dependency'], compile_only: bool, endl: str) -> str:
        msg_single = 'with dependency {}'
        msg_many = 'with dependencies {}'
        names = []
        for d in deps:
            if isinstance(d, dependencies.InternalDependency):
                FeatureNew.single_use('compiler method "dependencies" kwarg with internal dep', '0.57.0', self.subproject,
                                      location=self.current_node)
                continue
            if isinstance(d, dependencies.ExternalLibrary):
                if compile_only:
                    continue
                name = '-l' + d.name
            else:
                name = d.name
            names.append(name)
        if not names:
            return endl
        tpl = msg_many if len(names) > 1 else msg_single
        if endl is None:
            endl = ''
        return tpl.format(', '.join(names)) + endl

    @noPosargs
    @noKwargs
    def version_method(self, args: T.List['TYPE_var'], kwargs: 'TYPE_kwargs') -> str:
        return self.compiler.version

    @noPosargs
    @noKwargs
    def cmd_array_method(self, args: T.List['TYPE_var'], kwargs: 'TYPE_kwargs') -> T.List[str]:
        return self.compiler.exelist

    def _determine_args(self, nobuiltins: bool,
                        incdirs: T.List[build.IncludeDirs],
                        extra_args: T.List[str],
                        mode: CompileCheckMode = CompileCheckMode.LINK) -> T.List[str]:
        args: T.List[str] = []
        for i in incdirs:
            for idir in i.to_string_list(self.environment.get_source_dir()):
                args.extend(self.compiler.get_include_args(idir, False))
        if not nobuiltins:
            opts = self.environment.coredata.options
            args += self.compiler.get_option_compile_args(opts)
            if mode is CompileCheckMode.LINK:
                args.extend(self.compiler.get_option_link_args(opts))
        args.extend(extra_args)
        return args

    def _determine_dependencies(self, deps: T.List['dependencies.Dependency'], compile_only: bool = False, endl: str = ':') -> T.Tuple[T.List['dependencies.Dependency'], str]:
        deps = dependencies.get_leaf_external_dependencies(deps)
        return deps, self._dep_msg(deps, compile_only, endl)

    @typed_pos_args('compiler.alignment', str)
    @typed_kwargs(
        'compiler.alignment',
        _PREFIX_KW,
        _ARGS_KW,
        _DEPENDENCIES_KW,
    )
    def alignment_method(self, args: T.Tuple[str], kwargs: 'AlignmentKw') -> int:
        typename = args[0]
        deps, msg = self._determine_dependencies(kwargs['dependencies'], compile_only=self.compiler.is_cross)
        result, cached = self.compiler.alignment(typename, kwargs['prefix'], self.environment,
                                                 extra_args=kwargs['args'],
                                                 dependencies=deps)
        cached_msg = mlog.blue('(cached)') if cached else ''
        mlog.log('Checking for alignment of',
                 mlog.bold(typename, True), msg, mlog.bold(str(result)), cached_msg)
        return result

    @typed_pos_args('compiler.run', (str, mesonlib.File))
    @typed_kwargs('compiler.run', *_COMPILES_KWS)
    def run_method(self, args: T.Tuple['mesonlib.FileOrString'], kwargs: 'CompileKW') -> 'RunResult':
        code = args[0]
        if isinstance(code, mesonlib.File):
            self.interpreter.add_build_def_file(code)
            code = mesonlib.File.from_absolute_file(
                code.rel_to_builddir(self.environment.source_dir))
        testname = kwargs['name']
        extra_args = functools.partial(self._determine_args, kwargs['no_builtin_args'], kwargs['include_directories'], kwargs['args'])
        deps, msg = self._determine_dependencies(kwargs['dependencies'], compile_only=False, endl=None)
        result = self.compiler.run(code, self.environment, extra_args=extra_args,
                                   dependencies=deps)
        if testname:
            if not result.compiled:
                h = mlog.red('DID NOT COMPILE')
            elif result.returncode == 0:
                h = mlog.green('YES')
            else:
                h = mlog.red(f'NO ({result.returncode})')
            mlog.log('Checking if', mlog.bold(testname, True), msg, 'runs:', h)
        return result

    @noPosargs
    @noKwargs
    def get_id_method(self, args: T.List['TYPE_var'], kwargs: 'TYPE_kwargs') -> str:
        return self.compiler.get_id()

    @noPosargs
    @noKwargs
    @FeatureNew('compiler.get_linker_id', '0.53.0')
    def get_linker_id_method(self, args: T.List['TYPE_var'], kwargs: 'TYPE_kwargs') -> str:
        return self.compiler.get_linker_id()

    @noPosargs
    @noKwargs
    def symbols_have_underscore_prefix_method(self, args: T.List['TYPE_var'], kwargs: 'TYPE_kwargs') -> bool:
        '''
        Check if the compiler prefixes _ (underscore) to global C symbols
        See: https://en.wikipedia.org/wiki/Name_mangling#C
        '''
        return self.compiler.symbols_have_underscore_prefix(self.environment)

    @typed_pos_args('compiler.has_member', str, str)
    @typed_kwargs('compiler.has_member', *_COMMON_KWS)
    def has_member_method(self, args: T.Tuple[str, str], kwargs: 'CommonKW') -> bool:
        typename, membername = args
        extra_args = functools.partial(self._determine_args, kwargs['no_builtin_args'], kwargs['include_directories'], kwargs['args'])
        deps, msg = self._determine_dependencies(kwargs['dependencies'])
        had, cached = self.compiler.has_members(typename, [membername], kwargs['prefix'],
                                                self.environment,
                                                extra_args=extra_args,
                                                dependencies=deps)
        cached_msg = mlog.blue('(cached)') if cached else ''
        if had:
            hadtxt = mlog.green('YES')
        else:
            hadtxt = mlog.red('NO')
        mlog.log('Checking whether type', mlog.bold(typename, True),
                 'has member', mlog.bold(membername, True), msg, hadtxt, cached_msg)
        return had

    @typed_pos_args('compiler.has_members', str, varargs=str, min_varargs=1)
    @typed_kwargs('compiler.has_members', *_COMMON_KWS)
    def has_members_method(self, args: T.Tuple[str, T.List[str]], kwargs: 'CommonKW') -> bool:
        typename, membernames = args
        extra_args = functools.partial(self._determine_args, kwargs['no_builtin_args'], kwargs['include_directories'], kwargs['args'])
        deps, msg = self._determine_dependencies(kwargs['dependencies'])
        had, cached = self.compiler.has_members(typename, membernames, kwargs['prefix'],
                                                self.environment,
                                                extra_args=extra_args,
                                                dependencies=deps)
        cached_msg = mlog.blue('(cached)') if cached else ''
        if had:
            hadtxt = mlog.green('YES')
        else:
            hadtxt = mlog.red('NO')
        members = mlog.bold(', '.join([f'"{m}"' for m in membernames]))
        mlog.log('Checking whether type', mlog.bold(typename, True),
                 'has members', members, msg, hadtxt, cached_msg)
        return had

    @typed_pos_args('compiler.has_function', str)
    @typed_kwargs('compiler.has_function', *_COMMON_KWS)
    def has_function_method(self, args: T.Tuple[str], kwargs: 'CommonKW') -> bool:
        funcname = args[0]
        extra_args = self._determine_args(kwargs['no_builtin_args'], kwargs['include_directories'], kwargs['args'])
        deps, msg = self._determine_dependencies(kwargs['dependencies'], compile_only=False)
        had, cached = self.compiler.has_function(funcname, kwargs['prefix'], self.environment,
                                                 extra_args=extra_args,
                                                 dependencies=deps)
        cached_msg = mlog.blue('(cached)') if cached else ''
        if had:
            hadtxt = mlog.green('YES')
        else:
            hadtxt = mlog.red('NO')
        mlog.log('Checking for function', mlog.bold(funcname, True), msg, hadtxt, cached_msg)
        return had

    @typed_pos_args('compiler.has_type', str)
    @typed_kwargs('compiler.has_type', *_COMMON_KWS)
    def has_type_method(self, args: T.Tuple[str], kwargs: 'CommonKW') -> bool:
        typename = args[0]
        extra_args = functools.partial(self._determine_args, kwargs['no_builtin_args'], kwargs['include_directories'], kwargs['args'])
        deps, msg = self._determine_dependencies(kwargs['dependencies'])
        had, cached = self.compiler.has_type(typename, kwargs['prefix'], self.environment,
                                             extra_args=extra_args, dependencies=deps)
        cached_msg = mlog.blue('(cached)') if cached else ''
        if had:
            hadtxt = mlog.green('YES')
        else:
            hadtxt = mlog.red('NO')
        mlog.log('Checking for type', mlog.bold(typename, True), msg, hadtxt, cached_msg)
        return had

    @FeatureNew('compiler.compute_int', '0.40.0')
    @typed_pos_args('compiler.compute_int', str)
    @typed_kwargs(
        'compiler.compute_int',
        KwargInfo('low', (int, NoneType)),
        KwargInfo('high', (int, NoneType)),
        KwargInfo('guess', (int, NoneType)),
        *_COMMON_KWS,
    )
    def compute_int_method(self, args: T.Tuple[str], kwargs: 'ComputeIntKW') -> int:
        expression = args[0]
        extra_args = functools.partial(self._determine_args, kwargs['no_builtin_args'], kwargs['include_directories'], kwargs['args'])
        deps, msg = self._determine_dependencies(kwargs['dependencies'], compile_only=self.compiler.is_cross)
        res = self.compiler.compute_int(expression, kwargs['low'], kwargs['high'],
                                        kwargs['guess'], kwargs['prefix'],
                                        self.environment, extra_args=extra_args,
                                        dependencies=deps)
        mlog.log('Computing int of', mlog.bold(expression, True), msg, res)
        return res

    @typed_pos_args('compiler.sizeof', str)
    @typed_kwargs('compiler.sizeof', *_COMMON_KWS)
    def sizeof_method(self, args: T.Tuple[str], kwargs: 'CommonKW') -> int:
        element = args[0]
        extra_args = functools.partial(self._determine_args, kwargs['no_builtin_args'], kwargs['include_directories'], kwargs['args'])
        deps, msg = self._determine_dependencies(kwargs['dependencies'], compile_only=self.compiler.is_cross)
        esize, cached = self.compiler.sizeof(element, kwargs['prefix'], self.environment,
                                             extra_args=extra_args, dependencies=deps)
        cached_msg = mlog.blue('(cached)') if cached else ''
        mlog.log('Checking for size of',
                 mlog.bold(element, True), msg, mlog.bold(str(esize)), cached_msg)
        return esize

    @FeatureNew('compiler.get_define', '0.40.0')
    @typed_pos_args('compiler.get_define', str)
    @typed_kwargs('compiler.get_define', *_COMMON_KWS)
    def get_define_method(self, args: T.Tuple[str], kwargs: 'CommonKW') -> str:
        element = args[0]
        extra_args = functools.partial(self._determine_args, kwargs['no_builtin_args'], kwargs['include_directories'], kwargs['args'])
        deps, msg = self._determine_dependencies(kwargs['dependencies'])
        value, cached = self.compiler.get_define(element, kwargs['prefix'], self.environment,
                                                 extra_args=extra_args,
                                                 dependencies=deps)
        cached_msg = mlog.blue('(cached)') if cached else ''
        mlog.log('Fetching value of define', mlog.bold(element, True), msg, value, cached_msg)
        return value

    @typed_pos_args('compiler.compiles', (str, mesonlib.File))
    @typed_kwargs('compiler.compiles', *_COMPILES_KWS)
    def compiles_method(self, args: T.Tuple['mesonlib.FileOrString'], kwargs: 'CompileKW') -> bool:
        code = args[0]
        if isinstance(code, mesonlib.File):
            if code.is_built:
                FeatureNew.single_use('compiler.compiles with file created at setup time', '1.2.0', self.subproject,
                                      'It was broken and either errored or returned false.', self.current_node)
            self.interpreter.add_build_def_file(code)
            code = mesonlib.File.from_absolute_file(
                code.absolute_path(self.environment.source_dir, self.environment.build_dir))
        testname = kwargs['name']
        extra_args = functools.partial(self._determine_args, kwargs['no_builtin_args'], kwargs['include_directories'], kwargs['args'])
        deps, msg = self._determine_dependencies(kwargs['dependencies'], endl=None)
        result, cached = self.compiler.compiles(code, self.environment,
                                                extra_args=extra_args,
                                                dependencies=deps)
        if testname:
            if result:
                h = mlog.green('YES')
            else:
                h = mlog.red('NO')
            cached_msg = mlog.blue('(cached)') if cached else ''
            mlog.log('Checking if', mlog.bold(testname, True), msg, 'compiles:', h, cached_msg)
        return result

    @typed_pos_args('compiler.links', (str, mesonlib.File))
    @typed_kwargs('compiler.links', *_COMPILES_KWS)
    def links_method(self, args: T.Tuple['mesonlib.FileOrString'], kwargs: 'CompileKW') -> bool:
        code = args[0]
        compiler = None
        if isinstance(code, mesonlib.File):
            if code.is_built:
                FeatureNew.single_use('compiler.links with file created at setup time', '1.2.0', self.subproject,
                                      'It was broken and either errored or returned false.', self.current_node)
            self.interpreter.add_build_def_file(code)
            code = mesonlib.File.from_absolute_file(
                code.absolute_path(self.environment.source_dir, self.environment.build_dir))
            suffix = code.suffix
            if suffix not in self.compiler.file_suffixes:
                for_machine = self.compiler.for_machine
                clist = self.interpreter.coredata.compilers[for_machine]
                if suffix not in SUFFIX_TO_LANG:
                    # just pass it to the compiler driver
                    mlog.warning(f'Unknown suffix for test file {code}')
                elif SUFFIX_TO_LANG[suffix] not in clist:
                    mlog.warning(f'Passed {SUFFIX_TO_LANG[suffix]} source to links method, not specified for {for_machine.get_lower_case_name()} machine.')
                else:
                    compiler = clist[SUFFIX_TO_LANG[suffix]]

        testname = kwargs['name']
        extra_args = functools.partial(self._determine_args, kwargs['no_builtin_args'], kwargs['include_directories'], kwargs['args'])
        deps, msg = self._determine_dependencies(kwargs['dependencies'], compile_only=False)
        result, cached = self.compiler.links(code, self.environment,
                                             compiler=compiler,
                                             extra_args=extra_args,
                                             dependencies=deps)
        cached_msg = mlog.blue('(cached)') if cached else ''
        if testname:
            if result:
                h = mlog.green('YES')
            else:
                h = mlog.red('NO')
            mlog.log('Checking if', mlog.bold(testname, True), msg, 'links:', h, cached_msg)
        return result

    @FeatureNew('compiler.check_header', '0.47.0')
    @typed_pos_args('compiler.check_header', str)
    @typed_kwargs('compiler.check_header', *_HEADER_KWS)
    def check_header_method(self, args: T.Tuple[str], kwargs: 'HeaderKW') -> bool:
        hname = args[0]
        disabled, required, feature = extract_required_kwarg(kwargs, self.subproject, default=False)
        if disabled:
            mlog.log('Check usable header', mlog.bold(hname, True), 'skipped: feature', mlog.bold(feature), 'disabled')
            return False
        extra_args = functools.partial(self._determine_args, kwargs['no_builtin_args'], kwargs['include_directories'], kwargs['args'])
        deps, msg = self._determine_dependencies(kwargs['dependencies'])
        haz, cached = self.compiler.check_header(hname, kwargs['prefix'], self.environment,
                                                 extra_args=extra_args,
                                                 dependencies=deps)
        cached_msg = mlog.blue('(cached)') if cached else ''
        if required and not haz:
            raise InterpreterException(f'{self.compiler.get_display_language()} header {hname!r} not usable')
        elif haz:
            h = mlog.green('YES')
        else:
            h = mlog.red('NO')
        mlog.log('Check usable header', mlog.bold(hname, True), msg, h, cached_msg)
        return haz

    def _has_header_impl(self, hname: str, kwargs: 'HeaderKW') -> bool:
        disabled, required, feature = extract_required_kwarg(kwargs, self.subproject, default=False)
        if disabled:
            mlog.log('Has header', mlog.bold(hname, True), 'skipped: feature', mlog.bold(feature), 'disabled')
            return False
        extra_args = functools.partial(self._determine_args, kwargs['no_builtin_args'], kwargs['include_directories'], kwargs['args'])
        deps, msg = self._determine_dependencies(kwargs['dependencies'])
        haz, cached = self.compiler.has_header(hname, kwargs['prefix'], self.environment,
                                               extra_args=extra_args, dependencies=deps)
        cached_msg = mlog.blue('(cached)') if cached else ''
        if required and not haz:
            raise InterpreterException(f'{self.compiler.get_display_language()} header {hname!r} not found')
        elif haz:
            h = mlog.green('YES')
        else:
            h = mlog.red('NO')
        mlog.log('Has header', mlog.bold(hname, True), msg, h, cached_msg)
        return haz

    @typed_pos_args('compiler.has_header', str)
    @typed_kwargs('compiler.has_header', *_HEADER_KWS)
    def has_header_method(self, args: T.Tuple[str], kwargs: 'HeaderKW') -> bool:
        return self._has_header_impl(args[0], kwargs)

    @typed_pos_args('compiler.has_header_symbol', str, str)
    @typed_kwargs('compiler.has_header_symbol', *_HEADER_KWS)
    def has_header_symbol_method(self, args: T.Tuple[str, str], kwargs: 'HeaderKW') -> bool:
        hname, symbol = args
        disabled, required, feature = extract_required_kwarg(kwargs, self.subproject, default=False)
        if disabled:
            mlog.log('Header', mlog.bold(hname, True), 'has symbol', mlog.bold(symbol, True), 'skipped: feature', mlog.bold(feature), 'disabled')
            return False
        extra_args = functools.partial(self._determine_args, kwargs['no_builtin_args'], kwargs['include_directories'], kwargs['args'])
        deps, msg = self._determine_dependencies(kwargs['dependencies'])
        haz, cached = self.compiler.has_header_symbol(hname, symbol, kwargs['prefix'], self.environment,
                                                      extra_args=extra_args,
                                                      dependencies=deps)
        if required and not haz:
            raise InterpreterException(f'{self.compiler.get_display_language()} symbol {symbol} not found in header {hname}')
        elif haz:
            h = mlog.green('YES')
        else:
            h = mlog.red('NO')
        cached_msg = mlog.blue('(cached)') if cached else ''
        mlog.log('Header', mlog.bold(hname, True), 'has symbol', mlog.bold(symbol, True), msg, h, cached_msg)
        return haz

    def notfound_library(self, libname: str) -> 'dependencies.ExternalLibrary':
        lib = dependencies.ExternalLibrary(libname, None,
                                           self.environment,
                                           self.compiler.language,
                                           silent=True)
        return lib

    @disablerIfNotFound
    @typed_pos_args('compiler.find_library', str)
    @typed_kwargs(
        'compiler.find_library',
        KwargInfo('required', (bool, coredata.UserFeatureOption), default=True),
        KwargInfo('has_headers', ContainerTypeInfo(list, str), listify=True, default=[], since='0.50.0'),
        KwargInfo('static', (bool, NoneType), since='0.51.0'),
        KwargInfo('disabler', bool, default=False, since='0.49.0'),
        KwargInfo('dirs', ContainerTypeInfo(list, str), listify=True, default=[]),
        *(k.evolve(name=f'header_{k.name}') for k in _HEADER_KWS)
    )
    def find_library_method(self, args: T.Tuple[str], kwargs: 'FindLibraryKW') -> 'dependencies.ExternalLibrary':
        # TODO add dependencies support?
        libname = args[0]

        disabled, required, feature = extract_required_kwarg(kwargs, self.subproject)
        if disabled:
            mlog.log('Library', mlog.bold(libname), 'skipped: feature', mlog.bold(feature), 'disabled')
            return self.notfound_library(libname)

        # This could be done with a comprehension, but that confuses the type
        # checker, and having it check this seems valuable
        has_header_kwargs: 'HeaderKW' = {
            'required': required,
            'args': kwargs['header_args'],
            'dependencies': kwargs['header_dependencies'],
            'include_directories': kwargs['header_include_directories'],
            'prefix': kwargs['header_prefix'],
            'no_builtin_args': kwargs['header_no_builtin_args'],
        }
        for h in kwargs['has_headers']:
            if not self._has_header_impl(h, has_header_kwargs):
                return self.notfound_library(libname)

        search_dirs = extract_search_dirs(kwargs)

        prefer_static = self.environment.coredata.get_option(OptionKey('prefer_static'))
        if kwargs['static'] is True:
            libtype = mesonlib.LibType.STATIC
        elif kwargs['static'] is False:
            libtype = mesonlib.LibType.SHARED
        elif prefer_static:
            libtype = mesonlib.LibType.PREFER_STATIC
        else:
            libtype = mesonlib.LibType.PREFER_SHARED
        linkargs = self.compiler.find_library(libname, self.environment, search_dirs, libtype)
        if required and not linkargs:
            if libtype == mesonlib.LibType.PREFER_SHARED:
                libtype_s = 'shared or static'
            else:
                libtype_s = libtype.name.lower()
            raise InterpreterException('{} {} library {!r} not found'
                                       .format(self.compiler.get_display_language(),
                                               libtype_s, libname))
        lib = dependencies.ExternalLibrary(libname, linkargs, self.environment,
                                           self.compiler.language)
        return lib

    def _has_argument_impl(self, arguments: T.Union[str, T.List[str]],
                           mode: _TestMode = _TestMode.COMPILER) -> bool:
        """Shared implementation for methods checking compiler and linker arguments."""
        # This simplifies the callers
        if isinstance(arguments, str):
            arguments = [arguments]
        test = self.compiler.has_multi_link_arguments if mode is _TestMode.LINKER else self.compiler.has_multi_arguments
        result, cached = test(arguments, self.environment)
        cached_msg = mlog.blue('(cached)') if cached else ''
        mlog.log(
            'Compiler for',
            self.compiler.get_display_language(),
            'supports{}'.format(' link' if mode is _TestMode.LINKER else ''),
            'arguments {}:'.format(' '.join(arguments)),
            mlog.green('YES') if result else mlog.red('NO'),
            cached_msg)
        return result

    @noKwargs
    @typed_pos_args('compiler.has_argument', str)
    def has_argument_method(self, args: T.Tuple[str], kwargs: 'TYPE_kwargs') -> bool:
        return self._has_argument_impl([args[0]])

    @noKwargs
    @typed_pos_args('compiler.has_multi_arguments', varargs=str)
    @FeatureNew('compiler.has_multi_arguments', '0.37.0')
    def has_multi_arguments_method(self, args: T.Tuple[T.List[str]], kwargs: 'TYPE_kwargs') -> bool:
        return self._has_argument_impl(args[0])

    @FeatureNew('compiler.get_supported_arguments', '0.43.0')
    @typed_pos_args('compiler.get_supported_arguments', varargs=str)
    @typed_kwargs(
        'compiler.get_supported_arguments',
        KwargInfo('checked', str, default='off', since='0.59.0',
                  validator=in_set_validator({'warn', 'require', 'off'})),
    )
    def get_supported_arguments_method(self, args: T.Tuple[T.List[str]], kwargs: 'GetSupportedArgumentKw') -> T.List[str]:
        supported_args: T.List[str] = []
        checked = kwargs['checked']

        for arg in args[0]:
            if not self._has_argument_impl([arg]):
                msg = f'Compiler for {self.compiler.get_display_language()} does not support "{arg}"'
                if checked == 'warn':
                    mlog.warning(msg)
                elif checked == 'require':
                    raise mesonlib.MesonException(msg)
            else:
                supported_args.append(arg)
        return supported_args

    @noKwargs
    @typed_pos_args('compiler.first_supported_argument', varargs=str)
    def first_supported_argument_method(self, args: T.Tuple[T.List[str]], kwargs: 'TYPE_kwargs') -> T.List[str]:
        for arg in args[0]:
            if self._has_argument_impl([arg]):
                mlog.log('First supported argument:', mlog.bold(arg))
                return [arg]
        mlog.log('First supported argument:', mlog.red('None'))
        return []

    @FeatureNew('compiler.has_link_argument', '0.46.0')
    @noKwargs
    @typed_pos_args('compiler.has_link_argument', str)
    def has_link_argument_method(self, args: T.Tuple[str], kwargs: 'TYPE_kwargs') -> bool:
        return self._has_argument_impl([args[0]], mode=_TestMode.LINKER)

    @FeatureNew('compiler.has_multi_link_argument', '0.46.0')
    @noKwargs
    @typed_pos_args('compiler.has_multi_link_argument', varargs=str)
    def has_multi_link_arguments_method(self, args: T.Tuple[T.List[str]], kwargs: 'TYPE_kwargs') -> bool:
        return self._has_argument_impl(args[0], mode=_TestMode.LINKER)

    @FeatureNew('compiler.get_supported_link_arguments', '0.46.0')
    @noKwargs
    @typed_pos_args('compiler.get_supported_link_arguments', varargs=str)
    def get_supported_link_arguments_method(self, args: T.Tuple[T.List[str]], kwargs: 'TYPE_kwargs') -> T.List[str]:
        supported_args: T.List[str] = []
        for arg in args[0]:
            if self._has_argument_impl([arg], mode=_TestMode.LINKER):
                supported_args.append(arg)
        return supported_args

    @FeatureNew('compiler.first_supported_link_argument_method', '0.46.0')
    @noKwargs
    @typed_pos_args('compiler.first_supported_link_argument', varargs=str)
    def first_supported_link_argument_method(self, args: T.Tuple[T.List[str]], kwargs: 'TYPE_kwargs') -> T.List[str]:
        for arg in args[0]:
            if self._has_argument_impl([arg], mode=_TestMode.LINKER):
                mlog.log('First supported link argument:', mlog.bold(arg))
                return [arg]
        mlog.log('First supported link argument:', mlog.red('None'))
        return []

    def _has_function_attribute_impl(self, attr: str) -> bool:
        """Common helper for function attribute testing."""
        result, cached = self.compiler.has_func_attribute(attr, self.environment)
        cached_msg = mlog.blue('(cached)') if cached else ''
        h = mlog.green('YES') if result else mlog.red('NO')
        mlog.log(f'Compiler for {self.compiler.get_display_language()} supports function attribute {attr}:', h, cached_msg)
        return result

    @FeatureNew('compiler.has_function_attribute', '0.48.0')
    @noKwargs
    @typed_pos_args('compiler.has_function_attribute', str)
    def has_func_attribute_method(self, args: T.Tuple[str], kwargs: 'TYPE_kwargs') -> bool:
        return self._has_function_attribute_impl(args[0])

    @FeatureNew('compiler.get_supported_function_attributes', '0.48.0')
    @noKwargs
    @typed_pos_args('compiler.get_supported_function_attributes', varargs=str)
    def get_supported_function_attributes_method(self, args: T.Tuple[T.List[str]], kwargs: 'TYPE_kwargs') -> T.List[str]:
        return [a for a in args[0] if self._has_function_attribute_impl(a)]

    @FeatureNew('compiler.get_argument_syntax_method', '0.49.0')
    @noPosargs
    @noKwargs
    def get_argument_syntax_method(self, args: T.List['TYPE_var'], kwargs: 'TYPE_kwargs') -> str:
        return self.compiler.get_argument_syntax()

    @FeatureNew('compiler.preprocess', '0.64.0')
    @typed_pos_args('compiler.preprocess', varargs=(str, mesonlib.File, build.CustomTarget, build.CustomTargetIndex, build.GeneratedList), min_varargs=1)
    @typed_kwargs(
        'compiler.preprocess',
        KwargInfo('output', str, default='@PLAINNAME@.i'),
        KwargInfo('compile_args', ContainerTypeInfo(list, str), listify=True, default=[]),
        _INCLUDE_DIRS_KW,
        _DEPENDENCIES_KW.evolve(since='1.1.0'),
    )
    def preprocess_method(self, args: T.Tuple[T.List['mesonlib.FileOrString']], kwargs: 'PreprocessKW') -> T.List[build.CustomTargetIndex]:
        compiler = self.compiler.get_preprocessor()
        sources: 'SourceOutputs' = self.interpreter.source_strings_to_files(args[0])
        if any(isinstance(s, (build.CustomTarget, build.CustomTargetIndex, build.GeneratedList)) for s in sources):
            FeatureNew.single_use('compiler.preprocess with generated sources', '1.1.0', self.subproject,
                                  location=self.current_node)

        tg_counter = next(self.preprocess_uid[self.interpreter.subdir])
        if tg_counter > 0:
            FeatureNew.single_use('compiler.preprocess used multiple times', '1.1.0', self.subproject,
                                  location=self.current_node)
        tg_name = f'preprocessor_{tg_counter}'
        tg = build.CompileTarget(
            tg_name,
            self.interpreter.subdir,
            self.subproject,
            self.environment,
            sources,
            kwargs['output'],
            compiler,
            self.interpreter.backend,
            kwargs['compile_args'],
            kwargs['include_directories'],
            kwargs['dependencies'])
        self.interpreter.add_target(tg.name, tg)
        # Expose this target as list of its outputs, so user can pass them to
        # other targets, list outputs, etc.
        private_dir = os.path.relpath(self.interpreter.backend.get_target_private_dir(tg), self.interpreter.subdir)
        return [build.CustomTargetIndex(tg, os.path.join(private_dir, o)) for o in tg.outputs]
