# Copyright 2018 The Meson development team

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# This class contains the basic functionality needed to run any interpreter
# or an interpreter-based tool

from __future__ import annotations
import argparse
import copy
import os
import typing as T

from .. import compilers, environment, mesonlib, optinterpreter
from .. import coredata as cdata
from ..build import Executable, Jar, SharedLibrary, SharedModule, StaticLibrary
from ..compilers import detect_compiler_for
from ..interpreterbase import InvalidArguments
from ..mesonlib import MachineChoice, OptionKey
from ..mparser import BaseNode, ArithmeticNode, ArrayNode, ElementaryNode, IdNode, FunctionNode, StringNode
from .interpreter import AstInterpreter

if T.TYPE_CHECKING:
    from ..build import BuildTarget
    from ..interpreterbase import TYPE_nvar
    from .visitor import AstVisitor


# TODO: it would be nice to not have to duplicate this
BUILD_TARGET_FUNCTIONS = [
    'executable', 'jar', 'library', 'shared_library', 'shared_module',
    'static_library', 'both_libraries'
]

class IntrospectionHelper(argparse.Namespace):
    # mimic an argparse namespace
    def __init__(self, cross_file: str):
        super().__init__()
        self.cross_file = cross_file  # type: str
        self.native_file = None       # type: str
        self.cmd_line_options = {}    # type: T.Dict[str, str]

    def __eq__(self, other: object) -> bool:
        return NotImplemented

class IntrospectionInterpreter(AstInterpreter):
    # Interpreter to detect the options without a build directory
    # Most of the code is stolen from interpreter.Interpreter
    def __init__(self,
                 source_root: str,
                 subdir: str,
                 backend: str,
                 visitors: T.Optional[T.List[AstVisitor]] = None,
                 cross_file: T.Optional[str] = None,
                 subproject: str = '',
                 subproject_dir: str = 'subprojects',
                 env: T.Optional[environment.Environment] = None):
        visitors = visitors if visitors is not None else []
        super().__init__(source_root, subdir, subproject, visitors=visitors)

        options = IntrospectionHelper(cross_file)
        self.cross_file = cross_file
        if env is None:
            self.environment = environment.Environment(source_root, None, options)
        else:
            self.environment = env
        self.subproject_dir = subproject_dir
        self.coredata = self.environment.get_coredata()
        self.backend = backend
        self.default_options = {OptionKey('backend'): self.backend}
        self.project_data = {}    # type: T.Dict[str, T.Any]
        self.targets = []         # type: T.List[T.Dict[str, T.Any]]
        self.dependencies = []    # type: T.List[T.Dict[str, T.Any]]
        self.project_node = None  # type: BaseNode

        self.funcs.update({
            'add_languages': self.func_add_languages,
            'dependency': self.func_dependency,
            'executable': self.func_executable,
            'jar': self.func_jar,
            'library': self.func_library,
            'project': self.func_project,
            'shared_library': self.func_shared_lib,
            'shared_module': self.func_shared_module,
            'static_library': self.func_static_lib,
            'both_libraries': self.func_both_lib,
        })

    def func_project(self, node: BaseNode, args: T.List[TYPE_nvar], kwargs: T.Dict[str, TYPE_nvar]) -> None:
        if self.project_node:
            raise InvalidArguments('Second call to project()')
        self.project_node = node
        if len(args) < 1:
            raise InvalidArguments('Not enough arguments to project(). Needs at least the project name.')

        proj_name = args[0]
        proj_vers = kwargs.get('version', 'undefined')
        proj_langs = self.flatten_args(args[1:])
        if isinstance(proj_vers, ElementaryNode):
            proj_vers = proj_vers.value
        if not isinstance(proj_vers, str):
            proj_vers = 'undefined'
        self.project_data = {'descriptive_name': proj_name, 'version': proj_vers}

        optfile = os.path.join(self.source_root, self.subdir, 'meson.options')
        if not os.path.exists(optfile):
            optfile = os.path.join(self.source_root, self.subdir, 'meson_options.txt')
        if os.path.exists(optfile):
            oi = optinterpreter.OptionInterpreter(self.subproject)
            oi.process(optfile)
            self.coredata.update_project_options(oi.options)

        def_opts = self.flatten_args(kwargs.get('default_options', []))
        _project_default_options = mesonlib.stringlistify(def_opts)
        self.project_default_options = cdata.create_options_dict(_project_default_options, self.subproject)
        self.default_options.update(self.project_default_options)
        self.coredata.set_default_options(self.default_options, self.subproject, self.environment)

        if not self.is_subproject() and 'subproject_dir' in kwargs:
            spdirname = kwargs['subproject_dir']
            if isinstance(spdirname, StringNode):
                assert isinstance(spdirname.value, str)
                self.subproject_dir = spdirname.value
        if not self.is_subproject():
            self.project_data['subprojects'] = []
            subprojects_dir = os.path.join(self.source_root, self.subproject_dir)
            if os.path.isdir(subprojects_dir):
                for i in os.listdir(subprojects_dir):
                    if os.path.isdir(os.path.join(subprojects_dir, i)):
                        self.do_subproject(i)

        self.coredata.init_backend_options(self.backend)
        options = {k: v for k, v in self.environment.options.items() if k.is_backend()}

        self.coredata.set_options(options)
        self._add_languages(proj_langs, True, MachineChoice.HOST)
        self._add_languages(proj_langs, True, MachineChoice.BUILD)

    def do_subproject(self, dirname: str) -> None:
        subproject_dir_abs = os.path.join(self.environment.get_source_dir(), self.subproject_dir)
        subpr = os.path.join(subproject_dir_abs, dirname)
        try:
            subi = IntrospectionInterpreter(subpr, '', self.backend, cross_file=self.cross_file, subproject=dirname, subproject_dir=self.subproject_dir, env=self.environment, visitors=self.visitors)
            subi.analyze()
            subi.project_data['name'] = dirname
            self.project_data['subprojects'] += [subi.project_data]
        except (mesonlib.MesonException, RuntimeError):
            return

    def func_add_languages(self, node: BaseNode, args: T.List[TYPE_nvar], kwargs: T.Dict[str, TYPE_nvar]) -> None:
        kwargs = self.flatten_kwargs(kwargs)
        required = kwargs.get('required', True)
        if isinstance(required, cdata.UserFeatureOption):
            required = required.is_enabled()
        if 'native' in kwargs:
            native = kwargs.get('native', False)
            self._add_languages(args, required, MachineChoice.BUILD if native else MachineChoice.HOST)
        else:
            for for_machine in [MachineChoice.BUILD, MachineChoice.HOST]:
                self._add_languages(args, required, for_machine)

    def _add_languages(self, raw_langs: T.List[TYPE_nvar], required: bool, for_machine: MachineChoice) -> None:
        langs = []  # type: T.List[str]
        for l in self.flatten_args(raw_langs):
            if isinstance(l, str):
                langs.append(l)
            elif isinstance(l, StringNode):
                langs.append(l.value)

        for lang in sorted(langs, key=compilers.sort_clink):
            lang = lang.lower()
            if lang not in self.coredata.compilers[for_machine]:
                try:
                    comp = detect_compiler_for(self.environment, lang, for_machine, True)
                except mesonlib.MesonException:
                    # do we even care about introspecting this language?
                    if required:
                        raise
                    else:
                        continue
                if self.subproject:
                    options = {}
                    for k in comp.get_options():
                        v = copy.copy(self.coredata.options[k])
                        k = k.evolve(subproject=self.subproject)
                        options[k] = v
                    self.coredata.add_compiler_options(options, lang, for_machine, self.environment)

    def func_dependency(self, node: BaseNode, args: T.List[TYPE_nvar], kwargs: T.Dict[str, TYPE_nvar]) -> None:
        args = self.flatten_args(args)
        kwargs = self.flatten_kwargs(kwargs)
        if not args:
            return
        name = args[0]
        has_fallback = 'fallback' in kwargs
        required = kwargs.get('required', True)
        version = kwargs.get('version', [])
        if not isinstance(version, list):
            version = [version]
        if isinstance(required, ElementaryNode):
            required = required.value
        if not isinstance(required, bool):
            required = False
        self.dependencies += [{
            'name': name,
            'required': required,
            'version': version,
            'has_fallback': has_fallback,
            'conditional': node.condition_level > 0,
            'node': node
        }]

    def build_target(self, node: BaseNode, args: T.List[TYPE_nvar], kwargs_raw: T.Dict[str, TYPE_nvar], targetclass: T.Type[BuildTarget]) -> T.Optional[T.Dict[str, T.Any]]:
        args = self.flatten_args(args)
        if not args or not isinstance(args[0], str):
            return None
        name = args[0]
        srcqueue = [node]
        extra_queue = []

        # Process the sources BEFORE flattening the kwargs, to preserve the original nodes
        if 'sources' in kwargs_raw:
            srcqueue += mesonlib.listify(kwargs_raw['sources'])

        if 'extra_files' in kwargs_raw:
            extra_queue += mesonlib.listify(kwargs_raw['extra_files'])

        kwargs = self.flatten_kwargs(kwargs_raw, True)

        def traverse_nodes(inqueue: T.List[BaseNode]) -> T.List[BaseNode]:
            res = []  # type: T.List[BaseNode]
            while inqueue:
                curr = inqueue.pop(0)
                arg_node = None
                assert isinstance(curr, BaseNode)
                if isinstance(curr, FunctionNode):
                    arg_node = curr.args
                elif isinstance(curr, ArrayNode):
                    arg_node = curr.args
                elif isinstance(curr, IdNode):
                    # Try to resolve the ID and append the node to the queue
                    assert isinstance(curr.value, str)
                    var_name = curr.value
                    if var_name in self.assignments:
                        tmp_node = self.assignments[var_name]
                        if isinstance(tmp_node, (ArrayNode, IdNode, FunctionNode)):
                            inqueue += [tmp_node]
                elif isinstance(curr, ArithmeticNode):
                    inqueue += [curr.left, curr.right]
                if arg_node is None:
                    continue
                arg_nodes = arg_node.arguments.copy()
                # Pop the first element if the function is a build target function
                if isinstance(curr, FunctionNode) and curr.func_name in BUILD_TARGET_FUNCTIONS:
                    arg_nodes.pop(0)
                elementary_nodes = [x for x in arg_nodes if isinstance(x, (str, StringNode))]
                inqueue += [x for x in arg_nodes if isinstance(x, (FunctionNode, ArrayNode, IdNode, ArithmeticNode))]
                if elementary_nodes:
                    res += [curr]
            return res

        source_nodes = traverse_nodes(srcqueue)
        extraf_nodes = traverse_nodes(extra_queue)

        # Make sure nothing can crash when creating the build class
        kwargs_reduced = {k: v for k, v in kwargs.items() if k in targetclass.known_kwargs and k in {'install', 'build_by_default', 'build_always'}}
        kwargs_reduced = {k: v.value if isinstance(v, ElementaryNode) else v for k, v in kwargs_reduced.items()}
        kwargs_reduced = {k: v for k, v in kwargs_reduced.items() if not isinstance(v, BaseNode)}
        for_machine = MachineChoice.HOST
        objects = []        # type: T.List[T.Any]
        empty_sources = []  # type: T.List[T.Any]
        # Passing the unresolved sources list causes errors
        kwargs_reduced['_allow_no_sources'] = True
        target = targetclass(name, self.subdir, self.subproject, for_machine, empty_sources, [], objects,
                             self.environment, self.coredata.compilers[for_machine], kwargs_reduced)
        target.process_compilers_late()

        new_target = {
            'name': target.get_basename(),
            'id': target.get_id(),
            'type': target.get_typename(),
            'defined_in': os.path.normpath(os.path.join(self.source_root, self.subdir, environment.build_filename)),
            'subdir': self.subdir,
            'build_by_default': target.build_by_default,
            'installed': target.should_install(),
            'outputs': target.get_outputs(),
            'sources': source_nodes,
            'extra_files': extraf_nodes,
            'kwargs': kwargs,
            'node': node,
        }

        self.targets += [new_target]
        return new_target

    def build_library(self, node: BaseNode, args: T.List[TYPE_nvar], kwargs: T.Dict[str, TYPE_nvar]) -> T.Optional[T.Dict[str, T.Any]]:
        default_library = self.coredata.get_option(OptionKey('default_library'))
        if default_library == 'shared':
            return self.build_target(node, args, kwargs, SharedLibrary)
        elif default_library == 'static':
            return self.build_target(node, args, kwargs, StaticLibrary)
        elif default_library == 'both':
            return self.build_target(node, args, kwargs, SharedLibrary)
        return None

    def func_executable(self, node: BaseNode, args: T.List[TYPE_nvar], kwargs: T.Dict[str, TYPE_nvar]) -> T.Optional[T.Dict[str, T.Any]]:
        return self.build_target(node, args, kwargs, Executable)

    def func_static_lib(self, node: BaseNode, args: T.List[TYPE_nvar], kwargs: T.Dict[str, TYPE_nvar]) -> T.Optional[T.Dict[str, T.Any]]:
        return self.build_target(node, args, kwargs, StaticLibrary)

    def func_shared_lib(self, node: BaseNode, args: T.List[TYPE_nvar], kwargs: T.Dict[str, TYPE_nvar]) -> T.Optional[T.Dict[str, T.Any]]:
        return self.build_target(node, args, kwargs, SharedLibrary)

    def func_both_lib(self, node: BaseNode, args: T.List[TYPE_nvar], kwargs: T.Dict[str, TYPE_nvar]) -> T.Optional[T.Dict[str, T.Any]]:
        return self.build_target(node, args, kwargs, SharedLibrary)

    def func_shared_module(self, node: BaseNode, args: T.List[TYPE_nvar], kwargs: T.Dict[str, TYPE_nvar]) -> T.Optional[T.Dict[str, T.Any]]:
        return self.build_target(node, args, kwargs, SharedModule)

    def func_library(self, node: BaseNode, args: T.List[TYPE_nvar], kwargs: T.Dict[str, TYPE_nvar]) -> T.Optional[T.Dict[str, T.Any]]:
        return self.build_library(node, args, kwargs)

    def func_jar(self, node: BaseNode, args: T.List[TYPE_nvar], kwargs: T.Dict[str, TYPE_nvar]) -> T.Optional[T.Dict[str, T.Any]]:
        return self.build_target(node, args, kwargs, Jar)

    def func_build_target(self, node: BaseNode, args: T.List[TYPE_nvar], kwargs: T.Dict[str, TYPE_nvar]) -> T.Optional[T.Dict[str, T.Any]]:
        if 'target_type' not in kwargs:
            return None
        target_type = kwargs.pop('target_type')
        if isinstance(target_type, ElementaryNode):
            target_type = target_type.value
        if target_type == 'executable':
            return self.build_target(node, args, kwargs, Executable)
        elif target_type == 'shared_library':
            return self.build_target(node, args, kwargs, SharedLibrary)
        elif target_type == 'static_library':
            return self.build_target(node, args, kwargs, StaticLibrary)
        elif target_type == 'both_libraries':
            return self.build_target(node, args, kwargs, SharedLibrary)
        elif target_type == 'library':
            return self.build_library(node, args, kwargs)
        elif target_type == 'jar':
            return self.build_target(node, args, kwargs, Jar)
        return None

    def is_subproject(self) -> bool:
        return self.subproject != ''

    def analyze(self) -> None:
        self.load_root_meson_file()
        self.sanity_check_ast()
        self.parse_project()
        self.run()
