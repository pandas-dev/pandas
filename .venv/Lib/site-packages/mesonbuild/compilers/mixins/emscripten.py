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

"""Provides a mixin for shared code between C and C++ Emscripten compilers."""

import os.path
import typing as T

from ... import coredata
from ... import mesonlib
from ...mesonlib import OptionKey
from ...mesonlib import LibType

if T.TYPE_CHECKING:
    from ...environment import Environment
    from ...compilers.compilers import Compiler
    from ...dependencies import Dependency
else:
    # This is a bit clever, for mypy we pretend that these mixins descend from
    # Compiler, so we get all of the methods and attributes defined for us, but
    # for runtime we make them descend from object (which all classes normally
    # do). This gives up DRYer type checking, with no runtime impact
    Compiler = object


def wrap_js_includes(args: T.List[str]) -> T.List[str]:
    final_args = []
    for i in args:
        if i.endswith('.js') and not i.startswith('-'):
            final_args += ['--js-library', i]
        else:
            final_args += [i]
    return final_args

class EmscriptenMixin(Compiler):

    def _get_compile_output(self, dirname: str, mode: str) -> str:
        # In pre-processor mode, the output is sent to stdout and discarded
        if mode == 'preprocess':
            return None
        # Unlike sane toolchains, emcc infers the kind of output from its name.
        # This is the only reason why this method is overridden; compiler tests
        # do not work well with the default exe/obj suffices.
        if mode == 'link':
            suffix = 'js'
        else:
            suffix = 'o'
        return os.path.join(dirname, 'output.' + suffix)

    def thread_link_flags(self, env: 'Environment') -> T.List[str]:
        args = ['-pthread']
        count: int = env.coredata.options[OptionKey('thread_count', lang=self.language, machine=self.for_machine)].value
        if count:
            args.append(f'-sPTHREAD_POOL_SIZE={count}')
        return args

    def get_options(self) -> 'coredata.MutableKeyedOptionDictType':
        opts = super().get_options()
        key = OptionKey('thread_count', machine=self.for_machine, lang=self.language)
        opts.update({
            key: coredata.UserIntegerOption(
                'Number of threads to use in web assembly, set to 0 to disable',
                (0, None, 4),  # Default was picked at random
            ),
        })

        return opts

    @classmethod
    def native_args_to_unix(cls, args: T.List[str]) -> T.List[str]:
        return wrap_js_includes(super().native_args_to_unix(args))

    def get_dependency_link_args(self, dep: 'Dependency') -> T.List[str]:
        return wrap_js_includes(super().get_dependency_link_args(dep))

    def find_library(self, libname: str, env: 'Environment', extra_dirs: T.List[str],
                     libtype: LibType = LibType.PREFER_SHARED, lib_prefix_warning: bool = True) -> T.Optional[T.List[str]]:
        if not libname.endswith('.js'):
            return super().find_library(libname, env, extra_dirs, libtype, lib_prefix_warning)
        if os.path.isabs(libname):
            if os.path.exists(libname):
                return [libname]
        if len(extra_dirs) == 0:
            raise mesonlib.EnvironmentException('Looking up Emscripten JS libraries requires either an absolute path or specifying extra_dirs.')
        for d in extra_dirs:
            abs_path = os.path.join(d, libname)
            if os.path.exists(abs_path):
                return [abs_path]
        return None
