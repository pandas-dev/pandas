# Copyright 2013-2017 The Meson development team

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# This file contains the detection logic for external dependencies that
# are UI-related.
from __future__ import annotations

import os
import subprocess
import typing as T

from .. import mlog
from .. import mesonlib
from ..mesonlib import (
    Popen_safe, extract_as_list, version_compare_many
)
from ..environment import detect_cpu_family

from .base import DependencyException, DependencyMethods, DependencyTypeName, SystemDependency
from .configtool import ConfigToolDependency
from .detect import packages
from .factory import DependencyFactory

if T.TYPE_CHECKING:
    from ..environment import Environment


class GLDependencySystem(SystemDependency):
    def __init__(self, name: str, environment: 'Environment', kwargs: T.Dict[str, T.Any]) -> None:
        super().__init__(name, environment, kwargs)

        if self.env.machines[self.for_machine].is_darwin():
            self.is_found = True
            # FIXME: Use AppleFrameworks dependency
            self.link_args = ['-framework', 'OpenGL']
            # FIXME: Detect version using self.clib_compiler
            return
        elif self.env.machines[self.for_machine].is_windows():
            self.is_found = True
            # FIXME: Use self.clib_compiler.find_library()
            self.link_args = ['-lopengl32']
            # FIXME: Detect version using self.clib_compiler
            return
        else:
            links = self.clib_compiler.find_library('GL', environment, [])
            has_header = self.clib_compiler.has_header('GL/gl.h', '', environment)[0]
            if links and has_header:
                self.is_found = True
                self.link_args = links
            elif links:
                raise DependencyException('Found GL runtime library but no development header files')

class GnuStepDependency(ConfigToolDependency):

    tools = ['gnustep-config']
    tool_name = 'gnustep-config'

    def __init__(self, environment: 'Environment', kwargs: T.Dict[str, T.Any]) -> None:
        super().__init__('gnustep', environment, kwargs, language='objc')
        if not self.is_found:
            return
        self.modules = kwargs.get('modules', [])
        self.compile_args = self.filter_args(
            self.get_config_value(['--objc-flags'], 'compile_args'))
        self.link_args = self.weird_filter(self.get_config_value(
            ['--gui-libs' if 'gui' in self.modules else '--base-libs'],
            'link_args'))

    def find_config(self, versions: T.Optional[T.List[str]] = None, returncode: int = 0) -> T.Tuple[T.Optional[T.List[str]], T.Optional[str]]:
        tool = [self.tools[0]]
        try:
            p, out = Popen_safe(tool + ['--help'])[:2]
        except (FileNotFoundError, PermissionError):
            return (None, None)
        if p.returncode != returncode:
            return (None, None)
        self.config = tool
        found_version = self.detect_version()
        if versions and not version_compare_many(found_version, versions)[0]:
            return (None, found_version)

        return (tool, found_version)

    @staticmethod
    def weird_filter(elems: T.List[str]) -> T.List[str]:
        """When building packages, the output of the enclosing Make is
        sometimes mixed among the subprocess output. I have no idea why. As a
        hack filter out everything that is not a flag.
        """
        return [e for e in elems if e.startswith('-')]

    @staticmethod
    def filter_args(args: T.List[str]) -> T.List[str]:
        """gnustep-config returns a bunch of garbage args such as -O2 and so
        on. Drop everything that is not needed.
        """
        result = []
        for f in args:
            if f.startswith('-D') \
                    or f.startswith('-f') \
                    or f.startswith('-I') \
                    or f == '-pthread' \
                    or (f.startswith('-W') and not f == '-Wall'):
                result.append(f)
        return result

    def detect_version(self) -> str:
        gmake = self.get_config_value(['--variable=GNUMAKE'], 'variable')[0]
        makefile_dir = self.get_config_value(['--variable=GNUSTEP_MAKEFILES'], 'variable')[0]
        # This Makefile has the GNUStep version set
        base_make = os.path.join(makefile_dir, 'Additional', 'base.make')
        # Print the Makefile variable passed as the argument. For instance, if
        # you run the make target `print-SOME_VARIABLE`, this will print the
        # value of the variable `SOME_VARIABLE`.
        printver = "print-%:\n\t@echo '$($*)'"
        env = os.environ.copy()
        # See base.make to understand why this is set
        env['FOUNDATION_LIB'] = 'gnu'
        p, o, e = Popen_safe([gmake, '-f', '-', '-f', base_make,
                              'print-GNUSTEP_BASE_VERSION'],
                             env=env, write=printver, stdin=subprocess.PIPE)
        version = o.strip()
        if not version:
            mlog.debug("Couldn't detect GNUStep version, falling back to '1'")
            # Fallback to setting some 1.x version
            version = '1'
        return version

packages['gnustep'] = GnuStepDependency


class SDL2DependencyConfigTool(ConfigToolDependency):

    tools = ['sdl2-config']
    tool_name = 'sdl2-config'

    def __init__(self, name: str, environment: 'Environment', kwargs: T.Dict[str, T.Any]):
        super().__init__(name, environment, kwargs)
        if not self.is_found:
            return
        self.compile_args = self.get_config_value(['--cflags'], 'compile_args')
        self.link_args = self.get_config_value(['--libs'], 'link_args')


class WxDependency(ConfigToolDependency):

    tools = ['wx-config-3.0', 'wx-config-3.1', 'wx-config', 'wx-config-gtk3']
    tool_name = 'wx-config'

    def __init__(self, environment: 'Environment', kwargs: T.Dict[str, T.Any]):
        super().__init__('WxWidgets', environment, kwargs, language='cpp')
        if not self.is_found:
            return
        self.requested_modules = self.get_requested(kwargs)

        extra_args = []
        if self.static:
            extra_args.append('--static=yes')

            # Check to make sure static is going to work
            err = Popen_safe(self.config + extra_args)[2]
            if 'No config found to match' in err:
                mlog.debug('WxWidgets is missing static libraries.')
                self.is_found = False
                return

        # wx-config seems to have a cflags as well but since it requires C++,
        # this should be good, at least for now.
        self.compile_args = self.get_config_value(['--cxxflags'] + extra_args + self.requested_modules, 'compile_args')
        self.link_args = self.get_config_value(['--libs'] + extra_args + self.requested_modules, 'link_args')

    @staticmethod
    def get_requested(kwargs: T.Dict[str, T.Any]) -> T.List[str]:
        if 'modules' not in kwargs:
            return []
        candidates = extract_as_list(kwargs, 'modules')
        for c in candidates:
            if not isinstance(c, str):
                raise DependencyException('wxwidgets module argument is not a string')
        return candidates

packages['wxwidgets'] = WxDependency

class VulkanDependencySystem(SystemDependency):

    def __init__(self, name: str, environment: 'Environment', kwargs: T.Dict[str, T.Any], language: T.Optional[str] = None) -> None:
        super().__init__(name, environment, kwargs, language=language)

        try:
            self.vulkan_sdk = os.environ['VULKAN_SDK']
            if not os.path.isabs(self.vulkan_sdk):
                raise DependencyException('VULKAN_SDK must be an absolute path.')
        except KeyError:
            self.vulkan_sdk = None

        if self.vulkan_sdk:
            # TODO: this config might not work on some platforms, fix bugs as reported
            # we should at least detect other 64-bit platforms (e.g. armv8)
            lib_name = 'vulkan'
            lib_dir = 'lib'
            inc_dir = 'include'
            if mesonlib.is_windows():
                lib_name = 'vulkan-1'
                lib_dir = 'Lib32'
                inc_dir = 'Include'
                if detect_cpu_family(self.env.coredata.compilers.host) == 'x86_64':
                    lib_dir = 'Lib'

            # make sure header and lib are valid
            inc_path = os.path.join(self.vulkan_sdk, inc_dir)
            header = os.path.join(inc_path, 'vulkan', 'vulkan.h')
            lib_path = os.path.join(self.vulkan_sdk, lib_dir)
            find_lib = self.clib_compiler.find_library(lib_name, environment, [lib_path])

            if not find_lib:
                raise DependencyException('VULKAN_SDK point to invalid directory (no lib)')

            if not os.path.isfile(header):
                raise DependencyException('VULKAN_SDK point to invalid directory (no include)')

            # XXX: this is very odd, and may deserve being removed
            self.type_name = DependencyTypeName('vulkan_sdk')
            self.is_found = True
            self.compile_args.append('-I' + inc_path)
            self.link_args.append('-L' + lib_path)
            self.link_args.append('-l' + lib_name)

            # TODO: find a way to retrieve the version from the sdk?
            # Usually it is a part of the path to it (but does not have to be)
            return
        else:
            # simply try to guess it, usually works on linux
            libs = self.clib_compiler.find_library('vulkan', environment, [])
            if libs is not None and self.clib_compiler.has_header('vulkan/vulkan.h', '', environment, disable_cache=True)[0]:
                self.is_found = True
                for lib in libs:
                    self.link_args.append(lib)
                return

packages['gl'] = gl_factory = DependencyFactory(
    'gl',
    [DependencyMethods.PKGCONFIG, DependencyMethods.SYSTEM],
    system_class=GLDependencySystem,
)

packages['sdl2'] = sdl2_factory = DependencyFactory(
    'sdl2',
    [DependencyMethods.PKGCONFIG, DependencyMethods.CONFIG_TOOL, DependencyMethods.EXTRAFRAMEWORK, DependencyMethods.CMAKE],
    configtool_class=SDL2DependencyConfigTool,
    cmake_name='SDL2',
)

packages['vulkan'] = vulkan_factory = DependencyFactory(
    'vulkan',
    [DependencyMethods.PKGCONFIG, DependencyMethods.SYSTEM],
    system_class=VulkanDependencySystem,
)
