# Copyright 2019 The Meson development team

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
from __future__ import annotations

from mesonbuild.templates.sampleimpl import SampleImpl
import re


hello_d_template = '''module main;
import std.stdio;

enum PROJECT_NAME = "{project_name}";

int main(string[] args) {{
    if (args.length != 1){{
        writefln("%s takes no arguments.\\n", args[0]);
        return 1;
    }}
    writefln("This is project %s.\\n", PROJECT_NAME);
    return 0;
}}
'''

hello_d_meson_template = '''project('{project_name}', 'd',
    version : '{version}',
    default_options: ['warning_level=3'])

exe = executable('{exe_name}', '{source_name}',
  install : true)

test('basic', exe)
'''

lib_d_template = '''module {module_file};

/* This function will not be exported and is not
 * directly callable by users of this library.
 */
int internal_function() {{
    return 0;
}}

int {function_name}() {{
    return internal_function();
}}
'''

lib_d_test_template = '''module {module_file}_test;
import std.stdio;
import {module_file};


int main(string[] args) {{
    if (args.length != 1){{
        writefln("%s takes no arguments.\\n", args[0]);
        return 1;
    }}
    return {function_name}();
}}
'''

lib_d_meson_template = '''project('{project_name}', 'd',
  version : '{version}',
  default_options : ['warning_level=3'])

stlib = static_library('{lib_name}', '{source_file}',
  install : true,
  gnu_symbol_visibility : 'hidden',
)

test_exe = executable('{test_exe_name}', '{test_source_file}',
  link_with : stlib)
test('{test_name}', test_exe)

# Make this library usable as a Meson subproject.
{ltoken}_dep = declare_dependency(
  include_directories: include_directories('.'),
  link_with : stlib)

# Make this library usable from the Dlang
# build system.
dlang_mod = import('dlang')
if find_program('dub', required: false).found()
  dlang_mod.generate_dub_file(meson.project_name().to_lower(), meson.source_root(),
    name : meson.project_name(),
    license: meson.project_license(),
    sourceFiles : '{source_file}',
    description : 'Meson sample project.',
    version : '{version}',
  )
endif
'''


class DlangProject(SampleImpl):
    def __init__(self, options):
        super().__init__()
        self.name = options.name
        self.version = options.version

    def create_executable(self) -> None:
        lowercase_token = re.sub(r'[^a-z0-9]', '_', self.name.lower())
        source_name = lowercase_token + '.d'
        open(source_name, 'w', encoding='utf-8').write(hello_d_template.format(project_name=self.name))
        open('meson.build', 'w', encoding='utf-8').write(
            hello_d_meson_template.format(project_name=self.name,
                                          exe_name=lowercase_token,
                                          source_name=source_name,
                                          version=self.version))

    def create_library(self) -> None:
        lowercase_token = re.sub(r'[^a-z0-9]', '_', self.name.lower())
        uppercase_token = lowercase_token.upper()
        function_name = lowercase_token[0:3] + '_func'
        test_exe_name = lowercase_token + '_test'
        lib_m_name = lowercase_token
        lib_d_name = lowercase_token + '.d'
        test_d_name = lowercase_token + '_test.d'
        kwargs = {'utoken': uppercase_token,
                  'ltoken': lowercase_token,
                  'header_dir': lowercase_token,
                  'function_name': function_name,
                  'module_file': lib_m_name,
                  'source_file': lib_d_name,
                  'test_source_file': test_d_name,
                  'test_exe_name': test_exe_name,
                  'project_name': self.name,
                  'lib_name': lowercase_token,
                  'test_name': lowercase_token,
                  'version': self.version,
                  }
        open(lib_d_name, 'w', encoding='utf-8').write(lib_d_template.format(**kwargs))
        open(test_d_name, 'w', encoding='utf-8').write(lib_d_test_template.format(**kwargs))
        open('meson.build', 'w', encoding='utf-8').write(lib_d_meson_template.format(**kwargs))
