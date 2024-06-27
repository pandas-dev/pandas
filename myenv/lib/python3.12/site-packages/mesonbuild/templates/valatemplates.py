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


hello_vala_template = '''void main (string[] args) {{
    stdout.printf ("Hello {project_name}!\\n");
}}
'''

hello_vala_meson_template = '''project('{project_name}', ['c', 'vala'],
  version : '{version}')

dependencies = [
    dependency('glib-2.0'),
    dependency('gobject-2.0'),
]

exe = executable('{exe_name}', '{source_name}', dependencies : dependencies,
  install : true)

test('basic', exe)
'''


lib_vala_template = '''namespace {namespace} {{
    public int sum(int a, int b) {{
        return(a + b);
    }}

    public int square(int a) {{
        return(a * a);
    }}
}}
'''

lib_vala_test_template = '''using {namespace};

public void main() {{
    stdout.printf("\nTesting shlib");
    stdout.printf("\n\t2 + 3 is %d", sum(2, 3));
    stdout.printf("\n\t8 squared is %d\\n", square(8));
}}
'''

lib_vala_meson_template = '''project('{project_name}', ['c', 'vala'],
  version : '{version}')

dependencies = [
    dependency('glib-2.0'),
    dependency('gobject-2.0'),
]

# These arguments are only used to build the shared library
# not the executables that use the library.
shlib = shared_library('foo', '{source_file}',
               dependencies: dependencies,
               install: true,
               install_dir: [true, true, true])

test_exe = executable('{test_exe_name}', '{test_source_file}', dependencies : dependencies,
  link_with : shlib)
test('{test_name}', test_exe)

# Make this library usable as a Meson subproject.
{ltoken}_dep = declare_dependency(
  include_directories: include_directories('.'),
  link_with : shlib)
'''


class ValaProject(SampleImpl):
    def __init__(self, options):
        super().__init__()
        self.name = options.name
        self.version = options.version

    def create_executable(self) -> None:
        lowercase_token = re.sub(r'[^a-z0-9]', '_', self.name.lower())
        source_name = lowercase_token + '.vala'
        open(source_name, 'w', encoding='utf-8').write(hello_vala_template.format(project_name=self.name))
        open('meson.build', 'w', encoding='utf-8').write(
            hello_vala_meson_template.format(project_name=self.name,
                                             exe_name=lowercase_token,
                                             source_name=source_name,
                                             version=self.version))

    def create_library(self) -> None:
        lowercase_token = re.sub(r'[^a-z0-9]', '_', self.name.lower())
        uppercase_token = lowercase_token.upper()
        class_name = uppercase_token[0] + lowercase_token[1:]
        test_exe_name = lowercase_token + '_test'
        namespace = lowercase_token
        lib_vala_name = lowercase_token + '.vala'
        test_vala_name = lowercase_token + '_test.vala'
        kwargs = {'utoken': uppercase_token,
                  'ltoken': lowercase_token,
                  'header_dir': lowercase_token,
                  'class_name': class_name,
                  'namespace': namespace,
                  'source_file': lib_vala_name,
                  'test_source_file': test_vala_name,
                  'test_exe_name': test_exe_name,
                  'project_name': self.name,
                  'lib_name': lowercase_token,
                  'test_name': lowercase_token,
                  'version': self.version,
                  }
        open(lib_vala_name, 'w', encoding='utf-8').write(lib_vala_template.format(**kwargs))
        open(test_vala_name, 'w', encoding='utf-8').write(lib_vala_test_template.format(**kwargs))
        open('meson.build', 'w', encoding='utf-8').write(lib_vala_meson_template.format(**kwargs))
