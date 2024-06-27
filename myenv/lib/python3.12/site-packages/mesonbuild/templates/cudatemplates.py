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


hello_cuda_template = '''#include <iostream>

#define PROJECT_NAME "{project_name}"

int main(int argc, char **argv) {{
    if(argc != 1) {{
        std::cout << argv[0] << " takes no arguments.\\n";
        return 1;
    }}
    std::cout << "This is project " << PROJECT_NAME << ".\\n";
    return 0;
}}
'''

hello_cuda_meson_template = '''project('{project_name}', ['cuda', 'cpp'],
  version : '{version}',
  default_options : ['warning_level=3',
                     'cpp_std=c++14'])

exe = executable('{exe_name}', '{source_name}',
  install : true)

test('basic', exe)
'''

lib_h_template = '''#pragma once
#if defined _WIN32 || defined __CYGWIN__
  #ifdef BUILDING_{utoken}
    #define {utoken}_PUBLIC __declspec(dllexport)
  #else
    #define {utoken}_PUBLIC __declspec(dllimport)
  #endif
#else
  #ifdef BUILDING_{utoken}
      #define {utoken}_PUBLIC __attribute__ ((visibility ("default")))
  #else
      #define {utoken}_PUBLIC
  #endif
#endif

namespace {namespace} {{

class {utoken}_PUBLIC {class_name} {{

public:
  {class_name}();
  int get_number() const;

private:

  int number;

}};

}}

'''

lib_cuda_template = '''#include <{header_file}>

namespace {namespace} {{

{class_name}::{class_name}() {{
    number = 6;
}}

int {class_name}::get_number() const {{
  return number;
}}

}}
'''

lib_cuda_test_template = '''#include <{header_file}>
#include <iostream>

int main(int argc, char **argv) {{
    if(argc != 1) {{
        std::cout << argv[0] << " takes no arguments.\\n";
        return 1;
    }}
    {namespace}::{class_name} c;
    return c.get_number() != 6;
}}
'''

lib_cuda_meson_template = '''project('{project_name}', ['cuda', 'cpp'],
  version : '{version}',
  default_options : ['warning_level=3'])

# These arguments are only used to build the shared library
# not the executables that use the library.
lib_args = ['-DBUILDING_{utoken}']

shlib = shared_library('{lib_name}', '{source_file}',
  install : true,
  cpp_args : lib_args,
  gnu_symbol_visibility : 'hidden',
)

test_exe = executable('{test_exe_name}', '{test_source_file}',
  link_with : shlib)
test('{test_name}', test_exe)

# Make this library usable as a Meson subproject.
{ltoken}_dep = declare_dependency(
  include_directories: include_directories('.'),
  link_with : shlib)

# Make this library usable from the system's
# package manager.
install_headers('{header_file}', subdir : '{header_dir}')

pkg_mod = import('pkgconfig')
pkg_mod.generate(
  name : '{project_name}',
  filebase : '{ltoken}',
  description : 'Meson sample project.',
  subdirs : '{header_dir}',
  libraries : shlib,
  version : '{version}',
)
'''


class CudaProject(SampleImpl):
    def __init__(self, options):
        super().__init__()
        self.name = options.name
        self.version = options.version

    def create_executable(self) -> None:
        lowercase_token = re.sub(r'[^a-z0-9]', '_', self.name.lower())
        source_name = lowercase_token + '.cu'
        open(source_name, 'w', encoding='utf-8').write(hello_cuda_template.format(project_name=self.name))
        open('meson.build', 'w', encoding='utf-8').write(
            hello_cuda_meson_template.format(project_name=self.name,
                                             exe_name=lowercase_token,
                                             source_name=source_name,
                                             version=self.version))

    def create_library(self) -> None:
        lowercase_token = re.sub(r'[^a-z0-9]', '_', self.name.lower())
        uppercase_token = lowercase_token.upper()
        class_name = uppercase_token[0] + lowercase_token[1:]
        test_exe_name = lowercase_token + '_test'
        namespace = lowercase_token
        lib_h_name = lowercase_token + '.h'
        lib_cuda_name = lowercase_token + '.cu'
        test_cuda_name = lowercase_token + '_test.cu'
        kwargs = {'utoken': uppercase_token,
                  'ltoken': lowercase_token,
                  'header_dir': lowercase_token,
                  'class_name': class_name,
                  'namespace': namespace,
                  'header_file': lib_h_name,
                  'source_file': lib_cuda_name,
                  'test_source_file': test_cuda_name,
                  'test_exe_name': test_exe_name,
                  'project_name': self.name,
                  'lib_name': lowercase_token,
                  'test_name': lowercase_token,
                  'version': self.version,
                  }
        open(lib_h_name, 'w', encoding='utf-8').write(lib_h_template.format(**kwargs))
        open(lib_cuda_name, 'w', encoding='utf-8').write(lib_cuda_template.format(**kwargs))
        open(test_cuda_name, 'w', encoding='utf-8').write(lib_cuda_test_template.format(**kwargs))
        open('meson.build', 'w', encoding='utf-8').write(lib_cuda_meson_template.format(**kwargs))
