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


hello_cs_template = '''using System;

public class {class_name} {{
    const String PROJECT_NAME = "{project_name}";

    static int Main(String[] args) {{
      if (args.Length > 0) {{
          System.Console.WriteLine(String.Format("{project_name} takes no arguments.."));
          return 1;
      }}
      Console.WriteLine(String.Format("This is project {{0}}.", PROJECT_NAME));
      return 0;
    }}
}}

'''

hello_cs_meson_template = '''project('{project_name}', 'cs',
  version : '{version}',
  default_options : ['warning_level=3'])

exe = executable('{exe_name}', '{source_name}',
  install : true)

test('basic', exe)
'''

lib_cs_template = '''
public class {class_name} {{
    private const int number = 6;

    public int get_number() {{
      return number;
    }}
}}

'''

lib_cs_test_template = '''using System;

public class {class_test} {{
    static int Main(String[] args) {{
      if (args.Length > 0) {{
          System.Console.WriteLine("{project_name} takes no arguments..");
          return 1;
      }}
      {class_name} c = new {class_name}();
      Boolean result = true;
      return result.CompareTo(c.get_number() != 6);
    }}
}}

'''

lib_cs_meson_template = '''project('{project_name}', 'cs',
  version : '{version}',
  default_options : ['warning_level=3'])

stlib = shared_library('{lib_name}', '{source_file}',
  install : true,
)

test_exe = executable('{test_exe_name}', '{test_source_file}',
  link_with : stlib)
test('{test_name}', test_exe)

# Make this library usable as a Meson subproject.
{ltoken}_dep = declare_dependency(
  include_directories: include_directories('.'),
  link_with : stlib)

'''


class CSharpProject(SampleImpl):
    def __init__(self, options):
        super().__init__()
        self.name = options.name
        self.version = options.version

    def create_executable(self) -> None:
        lowercase_token = re.sub(r'[^a-z0-9]', '_', self.name.lower())
        uppercase_token = lowercase_token.upper()
        class_name = uppercase_token[0] + lowercase_token[1:]
        source_name = uppercase_token[0] + lowercase_token[1:] + '.cs'
        open(source_name, 'w', encoding='utf-8').write(
            hello_cs_template.format(project_name=self.name,
                                     class_name=class_name))
        open('meson.build', 'w', encoding='utf-8').write(
          hello_cs_meson_template.format(project_name=self.name,
                                         exe_name=self.name,
                                         source_name=source_name,
                                         version=self.version))

    def create_library(self) -> None:
        lowercase_token = re.sub(r'[^a-z0-9]', '_', self.name.lower())
        uppercase_token = lowercase_token.upper()
        class_name = uppercase_token[0] + lowercase_token[1:]
        class_test = uppercase_token[0] + lowercase_token[1:] + '_test'
        project_test = lowercase_token + '_test'
        lib_cs_name = uppercase_token[0] + lowercase_token[1:] + '.cs'
        test_cs_name = uppercase_token[0] + lowercase_token[1:] + '_test.cs'
        kwargs = {'utoken': uppercase_token,
                  'ltoken': lowercase_token,
                  'class_test': class_test,
                  'class_name': class_name,
                  'source_file': lib_cs_name,
                  'test_source_file': test_cs_name,
                  'test_exe_name': project_test,
                  'project_name': self.name,
                  'lib_name': lowercase_token,
                  'test_name': lowercase_token,
                  'version': self.version,
                  }
        open(lib_cs_name, 'w', encoding='utf-8').write(lib_cs_template.format(**kwargs))
        open(test_cs_name, 'w', encoding='utf-8').write(lib_cs_test_template.format(**kwargs))
        open('meson.build', 'w', encoding='utf-8').write(lib_cs_meson_template.format(**kwargs))
