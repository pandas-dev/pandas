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

import argparse

meson_executable_template = '''project('{project_name}', {language},
  version : '{version}',
  default_options : [{default_options}])

executable('{executable}',
           {sourcespec},{depspec}
           install : true)
'''


meson_jar_template = '''project('{project_name}', '{language}',
  version : '{version}',
  default_options : [{default_options}])

jar('{executable}',
    {sourcespec},{depspec}
    main_class: '{main_class}',
    install : true)
'''


def create_meson_build(options: argparse.Namespace) -> None:
    if options.type != 'executable':
        raise SystemExit('\nGenerating a meson.build file from existing sources is\n'
                         'supported only for project type "executable".\n'
                         'Run meson init in an empty directory to create a sample project.')
    default_options = ['warning_level=3']
    if options.language == 'cpp':
        # This shows how to set this very common option.
        default_options += ['cpp_std=c++14']
    # If we get a meson.build autoformatter one day, this code could
    # be simplified quite a bit.
    formatted_default_options = ', '.join(f"'{x}'" for x in default_options)
    sourcespec = ',\n           '.join(f"'{x}'" for x in options.srcfiles)
    depspec = ''
    if options.deps:
        depspec = '\n           dependencies : [\n              '
        depspec += ',\n              '.join(f"dependency('{x}')"
                                            for x in options.deps.split(','))
        depspec += '],'
    if options.language != 'java':
        language = f"'{options.language}'" if options.language != 'vala' else ['c', 'vala']
        content = meson_executable_template.format(project_name=options.name,
                                                   language=language,
                                                   version=options.version,
                                                   executable=options.executable,
                                                   sourcespec=sourcespec,
                                                   depspec=depspec,
                                                   default_options=formatted_default_options)
    else:
        content = meson_jar_template.format(project_name=options.name,
                                            language=options.language,
                                            version=options.version,
                                            executable=options.executable,
                                            main_class=options.name,
                                            sourcespec=sourcespec,
                                            depspec=depspec,
                                            default_options=formatted_default_options)
    open('meson.build', 'w', encoding='utf-8').write(content)
    print('Generated meson.build file:\n\n' + content)
