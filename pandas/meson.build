incdir_numpy = run_command(py,
  [
    '-c',
    'import os; os.chdir(".."); import numpy; print(numpy.get_include())'
  ],
  check: true
).stdout().strip()

inc_np = include_directories(incdir_numpy)
klib_include = include_directories('_libs/src/klib')
inc_datetime = include_directories('_libs/tslibs')

configure_file(
    input: '__init__.py',
    output: '__init__.py',
    copy: true
)
subdir('_libs')
subdir('io')
# TODO: this works but needs a patch in meson-python
# meson-python only looks in the meson-info/install_plan.json not in meson-info/installed.json
# However, subdirs are not listed in the install_plan only in installed.json
subdirs_list = [
    '_config',
    '_libs',
    '_testing',
    'api',
    'arrays',
    'compat',
    'core',
    'errors',
    'plotting',
    'tests',
    'tseries',
    'util'
]
foreach subdir: subdirs_list
    install_subdir(subdir, install_dir: py.get_install_dir(pure: false) / 'pandas')
endforeach
top_level_py_list = [
    '__init__.py',
    '_typing.py',
    '_version.py',
    'conftest.py',
    'testing.py'
]
foreach file: top_level_py_list
    py.install_sources(file,
                       pure: false,
                       subdir: 'pandas')
endforeach
