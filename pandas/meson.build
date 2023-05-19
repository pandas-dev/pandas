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

fs.copyfile('__init__.py')

subdir('_libs')

subdirs_list = [
    '_config',
    '_libs',
    '_testing',
    'api',
    'arrays',
    'compat',
    'core',
    'errors',
    'io',
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
