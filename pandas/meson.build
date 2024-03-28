# NumPy 2.0 can be detected with meson >= 0.14.0
np_dep = dependency('numpy', required : false)

if np_dep.found()
    inc_np = []
else
    incdir_numpy = run_command(py,
      [
        '-c',
        '''
import os
import numpy as np
try:
    # Check if include directory is inside the pandas dir
    # e.g. a venv created inside the pandas dir
    # If so, convert it to a relative path
    incdir = os.path.relpath(np.get_include())
except Exception:
    incdir = np.get_include()
print(incdir)
         '''
      ],
      check: true
    ).stdout().strip()

    inc_np = include_directories(incdir_numpy)
endif

inc_pd = include_directories('_libs/include')

fs.copyfile('__init__.py')

subdir('_libs')

subdirs_list = [
    '_config',
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
    install_subdir(subdir, install_dir: py.get_install_dir() / 'pandas')
endforeach

top_level_py_list = [
    '__init__.py',
    '_typing.py',
    '_version.py',
    'conftest.py',
    'testing.py'
]
py.install_sources(top_level_py_list, subdir: 'pandas')
