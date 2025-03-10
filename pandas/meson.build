incdir_numpy = run_command(
    py,
    [
        '-c',
        '''
import os

try:
    import numpy as np
    base_incdir = np.get_include()
except Exception:
    base_incdir = os.getenv('NUMPY_INC_DIR')

try:
    # Check if include directory is inside the pandas dir
    # e.g. a venv created inside the pandas dir
    # If so, convert it to a relative path
    incdir = os.path.relpath(base_incdir)
except Exception:
    incdir = base_incdir
print(incdir)
     ''',
    ],
    env: {'NUMPY_INC_DIR': get_option('numpy_inc_dir')},
    check: true,
).stdout().strip()

inc_np = include_directories(incdir_numpy)
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
    'util',
]
foreach subdir : subdirs_list
    install_subdir(subdir, install_dir: py.get_install_dir() / 'pandas')
endforeach

top_level_py_list = [
    '__init__.py',
    '_typing.py',
    '_version.py',
    'conftest.py',
    'testing.py',
]
py.install_sources(top_level_py_list, subdir: 'pandas')
