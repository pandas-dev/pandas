# Allow the NumPy include directory to be supplied via a meson property, which
# is needed when cross-compiling and the target's numpy cannot be imported on
# the build machine (GH#55305, GH#62883). Fall back to querying numpy directly.
incdir_numpy = meson.get_external_property('numpy-include-dir', 'not-given')
if incdir_numpy == 'not-given'
    incdir_numpy = run_command(
        py,
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
     ''',
        ],
        check: true,
    ).stdout().strip()
endif

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
    'conftest.py',
    'testing.py',
]
py.install_sources(top_level_py_list, subdir: 'pandas')
