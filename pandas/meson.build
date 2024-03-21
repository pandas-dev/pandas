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
inc_pd = include_directories('_libs/include')

fs.copyfile('__init__.py')

pd_dir = py.get_install_dir() / 'pandas'

# Generate version.py for sdist
meson.add_dist_script(
   ['_build_utils/gitversion.py', '--meson-dist', '--write',
     'pandas/_version.py']
)
if not fs.exists('_version.py')
  generate_version = custom_target(
    'generate-version',
    install: true,
    build_always_stale: true,
    build_by_default: true,
    output: '_version.py',
    input: '_build_utils/gitversion.py',
    command: [py, '@INPUT@', '--write', '@OUTPUT@'],
    install_dir: pd_dir,
    install_tag: 'python-runtime'
  )
else
  # When building from sdist, version.py exists and should be included
  py.install_sources(
    ['_version.py'],
    subdir : 'pandas'
  )
endif

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
    'conftest.py',
    'testing.py'
]
py.install_sources(top_level_py_list, subdir: 'pandas')
