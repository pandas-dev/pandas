# This file is adapted from https://github.com/scipy/scipy/blob/main/meson.build
project(
    'pandas',
    'c', 'cpp', 'cython',
    version: run_command(
        ['python', 'pandas/_build_utils/gitversion.py'],
        check: true).stdout().strip(),
    license: 'BSD-3',
    meson_version: '>=1.2.1',
    default_options: [
        'buildtype=release',
        'c_std=c11',
        'warning_level=2',
    ]
)

fs = import('fs')
py = import('python').find_installation(pure: false)
tempita = files('generate_pxi.py')

add_project_arguments('-DNPY_NO_DEPRECATED_API=0', language : 'c')
add_project_arguments('-DNPY_NO_DEPRECATED_API=0', language : 'cpp')

# Allow supporting older numpys than the version compiled against
# Set the define to the min supported version of numpy for pandas
# e.g. right now this is targeting numpy 1.21+
add_project_arguments('-DNPY_TARGET_VERSION=NPY_1_21_API_VERSION', language : 'c')
add_project_arguments('-DNPY_TARGET_VERSION=NPY_1_21_API_VERSION', language : 'cpp')

# Needed by pandas.test() when it looks for the pytest ini options
py.install_sources(
    'pyproject.toml',
    subdir: 'pandas'
)

subdir('pandas')
