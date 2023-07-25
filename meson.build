# This file is adapted from https://github.com/scipy/scipy/blob/main/meson.build
project(
    'pandas',
    'c', 'cpp', 'cython',
    version: run_command(['python', 'generate_version.py', '--print'], check: true).stdout().strip(),
    license: 'BSD-3',
    meson_version: '>=1.0.1',
    default_options: [
        # TODO: investigate, does meson try to compile against debug Python
        # when buildtype = debug, this seems to be causing problems on CI
        # where provided Python is not compiled in debug mode
        'buildtype=release',
        # TODO: Reactivate werror, some warnings on Windows
        #'werror=true',
        'c_std=c99'
    ]
)

py_mod = import('python')
fs = import('fs')
py = py_mod.find_installation('python')
py_dep = py.dependency()
tempita = files('generate_pxi.py')
versioneer = files('generate_version.py')


add_project_arguments('-DNPY_NO_DEPRECATED_API=0', language : 'c')
add_project_arguments('-DNPY_NO_DEPRECATED_API=0', language : 'cpp')

# Allow supporting older numpys than the version compiled against
# Set the define to the min supported version of numpy for pandas
# e.g. right now this is targeting numpy 1.21+
add_project_arguments('-DNPY_TARGET_VERSION=NPY_1_21_API_VERSION', language : 'c')
add_project_arguments('-DNPY_TARGET_VERSION=NPY_1_21_API_VERSION', language : 'cpp')


if fs.exists('_version_meson.py')
    py.install_sources('_version_meson.py', pure: false, subdir: 'pandas')
else
    custom_target('write_version_file',
        output: '_version_meson.py',
        command: [
            py, versioneer, '-o', '@OUTPUT@'
        ],
        build_by_default: true,
        build_always_stale: true,
        install: true,
        install_dir: py.get_install_dir(pure: false) / 'pandas'
    )
    meson.add_dist_script(py, versioneer, '-o', '_version_meson.py')
endif

# Needed by pandas.test() when it looks for the pytest ini options
py.install_sources('pyproject.toml', pure: false, subdir: 'pandas')
subdir('pandas')
