# This file is adapted from https://github.com/scipy/scipy/blob/main/meson.build
project(
    'pandas',
    'c', 'cpp', 'cython',
    version: '2.0.0.dev0',
    license: 'BSD-3',
    meson_version: '>=0.64',
    default_options: [
        # TODO: investigate, does meson try to compile against debug Python
        # when buildtype = debug, this seems to be causing problems on CI
        # where provided Python is not compiled in debug mode
        'buildtype=release',
        'werror=true',
        'c_std=c99'
    ]
)

add_project_arguments('-DNPY_NO_DEPRECATED_API=0', language : 'c')
add_project_arguments('-DNPY_NO_DEPRECATED_API=0', language : 'cpp')

# This is a cython bug
# TODO: Remove once cython/cython#4804 addressed
add_project_arguments('-DNDEBUG', language : 'c')
add_project_arguments('-DNDEBUG', language : 'cpp')

py_mod = import('python')
fs = import('fs')
py = py_mod.find_installation('python')
py_dep = py.dependency()
tempita = files('generate_pxi.py')
versioneer = files('generate_version.py')

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
