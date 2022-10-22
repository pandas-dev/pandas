# This file is adapted from https://github.com/scipy/scipy/blob/main/meson.build
project(
    'pandas',
    'c', 'cpp', 'cython',
    version: '1.6.0.dev0',
    license: 'BSD-3',
    # TODO: bump when meson 0.64.0 comes out,
    # we are relying on 0.64.0 features
    meson_version: '>=0.63',
    default_options: [
        'buildtype=release',
        'c_std=c99'
    ]
)

py_mod = import('python')
fs = import('fs')
py = py_mod.find_installation('python')
py_dep = py.dependency()
tempita = files('generate_pxi.py')
versioneer = files('generate_version.py')

if fs.exists('_version_meson.py')
    py.install_sources('_version_meson.py', subdir: 'pandas')
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

subdir('pandas')
