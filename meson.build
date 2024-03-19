# This file is adapted from https://github.com/scipy/scipy/blob/main/meson.build
project(
    'pandas',
    'c', 'cpp', 'cython',
    version: run_command(['generate_version.py', '--print'], check: true).stdout().strip(),
    license: 'BSD-3',
    meson_version: '>=1.4.0',
    default_options: [
        'buildtype=release',
        'c_std=c11',
    ]
)

fs = import('fs')
py = import('python').find_installation(pure: false)
tempita = files('generate_pxi.py')
versioneer = files('generate_version.py')

add_project_arguments('-DNPY_NO_DEPRECATED_API=0', language : 'c')
add_project_arguments('-DNPY_NO_DEPRECATED_API=0', language : 'cpp')

# Allow supporting older numpys than the version compiled against
# Set the define to the min supported version of numpy for pandas
# e.g. right now this is targeting numpy 1.21+
add_project_arguments('-DNPY_TARGET_VERSION=NPY_1_21_API_VERSION', language : 'c')
add_project_arguments('-DNPY_TARGET_VERSION=NPY_1_21_API_VERSION', language : 'cpp')

if meson.get_compiler('c').get_id() == 'clang'
    #  error: pragma diagnostic pop could not pop, no matching push
    add_project_arguments('-Wno-unknown-pragmas', language: 'c')
endif

if meson.get_compiler('c').get_id() == 'msvc'
    # MSVC throws a lot of C4244 warnings around loss of data with warning_level=2
    # On other compilers we would have to add -Wconversion to see these
    # these expose a lot of potential issues in our code base
    add_project_arguments(['/W2'], language: 'c')
    add_project_arguments(['/W2'], language: 'cpp')
else
    add_project_arguments(['-Wall', '-Wextra'], language: 'c')
    add_project_arguments(['-Wall', '-Wextra'], language: 'cpp')
endif


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
        install_dir: py.get_install_dir() / 'pandas'
    )
    meson.add_dist_script(py, versioneer, '-o', '_version_meson.py')
endif

# Needed by pandas.test() when it looks for the pytest ini options
py.install_sources(
    'pyproject.toml',
    subdir: 'pandas'
)

subdir('pandas')
