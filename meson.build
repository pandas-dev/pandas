# This file is adapted from https://github.com/scipy/scipy/blob/main/meson.build
project(
    'pandas',
    'c',
    'cpp',
    'cython',
    version: run_command(['generate_version.py', '--print'], check: true).stdout().strip(),
    license: 'BSD-3',
    meson_version: '>=1.2.1',
    default_options: ['buildtype=release', 'c_std=c11', 'warning_level=2'],
)

fs = import('fs')
py = import('python').find_installation(pure: false)
tempita = files('generate_pxi.py')
versioneer = files('generate_version.py')


add_project_arguments('-DNPY_NO_DEPRECATED_API=0', language: 'c')
add_project_arguments('-DNPY_NO_DEPRECATED_API=0', language: 'cpp')
# Enables settings __module__ on cdef classes
# https://github.com/cython/cython/issues/7231
add_project_arguments('-DCYTHON_USE_TYPE_SPECS=1', language: 'c')

# Allow supporting older numpys than the version compiled against
# Set the define to the min supported version of numpy for pandas
# e.g. right now this is targeting numpy 1.21+
add_project_arguments('-DNPY_TARGET_VERSION=NPY_1_21_API_VERSION', language: 'c')
add_project_arguments(
    '-DNPY_TARGET_VERSION=NPY_1_21_API_VERSION',
    language: 'cpp',
)


if fs.exists('_version_meson.py')
    py.install_sources('_version_meson.py', subdir: 'pandas')
else
    custom_target(
        'write_version_file',
        output: '_version_meson.py',
        command: [py, versioneer, '-o', '@OUTPUT@'],
        build_by_default: true,
        build_always_stale: true,
        install: true,
        install_dir: py.get_install_dir() / 'pandas',
    )
    meson.add_dist_script(py, versioneer, '-o', '_version_meson.py')
endif

cy = meson.get_compiler('cython')
if cy.version().version_compare('>=3.1.0')
    add_project_arguments('-Xfreethreading_compatible=true', language: 'cython')

    # Use shared utility code to reduce wheel sizes
    # copied from https://github.com/scikit-learn/scikit-learn/pull/31151/files
    cy = find_program(cy.cmd_array()[0])
    cython_shared_src = custom_target(
        install: false,
        output: '_cyutility.c',
        command: [
            cy,
            '-3',
            '-Xfreethreading_compatible=true',
            '--fast-fail',
            '--generate-shared=' + meson.current_build_dir() / '_cyutility.c',
        ],
    )

    py.extension_module(
        '_cyutility',
        cython_shared_src,
        subdir: 'pandas/_libs',
        install: true,
    )
endif

# Needed by pandas.test() when it looks for the pytest ini options
py.install_sources('pyproject.toml', subdir: 'pandas')

subdir('pandas')
