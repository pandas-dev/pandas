project(
    'xsimd',
    'cpp',
    meson_version: '>=0.58.0',
    license: 'BSD-3-Clause',
    version: '14.2.0',
)

xsimd_inc = include_directories('include')

xsimd_dep = declare_dependency(include_directories: xsimd_inc)
meson.override_dependency('xsimd', xsimd_dep)
