py.extension_module(
    'aggregations',
    ['aggregations.pyx'],
    include_directories: [inc_np, '../src'],
    dependencies: [py_dep],
    subdir: 'pandas/_libs/window',
    override_options : ['cython_language=cpp'],
    install: true
)

py.extension_module(
    'indexers',
    ['indexers.pyx'],
    include_directories: [inc_np],
    dependencies: [py_dep],
    subdir: 'pandas/_libs/window',
    install: true
)
