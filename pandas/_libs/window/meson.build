cy_args = ['-X always_allow_keywords=true']
# Use shared utility code to reduce wheel sizes
# copied from https://github.com/scikit-learn/scikit-learn/pull/31151/files
if cy.version().version_compare('>=3.1.0')
    cython_args += ['--shared=pandas._libs._cyutility']
endif

py.extension_module(
    'aggregations',
    ['aggregations.pyx'],
    cython_args: cy_args,
    include_directories: [inc_np, inc_pd],
    subdir: 'pandas/_libs/window',
    override_options: ['cython_language=cpp'],
    install: true,
)

py.extension_module(
    'indexers',
    ['indexers.pyx'],
    cython_args: cy_args,
    include_directories: [inc_np, inc_pd],
    subdir: 'pandas/_libs/window',
    install: true,
)

sources_to_install = ['__init__.py', 'aggregations.pyi', 'indexers.pyi']

foreach source : sources_to_install
    py.install_sources(source, subdir: 'pandas/_libs/window')
endforeach
