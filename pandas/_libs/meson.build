_algos_take_helper = custom_target('algos_take_helper_pxi',
    output: 'algos_take_helper.pxi',
    input: 'algos_take_helper.pxi.in',
    command: [
        py, tempita, '@INPUT@', '-o', '@OUTDIR@'
    ]
)
_algos_common_helper = custom_target('algos_common_helper_pxi',
    output: 'algos_common_helper.pxi',
    input: 'algos_common_helper.pxi.in',
    command: [
        py, tempita, '@INPUT@', '-o', '@OUTDIR@'
    ]
)
_khash_primitive_helper = custom_target('khash_primitive_helper_pxi',
    output: 'khash_for_primitive_helper.pxi',
    input: 'khash_for_primitive_helper.pxi.in',
    command: [
        py, tempita, '@INPUT@', '-o', '@OUTDIR@'
    ]
)
_hashtable_class_helper = custom_target('hashtable_class_helper_pxi',
    output: 'hashtable_class_helper.pxi',
    input: 'hashtable_class_helper.pxi.in',
    command: [
        py, tempita, '@INPUT@', '-o', '@OUTDIR@'
    ]
)
_hashtable_func_helper = custom_target('hashtable_func_helper_pxi',
    output: 'hashtable_func_helper.pxi',
    input: 'hashtable_func_helper.pxi.in',
    command: [
        py, tempita, '@INPUT@', '-o', '@OUTDIR@'
    ]
)
_index_class_helper = custom_target('index_class_helper_pxi',
    output: 'index_class_helper.pxi',
    input: 'index_class_helper.pxi.in',
    command: [
        py, tempita, '@INPUT@', '-o', '@OUTDIR@'
    ]
)
_sparse_op_helper = custom_target('sparse_op_helper_pxi',
    output: 'sparse_op_helper.pxi',
    input: 'sparse_op_helper.pxi.in',
    command: [
        py, tempita, '@INPUT@', '-o', '@OUTDIR@'
    ]
)
_intervaltree_helper = custom_target('intervaltree_helper_pxi',
    output: 'intervaltree.pxi',
    input: 'intervaltree.pxi.in',
    command: [
        py, tempita, '@INPUT@', '-o', '@OUTDIR@'
    ]
)

_algos_pxi_dep = declare_dependency(sources: [_algos_take_helper, _algos_common_helper])
_khash_pxi_dep = declare_dependency(sources: _khash_primitive_helper)
_hashtable_pxi_dep = declare_dependency(
    sources: [_hashtable_class_helper, _hashtable_func_helper]
)
_index_pxi_dep = declare_dependency(sources: _index_class_helper)
_intervaltree_pxi_dep = declare_dependency(sources: _intervaltree_helper)
_sparse_pxi_dep = declare_dependency(sources: _sparse_op_helper)


subdir('tslibs')

libs_sources = {
    # Dict of extension name -> dict of {sources, include_dirs, and deps}
    # numpy include dir is implicitly included
    'algos': {'sources': ['algos.pyx'],
              'deps': [_khash_pxi_dep, _algos_pxi_dep]},
    'arrays': {'sources': ['arrays.pyx']},
    'groupby': {'sources': ['groupby.pyx']},
    'hashing': {'sources': ['hashing.pyx']},
    'hashtable': {'sources': ['hashtable.pyx'],
                  'deps': [_khash_pxi_dep, _hashtable_pxi_dep]},
    'index': {'sources': ['index.pyx'],
              'deps': [_khash_pxi_dep, _index_pxi_dep]},
    'indexing': {'sources': ['indexing.pyx']},
    'internals': {'sources': ['internals.pyx']},
    'interval': {'sources': ['interval.pyx'],
                 'deps': [_khash_pxi_dep, _intervaltree_pxi_dep]},
    'join': {'sources': ['join.pyx'],
             'deps': [_khash_pxi_dep]},
    'lib': {'sources': ['lib.pyx', 'src/parser/tokenizer.c']},
    'missing': {'sources': ['missing.pyx']},
    'pandas_datetime': {'sources': ['src/vendored/numpy/datetime/np_datetime.c',
                                    'src/vendored/numpy/datetime/np_datetime_strings.c',
                                    'src/datetime/date_conversions.c',
                                    'src/datetime/pd_datetime.c']},
    'pandas_parser': {'sources': ['src/parser/tokenizer.c',
                                  'src/parser/io.c',
                                  'src/parser/pd_parser.c']},
    'parsers': {'sources': ['parsers.pyx', 'src/parser/tokenizer.c', 'src/parser/io.c'],
                'deps': [_khash_pxi_dep]},
    'json': {'sources': ['src/vendored/ujson/python/ujson.c',
                         'src/vendored/ujson/python/objToJSON.c',
                         'src/vendored/ujson/python/JSONtoObj.c',
                         'src/vendored/ujson/lib/ultrajsonenc.c',
                         'src/vendored/ujson/lib/ultrajsondec.c']},
    'ops': {'sources': ['ops.pyx']},
    'ops_dispatch': {'sources': ['ops_dispatch.pyx']},
    'properties': {'sources': ['properties.pyx']},
    'reshape': {'sources': ['reshape.pyx']},
    'sas': {'sources': ['sas.pyx']},
    'byteswap': {'sources': ['byteswap.pyx']},
    'sparse': {'sources': ['sparse.pyx'],
               'deps': [_sparse_pxi_dep]},
    'tslib': {'sources': ['tslib.pyx']},
    'testing': {'sources': ['testing.pyx']},
    'writers': {'sources': ['writers.pyx']}
}

cython_args = [
    '--include-dir',
    meson.current_build_dir(),
    '-X always_allow_keywords=true'
]
if get_option('buildtype') == 'debug'
    cython_args += ['--gdb']
endif

foreach ext_name, ext_dict : libs_sources
    py.extension_module(
        ext_name,
        ext_dict.get('sources'),
        cython_args: cython_args,
        include_directories: [inc_np, inc_pd],
        dependencies: ext_dict.get('deps', ''),
        subdir: 'pandas/_libs',
        install: true
    )
endforeach

# Basically just __init__.py and the .pyi files
sources_to_install = [
    '__init__.py',
    'algos.pyi',
    'arrays.pyi',
    'byteswap.pyi',
    'groupby.pyi',
    'hashing.pyi',
    'hashtable.pyi',
    'index.pyi',
    'indexing.pyi',
    'internals.pyi',
    'interval.pyi',
    'join.pyi',
    'json.pyi',
    'lib.pyi',
    'missing.pyi',
    'ops.pyi',
    'ops_dispatch.pyi',
    'parsers.pyi',
    'properties.pyi',
    'reshape.pyi',
    'sas.pyi',
    'sparse.pyi',
    'testing.pyi',
    'tslib.pyi',
    'writers.pyi'
]

foreach source: sources_to_install
    py.install_sources(
        source,
        subdir: 'pandas/_libs'
    )
endforeach

subdir('window')
