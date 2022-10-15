_algos_take_helper = custom_target('algos_take_helper_pxi',
    output: 'algos_take_helper.pxi',
    input: 'algos_take_helper.pxi.in',
    command: [
        py, tempita, '@INPUT@', '-o', '@OUTDIR@'
    ],
    # TODO: remove these two below lines
    # Weird bug in meson that only repros on my potato computer
    # (possibly b/c of low number of threads)?
    # The first custom_target is never built for some reason
    # so algos.pyx will error out later in the build.
    build_by_default: true,
    build_always_stale: true
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
_khash_primitive_helper_dep = declare_dependency(sources: _khash_primitive_helper)
# TODO: can this be removed, I wish meson copied .pyx source to the build dir automatically
# The reason we can't build the pyx files inplace and copy to build dir is because
# the generated pxi files cannot be written to the source directory.
# (Meson only supports out of tree builds)
cython_sources_list = [
    # List of cython sources e.g. .pyx, .pxd & __init__.py
    # Does NOT include .pxi.in
    '__init__.py',
    'algos.pxd',
    'algos.pyx',
    'arrays.pxd',
    'dtypes.pxd',
    'hashtable.pxd',
    'hashtable.pyx',
    'index.pyx',
    'indexing.pyx',
    'internals.pyx',
    'interval.pyx',
    'join.pyx',
    'khash.pxd',
    'lib.pxd',
    'missing.pxd',
    'parsers.pyx',
    'sparse.pyx',
    'testing.pyx',
    'tslib.pyx',
    'util.pxd',
]
cython_sources = {}

foreach source: cython_sources_list
    source_pyx = configure_file(
        input: source,
        output: source,
        copy: true
    )
    cython_sources += {source: source_pyx}
endforeach

subdir('tslibs')

libs_sources = {
    # Dict of extension name -> dict of {sources, include_dirs, and deps}
    # numpy include dir is implicitly included
    'algos': {'sources': [cython_sources['algos.pyx'], _algos_common_helper, _algos_take_helper, _khash_primitive_helper],
              'include_dirs': klib_include},
    'arrays': {'sources': ['arrays.pyx']},
    'groupby': {'sources': ['groupby.pyx']},
    'hashing': {'sources': ['hashing.pyx']},
    'hashtable': {'sources': [cython_sources['hashtable.pyx'], _khash_primitive_helper, _hashtable_class_helper, _hashtable_func_helper],
                  'include_dirs': klib_include},
    'index': {'sources': [cython_sources['index.pyx'], _index_class_helper],
              'include_dirs': [klib_include, 'tslibs']},
    'indexing': {'sources': ['indexing.pyx']},
    'internals': {'sources': ['internals.pyx']},
    'interval': {'sources': [cython_sources['interval.pyx'], _intervaltree_helper],
                'include_dirs': [klib_include, 'tslibs']},
    'join': {'sources': [cython_sources['join.pyx'], _khash_primitive_helper],
             'include_dirs': klib_include,
             'deps': _khash_primitive_helper_dep},
    'lib': {'sources': ['lib.pyx', 'src/parser/tokenizer.c'],
            'include_dirs': [klib_include, inc_datetime]},
    'missing': {'sources': ['missing.pyx'],
                'include_dirs': [inc_datetime]},
    'parsers': {'sources': [cython_sources['parsers.pyx'], 'src/parser/tokenizer.c', 'src/parser/io.c'],
                'include_dirs': [klib_include, 'src'],
                'deps': _khash_primitive_helper_dep},
    'json': {'sources': ['src/ujson/python/ujson.c',
                         'src/ujson/python/objToJSON.c',
                         'src/ujson/python/date_conversions.c',
                         'src/ujson/python/JSONtoObj.c',
                         'src/ujson/lib/ultrajsonenc.c',
                         'src/ujson/lib/ultrajsondec.c',
                         'tslibs/src/datetime/np_datetime.c',
                         'tslibs/src/datetime/np_datetime_strings.c'],
             'include_dirs': [inc_datetime, 'src/ujson/lib', 'src/ujson/python']},
    'reduction': {'sources': ['reduction.pyx']},
    'ops': {'sources': ['ops.pyx']},
    'ops_dispatch': {'sources': ['ops_dispatch.pyx']},
    'properties': {'sources': ['properties.pyx']},
    'reshape': {'sources': ['reshape.pyx']},
    'sparse': {'sources': [cython_sources['sparse.pyx'], _sparse_op_helper]},
    'tslib': {'sources': ['tslib.pyx', 'tslibs/src/datetime/np_datetime.c'],
              'include_dirs': inc_datetime},
    'testing': {'sources': ['testing.pyx']},
    'writers': {'sources': ['writers.pyx']}
}


foreach ext_name, ext_dict : libs_sources
    py.extension_module(
        ext_name,
        ext_dict.get('sources'),
        include_directories: [inc_np] + ext_dict.get('include_dirs', ''),
        dependencies: ext_dict.get('deps', ''),
        subdir: 'pandas/_libs',
        install: true
    )
endforeach

py.install_sources('__init__.py',
                    pure: false,
                    subdir: 'pandas/_libs')

subdir('window')
