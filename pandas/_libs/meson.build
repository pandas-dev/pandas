_algos_take_helper = custom_target(
    'algos_take_helper_pxi',
    output: 'algos_take_helper.pxi',
    input: 'algos_take_helper.pxi.in',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@'],
)
_algos_common_helper = custom_target(
    'algos_common_helper_pxi',
    output: 'algos_common_helper.pxi',
    input: 'algos_common_helper.pxi.in',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@'],
)
_khash_primitive_helper = custom_target(
    'khash_primitive_helper_pxi',
    output: 'khash_for_primitive_helper.pxi',
    input: 'khash_for_primitive_helper.pxi.in',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@'],
)
_hashtable_class_helper = custom_target(
    'hashtable_class_helper_pxi',
    output: 'hashtable_class_helper.pxi',
    input: 'hashtable_class_helper.pxi.in',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@'],
)
_hashtable_func_helper = custom_target(
    'hashtable_func_helper_pxi',
    output: 'hashtable_func_helper.pxi',
    input: 'hashtable_func_helper.pxi.in',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@'],
)
_index_class_helper = custom_target(
    'index_class_helper_pxi',
    output: 'index_class_helper.pxi',
    input: 'index_class_helper.pxi.in',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@'],
)
_sparse_op_helper = custom_target(
    'sparse_op_helper_pxi',
    output: 'sparse_op_helper.pxi',
    input: 'sparse_op_helper.pxi.in',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@'],
)
_intervaltree_helper = custom_target(
    'intervaltree_helper_pxi',
    output: 'intervaltree.pxi',
    input: 'intervaltree.pxi.in',
    command: [py, tempita, '@INPUT@', '-o', '@OUTDIR@'],
)
_khash_primitive_helper_dep = declare_dependency(
    sources: _khash_primitive_helper,
)

m_dep = cc.find_library('m', required: false)
fast_float = subproject('fast_float')
fast_float_dep = fast_float.get_variable('fast_float_dep')
xsimd_dep = dependency('xsimd', version: '>=14.2')

subdir('simd')

# Floating-point std::to_chars needs libc++ >= deployment-target ~13.3 on
# macOS, and __cpp_lib_to_chars is unreliable on Apple clang, so try-compile
# instead of checking the feature-test macro. Without it, double_repr.cpp
# falls back to PyOS_double_to_string.
cpp = meson.get_compiler('cpp')
if cpp.get_argument_syntax() == 'msvc'
    _fp_to_chars_std_args = ['/std:c++17']
else
    _fp_to_chars_std_args = ['-std=c++17']
endif
fp_to_chars_test = '''
#include <charconv>
int main() {
    char buf[64];
    auto res = std::to_chars(buf, buf + sizeof(buf), 1.5,
                             std::chars_format::scientific);
    auto res32 = std::to_chars(buf, buf + sizeof(buf), 1.5f,
                               std::chars_format::scientific);
    return (res.ec == std::errc() && res32.ec == std::errc()) ? 0 : 1;
}
'''
if cpp.links(
    fp_to_chars_test,
    name: 'floating-point std::to_chars',
    args: _fp_to_chars_std_args,
)
    double_repr_cpp_args = ['-DPD_HAVE_FP_TO_CHARS=1']
else
    double_repr_cpp_args = []
endif

subdir('tslibs')

libs_sources = {
    # Dict of extension name -> dict of {sources, include_dirs, and deps}
    # numpy include dir is implicitly included
    'algos': {
        'sources': ['algos.pyx', _algos_common_helper, _algos_take_helper],
        'deps': [_khash_primitive_helper_dep, m_dep],
    },
    'arrays': {'sources': ['arrays.pyx']},
    'groupby': {'sources': ['groupby.pyx'], 'deps': [m_dep]},
    'hashing': {'sources': ['hashing.pyx']},
    'hashtable': {
        'sources': [
            'hashtable.pyx',
            _hashtable_class_helper,
            _hashtable_func_helper,
        ],
        'deps': _khash_primitive_helper_dep,
    },
    'index': {
        'sources': ['index.pyx', _index_class_helper],
        'deps': _khash_primitive_helper_dep,
    },
    'indexing': {'sources': ['indexing.pyx']},
    'internals': {'sources': ['internals.pyx']},
    'interval': {
        'sources': ['interval.pyx', _intervaltree_helper],
        'deps': _khash_primitive_helper_dep,
    },
    'join': {
        'sources': ['join.pyx', _khash_primitive_helper],
        'deps': _khash_primitive_helper_dep,
    },
    'lib': {
        'sources': [
            'lib.pyx',
            'src/parser/tokenizer.c',
            'src/parser/fast_float_strtod.cpp',
            'src/parser/pd_strtoi.cpp',
        ],
        'deps': [fast_float_dep],
    },
    'missing': {'sources': ['missing.pyx']},
    'pandas_datetime': {
        'sources': [
            'src/vendored/numpy/datetime/np_datetime.c',
            'src/vendored/numpy/datetime/np_datetime_strings.c',
            'src/datetime/date_conversions.c',
            'src/datetime/pd_datetime.c',
        ],
    },
    'pandas_parser': {
        'sources': [
            'src/parser/tokenizer.c',
            'src/parser/fast_float_strtod.cpp',
            'src/parser/pd_strtoi.cpp',
            'src/parser/io.c',
            'src/parser/pd_parser.c',
        ],
        'deps': [fast_float_dep],
    },
    'parsers': {
        'sources': [
            'parsers.pyx',
            'src/parser/tokenizer.c',
            'src/parser/fast_float_strtod.cpp',
            'src/parser/pd_strtoi.cpp',
            'src/parser/io.c',
        ],
        'deps': [fast_float_dep, _khash_primitive_helper_dep],
    },
    '_ujson': {
        'sources': [
            'src/vendored/ujson/python/ujson.c',
            'src/vendored/ujson/python/objToJSON.c',
            'src/vendored/ujson/python/JSONtoObj.c',
            'src/vendored/ujson/lib/ultrajsonenc.c',
            'src/vendored/ujson/lib/ultrajsondec.c',
        ],
        'deps': [m_dep],
    },
    'ops': {'sources': ['ops.pyx']},
    'ops_dispatch': {'sources': ['ops_dispatch.pyx']},
    'properties': {'sources': ['properties.pyx']},
    'reshape': {'sources': ['reshape.pyx']},
    'sas': {'sources': ['sas.pyx']},
    'byteswap': {'sources': ['byteswap.pyx']},
    'sparse': {'sources': ['sparse.pyx', _sparse_op_helper], 'deps': [m_dep]},
    'tslib': {'sources': ['tslib.pyx']},
    'testing': {'sources': ['testing.pyx']},
    'writers': {
        'sources': ['writers.pyx', 'src/writer/double_repr.cpp'],
        'cpp_args': double_repr_cpp_args,
    },
}

cython_args = [
    '--include-dir',
    meson.current_build_dir(),
    '-X always_allow_keywords=true',
]
if get_option('buildtype') == 'debug'
    cython_args += ['--gdb']
endif

# Use shared utility code to reduce wheel sizes
# copied from https://github.com/scikit-learn/scikit-learn/pull/31151/files
if cy.version().version_compare('>=3.1.0')
    cython_args += ['--shared=pandas._libs._cyutility']
endif

foreach ext_name, ext_dict : libs_sources
    py.extension_module(
        ext_name,
        ext_dict.get('sources'),
        cython_args: cython_args,
        cpp_args: ext_dict.get('cpp_args', []),
        include_directories: [inc_np, inc_pd],
        dependencies: ext_dict.get('deps', []),
        subdir: 'pandas/_libs',
        install: true,
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
    '_ujson.pyi',
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
    'writers.pyi',
]

foreach source : sources_to_install
    py.install_sources(source, subdir: 'pandas/_libs')
endforeach

subdir('window')
