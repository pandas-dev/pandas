tslibs_sources = {
    # Dict of extension name -> dict of {sources, include_dirs, and deps}
    # numpy include dir is implicitly included
    'base': {'sources': ['base.pyx']},
    'ccalendar': {'sources': ['ccalendar.pyx']},
    'dtypes': {'sources': ['dtypes.pyx']},
    'conversion': {'sources': ['conversion.pyx']},
    'fields': {'sources': ['fields.pyx']},
    'nattype': {'sources': ['nattype.pyx']},
    'np_datetime': {'sources': ['np_datetime.pyx']},
    'offsets': {'sources': ['offsets.pyx']},
    'parsing': {'sources': ['parsing.pyx', '../src/parser/tokenizer.c']},
    'period': {'sources': ['period.pyx']},
    'strptime': {'sources': ['strptime.pyx']},
    'timedeltas': {'sources': ['timedeltas.pyx']},
    'timestamps': {'sources': ['timestamps.pyx']},
    'timezones': {'sources': ['timezones.pyx']},
    'tzconversion': {'sources': ['tzconversion.pyx']},
    'vectorized': {'sources': ['vectorized.pyx']},
}

cython_args = [
    '--include-dir',
    meson.current_build_dir(),
    '-X always_allow_keywords=true'
]
if get_option('buildtype') == 'debug'
    cython_args += ['--gdb']
endif

foreach ext_name, ext_dict : tslibs_sources
    py.extension_module(
        ext_name,
        ext_dict.get('sources'),
        cython_args: cython_args,
        include_directories: [inc_np, inc_pd],
        dependencies: ext_dict.get('deps', ''),
        subdir: 'pandas/_libs/tslibs',
        install: true
    )
endforeach

sources_to_install = [
    '__init__.py',
    'ccalendar.pyi',
    'conversion.pyi',
    'dtypes.pyi',
    'fields.pyi',
    'nattype.pyi',
    'np_datetime.pyi',
    'offsets.pyi',
    'parsing.pyi',
    'period.pyi',
    'strftime.py',
    'strptime.pyi',
    'timedeltas.pyi',
    'timestamps.pyi',
    'timezones.pyi',
    'tzconversion.pyi',
    'vectorized.pyi'
]

foreach source: sources_to_install
    py.install_sources(
        source,
        subdir: 'pandas/_libs/tslibs'
    )
endforeach
