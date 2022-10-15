# TODO: can this be removed, I wish meson copied .pyx source to the build dir automatically
tslibs_pxd_sources_list = [
    # List of cython sources e.g. .pyx, .pxd & __init__.py
    # Does NOT include .pxi.in
    '__init__.py',
    'base.pxd',
    'ccalendar.pxd',
    'conversion.pxd',
    'dtypes.pxd',
    'nattype.pxd',
    'np_datetime.pxd',
    'offsets.pxd',
    'parsing.pxd',
    'period.pxd',
    'timedeltas.pxd',
    'timestamps.pxd',
    'timezones.pxd',
    'tzconversion.pxd',
    'util.pxd',
]

foreach source: tslibs_pxd_sources_list
    source_pxd = configure_file(
        input: source,
        output: source,
        copy: true
    )
endforeach

tslibs_sources = {
    # Dict of extension name -> dict of {sources, include_dirs, and deps}
    # numpy include dir is implicitly included
    'base': {'sources': ['base.pyx']},
    'ccalendar': {'sources': ['ccalendar.pyx']},
    'dtypes': {'sources': ['dtypes.pyx']},
    'conversion': {'sources': ['conversion.pyx', 'src/datetime/np_datetime.c'],
                   'include_dirs': inc_datetime},
    'fields': {'sources': ['fields.pyx', 'src/datetime/np_datetime.c']},
    'nattype': {'sources': ['nattype.pyx']},
    'np_datetime': {'sources': ['np_datetime.pyx', 'src/datetime/np_datetime.c', 'src/datetime/np_datetime_strings.c'],
                    'include_dirs': inc_datetime},
    'offsets': {'sources': ['offsets.pyx', 'src/datetime/np_datetime.c'],
                'include_dirs': inc_datetime},
    'parsing': {'sources': ['parsing.pyx', '../src/parser/tokenizer.c'],
                'include_dirs': klib_include},
    'period': {'sources': ['period.pyx', 'src/datetime/np_datetime.c'],
                'include_dirs': inc_datetime},
    'strptime': {'sources': ['strptime.pyx', 'src/datetime/np_datetime.c'],
                'include_dirs': inc_datetime},
    'timedeltas': {'sources': ['timedeltas.pyx', 'src/datetime/np_datetime.c'],
                   'include_dirs': inc_datetime},
    'timestamps': {'sources': ['timestamps.pyx', 'src/datetime/np_datetime.c'],
                   'include_dirs': inc_datetime},
    'timezones': {'sources': ['timezones.pyx', 'src/datetime/np_datetime.c'],
                  'include_dirs': inc_datetime},
    'tzconversion': {'sources': ['tzconversion.pyx', 'src/datetime/np_datetime.c'],
                     'include_dirs': inc_datetime},
    'vectorized': {'sources': ['vectorized.pyx', 'src/datetime/np_datetime.c'],
                    'include_dirs': inc_datetime}
}

foreach ext_name, ext_dict : tslibs_sources
    py.extension_module(
        ext_name,
        ext_dict.get('sources'),
        include_directories: [inc_np] + ext_dict.get('include_dirs', ''),
        dependencies: ext_dict.get('deps', ''),
        subdir: 'pandas/_libs/tslibs',
        install: true
    )
endforeach

py.install_sources('__init__.py',
                    pure: false,
                    subdir: 'pandas/_libs/tslibs')
