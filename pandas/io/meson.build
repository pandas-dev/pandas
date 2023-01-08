subdirs_list = [
    # exclude sas, since it contains extension modules
    # and has its own meson.build
    'clipboard',
    'excel',
    'formats',
    'json',
    'parsers'
]
foreach subdir: subdirs_list
    install_subdir(subdir, install_dir: py.get_install_dir(pure: false) / 'pandas/io')
endforeach
top_level_py_list = [
    '__init__.py',
    '_util.py',
    'api.py',
    'clipboards.py',
    'common.py',
    'feather_format.py',
    'gbq.py',
    'html.py',
    'orc.py',
    'parquet.py',
    'pickle.py',
    'pytables.py',
    'spss.py',
    'sql.py',
    'stata.py',
    'xml.py'
]
foreach file: top_level_py_list
    py.install_sources(file,
                       pure: false,
                       subdir: 'pandas/io')
endforeach
subdir('sas')
