py.extension_module(
    '_sas',
    ['sas.pyx'],
    include_directories: [inc_np],
    dependencies: [py_dep],
    # The file is named sas.pyx but we want the
    # extension module to be named _sas
    cython_args: ['--module-name=pandas.io.sas._sas'],
    subdir: 'pandas/io/sas',
    install: true
)
py.extension_module(
    '_byteswap',
    ['byteswap.pyx'],
    include_directories: [inc_np],
    dependencies: [py_dep],
    # The file is named byteswap.pyx but we want the
    # extension module to be named _byteswap
    cython_args: ['--module-name=pandas.io.sas._byteswap'],
    subdir: 'pandas/io/sas',
    install: true
)
top_level_py_list = [
    '__init__.py',
    'sas7bdat.py',
    'sas_constants.py',
    'sas_xport.py',
    'sasreader.py'
]
foreach file: top_level_py_list
    py.install_sources(file,
                       pure: false,
                       subdir: 'pandas/io/sas')
endforeach
