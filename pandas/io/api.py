"""
Data IO api
"""

# flake8: noqa

from .clipboards import read_clipboard
from .excel import (
    ExcelFile,
    ExcelWriter,
    read_excel,
)
from .feather_format import read_feather
from .gbq import read_gbq
from .html import read_html
from .json import read_json
from .orc import read_orc
from .parquet import read_parquet
from .parsers import (
    read_csv,
    read_fwf,
    read_table,
)
from .pickle import (
    read_pickle,
    to_pickle,
)
from .pytables import (
    HDFStore,
    read_hdf,
)
from .sas import read_sas
from .spss import read_spss
from .sql import (
    read_sql,
    read_sql_query,
    read_sql_table,
)
from .stata import read_stata
