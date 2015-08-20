"""
Data IO api
"""

from pandas.io.parsers import read_csv, read_table, read_fwf
from pandas.io.clipboard import read_clipboard
from pandas.io.excel import ExcelFile, ExcelWriter, read_excel
from pandas.io.pytables import HDFStore, Term, get_store, read_hdf
from pandas.io.json import read_json
from pandas.io.html import read_html
from pandas.io.sql import read_sql, read_sql_table, read_sql_query
from pandas.io.sas import read_sas
from pandas.io.stata import read_stata
from pandas.io.pickle import read_pickle, to_pickle
from pandas.io.packers import read_msgpack, to_msgpack
from pandas.io.gbq import read_gbq
