"""
Data IO api
"""

from pandas.io.parsers import read_csv, read_table, read_fwf
from pandas.io.clipboard import read_clipboard
from pandas.io.excel import ExcelFile, ExcelWriter, read_excel
from pandas.io.pytables import HDFStore, Term, get_store, read_hdf
from pandas.io.html import read_html
from pandas.io.sql import read_sql
from pandas.io.stata import read_stata

__all__ = ['read_csv', 'read_table', 'read_clipboard', 'read_fwf',
        'to_clipboard', 'ExcelFile', 'ExcelWriter', 'read_excel', 'HDFStore', 'Term',
        'get_store', 'read_hdf', 'read_html', 'read_sql', 'read_stata']
