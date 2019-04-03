from datetime import time
from distutils.version import LooseVersion
from io import UnsupportedOperation
from urllib.request import urlopen

import numpy as np

from pandas.io.common import _is_url, get_filepath_or_buffer
from pandas.io.excel._base import _BaseExcelReader


class _XlrdReader(_BaseExcelReader):

    def __init__(self, filepath_or_buffer):
        """Reader using xlrd engine.

        Parameters
        ----------
        filepath_or_buffer : string, path object or Workbook
            Object to be parsed.
        """
        err_msg = "Install xlrd >= 1.0.0 for Excel support"

        try:
            import xlrd
        except ImportError:
            raise ImportError(err_msg)
        else:
            if xlrd.__VERSION__ < LooseVersion("1.0.0"):
                raise ImportError(err_msg +
                                  ". Current version " + xlrd.__VERSION__)

        from pandas.io.excel._base import ExcelFile
        # If filepath_or_buffer is a url, want to keep the data as bytes so
        # can't pass to get_filepath_or_buffer()
        if _is_url(filepath_or_buffer):
            filepath_or_buffer = urlopen(filepath_or_buffer)
        elif not isinstance(filepath_or_buffer, (ExcelFile, xlrd.Book)):
            filepath_or_buffer, _, _, _ = get_filepath_or_buffer(
                filepath_or_buffer)

        if isinstance(filepath_or_buffer, xlrd.Book):
            self.book = filepath_or_buffer
        elif hasattr(filepath_or_buffer, "read"):
            # N.B. xlrd.Book has a read attribute too
            if hasattr(filepath_or_buffer, 'seek'):
                try:
                    # GH 19779
                    filepath_or_buffer.seek(0)
                except UnsupportedOperation:
                    # HTTPResponse does not support seek()
                    # GH 20434
                    pass

            data = filepath_or_buffer.read()
            self.book = xlrd.open_workbook(file_contents=data)
        elif isinstance(filepath_or_buffer, str):
            self.book = xlrd.open_workbook(filepath_or_buffer)
        else:
            raise ValueError('Must explicitly set engine if not passing in'
                             ' buffer or path for io.')

    @property
    def sheet_names(self):
        return self.book.sheet_names()

    def get_sheet_by_name(self, name):
        return self.book.sheet_by_name(name)

    def get_sheet_by_index(self, index):
        return self.book.sheet_by_index(index)

    def get_sheet_data(self, sheet, convert_float):
        from xlrd import (xldate, XL_CELL_DATE,
                          XL_CELL_ERROR, XL_CELL_BOOLEAN,
                          XL_CELL_NUMBER)

        epoch1904 = self.book.datemode

        def _parse_cell(cell_contents, cell_typ):
            """converts the contents of the cell into a pandas
               appropriate object"""

            if cell_typ == XL_CELL_DATE:

                # Use the newer xlrd datetime handling.
                try:
                    cell_contents = xldate.xldate_as_datetime(
                        cell_contents, epoch1904)
                except OverflowError:
                    return cell_contents

                # Excel doesn't distinguish between dates and time,
                # so we treat dates on the epoch as times only.
                # Also, Excel supports 1900 and 1904 epochs.
                year = (cell_contents.timetuple())[0:3]
                if ((not epoch1904 and year == (1899, 12, 31)) or
                        (epoch1904 and year == (1904, 1, 1))):
                    cell_contents = time(cell_contents.hour,
                                         cell_contents.minute,
                                         cell_contents.second,
                                         cell_contents.microsecond)

            elif cell_typ == XL_CELL_ERROR:
                cell_contents = np.nan
            elif cell_typ == XL_CELL_BOOLEAN:
                cell_contents = bool(cell_contents)
            elif convert_float and cell_typ == XL_CELL_NUMBER:
                # GH5394 - Excel 'numbers' are always floats
                # it's a minimal perf hit and less surprising
                val = int(cell_contents)
                if val == cell_contents:
                    cell_contents = val
            return cell_contents

        data = []

        for i in range(sheet.nrows):
            row = [_parse_cell(value, typ)
                   for value, typ in zip(sheet.row_values(i),
                                         sheet.row_types(i))]
            data.append(row)

        return data
