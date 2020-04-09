from datetime import time
from typing import List, Optional, Sequence

import numpy as np

from pandas._typing import Scalar, Union
from pandas.compat._optional import import_optional_dependency

from pandas.io.excel._base import _BaseExcelReader
from pandas.io.parsers import _validate_integer


class _XlrdReader(_BaseExcelReader):
    def __init__(self, filepath_or_buffer):
        """
        Reader using xlrd engine.

        Parameters
        ----------
        filepath_or_buffer : string, path object or Workbook
            Object to be parsed.
        """
        err_msg = "Install xlrd >= 1.0.0 for Excel support"
        import_optional_dependency("xlrd", extra=err_msg)
        super().__init__(filepath_or_buffer)

    @property
    def _workbook_class(self):
        from xlrd import Book

        return Book

    def load_workbook(self, filepath_or_buffer):
        from xlrd import open_workbook

        if hasattr(filepath_or_buffer, "read"):
            data = filepath_or_buffer.read()
            return open_workbook(file_contents=data)
        else:
            return open_workbook(filepath_or_buffer)

    @property
    def sheet_names(self):
        return self.book.sheet_names()

    def get_sheet_by_name(self, name):
        return self.book.sheet_by_name(name)

    def get_sheet_by_index(self, index):
        return self.book.sheet_by_index(index)

    def get_sheet_data(
        self,
        sheet,
        convert_float,
        header: Optional[Union[int, Sequence[int]]],
        skiprows: Optional[Union[int, Sequence[int]]],
        nrows: Optional[int],
    ) -> List[List[Scalar]]:
        from xlrd import (
            xldate,
            XL_CELL_DATE,
            XL_CELL_ERROR,
            XL_CELL_BOOLEAN,
            XL_CELL_NUMBER,
        )

        epoch1904 = self.book.datemode

        def _parse_cell(cell_contents, cell_typ):
            """
            converts the contents of the cell into a pandas appropriate object
            """
            if cell_typ == XL_CELL_DATE:

                # Use the newer xlrd datetime handling.
                try:
                    cell_contents = xldate.xldate_as_datetime(cell_contents, epoch1904)
                except OverflowError:
                    return cell_contents

                # Excel doesn't distinguish between dates and time,
                # so we treat dates on the epoch as times only.
                # Also, Excel supports 1900 and 1904 epochs.
                year = (cell_contents.timetuple())[0:3]
                if (not epoch1904 and year == (1899, 12, 31)) or (
                    epoch1904 and year == (1904, 1, 1)
                ):
                    cell_contents = time(
                        cell_contents.hour,
                        cell_contents.minute,
                        cell_contents.second,
                        cell_contents.microsecond,
                    )

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

        data: List[List[Scalar]] = []

        if nrows is not None:
            _validate_integer("nrows", nrows)
        header = 0 if header is None else header
        skiprows = 0 if skiprows is None else skiprows
        if isinstance(header, list) or isinstance(skiprows, list):
            nrows = None
        for i in range(sheet.nrows):

            if nrows is not None:
                if header > 1:
                    header -= 1
                    data.append([])
                    continue
                elif skiprows > 0:
                    skiprows -= 1
                    data.append([])
                    continue
                if nrows >= 0:
                    nrows -= 1
                else:
                    break

            row = [
                _parse_cell(value, typ)
                for value, typ in zip(sheet.row_values(i), sheet.row_types(i))
            ]
            data.append(row)

        return data
