from typing import List

from pandas._typing import FilePathOrBuffer, Scalar, StorageOptions
from pandas.compat._optional import import_optional_dependency

from pandas.io.excel._base import BaseExcelReader


class PyxlsbReader(BaseExcelReader):
    def __init__(
        self,
        filepath_or_buffer: FilePathOrBuffer,
        storage_options: StorageOptions = None,
    ):
        """
        Reader using pyxlsb engine.

        Parameters
        ----------
        filepath_or_buffer : str, path object, or Workbook
            Object to be parsed.
        storage_options : dict, optional
            passed to fsspec for appropriate URLs (see ``_get_filepath_or_buffer``)
        """
        import_optional_dependency("pyxlsb")
        # This will call load_workbook on the filepath or buffer
        # And set the result to the book-attribute
        super().__init__(filepath_or_buffer, storage_options=storage_options)

    @property
    def _workbook_class(self):
        from pyxlsb import Workbook

        return Workbook

    def load_workbook(self, filepath_or_buffer: FilePathOrBuffer):
        from pyxlsb import open_workbook

        # TODO: hack in buffer capability
        # This might need some modifications to the Pyxlsb library
        # Actual work for opening it is in xlsbpackage.py, line 20-ish

        return open_workbook(filepath_or_buffer)

    @property
    def sheet_names(self) -> List[str]:
        return self.book.sheets

    def get_sheet_by_name(self, name: str):
        return self.book.get_sheet(name)

    def get_sheet_by_index(self, index: int):
        # pyxlsb sheets are indexed from 1 onwards
        # There's a fix for this in the source, but the pypi package doesn't have it
        return self.book.get_sheet(index + 1)

    def _convert_cell(self, cell, convert_float: bool) -> Scalar:
        # TODO: there is no way to distinguish between floats and datetimes in pyxlsb
        # This means that there is no way to read datetime types from an xlsb file yet
        if cell.v is None:
            return ""  # Prevents non-named columns from not showing up as Unnamed: i
        if isinstance(cell.v, float) and convert_float:
            val = int(cell.v)
            if val == cell.v:
                return val
            else:
                return float(cell.v)

        return cell.v

    def get_sheet_data(self, sheet, convert_float: bool) -> List[List[Scalar]]:
        return [
            [self._convert_cell(c, convert_float) for c in r]
            for r in sheet.rows(sparse=False)
        ]
