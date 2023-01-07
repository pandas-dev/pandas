from __future__ import annotations

from io import (
    BufferedReader,
    BytesIO,
)
from pathlib import PurePath
from tempfile import NamedTemporaryFile
from typing import Any

from pandas.compat._optional import import_optional_dependency

from pandas.io.excel._base import (
    BaseExcelReader,
    inspect_excel_format,
)


class __calamine__:
    pass


class CalamineExcelReader(BaseExcelReader):
    book: str
    _sheet_names: list[str] | None = None
    import_optional_dependency("python_calamine")

    @property
    def _workbook_class(self) -> type[__calamine__]:
        return __calamine__

    def load_workbook(
        self, filepath_or_buffer: str | PurePath | BufferedReader | BytesIO
    ) -> str:
        if isinstance(filepath_or_buffer, BufferedReader):
            filepath_or_buffer = filepath_or_buffer.name

        elif isinstance(filepath_or_buffer, BytesIO):
            ext = inspect_excel_format(filepath_or_buffer)
            with NamedTemporaryFile(suffix=f".{ext}", delete=False) as tmp_file:
                tmp_file.write(filepath_or_buffer.getvalue())
                filepath_or_buffer = tmp_file.name

        elif isinstance(filepath_or_buffer, PurePath):
            filepath_or_buffer = filepath_or_buffer.as_posix()

        assert isinstance(filepath_or_buffer, str)

        from python_calamine import get_sheet_names

        self._sheet_names = get_sheet_names(filepath_or_buffer)
        return filepath_or_buffer

    @property
    def sheet_names(self) -> list[str]:
        from python_calamine import get_sheet_names

        if self._sheet_names is None:
            self._sheet_names = get_sheet_names(self.book)
        return self._sheet_names

    def get_sheet_by_name(self, name: str) -> int:
        self.raise_if_bad_sheet_by_name(name)
        return self.sheet_names.index(name)

    def get_sheet_data(self, sheet: int, convert_float: bool) -> list[list[Any]]:
        from python_calamine import get_sheet_data

        return get_sheet_data(self.book, sheet)
