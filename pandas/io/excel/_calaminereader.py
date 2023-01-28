from __future__ import annotations

from datetime import (
    date,
    datetime,
    time,
)
from tempfile import NamedTemporaryFile
from typing import Union

from pandas._typing import (
    FilePath,
    ReadBuffer,
    Scalar,
    StorageOptions,
)
from pandas.compat._optional import import_optional_dependency

import pandas as pd

from pandas.io.common import stringify_path
from pandas.io.excel._base import (
    BaseExcelReader,
    inspect_excel_format,
)

ValueT = Union[int, float, str, bool, time, date, datetime]


class __calamine__:
    pass


class CalamineExcelReader(BaseExcelReader):
    book: str
    _sheet_names: list[str] | None = None

    def __init__(
        self,
        filepath_or_buffer: FilePath | ReadBuffer[bytes],
        storage_options: StorageOptions = None,
    ) -> None:
        import_optional_dependency("python_calamine")
        super().__init__(filepath_or_buffer, storage_options=storage_options)

    @property
    def _workbook_class(self) -> type[__calamine__]:
        return __calamine__

    def load_workbook(self, filepath_or_buffer) -> str:
        if hasattr(filepath_or_buffer, "read") and hasattr(filepath_or_buffer, "seek"):
            ext = inspect_excel_format(filepath_or_buffer)
            with NamedTemporaryFile(suffix=f".{ext}", delete=False) as tmp_file:
                filepath_or_buffer.seek(0)
                tmp_file.write(filepath_or_buffer.read())
                filepath_or_buffer = tmp_file.name
        else:
            filepath_or_buffer = stringify_path(filepath_or_buffer)

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

    def get_sheet_by_index(self, index: int) -> int:
        self.raise_if_bad_sheet_by_index(index)
        return index

    def get_sheet_data(
        self, sheet: int, file_rows_needed: int | None = None
    ) -> list[list[Scalar]]:
        def _convert_cell(value: ValueT) -> Scalar:
            if isinstance(value, float):
                val = int(value)
                if val == value:
                    return val
                else:
                    return value
            elif isinstance(value, date):
                return pd.Timestamp(value)
            elif isinstance(value, time):
                return value.isoformat()

            return value

        from python_calamine import get_sheet_data

        rows = get_sheet_data(self.book, sheet, skip_empty_area=False)
        data: list[list[Scalar]] = []

        for row in rows:
            data.append([_convert_cell(cell) for cell in row])
            if file_rows_needed is not None and len(data) >= file_rows_needed:
                break

        return data
