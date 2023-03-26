from __future__ import annotations

from datetime import (
    date,
    datetime,
    time,
)
from typing import (
    TYPE_CHECKING,
    Union,
)

from pandas.compat._optional import import_optional_dependency
from pandas.util._decorators import doc

import pandas as pd
from pandas.core.shared_docs import _shared_docs

from pandas.io.excel._base import BaseExcelReader

if TYPE_CHECKING:
    from pandas._typing import (
        FilePath,
        ReadBuffer,
        Scalar,
        StorageOptions,
    )

_ValueT = Union[int, float, str, bool, time, date, datetime]


class CalamineExcelReader(BaseExcelReader):
    @doc(storage_options=_shared_docs["storage_options"])
    def __init__(
        self,
        filepath_or_buffer: FilePath | ReadBuffer[bytes],
        storage_options: StorageOptions = None,
    ) -> None:
        """
        Reader using calamine engine (xlsx/xls/xlsb/ods).

        Parameters
        ----------
        filepath_or_buffer : str, path to be parsed or
            an open readable stream.
        {storage_options}
        """
        import_optional_dependency("python_calamine")
        super().__init__(filepath_or_buffer, storage_options=storage_options)

    @property
    def _workbook_class(self):
        from python_calamine import CalamineWorkbook

        return CalamineWorkbook

    def load_workbook(self, filepath_or_buffer: FilePath | ReadBuffer[bytes]):
        from python_calamine import load_workbook

        return load_workbook(filepath_or_buffer)  # type: ignore[arg-type]

    @property
    def sheet_names(self) -> list[str]:
        return self.book.sheet_names  # pyright: ignore

    def get_sheet_by_name(self, name: str):
        self.raise_if_bad_sheet_by_name(name)
        return self.book.get_sheet_by_name(name)  # pyright: ignore

    def get_sheet_by_index(self, index: int):
        self.raise_if_bad_sheet_by_index(index)
        return self.book.get_sheet_by_index(index)  # pyright: ignore

    def get_sheet_data(
        self, sheet, file_rows_needed: int | None = None
    ) -> list[list[Scalar]]:
        def _convert_cell(value: _ValueT) -> Scalar:
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

        rows: list[list[_ValueT]] = sheet.to_python(skip_empty_area=False)
        data: list[list[Scalar]] = []

        for row in rows:
            data.append([_convert_cell(cell) for cell in row])
            if file_rows_needed is not None and len(data) >= file_rows_needed:
                break

        return data
