from datetime import date, datetime, time, timedelta
from decimal import Decimal
from typing import Any
from typing_extensions import TypeAlias

from openpyxl.cell.rich_text import CellRichText
from openpyxl.worksheet.formula import ArrayFormula, DataTableFormula

from .cell import Cell as Cell, MergedCell as MergedCell, WriteOnlyCell as WriteOnlyCell
from .read_only import ReadOnlyCell as ReadOnlyCell

_TimeTypes: TypeAlias = datetime | date | time | timedelta
_CellGetValue: TypeAlias = (  # noqa: Y047 # Used in other modules
    # if numpy is installed also numpy bool and number types
    bool
    | float
    | Decimal
    | str
    | CellRichText
    | _TimeTypes
    | DataTableFormula
    | ArrayFormula
    | None
)
_AnyCellValue: TypeAlias = Any  # AnyOf _CellGetValue # noqa: Y047 # Used in other modules
_CellSetValue: TypeAlias = _CellGetValue | bytes  # noqa: Y047 # Used in other modules

_CellOrMergedCell: TypeAlias = Cell | MergedCell  # noqa: Y047 # Used in other modules
