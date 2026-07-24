from _typeshed import ConvertibleToFloat, ConvertibleToInt, Incomplete, Unused
from datetime import datetime
from typing import ClassVar, Final, Literal, overload
from typing_extensions import TypeAlias

from openpyxl.descriptors.base import (
    Alias,
    Bool,
    DateTime,
    Float,
    Integer,
    MinMax,
    NoneSet,
    Set,
    String,
    Typed,
    _ConvertibleToBool,
)
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.serialisable import Serialisable

_SortConditionSortBy: TypeAlias = Literal["value", "cellColor", "fontColor", "icon"]
_IconSet: TypeAlias = Literal[
    "3Arrows",
    "3ArrowsGray",
    "3Flags",
    "3TrafficLights1",
    "3TrafficLights2",
    "3Signs",
    "3Symbols",
    "3Symbols2",
    "4Arrows",
    "4ArrowsGray",
    "4RedToBlack",
    "4Rating",
    "4TrafficLights",
    "5Arrows",
    "5ArrowsGray",
    "5Rating",
    "5Quarters",
]
_SortStateSortMethod: TypeAlias = Literal["stroke", "pinYin"]
_CustomFilterOperator: TypeAlias = Literal[
    "equal", "lessThan", "lessThanOrEqual", "notEqual", "greaterThanOrEqual", "greaterThan"
]
_StringFilterOperator: TypeAlias = Literal["contains", "startswith", "endswith", "wildcard"]
_FiltersCalendarType: TypeAlias = Literal[
    "gregorian",
    "gregorianUs",
    "gregorianMeFrench",
    "gregorianArabic",
    "hijri",
    "hebrew",
    "taiwan",
    "japan",
    "thai",
    "korea",
    "saka",
    "gregorianXlitEnglish",
    "gregorianXlitFrench",
]
_DynamicFilterType: TypeAlias = Literal[
    "null",
    "aboveAverage",
    "belowAverage",
    "tomorrow",
    "today",
    "yesterday",
    "nextWeek",
    "thisWeek",
    "lastWeek",
    "nextMonth",
    "thisMonth",
    "lastMonth",
    "nextQuarter",
    "thisQuarter",
    "lastQuarter",
    "nextYear",
    "thisYear",
    "lastYear",
    "yearToDate",
    "Q1",
    "Q2",
    "Q3",
    "Q4",
    "M1",
    "M2",
    "M3",
    "M4",
    "M5",
    "M6",
    "M7",
    "M8",
    "M9",
    "M10",
    "M11",
    "M12",
]
_DateGroupItemDateTimeGrouping: TypeAlias = Literal["year", "month", "day", "hour", "minute", "second"]

class SortCondition(Serialisable):
    tagname: ClassVar[str]
    descending: Bool[Literal[True]]
    sortBy: NoneSet[_SortConditionSortBy]
    ref: Incomplete
    customList: String[Literal[True]]
    dxfId: Integer[Literal[True]]
    iconSet: NoneSet[_IconSet]
    iconId: Integer[Literal[True]]
    def __init__(
        self,
        ref=None,
        descending: _ConvertibleToBool | None = None,
        sortBy: _SortConditionSortBy | Literal["none"] | None = None,
        customList: str | None = None,
        dxfId: ConvertibleToInt | None = None,
        iconSet: _IconSet | Literal["none"] | None = None,
        iconId: ConvertibleToInt | None = None,
    ) -> None: ...

class SortState(Serialisable):
    tagname: ClassVar[str]
    columnSort: Bool[Literal[True]]
    caseSensitive: Bool[Literal[True]]
    sortMethod: NoneSet[_SortStateSortMethod]
    ref: Incomplete
    sortCondition: Incomplete
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        columnSort: _ConvertibleToBool | None = None,
        caseSensitive: _ConvertibleToBool | None = None,
        sortMethod: _SortStateSortMethod | Literal["none"] | None = None,
        ref=None,
        sortCondition=(),
        extLst: Unused = None,
    ) -> None: ...
    def __bool__(self) -> bool: ...

class IconFilter(Serialisable):
    tagname: ClassVar[str]
    iconSet: Set[_IconSet]
    iconId: Integer[Literal[True]]
    def __init__(self, iconSet: _IconSet, iconId: ConvertibleToInt | None = None) -> None: ...

class ColorFilter(Serialisable):
    tagname: ClassVar[str]
    dxfId: Integer[Literal[True]]
    cellColor: Bool[Literal[True]]
    def __init__(self, dxfId: ConvertibleToInt | None = None, cellColor: _ConvertibleToBool | None = None) -> None: ...

class DynamicFilter(Serialisable):
    tagname: ClassVar[str]
    type: Set[_DynamicFilterType]
    val: Float[Literal[True]]
    valIso: DateTime[Literal[True]]
    maxVal: Float[Literal[True]]
    maxValIso: DateTime[Literal[True]]
    def __init__(
        self,
        type: _DynamicFilterType,
        val: ConvertibleToFloat | None = None,
        valIso: datetime | str | None = None,
        maxVal: ConvertibleToFloat | None = None,
        maxValIso: datetime | str | None = None,
    ) -> None: ...

class CustomFilter(Serialisable):
    tagname: ClassVar[str]
    val: String[Literal[False]]
    operator: Set[_CustomFilterOperator]
    def __init__(self, operator: _CustomFilterOperator = "equal", val: str | None = None) -> None: ...
    def convert(self) -> BlankFilter | NumberFilter | StringFilter: ...

class BlankFilter(CustomFilter):
    def __init__(self, **kw: Unused) -> None: ...
    @property
    def operator(self) -> Literal["notEqual"]: ...  # type: ignore[override]
    @property
    def val(self) -> Literal[" "]: ...  # type: ignore[override]

class NumberFilter(CustomFilter):
    val: Float[Literal[False]]  # type: ignore[assignment]
    def __init__(self, operator: _CustomFilterOperator = "equal", val: ConvertibleToFloat | None = None) -> None: ...

string_format_mapping: Final[dict[_StringFilterOperator, str]]

class StringFilter(CustomFilter):
    operator: Set[_StringFilterOperator]  # type: ignore[assignment]
    val: String[Literal[False]]
    exclude: Bool[Literal[False]]
    def __init__(
        self, operator: _StringFilterOperator = "contains", val: str | None = None, exclude: _ConvertibleToBool = False
    ) -> None: ...

class CustomFilters(Serialisable):
    tagname: ClassVar[str]
    _and: Bool[Literal[True]]  # Not private. Avoids name clash
    customFilter: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, _and: _ConvertibleToBool | None = None, customFilter=()) -> None: ...

class Top10(Serialisable):
    tagname: ClassVar[str]
    top: Bool[Literal[True]]
    percent: Bool[Literal[True]]
    val: Float[Literal[False]]
    filterVal: Float[Literal[True]]
    @overload
    def __init__(
        self,
        top: _ConvertibleToBool | None = None,
        percent: _ConvertibleToBool | None = None,
        *,
        val: ConvertibleToFloat,
        filterVal: ConvertibleToFloat | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        top: _ConvertibleToBool | None,
        percent: _ConvertibleToBool | None,
        val: ConvertibleToFloat,
        filterVal: ConvertibleToFloat | None = None,
    ) -> None: ...

class DateGroupItem(Serialisable):
    tagname: ClassVar[str]
    year: Integer[Literal[False]]
    month: MinMax[float, Literal[True]]
    day: MinMax[float, Literal[True]]
    hour: MinMax[float, Literal[True]]
    minute: MinMax[float, Literal[True]]
    second: Integer[Literal[True]]
    dateTimeGrouping: Set[_DateGroupItemDateTimeGrouping]
    @overload
    def __init__(
        self,
        year: ConvertibleToInt,
        month: ConvertibleToFloat | None = None,
        day: ConvertibleToFloat | None = None,
        hour: ConvertibleToFloat | None = None,
        minute: ConvertibleToFloat | None = None,
        second: ConvertibleToInt | None = None,
        *,
        dateTimeGrouping: _DateGroupItemDateTimeGrouping,
    ) -> None: ...
    @overload
    def __init__(
        self,
        year: ConvertibleToInt,
        month: ConvertibleToFloat | None,
        day: ConvertibleToFloat | None,
        hour: ConvertibleToFloat | None,
        minute: ConvertibleToFloat | None,
        second: ConvertibleToInt | None,
        dateTimeGrouping: _DateGroupItemDateTimeGrouping,
    ) -> None: ...

class Filters(Serialisable):
    tagname: ClassVar[str]
    blank: Bool[Literal[True]]
    calendarType: NoneSet[_FiltersCalendarType]
    filter: Incomplete
    dateGroupItem: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        blank: _ConvertibleToBool | None = None,
        calendarType: _FiltersCalendarType | Literal["none"] | None = None,
        filter=(),
        dateGroupItem=(),
    ) -> None: ...

class FilterColumn(Serialisable):
    tagname: ClassVar[str]
    colId: Integer[Literal[False]]
    col_id: Alias
    hiddenButton: Bool[Literal[True]]
    showButton: Bool[Literal[True]]
    filters: Typed[Filters, Literal[True]]
    top10: Typed[Top10, Literal[True]]
    customFilters: Typed[CustomFilters, Literal[True]]
    dynamicFilter: Typed[DynamicFilter, Literal[True]]
    colorFilter: Typed[ColorFilter, Literal[True]]
    iconFilter: Typed[IconFilter, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        colId: ConvertibleToInt,
        hiddenButton: _ConvertibleToBool | None = False,
        showButton: _ConvertibleToBool | None = True,
        filters: Filters | None = None,
        top10: Top10 | None = None,
        customFilters: CustomFilters | None = None,
        dynamicFilter: DynamicFilter | None = None,
        colorFilter: ColorFilter | None = None,
        iconFilter: IconFilter | None = None,
        extLst: Unused = None,
        blank=None,
        vals=None,
    ) -> None: ...

class AutoFilter(Serialisable):
    tagname: ClassVar[str]
    ref: Incomplete
    filterColumn: Incomplete
    sortState: Typed[SortState, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, ref=None, filterColumn=(), sortState: SortState | None = None, extLst: Unused = None) -> None: ...
    def __bool__(self) -> bool: ...
    def add_filter_column(self, col_id, vals, blank: bool = False) -> None: ...
    def add_sort_condition(self, ref, descending: bool = False) -> None: ...
