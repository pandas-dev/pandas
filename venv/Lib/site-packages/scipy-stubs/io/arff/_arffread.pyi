import abc
import re
from collections.abc import Iterable, Iterator, Sequence
from csv import Dialect
from typing import Any, ClassVar, Final, Generic, Literal, Self
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp

from scipy.io._typing import FileLike

__all__ = ["ArffError", "MetaData", "ParseArffError", "loadarff"]

_T_co = TypeVar("_T_co", covariant=True, default=object)

###

r_meta: Final[re.Pattern[str]] = ...  # undocumented
r_comment: Final[re.Pattern[str]] = ...  # undocumented
r_empty: Final[re.Pattern[str]] = ...  # undocumented
r_headerline: Final[re.Pattern[str]] = ...  # undocumented
r_datameta: Final[re.Pattern[str]] = ...  # undocumented
r_relation: Final[re.Pattern[str]] = ...  # undocumented
r_attribute: Final[re.Pattern[str]] = ...  # undocumented
r_nominal: Final[re.Pattern[str]] = ...  # undocumented
r_date: Final[re.Pattern[str]] = ...  # undocumented
r_comattrval: Final[re.Pattern[str]] = ...  # undocumented
r_wcomattrval: Final[re.Pattern[str]] = ...  # undocumented

class ArffError(OSError): ...
class ParseArffError(ArffError): ...

class Attribute(Generic[_T_co], metaclass=abc.ABCMeta):  # undocumented
    type_name: ClassVar[str | None]
    dtype: Any
    range: Any

    name: Final[str]

    def __init__(self, /, name: str) -> None: ...
    @classmethod
    def parse_attribute(cls, name: str, attr_string: str) -> Self | None: ...
    def parse_data(self, /, data_str: str) -> _T_co: ...

class NominalAttribute(Attribute[str]):  # undocumented
    type_name: ClassVar[str | None] = "nominal"
    dtype: tuple[type[np.bytes_], ...]
    range: Sequence[str]

    values: Final[Sequence[str]]

    def __init__(self, /, name: str, values: Sequence[str]) -> None: ...

class NumericAttribute(Attribute[float]):  # undocumented
    type_name: ClassVar[str | None] = "numeric"
    dtype: type[np.float64]
    range: None

class StringAttribute(Attribute[None]):  # undocumented
    type_name: ClassVar[str | None] = "string"
    dtype: type[np.object_]
    range: None

class DateAttribute(Attribute[np.datetime64]):  # undocumented
    type_name: ClassVar[str | None] = "date"
    dtype: np.datetime64
    range: str

    date_format: Final[str]
    datetime_unit: Final[str]

    def __init__(self, /, name: str, date_format: str, datetime_unit: str) -> None: ...

class RelationalAttribute(Attribute[onp.Array1D[np.void]]):  # undocumented
    type_name: ClassVar[str | None] = "relational"
    dtype: type[np.object_]
    range: None

    attributes: Final[list[Attribute]]
    dialect: Dialect | None

class MetaData:
    name: Final[str]
    def __init__(self, /, rel: str, attr: Iterable[Attribute]) -> None: ...
    def __iter__(self, /) -> Iterator[str]: ...
    def __getitem__(self, /, key: str) -> tuple[str, str | Sequence[str] | None]: ...
    def names(self, /) -> list[str]: ...
    def types(self, /) -> list[str]: ...

def to_attribute(name: str, attr_string: str) -> Attribute: ...  # undocumented
def csv_sniffer_has_bug_last_field() -> Literal[False]: ...  # undocumented
def workaround_csv_sniffer_bug_last_field(
    sniff_line: str, dialect: Dialect, delimiters: Iterable[str]
) -> None: ...  # undocumented
def split_data_line(line: str, dialect: Dialect | None = None) -> tuple[list[str], Dialect]: ...  # undocumented
def tokenize_attribute(iterable: Iterable[int | str], attribute: str) -> tuple[Attribute, object]: ...  # undocumented
def tokenize_single_comma(val: str) -> tuple[str, str]: ...  # undocumented
def tokenize_single_wcomma(val: str) -> tuple[str, str]: ...  # undocumented
def read_relational_attribute(ofile: Iterator[str], relational_attribute: RelationalAttribute, i: str) -> str: ...  # undocumented
def read_header(ofile: Iterator[str]) -> tuple[str, list[Attribute]]: ...  # undocumented

#
def loadarff(f: FileLike[str]) -> tuple[onp.Array1D[np.void], MetaData]: ...
