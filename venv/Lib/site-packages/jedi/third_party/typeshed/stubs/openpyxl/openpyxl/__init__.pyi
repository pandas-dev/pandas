from _typeshed import StrPath, SupportsRead, SupportsWrite
from typing import IO, Literal, Protocol, type_check_only
from typing_extensions import TypeAlias

from openpyxl.compat.numbers import NUMPY as NUMPY
from openpyxl.reader.excel import load_workbook as load_workbook
from openpyxl.workbook import Workbook as Workbook
from openpyxl.xml import DEFUSEDXML as DEFUSEDXML, LXML as LXML

from ._constants import (
    __author__ as __author__,
    __author_email__ as __author_email__,
    __license__ as __license__,
    __maintainer_email__ as __maintainer_email__,
    __url__ as __url__,
    __version__ as __version__,
)

DEBUG: bool
open = load_workbook

# Utility types reused elsewhere
_VisibilityType: TypeAlias = Literal["visible", "hidden", "veryHidden"]  # noqa: Y047

# TODO: Use a proper protocol from ZipFile. See: #10880
# This alias is to minimize false-positives
_ZipFileFileProtocol: TypeAlias = StrPath | IO[bytes] | SupportsRead[bytes]  # noqa: Y047
_ZipFileFileWriteProtocol: TypeAlias = StrPath | IO[bytes] | SupportsWrite[bytes]  # noqa: Y047

@type_check_only
class _Decodable(Protocol):  # noqa: Y046
    def decode(self, encoding: str, /) -> str: ...
