from _typeshed import ConvertibleToFloat
from typing import ClassVar, Literal

from openpyxl.descriptors.base import Alias, Float, Typed
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.worksheet.header_footer import HeaderFooter
from openpyxl.worksheet.page import PrintPageSetup

class PageMargins(Serialisable):
    tagname: ClassVar[str]
    l: Float[Literal[False]]
    left: Alias
    r: Float[Literal[False]]
    right: Alias
    t: Float[Literal[False]]
    top: Alias
    b: Float[Literal[False]]
    bottom: Alias
    header: Float[Literal[False]]
    footer: Float[Literal[False]]
    def __init__(
        self,
        l: ConvertibleToFloat = 0.75,
        r: ConvertibleToFloat = 0.75,
        t: ConvertibleToFloat = 1,
        b: ConvertibleToFloat = 1,
        header: ConvertibleToFloat = 0.5,
        footer: ConvertibleToFloat = 0.5,
    ) -> None: ...

class PrintSettings(Serialisable):
    tagname: ClassVar[str]
    headerFooter: Typed[HeaderFooter, Literal[True]]
    pageMargins: Typed[PageMargins, Literal[True]]
    pageSetup: Typed[PrintPageSetup, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        headerFooter: HeaderFooter | None = None,
        pageMargins: PageMargins | None = None,
        pageSetup: PrintPageSetup | None = None,
    ) -> None: ...
