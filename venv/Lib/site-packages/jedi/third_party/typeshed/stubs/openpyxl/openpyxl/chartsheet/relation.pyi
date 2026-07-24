from _typeshed import ConvertibleToInt, Incomplete
from typing import ClassVar, Literal

from openpyxl.descriptors.base import Alias, Integer
from openpyxl.descriptors.serialisable import Serialisable

class SheetBackgroundPicture(Serialisable):
    tagname: ClassVar[str]
    id: Incomplete
    def __init__(self, id) -> None: ...

class DrawingHF(Serialisable):
    id: Incomplete
    lho: Integer[Literal[True]]
    leftHeaderOddPages: Alias
    lhe: Integer[Literal[True]]
    leftHeaderEvenPages: Alias
    lhf: Integer[Literal[True]]
    leftHeaderFirstPage: Alias
    cho: Integer[Literal[True]]
    centerHeaderOddPages: Alias
    che: Integer[Literal[True]]
    centerHeaderEvenPages: Alias
    chf: Integer[Literal[True]]
    centerHeaderFirstPage: Alias
    rho: Integer[Literal[True]]
    rightHeaderOddPages: Alias
    rhe: Integer[Literal[True]]
    rightHeaderEvenPages: Alias
    rhf: Integer[Literal[True]]
    rightHeaderFirstPage: Alias
    lfo: Integer[Literal[True]]
    leftFooterOddPages: Alias
    lfe: Integer[Literal[True]]
    leftFooterEvenPages: Alias
    lff: Integer[Literal[True]]
    leftFooterFirstPage: Alias
    cfo: Integer[Literal[True]]
    centerFooterOddPages: Alias
    cfe: Integer[Literal[True]]
    centerFooterEvenPages: Alias
    cff: Integer[Literal[True]]
    centerFooterFirstPage: Alias
    rfo: Integer[Literal[True]]
    rightFooterOddPages: Alias
    rfe: Integer[Literal[True]]
    rightFooterEvenPages: Alias
    rff: Integer[Literal[True]]
    rightFooterFirstPage: Alias
    def __init__(
        self,
        id=None,
        lho: ConvertibleToInt | None = None,
        lhe: ConvertibleToInt | None = None,
        lhf: ConvertibleToInt | None = None,
        cho: ConvertibleToInt | None = None,
        che: ConvertibleToInt | None = None,
        chf: ConvertibleToInt | None = None,
        rho: ConvertibleToInt | None = None,
        rhe: ConvertibleToInt | None = None,
        rhf: ConvertibleToInt | None = None,
        lfo: ConvertibleToInt | None = None,
        lfe: ConvertibleToInt | None = None,
        lff: ConvertibleToInt | None = None,
        cfo: ConvertibleToInt | None = None,
        cfe: ConvertibleToInt | None = None,
        cff: ConvertibleToInt | None = None,
        rfo: ConvertibleToInt | None = None,
        rfe: ConvertibleToInt | None = None,
        rff: ConvertibleToInt | None = None,
    ) -> None: ...
