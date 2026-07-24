from _typeshed import Incomplete
from collections.abc import Callable
from typing import Literal

def register_reset(func: Callable[[], Callable[[], object] | None]) -> None: ...
def _reset() -> None: ...

allowTableBoundsErrors: int
shapeChecking: int
defaultEncoding: str
defaultGraphicsFontName: str
pageCompression: int
useA85: int
defaultPageSize: tuple[float, float]
defaultImageCaching: int
warnOnMissingFontGlyphs: int
verbose: int
showBoundary: int
emptyTableAction: str
invariant: int
eps_preview_transparent: Incomplete
eps_preview: int
eps_ttf_embed: int
eps_ttf_embed_uid: int
overlapAttachedSpace: int
longTableOptimize: int
autoConvertEncoding: int
_FUZZ: float
wrapA85: int
fsEncodings: tuple[Literal["utf8"], Literal["cp1252"], Literal["cp430"]]
odbc_driver: str
platypus_link_underline: int
canvas_basefontname: str
allowShortTableRows: int
imageReaderFlags: int
paraFontSizeHeightOffset: int
canvas_baseColor: Incomplete
ignoreContainerActions: int
ttfAsciiReadable: int
pdfMultiLine: int
pdfComments: int
debug: int
listWrapOnFakeWidth: int
underlineWidth: str
underlineOffset: str
underlineGap: str
strikeWidth: str
strikeOffset: str
strikeGap: str
decimalSymbol: str
errorOnDuplicatePageLabelPage: int
autoGenerateMissingTTFName: int
allowTTFSubsetting: list[str]
spaceShrinkage: float
hyphenationLang: str
uriWasteReduce: int
embeddedHyphenation: int
hyphenationMinWordLength: int
reserveTTFNotdef: int
documentLang: Incomplete
encryptionStrength: int
trustedHosts: list[str] | None
trustedSchemes: list[str]
renderPMBackend: str
xmlParser: str
textPaths: str
toColorCanUse: str
defCWRF: float
unShapedFontGlob: list[str] | None
T1SearchPath: list[str]
TTFSearchPath: list[str]
CMapSearchPath: list[str]
