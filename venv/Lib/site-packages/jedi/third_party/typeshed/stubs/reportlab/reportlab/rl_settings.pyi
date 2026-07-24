from _typeshed import Incomplete
from typing import Final, Literal

__version__: Final[str]

# NOTE: All the attributes in this module are Final
#       rl_config is the place to make changes
allowTableBoundsErrors: Final[int]
shapeChecking: Final[int]
defaultEncoding: Final[str]
defaultGraphicsFontName: Final[str]
pageCompression: Final[int]
useA85: Final[int]
defaultPageSize: Final[str]
defaultImageCaching: Final[int]
warnOnMissingFontGlyphs: Final[int]
verbose: Final[int]
showBoundary: Final[int]
emptyTableAction: Final[str]
invariant: Final[int]
eps_preview_transparent: Final[Incomplete]
eps_preview: Final[int]
eps_ttf_embed: Final[int]
eps_ttf_embed_uid: Final[int]
overlapAttachedSpace: Final[int]
longTableOptimize: Final[int]
autoConvertEncoding: Final[int]
_FUZZ: Final[float]
wrapA85: Final[int]
fsEncodings: Final[tuple[Literal["utf8"], Literal["cp1252"], Literal["cp430"]]]
odbc_driver: Final[str]
platypus_link_underline: Final[int]
canvas_basefontname: Final[str]
allowShortTableRows: Final[int]
imageReaderFlags: Final[int]
paraFontSizeHeightOffset: Final[int]
canvas_baseColor: Final[Incomplete]
ignoreContainerActions: Final[int]
ttfAsciiReadable: Final[int]
pdfMultiLine: Final[int]
pdfComments: Final[int]
debug: Final[int]
listWrapOnFakeWidth: Final[int]
underlineWidth: Final[str]
underlineOffset: Final[str]
underlineGap: Final[str]
strikeWidth: Final[str]
strikeOffset: Final[str]
strikeGap: Final[str]
decimalSymbol: Final[str]
errorOnDuplicatePageLabelPage: Final[int]
autoGenerateMissingTTFName: Final[int]
allowTTFSubsetting: Final[list[str]]
spaceShrinkage: Final[float]
hyphenationLang: Final[str]
uriWasteReduce: Final[int]
embeddedHyphenation: Final[int]
hyphenationMinWordLength: Final[int]
reserveTTFNotdef: Final[int]
documentLang: Final[Incomplete]
encryptionStrength: Final[int]
trustedHosts: Final[Incomplete]
trustedSchemes: Final[list[str]]
renderPMBackend: Final[str]
xmlParser: Final[str]
textPaths: Final[str]
toColorCanUse: Final[str]
defCWRF: Final[float]
unShapedFontGlob: list[str] | None
T1SearchPath: Final[tuple[str, ...]]
TTFSearchPath: Final[tuple[str, ...]]
CMapSearchPath: Final[tuple[str, ...]]

__all__ = (
    "allowTableBoundsErrors",
    "shapeChecking",
    "defaultEncoding",
    "defaultGraphicsFontName",
    "pageCompression",
    "useA85",
    "defaultPageSize",
    "defaultImageCaching",
    "warnOnMissingFontGlyphs",
    "verbose",
    "showBoundary",
    "emptyTableAction",
    "invariant",
    "eps_preview_transparent",
    "eps_preview",
    "eps_ttf_embed",
    "eps_ttf_embed_uid",
    "overlapAttachedSpace",
    "longTableOptimize",
    "autoConvertEncoding",
    "_FUZZ",
    "wrapA85",
    "fsEncodings",
    "odbc_driver",
    "platypus_link_underline",
    "canvas_basefontname",
    "allowShortTableRows",
    "imageReaderFlags",
    "paraFontSizeHeightOffset",
    "canvas_baseColor",
    "ignoreContainerActions",
    "ttfAsciiReadable",
    "pdfMultiLine",
    "pdfComments",
    "debug",
    "listWrapOnFakeWidth",
    "T1SearchPath",
    "TTFSearchPath",
    "CMapSearchPath",
    "decimalSymbol",
    "errorOnDuplicatePageLabelPage",
    "autoGenerateMissingTTFName",
    "allowTTFSubsetting",
    "spaceShrinkage",
    "underlineWidth",
    "underlineOffset",
    "underlineGap",
    "strikeWidth",
    "strikeOffset",
    "strikeGap",
    "hyphenationLang",
    "uriWasteReduce",
    "embeddedHyphenation",
    "hyphenationMinWordLength",
    "reserveTTFNotdef",
    "documentLang",
    "encryptionStrength",
    "trustedHosts",
    "trustedSchemes",
    "renderPMBackend",
    "xmlParser",
    "textPaths",
    "toColorCanUse",
    "defCWRF",
    "unShapedFontGlob",
)
