from typing import Final

MIN_ROW: Final = 0
MIN_COLUMN: Final = 0
MAX_COLUMN: Final = 16384
MAX_ROW: Final = 1048576

PACKAGE_PROPS: Final = "docProps"
PACKAGE_XL: Final = "xl"
PACKAGE_RELS: Final = "_rels"
PACKAGE_THEME: Final = "xl/theme"
PACKAGE_WORKSHEETS: Final = "xl/worksheets"
PACKAGE_CHARTSHEETS: Final = "xl/chartsheets"
PACKAGE_DRAWINGS: Final = "xl/drawings"
PACKAGE_CHARTS: Final = "xl/charts"
PACKAGE_IMAGES: Final = "xl/media"
PACKAGE_WORKSHEET_RELS: Final = "xl/worksheets/_rels"
PACKAGE_CHARTSHEETS_RELS: Final = "xl/chartsheets/_rels"
PACKAGE_PIVOT_TABLE: Final = "xl/pivotTables"
PACKAGE_PIVOT_CACHE: Final = "xl/pivotCache"

ARC_CONTENT_TYPES: Final = "[Content_Types].xml"
ARC_ROOT_RELS: Final = "_rels/.rels"
ARC_WORKBOOK_RELS: Final = "xl/_rels/workbook.xml.rels"
ARC_CORE: Final = "docProps/core.xml"
ARC_APP: Final = "docProps/app.xml"
ARC_CUSTOM: Final = "docProps/custom.xml"
ARC_WORKBOOK: Final = "xl/workbook.xml"
ARC_STYLE: Final = "xl/styles.xml"
ARC_THEME: Final = "xl/theme/theme1.xml"
ARC_SHARED_STRINGS: Final = "xl/sharedStrings.xml"
ARC_CUSTOM_UI: Final = "customUI/customUI.xml"

DCORE_NS: Final = "http://purl.org/dc/elements/1.1/"
DCTERMS_NS: Final = "http://purl.org/dc/terms/"
DCTERMS_PREFIX: Final = "dcterms"

DOC_NS: Final[str]
REL_NS: Final[str]
COMMENTS_NS: Final[str]
IMAGE_NS: Final[str]
VML_NS: Final[str]
VTYPES_NS: Final[str]
XPROPS_NS: Final[str]
CUSTPROPS_NS: Final[str]
EXTERNAL_LINK_NS: Final[str]

CPROPS_FMTID: Final = "{D5CDD505-2E9C-101B-9397-08002B2CF9AE}"

PKG_NS: Final = "http://schemas.openxmlformats.org/package/2006/"
PKG_REL_NS: Final[str]
COREPROPS_NS: Final[str]
CONTYPES_NS: Final[str]

XSI_NS: Final = "http://www.w3.org/2001/XMLSchema-instance"
XML_NS: Final = "http://www.w3.org/XML/1998/namespace"
SHEET_MAIN_NS: Final[str]

CHART_NS: Final[str]
DRAWING_NS: Final[str]
SHEET_DRAWING_NS: Final[str]
CHART_DRAWING_NS: Final[str]

CUSTOMUI_NS: Final[str]

NAMESPACES: Final[dict[str, str]]

WORKBOOK_MACRO: Final = "application/vnd.ms-excel.%s.macroEnabled.main+xml"
WORKBOOK: Final[str]
SPREADSHEET: Final[str]
SHARED_STRINGS: Final[str]
EXTERNAL_LINK: Final[str]
WORKSHEET_TYPE: Final[str]
COMMENTS_TYPE: Final[str]
STYLES_TYPE: Final[str]
CHARTSHEET_TYPE: Final[str]
DRAWING_TYPE: Final[str]
CHART_TYPE: Final[str]
CHARTSHAPE_TYPE: Final[str]
THEME_TYPE: Final[str]
CPROPS_TYPE: Final[str]
XLTM: Final[str]
XLSM: Final[str]
XLTX: Final[str]
XLSX: Final[str]

EXT_TYPES: Final[dict[str, str]]

CTRL: Final = "application/vnd.ms-excel.controlproperties+xml"
ACTIVEX: Final = "application/vnd.ms-office.activeX+xml"
VBA: Final = "application/vnd.ms-office.vbaProject"
