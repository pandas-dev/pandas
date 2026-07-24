from _typeshed import Incomplete, SupportsKeysAndGetItem
from collections.abc import Callable, Iterable, Mapping
from html.parser import HTMLParser
from logging import Logger
from typing import ClassVar, Final, Literal
from typing_extensions import TypeAlias

from fpdf import FPDF

from .enums import Align, TextEmphasis
from .fonts import FontFace
from .table import Row, Table

__author__: Final[str]
__copyright__: Final[str]

_OLType: TypeAlias = Literal["1", "a", "A", "I", "i"]

LOGGER: Logger
MESSAGE_WAITING_WIN1252: Final = "\x95"
BULLET_UNICODE: Final = "•"
DEGREE_SIGN_WIN1252: Final = "\xb0"
RING_OPERATOR_UNICODE: Final = "∘"
HEADING_TAGS: Final[tuple[str, ...]]
DEFAULT_TAG_STYLES: Final[dict[str, FontFace]]
INLINE_TAGS: Final[tuple[str, ...]]
BLOCK_TAGS: Final[tuple[str, ...]]

COLOR_DICT: Final[dict[str, str]]

def color_as_decimal(color: str | None = "#000000") -> tuple[int, int, int] | None: ...
def parse_css_style(style_attr: str) -> dict[str, str]: ...

class HTML2FPDF(HTMLParser):
    HTML_UNCLOSED_TAGS: ClassVar[tuple[str, ...]]
    TABLE_LINE_HEIGHT: ClassVar[float]

    pdf: FPDF
    image_map: Callable[[str], str]
    ul_bullet_char: str
    li_prefix_color: tuple[int, int, int]
    warn_on_tags_not_matching: bool

    font_family: str
    font_size_pt: float
    font_emphasis: TextEmphasis
    font_color: tuple[int, int, int]

    style_stack: list[FontFace]
    h: float
    follows_trailing_space: bool
    follows_heading: bool
    href: str
    align: float | Align | None
    indent: int
    line_height_stack: list[Incomplete]
    ol_type: dict[int, _OLType]
    bullet: list[Incomplete]
    heading_level: Incomplete | None
    render_title_tag: bool
    table_line_separators: bool
    table: Table | None
    table_row: Row | None
    tr: dict[str, str] | None
    td_th: dict[str, str] | None
    tag_indents: dict[str, int]
    tag_styles: dict[str, FontFace]

    def __init__(
        self,
        pdf: FPDF,
        image_map: Callable[[str], str] | None = None,
        li_tag_indent: int | None = None,
        dd_tag_indent: int | None = None,
        table_line_separators: bool = False,
        ul_bullet_char: str = "disc",
        li_prefix_color: tuple[int, int, int] = (190, 0, 0),
        heading_sizes: SupportsKeysAndGetItem[str, int] | Iterable[tuple[str, int]] | None = None,
        pre_code_font: str | None = None,
        warn_on_tags_not_matching: bool = True,
        tag_indents: dict[str, int] | None = None,
        tag_styles: Mapping[str, FontFace] | None = None,
        font_family: str = "times",
        render_title_tag: bool = False,
    ) -> None: ...
    def handle_data(self, data) -> None: ...
    def handle_starttag(self, tag, attrs) -> None: ...
    def handle_endtag(self, tag) -> None: ...
    def put_link(self, text) -> None: ...
    def render_toc(self, pdf, outline) -> None: ...
    def error(self, message: str) -> None: ...

def ul_prefix(ul_type: str, is_ttf_font: bool | None) -> str: ...
def ol_prefix(ol_type: _OLType, index: int) -> str: ...

class HTMLMixin:
    def __init__(self, *args, **kwargs) -> None: ...
