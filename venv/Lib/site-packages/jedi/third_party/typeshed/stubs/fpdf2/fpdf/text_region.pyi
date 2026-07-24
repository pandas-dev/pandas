from _typeshed import Incomplete
from collections.abc import Sequence
from typing import NamedTuple
from typing_extensions import Self

from .enums import Align, WrapMode
from .image_datastructures import RasterImageInfo, VectorImageInfo, _TextAlign
from .line_break import Fragment, TextLine

class Extents(NamedTuple):
    left: float
    right: float

class TextRegionMixin:
    def __init__(self, *args, **kwargs) -> None: ...
    def register_text_region(self, region) -> None: ...
    def is_current_text_region(self, region) -> bool: ...
    def clear_text_region(self) -> None: ...

class LineWrapper(NamedTuple):
    line: Sequence[Incomplete]
    paragraph: Paragraph
    first_line: bool = False
    last_line: bool = False

class Bullet:
    fragments: Sequence[Fragment]
    text_line: TextLine
    r_margin: float
    rendered_flag: bool
    def __init__(self, bullet_fragments: Sequence[Fragment], text_line: TextLine, bullet_r_margin: float) -> None: ...
    def get_fragments_width(self) -> float: ...

class Paragraph:
    pdf: Incomplete
    text_align: Align
    line_height: Incomplete
    top_margin: Incomplete
    bottom_margin: Incomplete
    indent: float
    skip_leading_spaces: bool
    wrapmode: Incomplete
    bullet: Bullet | None
    first_line_indent: float

    def __init__(
        self,
        region,
        text_align: _TextAlign | None = None,
        line_height=None,
        top_margin: float = 0,
        bottom_margin: float = 0,
        indent: float = 0,
        bullet_r_margin: float | None = None,
        bullet_string: str = "",
        skip_leading_spaces: bool = False,
        wrapmode: WrapMode | None = None,
        first_line_indent: float = 0,
    ) -> None: ...
    def __enter__(self): ...
    def __exit__(self, exc_type, exc_value, traceback) -> None: ...
    def write(self, text: str, link=None): ...
    def generate_bullet_frags_and_tl(
        self, bullet_string: str, bullet_r_margin: float
    ) -> tuple[tuple[Fragment, ...], TextLine] | None: ...
    def ln(self, h: float | None = None) -> None: ...
    def build_lines(self, print_sh: bool) -> list[LineWrapper]: ...

class ImageParagraph:
    region: Incomplete
    name: Incomplete
    align: Align | None
    width: float | None
    height: float | None
    fill_width: bool
    keep_aspect_ratio: bool
    top_margin: float
    bottom_margin: float
    link: Incomplete | None
    title: Incomplete | None
    alt_text: Incomplete | None
    img: Incomplete | None
    info: Incomplete | None

    def __init__(
        self,
        region,
        name,
        align: _TextAlign | None = None,
        width: float | None = None,
        height: float | None = None,
        fill_width: bool = False,
        keep_aspect_ratio: bool = False,
        top_margin: float = 0,
        bottom_margin: float = 0,
        link=None,
        title=None,
        alt_text=None,
    ) -> None: ...
    def build_line(self) -> Self: ...
    def render(self, col_left: float, col_width: float, max_height: float) -> VectorImageInfo | RasterImageInfo: ...

class ParagraphCollectorMixin:
    pdf: Incomplete
    text_align: Align
    line_height: Incomplete
    print_sh: Incomplete
    wrapmode: Incomplete
    skip_leading_spaces: Incomplete
    def __init__(
        self,
        pdf,
        *args,
        text: str | None = None,
        text_align: _TextAlign = "LEFT",
        line_height: float = 1.0,
        print_sh: bool = False,
        skip_leading_spaces: bool = False,
        wrapmode: WrapMode | None = None,
        img=None,
        img_fill_width: bool = False,
        **kwargs,
    ) -> None: ...
    def __enter__(self): ...
    def __exit__(self, exc_type, exc_value, traceback) -> None: ...
    def write(self, text: str, link=None): ...
    def ln(self, h: float | None = None) -> None: ...
    def paragraph(
        self,
        text_align: _TextAlign | None = None,
        line_height=None,
        skip_leading_spaces: bool = False,
        top_margin: int = 0,
        bottom_margin: int = 0,
        indent: int = 0,
        bullet_string: str = "",
        bullet_r_margin: float | None = None,
        wrapmode: WrapMode | None = None,
        first_line_indent: float = 0,
    ) -> Paragraph: ...
    def end_paragraph(self) -> None: ...
    def image(
        self,
        name,
        align: _TextAlign | None = None,
        width: float | None = None,
        height: float | None = None,
        fill_width: bool = False,
        keep_aspect_ratio: bool = False,
        top_margin: float = 0,
        bottom_margin: float = 0,
        link=None,
        title=None,
        alt_text=None,
    ) -> None: ...

class TextRegion(ParagraphCollectorMixin):
    def current_x_extents(self, y, height) -> None: ...
    def collect_lines(self): ...
    def render(self) -> None: ...
    def get_width(self, height): ...

class TextColumnarMixin:
    l_margin: Incomplete
    r_margin: Incomplete
    def __init__(self, pdf, *args, l_margin=None, r_margin=None, **kwargs) -> None: ...

class TextColumns(TextRegion, TextColumnarMixin):
    balance: Incomplete
    def __init__(self, pdf, *args, ncols: int = 1, gutter: float = 10, balance: bool = False, **kwargs) -> None: ...
    def __enter__(self) -> Self: ...
    def new_column(self) -> None: ...
    def render(self) -> None: ...
    def current_x_extents(self, y, height): ...
