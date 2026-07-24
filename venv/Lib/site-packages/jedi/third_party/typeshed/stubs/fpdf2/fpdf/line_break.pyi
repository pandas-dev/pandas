from _typeshed import Incomplete
from collections.abc import Callable, Sequence
from typing import Final, NamedTuple
from uuid import UUID

from .enums import Align, TextDirection, WrapMode

SOFT_HYPHEN: Final[str]
HYPHEN: Final[str]
SPACE: Final[str]
BREAKING_SPACE_SYMBOLS: Final[list[str]]
BREAKING_SPACE_SYMBOLS_STR: Final[str]
NBSP: Final[str]
NEWLINE: Final[str]
FORM_FEED: Final[str]

class Fragment:
    characters: list[str]
    graphics_state: dict[str, Incomplete]
    k: float
    url: str | None
    def __init__(
        self, characters: list[str] | str, graphics_state: dict[str, Incomplete], k: float, link: str | int | None = None
    ) -> None: ...
    @property
    def font(self): ...
    @font.setter
    def font(self, v) -> None: ...
    @property
    def is_ttf_font(self): ...
    @property
    def font_style(self): ...
    @property
    def font_family(self): ...
    @property
    def font_size_pt(self): ...
    @property
    def font_size(self): ...
    @property
    def font_stretching(self): ...
    @property
    def char_spacing(self): ...
    @property
    def text_mode(self): ...
    @property
    def underline(self) -> bool: ...
    @property
    def strikethrough(self) -> bool: ...
    @property
    def draw_color(self): ...
    @property
    def fill_color(self): ...
    @property
    def text_color(self): ...
    @property
    def line_width(self): ...
    @property
    def char_vpos(self): ...
    @property
    def lift(self): ...
    @property
    def string(self) -> str: ...
    @property
    def width(self) -> float: ...
    @property
    def text_shaping_parameters(self): ...
    @property
    def paragraph_direction(self) -> TextDirection: ...
    @property
    def fragment_direction(self) -> TextDirection: ...
    def trim(self, index: int) -> None: ...
    def __eq__(self, other: Fragment) -> bool: ...  # type: ignore[override]
    def get_width(self, start: int = 0, end: int | None = None, chars: str | None = None, initial_cs: bool = True) -> float: ...
    def has_same_style(self, other: Fragment) -> bool: ...
    def get_character_width(self, character: str, print_sh: bool = False, initial_cs: bool = True): ...
    def render_pdf_text(self, frag_ws, current_ws, word_spacing, adjust_x, adjust_y, h): ...
    def render_pdf_text_ttf(self, frag_ws, word_spacing): ...
    def render_with_text_shaping(self, pos_x: float, pos_y: float, h: float, word_spacing: float) -> str: ...
    def render_pdf_text_core(self, frag_ws, current_ws): ...

class TotalPagesSubstitutionFragment(Fragment):
    uuid: UUID
    def get_placeholder_string(self) -> str: ...
    def render_text_substitution(self, replacement_text: str) -> str: ...

class TextLine(NamedTuple):
    fragments: tuple[Fragment, ...]
    text_width: float
    number_of_spaces: int
    align: Align
    height: float
    max_width: float
    trailing_nl: bool = False
    trailing_form_feed: bool = False
    indent: float = 0
    def get_ordered_fragments(self) -> tuple[Fragment, ...]: ...

class SpaceHint(NamedTuple):
    original_fragment_index: int
    original_character_index: int
    current_line_fragment_index: int
    current_line_character_index: int
    line_width: float
    number_of_spaces: int

class HyphenHint(NamedTuple):
    original_fragment_index: int
    original_character_index: int
    current_line_fragment_index: int
    current_line_character_index: int
    line_width: float
    number_of_spaces: int
    curchar: str
    curchar_width: float
    graphics_state: dict[str, Incomplete]
    k: float

class CurrentLine:
    max_width: float
    print_sh: bool
    indent: float
    fragments: list[Fragment]
    height: int
    number_of_spaces: int
    space_break_hint: Incomplete
    hyphen_break_hint: Incomplete
    def __init__(self, max_width: float, print_sh: bool = False, indent: float = 0) -> None: ...
    @property
    def width(self) -> float: ...
    def add_character(
        self,
        character: str,
        character_width: float,
        original_fragment: Fragment,
        original_fragment_index: int,
        original_character_index: int,
        height: float,
        url: str | None = None,
    ) -> None: ...
    def trim_trailing_spaces(self) -> None: ...
    def manual_break(self, align: Align, trailing_nl: bool = False, trailing_form_feed: bool = False) -> TextLine: ...
    def automatic_break_possible(self) -> bool: ...
    def automatic_break(self, align: Align) -> tuple[Incomplete, Incomplete, TextLine]: ...

class MultiLineBreak:
    fragments: Sequence[Fragment]
    get_width: float
    margins: Sequence[float]
    align: Align
    print_sh: bool
    wrapmode: WrapMode
    line_height: float
    skip_leading_spaces: bool
    fragment_index: int
    character_index: int
    idx_last_forced_break: int | None
    first_line_indent: float
    def __init__(
        self,
        fragments: Sequence[Fragment],
        max_width: float | Callable[[float], float],
        margins: Sequence[float],
        align: Align = ...,
        print_sh: bool = False,
        wrapmode: WrapMode = ...,
        line_height: float = 1.0,
        skip_leading_spaces: bool = False,
        first_line_indent: float = 0,
    ) -> None: ...
    def get_line(self) -> TextLine: ...
