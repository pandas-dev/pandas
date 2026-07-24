from _typeshed import ConvertibleToInt
from collections.abc import Sequence
from typing import Any, Generic, NamedTuple, Protocol, TypeVar, overload, type_check_only
from typing_extensions import TypeAlias

from ._types import ErrorCorrect, MaskPattern
from .image.base import BaseImage
from .image.pil import PilImage
from .image.pure import PyPNGImage
from .util import QRData

ModulesType: TypeAlias = list[list[bool | None]]
precomputed_qr_blanks: dict[int, ModulesType]

_DefaultImage: TypeAlias = PilImage | PyPNGImage  # PilImage if Pillow is installed, PyPNGImage otherwise
_AnySeq = TypeVar("_AnySeq", bound=Sequence[Any])

@overload
def make(
    data: QRData | bytes | str,
    *,
    version: ConvertibleToInt | None = None,
    error_correction: ErrorCorrect = 0,
    box_size: ConvertibleToInt = 10,
    border: ConvertibleToInt = 4,
    image_factory: None = None,
    mask_pattern: MaskPattern | None = None,
) -> _DefaultImage: ...
@overload
def make(
    data: QRData | bytes | str,
    *,
    version: ConvertibleToInt | None = None,
    error_correction: ErrorCorrect = 0,
    box_size: ConvertibleToInt = 10,
    border: ConvertibleToInt = 4,
    image_factory: type[GenericImage],
    mask_pattern: MaskPattern | None = None,
) -> GenericImage: ...
def copy_2d_array(x: Sequence[_AnySeq]) -> list[_AnySeq]: ...

class ActiveWithNeighbors(NamedTuple):
    NW: bool
    N: bool
    NE: bool
    W: bool
    me: bool
    E: bool
    SW: bool
    S: bool
    SE: bool
    def __bool__(self) -> bool: ...

GenericImage = TypeVar("GenericImage", bound=BaseImage)  # noqa: Y001
GenericImageLocal = TypeVar("GenericImageLocal", bound=BaseImage)  # noqa: Y001

@type_check_only
class _TTYWriter(Protocol):
    def isatty(self) -> bool: ...
    def write(self, s: str, /) -> object: ...
    def flush(self) -> object: ...

class QRCode(Generic[GenericImage]):
    modules: ModulesType
    error_correction: ErrorCorrect
    box_size: int
    border: int
    image_factory: type[GenericImage] | None
    @overload
    def __init__(
        self,
        version: ConvertibleToInt | None,
        error_correction: ErrorCorrect,
        box_size: ConvertibleToInt,
        border: ConvertibleToInt,
        image_factory: type[GenericImage],
        mask_pattern: MaskPattern | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        version: ConvertibleToInt | None = None,
        error_correction: ErrorCorrect = 0,
        box_size: ConvertibleToInt = 10,
        border: ConvertibleToInt = 4,
        *,
        image_factory: type[GenericImage],
        mask_pattern: MaskPattern | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self: QRCode[_DefaultImage],
        version: ConvertibleToInt | None = None,
        error_correction: ErrorCorrect = 0,
        box_size: ConvertibleToInt = 10,
        border: ConvertibleToInt = 4,
        image_factory: None = None,
        mask_pattern: MaskPattern | None = None,
    ) -> None: ...
    @property
    def version(self) -> int: ...
    @version.setter
    def version(self, value: ConvertibleToInt | None) -> None: ...
    @property
    def mask_pattern(self) -> MaskPattern | None: ...
    @mask_pattern.setter
    def mask_pattern(self, pattern: MaskPattern | None) -> None: ...
    modules_count: int
    data_cache: list[int]
    data_list: list[QRData]
    def clear(self) -> None: ...
    def add_data(self, data: QRData | bytes | str, optimize: int = 20) -> None: ...
    def make(self, fit: bool = True) -> None: ...
    def makeImpl(self, test: bool, mask_pattern: MaskPattern) -> None: ...
    def setup_position_probe_pattern(self, row: int, col: int) -> None: ...
    def best_fit(self, start: int | None = None) -> int: ...
    def best_mask_pattern(self) -> int: ...
    def print_tty(self, out: _TTYWriter | None = None) -> None: ...
    def print_ascii(self, out: _TTYWriter | None = None, tty: bool = False, invert: bool = False) -> None: ...
    # kwargs are passed on to the specific image factory used, and in turn passed through to
    # their make_image method.
    @overload
    def make_image(self, image_factory: None = None, **kwargs: Any) -> GenericImage: ...
    @overload
    def make_image(self, image_factory: type[GenericImageLocal], **kwargs: Any) -> GenericImageLocal: ...
    def is_constrained(self, row: int, col: int) -> bool: ...
    def setup_timing_pattern(self) -> None: ...
    def setup_position_adjust_pattern(self) -> None: ...
    def setup_type_number(self, test: bool) -> None: ...
    def setup_type_info(self, test: bool, mask_pattern: MaskPattern) -> None: ...
    def map_data(self, data: Sequence[int], mask_pattern: MaskPattern) -> None: ...
    def get_matrix(self) -> list[list[bool]]: ...
    def active_with_neighbors(self, row: int, col: int) -> ActiveWithNeighbors: ...
