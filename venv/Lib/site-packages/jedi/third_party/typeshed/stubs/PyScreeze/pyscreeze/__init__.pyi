import sys
from _typeshed import ConvertibleToFloat, Incomplete, StrOrBytesPath, Unused
from collections.abc import Callable, Generator
from typing import Final, NamedTuple, TypeVar, overload
from typing_extensions import ParamSpec, TypeAlias

from PIL import Image

_P = ParamSpec("_P")
_R = TypeVar("_R")
# cv2.typing.MatLike: is an alias for `numpy.ndarray | cv2.mat_wrapper.Mat`, Mat extends ndarray.
# But can't import either, because pyscreeze does not declare them as dependencies, stub_uploader won't let it.
_MatLike: TypeAlias = Incomplete

PILLOW_VERSION: Final[tuple[int, int, int]]
RUNNING_PYTHON_2: Final = False
SCROT_EXISTS: Final[bool]
GNOMESCREENSHOT_EXISTS: Final[bool]

if sys.platform == "linux":
    RUNNING_X11: Final[bool]
    RUNNING_WAYLAND: Final[bool]

# Meant to be overridable as a setting
GRAYSCALE_DEFAULT: bool
# Meant to be overridable for backward-compatibility
USE_IMAGE_NOT_FOUND_EXCEPTION: bool

class Box(NamedTuple):
    left: int
    top: int
    width: int
    height: int

class Point(NamedTuple):
    x: int
    y: int

class RGB(NamedTuple):
    red: int
    green: int
    blue: int

class PyScreezeException(Exception): ...
class ImageNotFoundException(PyScreezeException): ...

def requiresPyGetWindow(wrappedFunction: Callable[_P, _R]) -> Callable[_P, _R]: ...

# _locateAll_opencv
@overload
def locate(
    needleImage: str | Image.Image | _MatLike,
    haystackImage: str | Image.Image | _MatLike,
    *,
    grayscale: bool | None = None,
    limit: Unused = 1,
    region: tuple[int, int, int, int] | None = None,
    step: int = 1,
    confidence: ConvertibleToFloat = 0.999,
) -> Box | None: ...

# _locateAll_pillow
@overload
def locate(
    needleImage: str | Image.Image,
    haystackImage: str | Image.Image,
    *,
    grayscale: bool | None = None,
    limit: Unused = 1,
    region: tuple[int, int, int, int] | None = None,
    step: int = 1,
    confidence: None = None,
) -> Box | None: ...

# _locateAll_opencv
@overload
def locateOnScreen(
    image: str | Image.Image | _MatLike,
    minSearchTime: float = 0,
    *,
    grayscale: bool | None = None,
    limit: Unused = 1,
    region: tuple[int, int, int, int] | None = None,
    step: int = 1,
    confidence: ConvertibleToFloat = 0.999,
) -> Box | None: ...

# _locateAll_pillow
@overload
def locateOnScreen(
    image: str | Image.Image,
    minSearchTime: float = 0,
    *,
    grayscale: bool | None = None,
    limit: Unused = 1,
    region: tuple[int, int, int, int] | None = None,
    step: int = 1,
    confidence: None = None,
) -> Box | None: ...

# _locateAll_opencv
@overload
def locateAllOnScreen(
    image: str | Image.Image | _MatLike,
    *,
    grayscale: bool | None = None,
    limit: int = 1000,
    region: tuple[int, int, int, int] | None = None,
    step: int = 1,
    confidence: ConvertibleToFloat = 0.999,
) -> Generator[Box]: ...

# _locateAll_pillow
@overload
def locateAllOnScreen(
    image: str | Image.Image,
    *,
    grayscale: bool | None = None,
    limit: int | None = None,
    region: tuple[int, int, int, int] | None = None,
    step: int = 1,
    confidence: None = None,
) -> Generator[Box]: ...

# _locateAll_opencv
@overload
def locateCenterOnScreen(
    image: str | Image.Image | _MatLike,
    *,
    minSearchTime: float = 0,
    grayscale: bool | None = None,
    limit: Unused = 1,
    region: tuple[int, int, int, int] | None = None,
    step: int = 1,
    confidence: ConvertibleToFloat = 0.999,
) -> Point | None: ...

# _locateAll_pillow
@overload
def locateCenterOnScreen(
    image: str | Image.Image,
    *,
    minSearchTime: float = 0,
    grayscale: bool | None = None,
    limit: Unused = 1,
    region: tuple[int, int, int, int] | None = None,
    step: int = 1,
    confidence: None = None,
) -> Point | None: ...
def locateOnScreenNear(image: str | Image.Image | _MatLike, x: int, y: int) -> Box: ...
def locateCenterOnScreenNear(image: str | Image.Image | _MatLike, x: int, y: int) -> Point | None: ...

# _locateAll_opencv
@overload
def locateOnWindow(
    image: str | Image.Image | _MatLike,
    title: str,
    *,
    grayscale: bool | None = None,
    limit: Unused = 1,
    step: int = 1,
    confidence: ConvertibleToFloat = 0.999,
) -> Box | None: ...

# _locateAll_pillow
@overload
def locateOnWindow(
    image: str | Image.Image,
    title: str,
    *,
    grayscale: bool | None = None,
    limit: Unused = 1,
    step: int = 1,
    confidence: None = None,
) -> Box | None: ...
def showRegionOnScreen(
    region: tuple[int, int, int, int], outlineColor: str = "red", filename: str = "_showRegionOnScreen.png"
) -> None: ...
def center(coords: tuple[int, int, int, int]) -> Point: ...
def pixelMatchesColor(
    x: int, y: int, expectedRGBColor: tuple[int, int, int] | tuple[int, int, int, int], tolerance: int = 0
) -> bool: ...
def pixel(x: int, y: int) -> tuple[int, int, int]: ...

if sys.platform == "win32":
    def screenshot(
        imageFilename: StrOrBytesPath | None = None, region: tuple[int, int, int, int] | None = None, allScreens: bool = False
    ) -> Image.Image: ...

else:
    def screenshot(
        imageFilename: StrOrBytesPath | None = None, region: tuple[int, int, int, int] | None = None
    ) -> Image.Image: ...

# _locateAll_opencv
@overload
def locateAll(
    needleImage: str | Image.Image | _MatLike,
    haystackImage: str | Image.Image | _MatLike,
    grayscale: bool | None = None,
    limit: int = 1000,
    region: tuple[int, int, int, int] | None = None,
    step: int = 1,
    confidence: ConvertibleToFloat = 0.999,
) -> Generator[Box]: ...

# _locateAll_pillow
@overload
def locateAll(
    needleImage: str | Image.Image,
    haystackImage: str | Image.Image,
    grayscale: bool | None = None,
    limit: int | None = None,
    region: tuple[int, int, int, int] | None = None,
    step: int = 1,
    confidence: None = None,
) -> Generator[Box]: ...
