from typing import Any
from typing_extensions import TypeAlias

_Image: TypeAlias = Any  # actually PIL.Image, but not worth adding a dependency

def shift_hue(image: _Image, hue: float) -> _Image: ...
def make_transparent(image: _Image) -> _Image: ...
