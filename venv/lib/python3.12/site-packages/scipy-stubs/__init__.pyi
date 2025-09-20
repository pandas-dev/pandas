from collections.abc import Sequence
from typing import Final, Literal, TypeAlias

from numpy import __version__ as __numpy_version__  # noqa: ICN003

from . import (
    cluster,
    constants,
    datasets,
    differentiate,
    fft,
    fftpack,
    integrate,
    interpolate,
    io,
    linalg,
    ndimage,
    odr,
    optimize,
    signal,
    sparse,
    spatial,
    special,
    stats,
)
from .__config__ import show as show_config
from ._lib._ccallback import LowLevelCallable
from ._lib._testutils import PytestTester
from .version import version as __version__

__all__ = [
    "LowLevelCallable",
    "__version__",
    "cluster",
    "constants",
    "datasets",
    "differentiate",
    "fft",
    "fftpack",
    "integrate",
    "interpolate",
    "io",
    "linalg",
    "ndimage",
    "odr",
    "optimize",
    "show_config",
    "signal",
    "sparse",
    "spatial",
    "special",
    "stats",
    "test",
]

###

_SubModule: TypeAlias = Literal[
    "cluster",
    "constants",
    "datasets",
    "differentiate",
    "fft",
    "fftpack",
    "integrate",
    "interpolate",
    "io",
    "linalg",
    "ndimage",
    "odr",
    "optimize",
    "signal",
    "sparse",
    "spatial",
    "special",
    "stats",
]

###

np_minversion: Final = "1.25.2"  # undocumented
np_maxversion: Final = "2.6.0"  # undocumented
test: Final[PytestTester] = ...  # undocumented
submodules: Final[Sequence[_SubModule]] = ...  # undocumented
