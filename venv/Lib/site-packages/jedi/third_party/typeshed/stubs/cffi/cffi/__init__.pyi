from typing import Final

from .api import FFI as FFI
from .error import (
    CDefError as CDefError,
    FFIError as FFIError,
    PkgConfigError as PkgConfigError,
    VerificationError as VerificationError,
    VerificationMissing as VerificationMissing,
)

__all__ = ["FFI", "VerificationError", "VerificationMissing", "CDefError", "FFIError"]
__version__: Final[str]
__version_info__: Final[tuple[int, int, int]]
__version_verifier_modules__: Final[str]
