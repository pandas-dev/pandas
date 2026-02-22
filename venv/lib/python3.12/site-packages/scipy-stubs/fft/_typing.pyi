# NOTE(scipy-stubs): This private module should not be used outside of scipy-stubs

from typing import Literal, TypeAlias

__all__ = "DCTType", "NormalizationMode"

DCTType: TypeAlias = Literal[1, 2, 3, 4]
NormalizationMode: TypeAlias = Literal["backward", "ortho", "forward"]
