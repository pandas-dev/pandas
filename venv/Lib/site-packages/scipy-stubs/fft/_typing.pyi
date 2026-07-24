# NOTE(scipy-stubs): This private module should not be used outside of scipy-stubs

from typing import Literal

__all__ = "DCTType", "NormalizationMode"

type DCTType = Literal[1, 2, 3, 4]
type NormalizationMode = Literal["backward", "ortho", "forward"]
