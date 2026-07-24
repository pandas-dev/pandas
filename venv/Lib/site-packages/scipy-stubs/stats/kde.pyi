# This module is not meant for public use and will be removed in SciPy v2.0.0.
from typing_extensions import deprecated

from ._kde import gaussian_kde as _gaussian_kde

__all__ = ["gaussian_kde"]

@deprecated("will be removed in SciPy v2.0.0")
class gaussian_kde(_gaussian_kde): ...
