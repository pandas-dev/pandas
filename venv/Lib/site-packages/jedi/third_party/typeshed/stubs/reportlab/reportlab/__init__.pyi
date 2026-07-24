from _typeshed import SupportsRichComparison
from typing import Final, Literal, TypeVar

_SupportsRichComparisonT = TypeVar("_SupportsRichComparisonT", bound=SupportsRichComparison)

Version: Final[str]
__version__: Final[str]
__date__: Final[str]
__min_python_version__: Final[tuple[int, int]]

def cmp(a: _SupportsRichComparisonT, b: _SupportsRichComparisonT) -> Literal[-1, 0, 1]: ...
