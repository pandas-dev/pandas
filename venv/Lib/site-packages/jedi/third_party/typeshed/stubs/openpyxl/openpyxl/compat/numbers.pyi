from decimal import Decimal
from typing import Final
from typing_extensions import TypeAlias

# NOTE: Can't specify numpy as a dependency because openpyxl doesn't declare it as one
# import numpy
# import numpy._typing
# _NBitBase: TypeAlias = numpy._typing.NBitBase
# _NumericTypes: TypeAlias = int | float | Decimal | numpy.bool_ | numpy.floating[_NBitBase] | numpy.integer[_NBitBase]

_NumericTypes: TypeAlias = int | float | Decimal
NUMERIC_TYPES: Final[tuple[type[_NumericTypes], ...]]

NUMPY: Final[bool]
