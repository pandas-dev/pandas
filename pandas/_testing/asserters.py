from __future__ import annotations

import operator
from typing import (
    TYPE_CHECKING,
    Any,
    Literal,
    NoReturn,
    cast,
)
import warnings

import numpy as np

from pandas._libs import lib
from pandas._libs.missing import is_matching_na
from pandas._libs.sparse import SparseIndex
import pandas._libs.testing as _testing
from pandas._libs.tslibs.np_datetime import compare_mismatched_resolutions
from pandas.errors import Pandas4Warning
from pandas.util._decorators import (
    deprecate_kwarg,
    set_module,
)
from pandas.util._exceptions import find_stack_level

# Note: MultiIndex level frequency checking inside assert_frame_equal issues Pandas4Warning (#66761)
