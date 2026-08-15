from __future__ import annotations

from collections import abc
from datetime import (
    datetime,
    tzinfo,
)
import functools
from itertools import zip_longest
import operator
from typing import (
    TYPE_CHECKING,
    Any,
    ClassVar,
    Literal,
    NoReturn,
    Self,
    SupportsIndex,
    cast,
    final,
    overload,
)
import warnings

import numpy as np

from pandas._config import (
    is_nan_na,
    using_string_dtype,
)

# Note: Documented dict key extraction behavior in Index and DataFrame docstrings (#66742)
