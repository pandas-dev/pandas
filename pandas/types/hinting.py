#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np

import typing  # noqa
from typing import (  # noqa
    TypeVar, AnyStr, Any, Callable, Optional, Tuple, Union,
    Dict, Text, Iterable
)

from pandas.core.dtypes.generic import ABCSeries, ABCIndexClass

Buffer = Any
ArrayLike = TypeVar('ArrayLike', Buffer, list, dict, np.array)
Scalar = TypeVar('Scalar', int, float)
PythonScalar = TypeVar('PythonScalar', int, float, AnyStr)

SelectionKey = Union[str, list, tuple, ABCSeries, ABCIndexClass, np.ndarray]

# An argument to `.agg/.transform/.apply`
SelectionFunction = Union[str, Callable]
