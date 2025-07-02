"""
All of pandas' ExtensionArrays.

See :ref:`extending.extension-types` for more.
"""

# Explanations for the ExtensionArrays listed below:
#
# ArrowExtensionArray: Wraps an Apache Arrow array for high-performance data handling.
#
# ArrowStringArray:    A specific version for string data backed by Apache Arrow.
#
# BooleanArray:        Stores boolean data (True, False) but with the
#                      ability to hold missing values (NA).
#
# Categorical:         For data that belongs to a fixed, finite set of
#                      categories. Very memory efficient.
#
# DatetimeArray:       Handles timezone-aware or timezone-naive date and time data.
#
# FloatingArray:       For floating-point (decimal) numbers, with support
#                      for missing values.
#
# IntegerArray:        For integer (whole) numbers, with support for missing values.
#
# IntervalArray:       Stores data representing intervals or ranges (e.g., 1-5, 6-10).
#
# NumpyExtensionArray: A wrapper around a standard NumPy array, for compatibility
#                      within the extension system.
#
# PeriodArray:         For data representing regular time periods (e.g., months,
#                      quarters, years).
#
# SparseArray:         Memory-efficient array for data that is mostly zero or NA.
#
# StringArray:         Dedicated array for handling string data, with support
#                      for missing values.
#
# TimedeltaArray:      For data representing durations or differences in time
#                      (e.g., 2 days, 5 hours).

from pandas.core.arrays import (
    ArrowExtensionArray,
    ArrowStringArray,
    BooleanArray,
    Categorical,
    DatetimeArray,
    FloatingArray,
    IntegerArray,
    IntervalArray,
    NumpyExtensionArray,
    PeriodArray,
    SparseArray,
    StringArray,
    TimedeltaArray,
)

__all__ = [
    "ArrowExtensionArray",
    "ArrowStringArray",
    "BooleanArray",
    "Categorical",
    "DatetimeArray",
    "FloatingArray",
    "IntegerArray",
    "IntervalArray",
    "NumpyExtensionArray",
    "PeriodArray",
    "SparseArray",
    "StringArray",
    "TimedeltaArray",
]
