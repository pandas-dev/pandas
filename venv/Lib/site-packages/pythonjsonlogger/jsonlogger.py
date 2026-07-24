"""Stub module retained for compatibility.

It retains access to old names whilst sending deprecation warnings.
"""

# pylint: disable=wrong-import-position,unused-import

import warnings

## Throw warning
warnings.warn(
    "pythonjsonlogger.jsonlogger has been moved to pythonjsonlogger.json",
    DeprecationWarning,
)

## Import names
from .json import JsonFormatter, JsonEncoder
from .core import RESERVED_ATTRS
