# Copyright 2017 Virgil Dupras

# This software is licensed under the "BSD" License as described in the "LICENSE" file,
# which should be included with this package. The terms are also available at
# http://www.hardcoded.net/licenses/bsd_license

import sys
import os

PY3 = sys.version_info[0] >= 3
if PY3:
    text_type = str
    binary_type = bytes
    if os.supports_bytes_environ:
        # environb will be unset under Windows, but then again we're not supposed to use it.
        environb = os.environb
else:
    text_type = unicode  # noqa: F821
    binary_type = str
    environb = os.environ

try:
    from collections.abc import Iterable as iterable_type
except ImportError:
    from collections import Iterable as iterable_type  # noqa: F401
