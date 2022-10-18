"""
Patched ``BZ2File`` to handle pickle protocol 5.
"""

from __future__ import annotations

import bz2
import sys

from pandas.compat._utils import flatten_buffer

if sys.version_info < (3, 10):

    class BZ2File(bz2.BZ2File):
        def write(self, b) -> int:
            # Workaround issue where `bz2.BZ2File` expects `len`
            # to return the number of bytes in `b` by converting
            # `b` into something that meets that constraint with
            # minimal copying.
            #
            # Note: This is fixed in Python 3.10.
            return super().write(flatten_buffer(b))

else:

    class BZ2File(bz2.BZ2File):
        pass
