"""
Patched ``LZMAFile`` to handle pickle protocol 5.
"""

from __future__ import annotations

import lzma
import sys

from pandas.compat._utils import flatten_buffer


if sys.version_info < (3, 10):

    class LZMAFile(lzma.LZMAFile):
        def write(self, b) -> int:
            # Workaround issue where `lzma.LZMAFile` expects `len`
            # to return the number of bytes in `b` by converting
            # `b` into something that meets that constraint with
            # minimal copying.
            #
            # Note: This is fixed in Python 3.10.
            return super().write(flatten_buffer(b))

else:

    class LZMAFile(lzma.LZMAFile):
        pass
