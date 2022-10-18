"""
Patched ``LZMAFile`` to handle pickle protocol 5.
"""

from __future__ import annotations

import lzma

from pandas.compat._utils import flatten_buffer


class LZMAFile(lzma.LZMAFile):
    def write(self, b) -> int:
        # Workaround issue where `lzma.LZMAFile` expects `len`
        # to return the number of bytes in `b` by converting
        # `b` into something that meets that constraint with
        # minimal copying.
        return super().write(flatten_buffer(b))
