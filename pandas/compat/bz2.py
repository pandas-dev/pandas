"""
Patched ``BZ2File`` to handle pickle protocol 5.
"""

from __future__ import annotations

import bz2

from pandas.compat.pickle_compat import flatten_buffer


class BZ2File(bz2.BZ2File):
    def write(self, b) -> int:
        # Workaround issue where `bz2.BZ2File` expects `len`
        # to return the number of bytes in `b` by converting
        # `b` into something that meets that constraint with
        # minimal copying.
        return super().write(flatten_buffer(b))
