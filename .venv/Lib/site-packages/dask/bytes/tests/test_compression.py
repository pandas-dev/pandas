from __future__ import annotations

from io import BytesIO

import pytest
from fsspec.compression import compr

from dask.bytes.utils import compress


@pytest.mark.parametrize("fmt,File", compr.items())
def test_files(fmt, File):
    if fmt not in compress:
        pytest.skip("compression function not provided")
    if fmt is None:
        return
    data = b"1234" * 1000
    compressed = compress[fmt](data)

    b = BytesIO(compressed)
    g = File(b, mode="rb")
    data2 = g.read()
    g.close()
    assert data == data2
