from __future__ import annotations

import bz2
import gzip
import io
import lzma
import zipfile


def zip_compress(data):
    """Write data into zipfile and return the bytes"""
    out = io.BytesIO()
    with zipfile.ZipFile(file=out, mode="w") as z:
        with z.open("myfile", "w") as zf:
            zf.write(data)
    out.seek(0)
    return out.read()


compress = {
    "gzip": gzip.compress,
    "bz2": bz2.compress,
    None: lambda x: x,
    "xz": lzma.compress,
    "zip": zip_compress,
}
