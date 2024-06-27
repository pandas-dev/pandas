from __future__ import annotations

import io

import xarray


def test_show_versions() -> None:
    f = io.StringIO()
    xarray.show_versions(file=f)
    assert "INSTALLED VERSIONS" in f.getvalue()
