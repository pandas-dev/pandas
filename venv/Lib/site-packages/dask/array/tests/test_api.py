from __future__ import annotations

from types import ModuleType

import pytest

pytest.importorskip("numpy")

DA_EXPORTED_SUBMODULES = {"backends", "fft", "lib", "linalg", "ma", "overlap", "random"}


def test_api():
    """Tests that `dask.array.__all__` is correct"""
    import dask.array as da

    member_dict = vars(da)
    members = set(member_dict)
    # unexported submodules
    members -= {"tests"}
    members -= {
        m
        for m, mod in member_dict.items()
        if m not in DA_EXPORTED_SUBMODULES
        if isinstance(mod, ModuleType) and mod.__package__ == "dask.array"
    }
    # imported utility modules
    members -= {"annotations", "importlib", "warnings"}
    # private utilities and `__dunder__` members
    members -= {"ARRAY_EXPR_ENABLED"}
    members -= {m for m in members if m.startswith("_")}

    assert set(da.__all__) == members
