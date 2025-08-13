from __future__ import annotations

from types import ModuleType

import pytest

pytest.importorskip("pandas")

DA_EXPORTED_SUBMODULES = {"backends", "dispatch", "demo"}


def test_api():
    """Tests that `dask.array.__all__` is correct"""
    import dask.dataframe as ddf

    member_dict = vars(ddf)
    members = set(member_dict)
    # unexported submodules
    members -= {"tests", "dask"}
    members -= {
        m
        for m, mod in member_dict.items()
        if m not in DA_EXPORTED_SUBMODULES
        if isinstance(mod, ModuleType)
        and mod.__package__
        and mod.__package__.startswith("dask.dataframe")
    }
    # imported utility modules
    members -= {"annotations"}
    # private utilities and `__dunder__` members
    members -= {m for m in members if m.startswith("_")}

    assert set(ddf.__all__) == members
