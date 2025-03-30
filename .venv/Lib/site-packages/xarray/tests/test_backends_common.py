from __future__ import annotations

import pytest

from xarray.backends.common import robust_getitem


class DummyFailure(Exception):
    pass


class DummyArray:
    def __init__(self, failures):
        self.failures = failures

    def __getitem__(self, key):
        if self.failures:
            self.failures -= 1
            raise DummyFailure
        return "success"


def test_robust_getitem() -> None:
    array = DummyArray(failures=2)
    with pytest.raises(DummyFailure):
        array[...]
    result = robust_getitem(array, ..., catch=DummyFailure, initial_delay=1)
    assert result == "success"

    array = DummyArray(failures=3)
    with pytest.raises(DummyFailure):
        robust_getitem(array, ..., catch=DummyFailure, initial_delay=1, max_retries=2)
