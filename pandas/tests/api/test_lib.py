# -*- coding: utf-8 -*-

import pytest
import pandas  # noqa


@pytest.mark.parametrize('f', ['infer_dtype'])
def test_importable(f):
    from pandas.api import lib
    e = getattr(lib, f)
    assert e is not None
