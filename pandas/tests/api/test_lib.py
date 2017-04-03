# -*- coding: utf-8 -*-

from warnings import catch_warnings
import pandas  # noqa


def test_moved_infer_dtype():
    with catch_warnings(record=True):
        e = pandas.lib.infer_dtype('foo')
        assert e is not None
