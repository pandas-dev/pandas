# -*- coding: utf-8 -*-

import inspect

import pytest


def check_accessible(x, depth=0):
    if depth == 0:
        return

    for _, member in inspect.getmembers(x):
        check_accessible(member, depth - 1)


@pytest.mark.filterwarnings("ignore::FutureWarning")
def test_libs():
    import pandas._libs as libs
    check_accessible(libs, 3)
