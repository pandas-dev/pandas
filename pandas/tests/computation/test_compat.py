import importlib

import pytest

from pandas.compat._optional import VERSIONS

import pandas as pd
import pandas._testing as tm
from pandas.core.computation import (
    check,
    expr,
)
from pandas.core.computation.engines import ENGINES
from pandas.util.version import Version


def test_compat():
    # test we have compat with our version of numexpr

    ne = pytest.importorskip("numexpr")

    ver = ne.__version__
    if Version(ver) < Version(VERSIONS["numexpr"]):
        assert not check.NUMEXPR_INSTALLED
    elif Version(ver).base_version in check.BLOCKED_NUMEXPR_VERSIONS:
        assert not check.NUMEXPR_INSTALLED
    else:
        assert check.NUMEXPR_INSTALLED


def test_blocked_numexpr_version_not_used(monkeypatch):
    # GH#65408 numexpr releases known to silently return wrong results are not
    #  used even though they satisfy the minimum version
    ne = pytest.importorskip("numexpr")
    if not check.BLOCKED_NUMEXPR_VERSIONS:
        pytest.skip("no numexpr versions are currently blocked")

    monkeypatch.setattr(ne, "__version__", min(check.BLOCKED_NUMEXPR_VERSIONS))
    try:
        with tm.assert_produces_warning(
            UserWarning,
            match="can silently return incorrect",
            check_stacklevel=False,
        ):
            importlib.reload(check)
        assert not check.NUMEXPR_INSTALLED
    finally:
        monkeypatch.undo()
        importlib.reload(check)


@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize("parser", expr.PARSERS)
def test_invalid_numexpr_version(engine, parser):
    if engine == "numexpr" and not check.NUMEXPR_INSTALLED:
        pytest.skip("numexpr not installed or an unsupported version")
    a, b = 1, 2  # noqa: F841
    res = pd.eval("a + b", engine=engine, parser=parser)
    assert res == 3
