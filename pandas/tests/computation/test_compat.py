import importlib
import operator

import numpy as np
import pytest

from pandas.compat._optional import VERSIONS

import pandas as pd
import pandas._testing as tm
from pandas.core.computation import (
    check,
    expr,
    expressions,
)
from pandas.core.computation.engines import ENGINES
from pandas.util.version import Version


def test_compat():
    # test we have compat with our version of numexpr

    ne = pytest.importorskip("numexpr")

    ver = ne.__version__
    if (
        Version(ver) < Version(VERSIONS["numexpr"])
        or Version(ver).base_version in check.BLOCKED_NUMEXPR_VERSIONS
    ):
        assert not check.NUMEXPR_INSTALLED
    else:
        assert check.NUMEXPR_INSTALLED


@pytest.fixture
def blocked_numexpr(monkeypatch):
    """
    Behave as though the installed numexpr were a blocked version, and yield it.
    """
    if not check.BLOCKED_NUMEXPR_VERSIONS:
        pytest.skip("no numexpr versions are currently blocked")

    version = min(check.BLOCKED_NUMEXPR_VERSIONS)
    monkeypatch.setattr(check, "NUMEXPR_INSTALLED", False)
    monkeypatch.setattr(check, "NUMEXPR_BLOCKED_VERSION", version)
    monkeypatch.setattr(check, "_warned_blocked", False)
    monkeypatch.setattr(expressions, "NUMEXPR_INSTALLED", False)
    monkeypatch.setattr(expressions, "NUMEXPR_BLOCKED_VERSION", version)
    monkeypatch.setattr(expressions, "USE_NUMEXPR", False)

    yield version

    # restore the _evaluate/_where dispatch a test may have re-chosen
    monkeypatch.undo()
    expressions.set_use_numexpr(pd.get_option("compute.use_numexpr"))


def test_blocked_numexpr_version_not_used(monkeypatch):
    # GH#65408 numexpr releases known to silently return wrong results are not
    #  used even though they satisfy the minimum version
    ne = pytest.importorskip("numexpr")
    if not check.BLOCKED_NUMEXPR_VERSIONS:
        pytest.skip("no numexpr versions are currently blocked")

    monkeypatch.setattr(ne, "__version__", min(check.BLOCKED_NUMEXPR_VERSIONS))
    try:
        # GH#66956 importing pandas does not warn; the warning is deferred to
        #  the first operation that would have used numexpr
        with tm.assert_produces_warning(None):
            importlib.reload(check)
        assert not check.NUMEXPR_INSTALLED
        assert check.NUMEXPR_BLOCKED_VERSION == min(check.BLOCKED_NUMEXPR_VERSIONS)
    finally:
        monkeypatch.undo()
        importlib.reload(check)


def test_blocked_numexpr_warns_once(blocked_numexpr):
    # GH#66956
    with tm.assert_produces_warning(UserWarning, match="can silently return incorrect"):
        check.warn_numexpr_blocked()

    with tm.assert_produces_warning(None):
        check.warn_numexpr_blocked()


def test_blocked_numexpr_warns_on_arithmetic(blocked_numexpr, monkeypatch):
    # GH#66956 an operation large enough to have been offloaded warns, a
    #  smaller one does not
    monkeypatch.setattr(expressions, "_MIN_ELEMENTS", 10)
    expressions.set_use_numexpr(True)
    assert expressions._evaluate is expressions._evaluate_blocked

    small = np.ones(5)
    with tm.assert_produces_warning(None):
        expressions.evaluate(operator.add, small, small)

    arr = np.ones(11)
    with tm.assert_produces_warning(UserWarning, match="can silently return incorrect"):
        expressions.evaluate(operator.add, arr, arr)


def test_blocked_numexpr_warns_on_where(blocked_numexpr, monkeypatch):
    # GH#66956
    monkeypatch.setattr(expressions, "_MIN_ELEMENTS", 10)
    expressions.set_use_numexpr(True)
    assert expressions._where is expressions._where_blocked

    arr = np.ones(11)
    cond = np.ones(11, dtype=bool)
    with tm.assert_produces_warning(UserWarning, match="can silently return incorrect"):
        expressions.where(cond, arr, arr)


def test_blocked_numexpr_warns_on_eval(blocked_numexpr):
    # GH#66956 eval falls back to the python engine and reports why
    a, b = 1, 2  # noqa: F841
    with tm.assert_produces_warning(UserWarning, match="can silently return incorrect"):
        res = pd.eval("a + b")
    assert res == 3


def test_blocked_numexpr_engine_raises_naming_the_version(blocked_numexpr):
    # GH#66956 this path raises instead of warning, so the message has to say
    #  which numexpr is refused and why
    a, b = 1, 2  # noqa: F841
    msg = f"numexpr {blocked_numexpr} is installed, but can silently return"
    with pytest.raises(ImportError, match=msg):
        pd.eval("a + b", engine="numexpr")


def test_blocked_numexpr_no_warning_if_disabled(blocked_numexpr):
    # GH#66956 a user who has turned numexpr off would not have used it either
    a, b = 1, 2  # noqa: F841
    with pd.option_context("compute.use_numexpr", False):
        with tm.assert_produces_warning(None):
            assert pd.eval("a + b") == 3


@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize("parser", expr.PARSERS)
def test_invalid_numexpr_version(engine, parser):
    if engine == "numexpr" and not check.NUMEXPR_INSTALLED:
        pytest.skip("numexpr not installed or an unsupported version")
    a, b = 1, 2  # noqa: F841
    res = pd.eval("a + b", engine=engine, parser=parser)
    assert res == 3
