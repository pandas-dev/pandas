from __future__ import annotations

import pytest

try:
    from dask.dataframe.utils import pyarrow_strings_enabled

    convert_string = pyarrow_strings_enabled()
except (ImportError, RuntimeError):
    convert_string = False

skip_with_pyarrow_strings = pytest.mark.skipif(
    convert_string,
    reason="No need to run with pyarrow strings",
)

xfail_with_pyarrow_strings = pytest.mark.xfail(
    convert_string,
    reason="Known failure with pyarrow strings",
)


def pytest_collection_modifyitems(config, items):
    for item in items:
        if "skip_with_pyarrow_strings" in item.keywords:
            item.add_marker(skip_with_pyarrow_strings)
        if "xfail_with_pyarrow_strings" in item.keywords:
            item.add_marker(xfail_with_pyarrow_strings)


@pytest.fixture(params=["disk", "tasks"])
def shuffle_method(request):
    import dask

    with dask.config.set({"dataframe.shuffle.method": request.param}):
        yield request.param


@pytest.fixture
def xfail(request: pytest.FixtureRequest):
    """
    This fixture returns a function that, when called inside the test,
    XFAILs the currently running test.

    You should prefer using `@pytest.mark.xfail` when possible. However, when a failure
    condition is too complicated to test statically (e.g. when it depends on the
    intersection of multiple pytest parameters), this is preferable to just calling
    `pytest.xfail`. The difference is that `pytest.xfail` will immediately end the test,
    while this fixture allows the test to execute, so that it may result in a XPASS.

    xref https://github.com/pandas-dev/pandas/issues/38902

    Usage
    -----
    ::
        @pytest.mark.parametrize("param1", [1, 2, 3])
        @pytest.mark.parametrize("param2", [4, 5, 6])
        def test1(xfail):
            if param1 == 2 and param2 == 6:
                xfail("This combination of parameters is known to fail")
            # Test continues, including in the XFAIL'ed use case

    Parameters
    ----------
    reason : str
        Reason for the expected failure.
    strict: bool, optional
        If True, the test will be marked as failed if it passes (XPASS).
        If False, the test will be marked as passed both if it fails (XFAIL)
        or passes (XPASS).
        Default: ``xfail_strict`` value in ``pyproject.toml``, or False if absent.
    """

    def _(reason: str, strict: bool | None = None) -> None:
        if strict is not None:
            marker = pytest.mark.xfail(reason=reason, strict=strict)
        else:
            marker = pytest.mark.xfail(reason=reason)
        request.node.add_marker(marker)

    return _
