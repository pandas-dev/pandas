from __future__ import annotations

import pytest

import dask


@pytest.mark.gpu
@pytest.mark.parametrize("backend", ["pandas", "cudf"])
def test_CreationDispatch_error_informative_message(backend):
    # Check that an informative error is emitted when a backend dispatch
    # method fails
    pytest.importorskip(backend)
    dd = pytest.importorskip("dask.dataframe")
    if dd._dask_expr_enabled():
        pytest.skip("not properly supported yet")
    data = {"a": [1, 2, 3, 4], "B": [10, 11, 12, 13]}
    with dask.config.set({"dataframe.backend": backend}):
        with pytest.raises(TypeError) as excinfo:
            dd.from_dict(data, npartitions=2, unsupported_kwarg=True)

        msg = str(excinfo.value)
        assert "error occurred while calling the from_dict method" in msg
        assert backend in msg
