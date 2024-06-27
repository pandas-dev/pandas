from __future__ import annotations


def test_basic():
    try:
        import dask_ml  # noqa: F401
    except ImportError:
        try:
            from dask.ml.model_selection import GridSearchCV  # noqa: F401
        except ImportError as e:
            assert "conda install dask-ml" in str(e)
        else:
            assert False
    else:
        from dask.ml.model_selection import GridSearchCV  # noqa: F401
