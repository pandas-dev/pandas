from __future__ import annotations


def __getattr__(value):
    try:
        import dask_ml
    except ImportError as e:
        msg = (
            "Dask-ML is not installed.\n\n"
            "Please either conda or pip install dask-ml:\n\n"
            "  conda install dask-ml                      # either conda install\n"
            "  python -m pip install dask-ml --upgrade    # or pip install"
        )
        raise ImportError(msg) from e
    return getattr(dask_ml, value)
