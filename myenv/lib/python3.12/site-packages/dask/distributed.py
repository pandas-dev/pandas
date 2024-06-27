# flake8: noqa

from __future__ import annotations

_import_error_message = (
    "dask.distributed is not installed.\n\n"
    "Please either conda or pip install distributed:\n\n"
    "  conda install dask distributed             # either conda install\n"
    '  python -m pip install "dask[distributed]" --upgrade    # or pip install'
)

try:
    from distributed import *
except ImportError as e:
    if e.msg == "No module named 'distributed'":
        raise ImportError(_import_error_message) from e
    else:
        raise


def __getattr__(value):
    try:
        import distributed
    except ImportError as e:
        raise ImportError(_import_error_message) from e
    return getattr(distributed, value)
