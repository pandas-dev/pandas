from __future__ import annotations

import contextlib
import warnings

import numpy as np

from dask.array._array_expr._expr import ArrayExpr
from dask.array.utils import meta_from_array
from dask.utils import has_keyword, is_arraylike


def compute_meta(func, _dtype, *args, **kwargs):
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)

        args_meta = [
            (
                x._meta
                if isinstance(x, ArrayExpr)
                else meta_from_array(x) if is_arraylike(x) else x
            )
            for x in args
        ]
        kwargs_meta = {
            k: (
                v._meta
                if isinstance(v, ArrayExpr)
                else meta_from_array(v) if is_arraylike(v) else v
            )
            for k, v in kwargs.items()
        }

        # todo: look for alternative to this, causes issues when using map_blocks()
        # with np.vectorize, such as dask.array.routines._isnonzero_vec().
        if isinstance(func, np.vectorize):
            meta = func(*args_meta)
        else:
            try:
                # some reduction functions need to know they are computing meta
                if has_keyword(func, "computing_meta"):
                    kwargs_meta["computing_meta"] = True
                meta = func(*args_meta, **kwargs_meta)
            except TypeError as e:
                if any(
                    s in str(e)
                    for s in [
                        "unexpected keyword argument",
                        "is an invalid keyword for",
                        "Did not understand the following kwargs",
                    ]
                ):
                    raise
                else:
                    return None
            except ValueError as e:
                # min/max functions have no identity, just use the same input type when there's only one
                if len(
                    args_meta
                ) == 1 and "zero-size array to reduction operation" in str(e):
                    meta = args_meta[0]
                else:
                    return None
            except Exception:
                return None

        if _dtype and getattr(meta, "dtype", None) != _dtype:
            with contextlib.suppress(AttributeError):
                meta = meta.astype(_dtype)

        if np.isscalar(meta):
            meta = np.array(meta)

        return meta
