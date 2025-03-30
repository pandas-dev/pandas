from __future__ import annotations

import numpy as np

from dask._dispatch import get_collection_type

try:
    import sparse

    sparse_installed = True
except ImportError:
    sparse_installed = False


try:
    import scipy.sparse as sp

    scipy_installed = True
except ImportError:
    scipy_installed = False


def create_array_collection(expr):
    # This is hacky and an abstraction leak, but utilizing get_collection_type
    # to infer that we want to create an array is the only way that is guaranteed
    # to be a general solution.
    # We can get rid of this when we have an Array expression
    import dask.array as da
    from dask.highlevelgraph import HighLevelGraph
    from dask.layers import Blockwise

    if da._array_expr_enabled():
        from dask.array._array_expr._expr import ArrayExpr

        if isinstance(expr, ArrayExpr):
            from dask.array._array_expr._collection import Array

            return Array(expr)

    result = expr.optimize()
    dsk = result.__dask_graph__()
    name = result._name
    meta = result._meta
    divisions = result.divisions

    chunks = ((np.nan,) * (len(divisions) - 1),) + tuple((d,) for d in meta.shape[1:])
    if len(chunks) > 1:
        if isinstance(dsk, HighLevelGraph):
            layer = dsk.layers[name]
        else:
            # dask-expr provides a dict only
            layer = dsk

        new_keys = []
        if isinstance(layer, Blockwise):
            layer.new_axes["j"] = chunks[1][0]
            layer.output_indices = layer.output_indices + ("j",)
        else:
            from dask._task_spec import Alias, Task

            suffix = (0,) * (len(chunks) - 1)
            for i in range(len(chunks[0])):
                task = layer.get((name, i))
                new_key = (name, i) + suffix
                if isinstance(task, Task):
                    task = Alias(new_key, task.key)
                layer[new_key] = task
                new_keys.append(new_key)
    else:
        new_keys = [(name, 0)]
    if da._array_expr_enabled():
        from dask.array._array_expr._collection import from_graph

        return from_graph(dsk, meta, chunks, set(new_keys), name)
    else:
        return da.Array(dsk, name=name, chunks=chunks, dtype=meta.dtype)


@get_collection_type.register(np.ndarray)
def get_collection_type_array(_):
    import dask.array as da

    if da._array_expr_enabled():
        from dask.array._array_expr._collection import Array

        return Array
    return create_array_collection


if sparse_installed:  # type: ignore[misc]

    @get_collection_type.register(sparse.COO)
    def get_collection_type_array(_):
        return create_array_collection


if scipy_installed:  # type: ignore[misc]

    @get_collection_type.register(sp.csr_matrix)
    def get_collection_type_array(_):
        return create_array_collection


if scipy_installed and hasattr(sp, "sparray"):  # type: ignore[misc]

    @get_collection_type.register(sp.csr_array)
    def get_collection_type_array(_):
        return create_array_collection


@get_collection_type.register(object)
def get_collection_type_object(_):

    return create_scalar_collection


def create_scalar_collection(expr):
    from dask.array._array_expr._expr import ArrayExpr

    if isinstance(expr, ArrayExpr):
        # This mirrors the legacy implementation in dask.array
        from dask.array._array_expr._collection import Array

        return Array(expr)

    from dask.dataframe.dask_expr._collection import Scalar

    return Scalar(expr)
