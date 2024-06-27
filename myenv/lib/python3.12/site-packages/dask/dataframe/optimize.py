""" Dataframe optimizations """
from __future__ import annotations

import operator

import numpy as np

from dask import config, core
from dask.blockwise import Blockwise, fuse_roots, optimize_blockwise
from dask.highlevelgraph import HighLevelGraph
from dask.optimization import cull, fuse
from dask.utils import ensure_dict


def optimize(dsk, keys, **kwargs):
    if not isinstance(keys, (list, set)):
        keys = [keys]
    keys = list(core.flatten(keys))

    if not isinstance(dsk, HighLevelGraph):
        dsk = HighLevelGraph.from_collections(id(dsk), dsk, dependencies=())
    else:
        # Perform Blockwise optimizations for HLG input
        dsk = optimize_dataframe_getitem(dsk, keys=keys)
        dsk = optimize_blockwise(dsk, keys=keys)
        dsk = fuse_roots(dsk, keys=keys)
    dsk = dsk.cull(set(keys))

    # Do not perform low-level fusion unless the user has
    # specified True explicitly. The configuration will
    # be None by default.
    if not config.get("optimization.fuse.active"):
        return dsk

    dependencies = dsk.get_all_dependencies()
    dsk = ensure_dict(dsk)

    fuse_subgraphs = config.get("optimization.fuse.subgraphs")
    if fuse_subgraphs is None:
        fuse_subgraphs = True
    dsk, _ = fuse(
        dsk,
        keys,
        dependencies=dependencies,
        fuse_subgraphs=fuse_subgraphs,
    )
    dsk, _ = cull(dsk, keys)
    return dsk


def optimize_dataframe_getitem(dsk, keys):
    # This optimization looks for all `DataFrameIOLayer` instances,
    # and calls `project_columns` on any IO layers that precede
    # a (qualified) `getitem` operation.

    from dask.layers import DataFrameIOLayer

    # Construct a list containing the names of all
    # DataFrameIOLayer layers in the graph
    io_layers = [k for k, v in dsk.layers.items() if isinstance(v, DataFrameIOLayer)]

    def _is_selection(layer):
        # Utility to check if layer is a getitem selection

        # Must be Blockwise
        if not isinstance(layer, Blockwise):
            return False

        # Callable must be `getitem`
        if layer.dsk[layer.output][0] != operator.getitem:
            return False

        return True

    def _kind(layer):
        # Utility to check type of getitem selection

        # Selection is second indice
        key, ind = layer.indices[1]
        if ind is None:
            if isinstance(key, (tuple, str, list, np.ndarray)) or np.isscalar(key):
                return "column-selection"
        return "row-selection"

    # Loop over each DataFrameIOLayer layer
    layers = dsk.layers.copy()
    dependencies = dsk.dependencies.copy()
    for io_layer_name in io_layers:
        columns = set()

        # Bail on the optimization if the  IO layer is in the
        # requested keys, because we cannot change the name
        # anymore. These keys are structured like
        # [('getitem-<token>', 0), ...] so we check for the
        # first item of the tuple.
        # (See https://github.com/dask/dask/issues/5893)
        if any(
            layers[io_layer_name].name == x[0] for x in keys if isinstance(x, tuple)
        ):
            continue

        # Inspect dependents of the current IO layer
        deps = dsk.dependents[io_layer_name]

        # This optimization currently supports two variations
        # of `io_layer_name` dependents (`deps`):
        #
        #  - CASE A
        #    - 1 Dependent: A column-based getitem layer
        #    - This corresponds to the simple case that a
        #      column selection directly follows the IO layer
        #
        #  - CASE B
        #    - >1 Dependents: An arbitrary number of column-
        #      based getitem layers, and a single row selection
        #    - Usually corresponds to a filter operation that
        #      may (or may not) precede a column selection
        #    - This pattern is typical in Dask-SQL SELECT
        #      queries that include a WHERE statement

        # Bail if dependent layer type(s) do not agree with
        # case A or case B

        if not all(_is_selection(dsk.layers[k]) for k in deps) or {
            _kind(dsk.layers[k]) for k in deps
        } not in (
            {"column-selection"},
            {"column-selection", "row-selection"},
        ):
            continue

        # Split the column- and row-selection layers.
        # For case A, we will simply use information
        # from col_select_layers to perform column
        # projection in the root IO layer. For case B,
        # these layers are not the final column selection
        # (they are only part of a filtering operation).
        row_select_layers = {k for k in deps if _kind(dsk.layers[k]) == "row-selection"}
        col_select_layers = deps - row_select_layers

        # Can only handle single row-selection dependent (case B)
        if len(row_select_layers) > 1:
            continue

        # Define utility to walk the dependency graph
        # and check that the graph terminates with
        # the `success` key
        def _walk_deps(dependents, key, success):
            if key == success:
                return True
            deps = dependents[key]
            if deps:  # noqa: B023
                return all(
                    _walk_deps(dependents, dep, success) for dep in deps  # noqa: B023
                )
            else:
                return False

        # If this is not case A, we now need to check if
        # we are dealing with case B (and should bail on
        # the optimization if not).
        #
        # For case B, we should be able to start at
        # col_select_layer, and follow the graph to
        # row_select_layer. The subgraph between these
        # layers must depend ONLY on col_select_layer,
        # and be consumed ONLY by row_select_layer.
        # If these conditions are met, then a column-
        # selection layer directly following
        # row_select_layer can be used for projection.
        if row_select_layers:
            # Before walking the subgraph, check that there
            # is a column-selection layer directly following
            # row_select_layer. Otherwise, we can bail now.
            row_select_layer = row_select_layers.pop()
            if len(dsk.dependents[row_select_layer]) != 1:
                continue  # Too many/few row_select_layer dependents
            _layer = dsk.layers[list(dsk.dependents[row_select_layer])[0]]
            if _is_selection(_layer) and _kind(_layer) == "column-selection":
                # Include this column selection in our list of columns
                selection = _layer.indices[1][0]
                columns |= set(
                    selection if isinstance(selection, list) else [selection]
                )
            else:
                continue  # row_select_layer dependent not column selection

            # Walk the subgraph to check that all dependencies flow
            # from col_select_layers to the same col_select_layer
            if not all(
                _walk_deps(dsk.dependents, col_select_layer, col_select_layer)
                for col_select_layer in col_select_layers
            ):
                continue

        # Update columns with selections in col_select_layers
        for col_select_layer in col_select_layers:
            selection = dsk.layers[col_select_layer].indices[1][0]
            columns |= set(selection if isinstance(selection, list) else [selection])

        # If we got here, column projection is supported.
        # Add deps to update_blocks
        update_blocks = {dep: dsk.layers[dep] for dep in deps}

        # Project columns and update blocks
        old = layers[io_layer_name]
        new = old.project_columns(columns)
        if new.name != old.name:
            assert len(update_blocks)
            for block_key, block in update_blocks.items():
                # (('read-parquet-old', (.,)), ( ... )) ->
                # (('read-parquet-new', (.,)), ( ... ))
                new_indices = ((new.name, block.indices[0][1]), block.indices[1])
                numblocks = {new.name: block.numblocks[old.name]}
                new_block = Blockwise(
                    block.output,
                    block.output_indices,
                    block.dsk,
                    new_indices,
                    numblocks,
                    block.concatenate,
                    block.new_axes,
                )
                layers[block_key] = new_block
                dependencies[block_key] = {new.name}
            dependencies[new.name] = dependencies.pop(io_layer_name)

        layers[new.name] = new
        if new.name != old.name:
            del layers[old.name]

    new_hlg = HighLevelGraph(layers, dependencies)
    return new_hlg
