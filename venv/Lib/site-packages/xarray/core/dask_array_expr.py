"""Optional xarray integration for ``dask_array`` expression-backed arrays.

xarray stays expression-free. ``dask_array`` owns the lazy array expression
objects, while xarray owns Dataset/DataArray semantics such as coordinates,
indexes, attrs, template validation, and rebuilding final xarray objects.

Most xarray operations only need the normal chunk-manager methods. The special
case here is ``xarray.map_blocks``: it can return multiple output variables from
one user function call per input block. The helper below converts xarray's block
metadata into a private ``dask_array`` multi-output map expression. Each output
variable is still a normal ``dask_array.Array`` child expression, so Dask can
group the children with the composite-collection protocol and ``dask_array`` can
optimize, cull, persist, and compute those arrays.

If a Dataset mixes ``dask_array.Array`` with legacy ``dask.array.Array`` objects,
this path raises before constructing a graph. The generic Dataset expression
protocol instead returns ``None`` for mixed datasets so Dask can use xarray's
existing HighLevelGraph fallback.
"""

from __future__ import annotations

import itertools
from collections.abc import Callable, Hashable, Mapping, Sequence
from importlib import import_module
from typing import Any

from xarray.core.coordinates import Coordinates
from xarray.core.dataarray import DataArray
from xarray.core.dataset import Dataset
from xarray.core.utils import is_dask_collection


def is_dask_array_expr_array(data: Any) -> bool:
    try:
        dask_array = import_module("dask_array")
    except ImportError:
        return False

    array_type = getattr(dask_array, "Array", None)
    return array_type is not None and isinstance(data, array_type)


def collect_dask_array_expr_chunked_data(
    xarray_objs: Sequence[Dataset],
) -> tuple[bool, list[Any]]:
    chunked_data = [
        variable.data
        for arg in xarray_objs
        for variable in arg.variables.values()
        if is_dask_collection(variable.data)
    ]
    has_dask_array_expr = any(is_dask_array_expr_array(data) for data in chunked_data)
    has_other_chunked = any(not is_dask_array_expr_array(data) for data in chunked_data)

    if has_dask_array_expr and has_other_chunked:
        raise TypeError(
            "xarray.map_blocks cannot mix dask_array.Array with legacy or other "
            "Dask-backed arrays. Convert inputs to one array backend first."
        )

    return has_dask_array_expr, chunked_data


def _execute_map_blocks_multi_output(block_spec: Mapping[str, Any], *blocks: Any):
    args = []
    for arg_spec in block_spec["args"]:
        if arg_spec[0] == "literal":
            args.append(arg_spec[1])
            continue

        _, data_vars, coords, attrs = arg_spec

        def build_variables(variable_specs):
            variables = []
            for name, dims, data, var_attrs in variable_specs:
                if data[0] == "block":
                    data = blocks[data[1]]
                else:
                    data = data[1]
                variables.append((name, (dims, data, var_attrs)))
            return dict(variables)

        args.append(Dataset(build_variables(data_vars), build_variables(coords), attrs))

    result = block_spec["wrapper"](
        block_spec["func"],
        args,
        block_spec["kwargs"],
        block_spec["arg_is_array"],
        block_spec["expected"],
        block_spec["expected_indexes"],
    )
    return {
        output_index: result[variable_name]
        for output_index, variable_name in enumerate(block_spec["output_variables"])
    }


def map_blocks_with_dask_array_expr(
    *,
    func: Callable[..., Any],
    npargs: Sequence[Any],
    kwargs: Mapping[str, Any],
    is_xarray: Sequence[bool],
    is_array: Sequence[bool],
    input_chunks: Mapping[Hashable, tuple[int, ...]],
    output_chunks: Mapping[Hashable, tuple[int, ...]],
    coordinates: Coordinates,
    template: Dataset,
    result_is_array: bool,
    template_name: Hashable | None,
    token: str,
    ichunk: Mapping[Hashable, range],
    input_chunk_bounds: Mapping[Hashable, Any],
    output_chunk_bounds: Mapping[Hashable, Any],
    computed_variables: Sequence[Hashable],
    new_indexes: set[Hashable],
    modified_indexes: set[Hashable],
    chunked_data: Sequence[Any],
    wrapper: Callable[..., Any],
    get_chunk_slicer: Callable[[Hashable, Mapping[Any, Any], Mapping[Any, Any]], slice],
    dataset_to_dataarray: Callable[[Dataset], DataArray],
) -> DataArray | Dataset:
    missing_chunked_dims = {
        dim
        for dim, chunks in input_chunks.items()
        if len(chunks) > 1 and dim not in output_chunks
    }
    if missing_chunked_dims:
        raise NotImplementedError(
            "dask_array-backed xarray.map_blocks does not yet support "
            "dropping multi-chunk dimensions. Rechunk these dimensions to "
            f"one chunk first: {sorted(missing_chunked_dims, key=repr)!r}."
        )

    from xarray.namedarray.parallelcompat import get_chunked_array_type

    chunkmanager = get_chunked_array_type(*chunked_data)
    map_blocks_multi_output = getattr(chunkmanager, "map_blocks_multi_output", None)
    if map_blocks_multi_output is None:
        raise NotImplementedError(
            "The dask_array chunk manager does not support map_blocks_multi_output."
        )

    input_exprs: list[Any] = []
    input_indices: list[Any] = []
    arg_templates: list[Any] = []
    for isxr, arg in zip(is_xarray, npargs, strict=True):
        if not isxr:
            if is_dask_collection(arg):
                raise TypeError(
                    "dask_array-backed xarray.map_blocks only supports Dask "
                    "collections inside xarray arguments."
                )
            arg_templates.append(("literal", arg))
            continue

        variable_templates: list[Any] = []
        for name, variable in arg.variables.items():
            is_coord = name in arg._coord_names
            if is_dask_collection(variable.data):
                data: Any = variable.data
                input_exprs.append(data.expr)
                input_indices.append(variable.dims)
                variable_templates.append(
                    (
                        "chunked",
                        name,
                        variable.dims,
                        variable.attrs,
                        len(input_exprs) - 1,
                        is_coord,
                        None,
                    )
                )
            else:
                variable_templates.append(
                    (
                        "static",
                        name,
                        variable.dims,
                        variable.attrs,
                        None,
                        is_coord,
                        variable,
                    )
                )
        arg_templates.append(("xarray", arg.attrs, variable_templates))

    def build_block_specs() -> dict[tuple[Any, ...], dict[str, Any]]:
        specs: dict[tuple[Any, ...], dict[str, Any]] = {}
        for chunk_tuple in itertools.product(*ichunk.values()):
            chunk_index = dict(zip(ichunk.keys(), chunk_tuple, strict=True))
            arg_specs: list[Any] = []
            for arg_template in arg_templates:
                if arg_template[0] == "literal":
                    arg_specs.append(arg_template)
                    continue

                _, dataset_attrs, variable_templates = arg_template
                data_vars: list[Any] = []
                coords: list[Any] = []
                for (
                    kind,
                    name,
                    dims,
                    variable_attrs,
                    input_position,
                    is_coord,
                    variable,
                ) in variable_templates:
                    if kind == "chunked":
                        data = ("block", input_position)
                    else:
                        assert variable is not None
                        subsetter = {
                            dim: get_chunk_slicer(dim, chunk_index, input_chunk_bounds)
                            for dim in variable.dims
                        }
                        data = ("static", variable.isel(subsetter)._data)

                    target = coords if is_coord else data_vars
                    target.append((name, dims, data, variable_attrs))

                arg_specs.append(("xarray", data_vars, coords, dataset_attrs))

            indexes = {
                dim: coordinates.xindexes[dim][
                    get_chunk_slicer(dim, chunk_index, output_chunk_bounds)
                ]
                for dim in (new_indexes | modified_indexes)
            }
            expected = {
                "shapes": {
                    k: output_chunks[k][v]
                    for k, v in chunk_index.items()
                    if k in output_chunks
                },
                "data_vars": set(template.data_vars.keys()),
                "coords": set(template.coords.keys()),
            }
            specs[chunk_tuple] = {
                "wrapper": wrapper,
                "func": func,
                "args": arg_specs,
                "kwargs": kwargs,
                "arg_is_array": is_array,
                "expected": expected,
                "expected_indexes": indexes,
                "output_variables": tuple(computed_variables),
            }
        return specs

    outputs = []
    for output_index, name in enumerate(computed_variables):
        variable = template.variables[name]
        var_chunks = []
        for dim in variable.dims:
            if dim in output_chunks:
                var_chunks.append(output_chunks[dim])
            elif dim in template.dims:
                var_chunks.append((template.sizes[dim],))

        outputs.append(
            {
                "key": output_index,
                "indices": variable.dims,
                "chunks": tuple(var_chunks),
                "dtype": variable.dtype,
                "name": f"{name}-{token}",
            }
        )

    mapped_arrays = map_blocks_multi_output(
        _execute_map_blocks_multi_output,
        input_exprs,
        input_indices,
        tuple(input_chunks),
        build_block_specs(),
        outputs,
        token=token,
    )

    result = Dataset(coords=coordinates, attrs=template.attrs)
    for index in result._indexes:
        result[index].attrs = template[index].attrs
        result[index].encoding = template[index].encoding

    for name, data in zip(computed_variables, mapped_arrays, strict=True):
        result[name] = (template[name].dims, data, template[name].attrs)
        result[name].encoding = template[name].encoding

    result = result.set_coords(template._coord_names)

    if result_is_array:
        da = dataset_to_dataarray(result)
        da.name = template_name
        return da
    return result
