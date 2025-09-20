from __future__ import annotations

import functools
import operator
import pickle

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose, assert_array_equal

import xarray as xr
from xarray.computation.apply_ufunc import (
    _UFuncSignature,
    apply_ufunc,
    broadcast_compat_data,
    collect_dict_values,
    join_dict_keys,
    ordered_set_intersection,
    ordered_set_union,
    unified_dim_sizes,
)
from xarray.core.utils import result_name
from xarray.structure.alignment import broadcast
from xarray.tests import (
    has_dask,
    raise_if_dask_computes,
    requires_cftime,
    requires_dask,
)


def assert_identical(a, b):
    """A version of this function which accepts numpy arrays"""
    __tracebackhide__ = True
    from xarray.testing import assert_identical as assert_identical_

    if hasattr(a, "identical"):
        assert_identical_(a, b)
    else:
        assert_array_equal(a, b)


def test_signature_properties() -> None:
    sig = _UFuncSignature([["x"], ["x", "y"]], [["z"]])
    assert sig.input_core_dims == (("x",), ("x", "y"))
    assert sig.output_core_dims == (("z",),)
    assert sig.all_input_core_dims == frozenset(["x", "y"])
    assert sig.all_output_core_dims == frozenset(["z"])
    assert sig.num_inputs == 2
    assert sig.num_outputs == 1
    assert str(sig) == "(x),(x,y)->(z)"
    assert sig.to_gufunc_string() == "(dim0),(dim0,dim1)->(dim2)"
    assert (
        sig.to_gufunc_string(exclude_dims=set("x")) == "(dim0_0),(dim0_1,dim1)->(dim2)"
    )
    # dimension names matter
    assert _UFuncSignature([["x"]]) != _UFuncSignature([["y"]])


def test_result_name() -> None:
    class Named:
        def __init__(self, name=None):
            self.name = name

    assert result_name([1, 2]) is None
    assert result_name([Named()]) is None
    assert result_name([Named("foo"), 2]) == "foo"
    assert result_name([Named("foo"), Named("bar")]) is None
    assert result_name([Named("foo"), Named()]) is None


def test_ordered_set_union() -> None:
    assert list(ordered_set_union([[1, 2]])) == [1, 2]
    assert list(ordered_set_union([[1, 2], [2, 1]])) == [1, 2]
    assert list(ordered_set_union([[0], [1, 2], [1, 3]])) == [0, 1, 2, 3]


def test_ordered_set_intersection() -> None:
    assert list(ordered_set_intersection([[1, 2]])) == [1, 2]
    assert list(ordered_set_intersection([[1, 2], [2, 1]])) == [1, 2]
    assert list(ordered_set_intersection([[1, 2], [1, 3]])) == [1]
    assert list(ordered_set_intersection([[1, 2], [2]])) == [2]


def test_join_dict_keys() -> None:
    dicts = [dict.fromkeys(keys) for keys in [["x", "y"], ["y", "z"]]]
    assert list(join_dict_keys(dicts, "left")) == ["x", "y"]
    assert list(join_dict_keys(dicts, "right")) == ["y", "z"]
    assert list(join_dict_keys(dicts, "inner")) == ["y"]
    assert list(join_dict_keys(dicts, "outer")) == ["x", "y", "z"]
    with pytest.raises(ValueError):
        join_dict_keys(dicts, "exact")
    with pytest.raises(KeyError):
        join_dict_keys(dicts, "foobar")


def test_collect_dict_values() -> None:
    dicts = [{"x": 1, "y": 2, "z": 3}, {"z": 4}, 5]
    expected = [[1, 0, 5], [2, 0, 5], [3, 4, 5]]
    collected = collect_dict_values(dicts, ["x", "y", "z"], fill_value=0)
    assert collected == expected


def identity(x):
    return x


def test_apply_identity() -> None:
    array = np.arange(10)
    variable = xr.Variable("x", array)
    data_array = xr.DataArray(variable, [("x", -array)])
    dataset = xr.Dataset({"y": variable}, {"x": -array})

    apply_identity = functools.partial(apply_ufunc, identity)

    assert_identical(array, apply_identity(array))
    assert_identical(variable, apply_identity(variable))
    assert_identical(data_array, apply_identity(data_array))
    assert_identical(data_array, apply_identity(data_array.groupby("x")))
    assert_identical(dataset, apply_identity(dataset))
    assert_identical(dataset, apply_identity(dataset.groupby("x")))


def add(a, b):
    return apply_ufunc(operator.add, a, b)


def test_apply_two_inputs() -> None:
    array = np.array([1, 2, 3])
    variable = xr.Variable("x", array)
    data_array = xr.DataArray(variable, [("x", -array)])
    dataset = xr.Dataset({"y": variable}, {"x": -array})

    zero_array = np.zeros_like(array)
    zero_variable = xr.Variable("x", zero_array)
    zero_data_array = xr.DataArray(zero_variable, [("x", -array)])
    zero_dataset = xr.Dataset({"y": zero_variable}, {"x": -array})

    assert_identical(array, add(array, zero_array))
    assert_identical(array, add(zero_array, array))

    assert_identical(variable, add(variable, zero_array))
    assert_identical(variable, add(variable, zero_variable))
    assert_identical(variable, add(zero_array, variable))
    assert_identical(variable, add(zero_variable, variable))

    assert_identical(data_array, add(data_array, zero_array))
    assert_identical(data_array, add(data_array, zero_variable))
    assert_identical(data_array, add(data_array, zero_data_array))
    assert_identical(data_array, add(zero_array, data_array))
    assert_identical(data_array, add(zero_variable, data_array))
    assert_identical(data_array, add(zero_data_array, data_array))

    assert_identical(dataset, add(dataset, zero_array))
    assert_identical(dataset, add(dataset, zero_variable))
    assert_identical(dataset, add(dataset, zero_data_array))
    assert_identical(dataset, add(dataset, zero_dataset))
    assert_identical(dataset, add(zero_array, dataset))
    assert_identical(dataset, add(zero_variable, dataset))
    assert_identical(dataset, add(zero_data_array, dataset))
    assert_identical(dataset, add(zero_dataset, dataset))

    assert_identical(data_array, add(data_array.groupby("x"), zero_data_array))
    assert_identical(data_array, add(zero_data_array, data_array.groupby("x")))

    assert_identical(dataset, add(data_array.groupby("x"), zero_dataset))
    assert_identical(dataset, add(zero_dataset, data_array.groupby("x")))

    assert_identical(dataset, add(dataset.groupby("x"), zero_data_array))
    assert_identical(dataset, add(dataset.groupby("x"), zero_dataset))
    assert_identical(dataset, add(zero_data_array, dataset.groupby("x")))
    assert_identical(dataset, add(zero_dataset, dataset.groupby("x")))


def test_apply_1d_and_0d() -> None:
    array = np.array([1, 2, 3])
    variable = xr.Variable("x", array)
    data_array = xr.DataArray(variable, [("x", -array)])
    dataset = xr.Dataset({"y": variable}, {"x": -array})

    zero_array = 0
    zero_variable = xr.Variable((), zero_array)
    zero_data_array = xr.DataArray(zero_variable)
    zero_dataset = xr.Dataset({"y": zero_variable})

    assert_identical(array, add(array, zero_array))
    assert_identical(array, add(zero_array, array))

    assert_identical(variable, add(variable, zero_array))
    assert_identical(variable, add(variable, zero_variable))
    assert_identical(variable, add(zero_array, variable))
    assert_identical(variable, add(zero_variable, variable))

    assert_identical(data_array, add(data_array, zero_array))
    assert_identical(data_array, add(data_array, zero_variable))
    assert_identical(data_array, add(data_array, zero_data_array))
    assert_identical(data_array, add(zero_array, data_array))
    assert_identical(data_array, add(zero_variable, data_array))
    assert_identical(data_array, add(zero_data_array, data_array))

    assert_identical(dataset, add(dataset, zero_array))
    assert_identical(dataset, add(dataset, zero_variable))
    assert_identical(dataset, add(dataset, zero_data_array))
    assert_identical(dataset, add(dataset, zero_dataset))
    assert_identical(dataset, add(zero_array, dataset))
    assert_identical(dataset, add(zero_variable, dataset))
    assert_identical(dataset, add(zero_data_array, dataset))
    assert_identical(dataset, add(zero_dataset, dataset))

    assert_identical(data_array, add(data_array.groupby("x"), zero_data_array))
    assert_identical(data_array, add(zero_data_array, data_array.groupby("x")))

    assert_identical(dataset, add(data_array.groupby("x"), zero_dataset))
    assert_identical(dataset, add(zero_dataset, data_array.groupby("x")))

    assert_identical(dataset, add(dataset.groupby("x"), zero_data_array))
    assert_identical(dataset, add(dataset.groupby("x"), zero_dataset))
    assert_identical(dataset, add(zero_data_array, dataset.groupby("x")))
    assert_identical(dataset, add(zero_dataset, dataset.groupby("x")))


def test_apply_two_outputs() -> None:
    array = np.arange(5)
    variable = xr.Variable("x", array)
    data_array = xr.DataArray(variable, [("x", -array)])
    dataset = xr.Dataset({"y": variable}, {"x": -array})

    def twice(obj):
        def func(x):
            return (x, x)

        return apply_ufunc(func, obj, output_core_dims=[[], []])

    out0, out1 = twice(array)
    assert_identical(out0, array)
    assert_identical(out1, array)

    out0, out1 = twice(variable)
    assert_identical(out0, variable)
    assert_identical(out1, variable)

    out0, out1 = twice(data_array)
    assert_identical(out0, data_array)
    assert_identical(out1, data_array)

    out0, out1 = twice(dataset)
    assert_identical(out0, dataset)
    assert_identical(out1, dataset)

    out0, out1 = twice(data_array.groupby("x"))
    assert_identical(out0, data_array)
    assert_identical(out1, data_array)

    out0, out1 = twice(dataset.groupby("x"))
    assert_identical(out0, dataset)
    assert_identical(out1, dataset)


def test_apply_missing_dims() -> None:
    ## Single arg

    def add_one(a, core_dims, on_missing_core_dim):
        return apply_ufunc(
            lambda x: x + 1,
            a,
            input_core_dims=core_dims,
            output_core_dims=core_dims,
            on_missing_core_dim=on_missing_core_dim,
        )

    array = np.arange(6).reshape(2, 3)
    variable = xr.Variable(["x", "y"], array)
    variable_no_y = xr.Variable(["x", "z"], array)

    ds = xr.Dataset({"x_y": variable, "x_z": variable_no_y})

    # Check the standard stuff works OK
    assert_identical(
        add_one(ds[["x_y"]], core_dims=[["y"]], on_missing_core_dim="raise"),
        ds[["x_y"]] + 1,
    )

    # `raise` — should raise on a missing dim
    with pytest.raises(ValueError):
        add_one(ds, core_dims=[["y"]], on_missing_core_dim="raise")

    # `drop` — should drop the var with the missing dim
    assert_identical(
        add_one(ds, core_dims=[["y"]], on_missing_core_dim="drop"),
        (ds + 1).drop_vars("x_z"),
    )

    # `copy` — should not add one to the missing with `copy`
    copy_result = add_one(ds, core_dims=[["y"]], on_missing_core_dim="copy")
    assert_identical(copy_result["x_y"], (ds + 1)["x_y"])
    assert_identical(copy_result["x_z"], ds["x_z"])

    ## Multiple args

    def sum_add(a, b, core_dims, on_missing_core_dim):
        return apply_ufunc(
            lambda a, b, axis=None: a.sum(axis) + b.sum(axis),
            a,
            b,
            input_core_dims=core_dims,
            on_missing_core_dim=on_missing_core_dim,
        )

    # Check the standard stuff works OK
    assert_identical(
        sum_add(
            ds[["x_y"]],
            ds[["x_y"]],
            core_dims=[["x", "y"], ["x", "y"]],
            on_missing_core_dim="raise",
        ),
        ds[["x_y"]].sum() * 2,
    )

    # `raise` — should raise on a missing dim
    with pytest.raises(
        ValueError,
        match=r".*Missing core dims \{'y'\} from arg number 1 on a variable named `x_z`:\n.*<xarray.Variable \(x: 2, z: ",
    ):
        sum_add(
            ds[["x_z"]],
            ds[["x_z"]],
            core_dims=[["x", "y"], ["x", "z"]],
            on_missing_core_dim="raise",
        )

    # `raise` on a missing dim on a non-first arg
    with pytest.raises(
        ValueError,
        match=r".*Missing core dims \{'y'\} from arg number 2 on a variable named `x_z`:\n.*<xarray.Variable \(x: 2, z: ",
    ):
        sum_add(
            ds[["x_z"]],
            ds[["x_z"]],
            core_dims=[["x", "z"], ["x", "y"]],
            on_missing_core_dim="raise",
        )

    # `drop` — should drop the var with the missing dim
    assert_identical(
        sum_add(
            ds[["x_z"]],
            ds[["x_z"]],
            core_dims=[["x", "y"], ["x", "y"]],
            on_missing_core_dim="drop",
        ),
        ds[[]],
    )

    # `copy` — should drop the var with the missing dim
    assert_identical(
        sum_add(
            ds[["x_z"]],
            ds[["x_z"]],
            core_dims=[["x", "y"], ["x", "y"]],
            on_missing_core_dim="copy",
        ),
        ds[["x_z"]],
    )

    ## Multiple vars per arg
    assert_identical(
        sum_add(
            ds[["x_y", "x_y"]],
            ds[["x_y", "x_y"]],
            core_dims=[["x", "y"], ["x", "y"]],
            on_missing_core_dim="raise",
        ),
        ds[["x_y", "x_y"]].sum() * 2,
    )

    assert_identical(
        sum_add(
            ds,
            ds,
            core_dims=[["x", "y"], ["x", "y"]],
            on_missing_core_dim="drop",
        ),
        ds[["x_y"]].sum() * 2,
    )

    assert_identical(
        sum_add(
            # The first one has the wrong dims — using `z` for the `x_t` var
            ds.assign(x_t=ds["x_z"])[["x_y", "x_t"]],
            ds.assign(x_t=ds["x_y"])[["x_y", "x_t"]],
            core_dims=[["x", "y"], ["x", "y"]],
            on_missing_core_dim="drop",
        ),
        ds[["x_y"]].sum() * 2,
    )

    assert_identical(
        sum_add(
            ds.assign(x_t=ds["x_y"])[["x_y", "x_t"]],
            # The second one has the wrong dims — using `z` for the `x_t` var (seems
            # duplicative but this was a bug in the initial impl...)
            ds.assign(x_t=ds["x_z"])[["x_y", "x_t"]],
            core_dims=[["x", "y"], ["x", "y"]],
            on_missing_core_dim="drop",
        ),
        ds[["x_y"]].sum() * 2,
    )

    assert_identical(
        sum_add(
            ds.assign(x_t=ds["x_y"])[["x_y", "x_t"]],
            ds.assign(x_t=ds["x_z"])[["x_y", "x_t"]],
            core_dims=[["x", "y"], ["x", "y"]],
            on_missing_core_dim="copy",
        ),
        ds.drop_vars("x_z").assign(x_y=30, x_t=ds["x_y"]),
    )


@requires_dask
def test_apply_dask_parallelized_two_outputs() -> None:
    data_array = xr.DataArray([[0, 1, 2], [1, 2, 3]], dims=("x", "y"))

    def twice(obj):
        def func(x):
            return (x, x)

        return apply_ufunc(func, obj, output_core_dims=[[], []], dask="parallelized")

    out0, out1 = twice(data_array.chunk({"x": 1}))
    assert_identical(data_array, out0)
    assert_identical(data_array, out1)


def test_apply_input_core_dimension() -> None:
    def first_element(obj, dim):
        def func(x):
            return x[..., 0]

        return apply_ufunc(func, obj, input_core_dims=[[dim]])

    array = np.array([[1, 2], [3, 4]])
    variable = xr.Variable(["x", "y"], array)
    data_array = xr.DataArray(variable, {"x": ["a", "b"], "y": [-1, -2]})
    dataset = xr.Dataset({"data": data_array})

    expected_variable_x = xr.Variable(["y"], [1, 2])
    expected_data_array_x = xr.DataArray(expected_variable_x, {"y": [-1, -2]})
    expected_dataset_x = xr.Dataset({"data": expected_data_array_x})

    expected_variable_y = xr.Variable(["x"], [1, 3])
    expected_data_array_y = xr.DataArray(expected_variable_y, {"x": ["a", "b"]})
    expected_dataset_y = xr.Dataset({"data": expected_data_array_y})

    assert_identical(expected_variable_x, first_element(variable, "x"))
    assert_identical(expected_variable_y, first_element(variable, "y"))

    assert_identical(expected_data_array_x, first_element(data_array, "x"))
    assert_identical(expected_data_array_y, first_element(data_array, "y"))

    assert_identical(expected_dataset_x, first_element(dataset, "x"))
    assert_identical(expected_dataset_y, first_element(dataset, "y"))

    assert_identical(expected_data_array_x, first_element(data_array.groupby("y"), "x"))
    assert_identical(expected_dataset_x, first_element(dataset.groupby("y"), "x"))

    def multiply(*args):
        val = args[0]
        for arg in args[1:]:
            val = val * arg
        return val

    # regression test for GH:2341
    with pytest.raises(ValueError):
        apply_ufunc(
            multiply,
            data_array,
            data_array["y"].values,
            input_core_dims=[["y"]],
            output_core_dims=[["y"]],
        )
    expected = xr.DataArray(
        multiply(data_array, data_array["y"]), dims=["x", "y"], coords=data_array.coords
    )
    actual = apply_ufunc(
        multiply,
        data_array,
        data_array["y"].values,
        input_core_dims=[["y"], []],
        output_core_dims=[["y"]],
    )
    assert_identical(expected, actual)


def test_apply_output_core_dimension() -> None:
    def stack_negative(obj):
        def func(x):
            return np.stack([x, -x], axis=-1)

        result = apply_ufunc(func, obj, output_core_dims=[["sign"]])
        if isinstance(result, xr.Dataset | xr.DataArray):
            result.coords["sign"] = [1, -1]
        return result

    array = np.array([[1, 2], [3, 4]])
    variable = xr.Variable(["x", "y"], array)
    data_array = xr.DataArray(variable, {"x": ["a", "b"], "y": [-1, -2]})
    dataset = xr.Dataset({"data": data_array})

    stacked_array = np.array([[[1, -1], [2, -2]], [[3, -3], [4, -4]]])
    stacked_variable = xr.Variable(["x", "y", "sign"], stacked_array)
    stacked_coords = {"x": ["a", "b"], "y": [-1, -2], "sign": [1, -1]}
    stacked_data_array = xr.DataArray(stacked_variable, stacked_coords)
    stacked_dataset = xr.Dataset({"data": stacked_data_array})

    assert_identical(stacked_array, stack_negative(array))
    assert_identical(stacked_variable, stack_negative(variable))
    assert_identical(stacked_data_array, stack_negative(data_array))
    assert_identical(stacked_dataset, stack_negative(dataset))
    assert_identical(stacked_data_array, stack_negative(data_array.groupby("x")))
    assert_identical(stacked_dataset, stack_negative(dataset.groupby("x")))

    def original_and_stack_negative(obj):
        def func(x):
            return (x, np.stack([x, -x], axis=-1))

        result = apply_ufunc(func, obj, output_core_dims=[[], ["sign"]])
        if isinstance(result[1], xr.Dataset | xr.DataArray):
            result[1].coords["sign"] = [1, -1]
        return result

    out0, out1 = original_and_stack_negative(array)
    assert_identical(array, out0)
    assert_identical(stacked_array, out1)

    out0, out1 = original_and_stack_negative(variable)
    assert_identical(variable, out0)
    assert_identical(stacked_variable, out1)

    out0, out1 = original_and_stack_negative(data_array)
    assert_identical(data_array, out0)
    assert_identical(stacked_data_array, out1)

    out0, out1 = original_and_stack_negative(dataset)
    assert_identical(dataset, out0)
    assert_identical(stacked_dataset, out1)

    out0, out1 = original_and_stack_negative(data_array.groupby("x"))
    assert_identical(data_array, out0)
    assert_identical(stacked_data_array, out1)

    out0, out1 = original_and_stack_negative(dataset.groupby("x"))
    assert_identical(dataset, out0)
    assert_identical(stacked_dataset, out1)


def test_apply_exclude() -> None:
    def concatenate(objects, dim="x"):
        def func(*x):
            return np.concatenate(x, axis=-1)

        result = apply_ufunc(
            func,
            *objects,
            input_core_dims=[[dim]] * len(objects),
            output_core_dims=[[dim]],
            exclude_dims={dim},
        )
        if isinstance(result, xr.Dataset | xr.DataArray):
            # note: this will fail if dim is not a coordinate on any input
            new_coord = np.concatenate([obj.coords[dim] for obj in objects])
            result.coords[dim] = new_coord
        return result

    arrays = [np.array([1]), np.array([2, 3])]
    variables = [xr.Variable("x", a) for a in arrays]
    data_arrays = [
        xr.DataArray(v, {"x": c, "y": ("x", range(len(c)))})
        for v, c in zip(variables, [["a"], ["b", "c"]], strict=True)
    ]
    datasets = [xr.Dataset({"data": data_array}) for data_array in data_arrays]

    expected_array = np.array([1, 2, 3])
    expected_variable = xr.Variable("x", expected_array)
    expected_data_array = xr.DataArray(expected_variable, [("x", list("abc"))])
    expected_dataset = xr.Dataset({"data": expected_data_array})

    assert_identical(expected_array, concatenate(arrays))
    assert_identical(expected_variable, concatenate(variables))
    assert_identical(expected_data_array, concatenate(data_arrays))
    assert_identical(expected_dataset, concatenate(datasets))

    # must also be a core dimension
    with pytest.raises(ValueError):
        apply_ufunc(identity, variables[0], exclude_dims={"x"})


def test_apply_groupby_add() -> None:
    array = np.arange(5)
    variable = xr.Variable("x", array)
    coords = {"x": -array, "y": ("x", [0, 0, 1, 1, 2])}
    data_array = xr.DataArray(variable, coords, dims="x")
    dataset = xr.Dataset({"z": variable}, coords)

    other_variable = xr.Variable("y", [0, 10])
    other_data_array = xr.DataArray(other_variable, dims="y")
    other_dataset = xr.Dataset({"z": other_variable})

    expected_variable = xr.Variable("x", [0, 1, 12, 13, np.nan])
    expected_data_array = xr.DataArray(expected_variable, coords, dims="x")
    expected_dataset = xr.Dataset({"z": expected_variable}, coords)

    assert_identical(
        expected_data_array, add(data_array.groupby("y"), other_data_array)
    )
    assert_identical(expected_dataset, add(data_array.groupby("y"), other_dataset))
    assert_identical(expected_dataset, add(dataset.groupby("y"), other_data_array))
    assert_identical(expected_dataset, add(dataset.groupby("y"), other_dataset))

    # cannot be performed with xarray.Variable objects that share a dimension
    with pytest.raises(ValueError):
        add(data_array.groupby("y"), other_variable)

    # if they are all grouped the same way
    with pytest.raises(ValueError):
        add(data_array.groupby("y"), data_array[:4].groupby("y"))
    with pytest.raises(ValueError):
        add(data_array.groupby("y"), data_array[1:].groupby("y"))
    with pytest.raises(ValueError):
        add(data_array.groupby("y"), other_data_array.groupby("y"))
    with pytest.raises(ValueError):
        add(data_array.groupby("y"), data_array.groupby("x"))


@pytest.mark.filterwarnings("ignore:Duplicate dimension names present")
def test_unified_dim_sizes() -> None:
    assert unified_dim_sizes([xr.Variable((), 0)]) == {}
    assert unified_dim_sizes([xr.Variable("x", [1]), xr.Variable("x", [1])]) == {"x": 1}
    assert unified_dim_sizes([xr.Variable("x", [1]), xr.Variable("y", [1, 2])]) == {
        "x": 1,
        "y": 2,
    }
    assert unified_dim_sizes(
        [xr.Variable(("x", "z"), [[1]]), xr.Variable(("y", "z"), [[1, 2], [3, 4]])],
        exclude_dims={"z"},
    ) == {"x": 1, "y": 2}

    with pytest.raises(ValueError, match="broadcasting cannot handle"):
        with pytest.warns(UserWarning, match="Duplicate dimension names"):
            unified_dim_sizes([xr.Variable(("x", "x"), [[1]])])

    # mismatched lengths
    with pytest.raises(ValueError):
        unified_dim_sizes([xr.Variable("x", [1]), xr.Variable("x", [1, 2])])


def test_broadcast_compat_data_1d() -> None:
    data = np.arange(5)
    var = xr.Variable("x", data)

    assert_identical(data, broadcast_compat_data(var, ("x",), ()))
    assert_identical(data, broadcast_compat_data(var, (), ("x",)))
    assert_identical(data[:], broadcast_compat_data(var, ("w",), ("x",)))
    assert_identical(data[:, None], broadcast_compat_data(var, ("w", "x", "y"), ()))

    with pytest.raises(ValueError):
        broadcast_compat_data(var, ("x",), ("w",))

    with pytest.raises(ValueError):
        broadcast_compat_data(var, (), ())


def test_broadcast_compat_data_2d() -> None:
    data = np.arange(12).reshape(3, 4)
    var = xr.Variable(["x", "y"], data)

    assert_identical(data, broadcast_compat_data(var, ("x", "y"), ()))
    assert_identical(data, broadcast_compat_data(var, ("x",), ("y",)))
    assert_identical(data, broadcast_compat_data(var, (), ("x", "y")))
    assert_identical(data.T, broadcast_compat_data(var, ("y", "x"), ()))
    assert_identical(data.T, broadcast_compat_data(var, ("y",), ("x",)))
    assert_identical(data, broadcast_compat_data(var, ("w", "x"), ("y",)))
    assert_identical(data, broadcast_compat_data(var, ("w",), ("x", "y")))
    assert_identical(data.T, broadcast_compat_data(var, ("w",), ("y", "x")))
    assert_identical(
        data[:, :, None], broadcast_compat_data(var, ("w", "x", "y", "z"), ())
    )
    assert_identical(
        data[None, :, :].T, broadcast_compat_data(var, ("w", "y", "x", "z"), ())
    )


def test_keep_attrs() -> None:
    def add(a, b, keep_attrs):
        if keep_attrs:
            return apply_ufunc(operator.add, a, b, keep_attrs=keep_attrs)
        else:
            return apply_ufunc(operator.add, a, b)

    a = xr.DataArray([0, 1], [("x", [0, 1])])
    a.attrs["attr"] = "da"
    a["x"].attrs["attr"] = "da_coord"
    b = xr.DataArray([1, 2], [("x", [0, 1])])

    actual = add(a, b, keep_attrs=False)
    assert not actual.attrs
    actual = add(a, b, keep_attrs=True)
    assert_identical(actual.attrs, a.attrs)
    assert_identical(actual["x"].attrs, a["x"].attrs)

    actual = add(a.variable, b.variable, keep_attrs=False)
    assert not actual.attrs
    actual = add(a.variable, b.variable, keep_attrs=True)
    assert_identical(actual.attrs, a.attrs)

    ds_a = xr.Dataset({"x": [0, 1]})
    ds_a.attrs["attr"] = "ds"
    ds_a.x.attrs["attr"] = "da"
    ds_b = xr.Dataset({"x": [0, 1]})

    actual = add(ds_a, ds_b, keep_attrs=False)
    assert not actual.attrs
    actual = add(ds_a, ds_b, keep_attrs=True)
    assert_identical(actual.attrs, ds_a.attrs)
    assert_identical(actual.x.attrs, ds_a.x.attrs)


@pytest.mark.parametrize(
    ["strategy", "attrs", "expected", "error"],
    (
        pytest.param(
            None,
            [{"a": 1}, {"a": 2}, {"a": 3}],
            {},
            False,
            id="default",
        ),
        pytest.param(
            False,
            [{"a": 1}, {"a": 2}, {"a": 3}],
            {},
            False,
            id="False",
        ),
        pytest.param(
            True,
            [{"a": 1}, {"a": 2}, {"a": 3}],
            {"a": 1},
            False,
            id="True",
        ),
        pytest.param(
            "override",
            [{"a": 1}, {"a": 2}, {"a": 3}],
            {"a": 1},
            False,
            id="override",
        ),
        pytest.param(
            "drop",
            [{"a": 1}, {"a": 2}, {"a": 3}],
            {},
            False,
            id="drop",
        ),
        pytest.param(
            "drop_conflicts",
            [{"a": 1, "b": 2}, {"b": 1, "c": 3}, {"c": 3, "d": 4}],
            {"a": 1, "c": 3, "d": 4},
            False,
            id="drop_conflicts",
        ),
        pytest.param(
            "no_conflicts",
            [{"a": 1}, {"b": 2}, {"b": 3}],
            None,
            True,
            id="no_conflicts",
        ),
    ),
)
def test_keep_attrs_strategies_variable(strategy, attrs, expected, error) -> None:
    a = xr.Variable("x", [0, 1], attrs=attrs[0])
    b = xr.Variable("x", [0, 1], attrs=attrs[1])
    c = xr.Variable("x", [0, 1], attrs=attrs[2])

    if error:
        with pytest.raises(xr.MergeError):
            apply_ufunc(lambda *args: sum(args), a, b, c, keep_attrs=strategy)
    else:
        expected = xr.Variable("x", [0, 3], attrs=expected)
        actual = apply_ufunc(lambda *args: sum(args), a, b, c, keep_attrs=strategy)

        assert_identical(actual, expected)


@pytest.mark.parametrize(
    ["strategy", "attrs", "expected", "error"],
    (
        pytest.param(
            None,
            [{"a": 1}, {"a": 2}, {"a": 3}],
            {},
            False,
            id="default",
        ),
        pytest.param(
            False,
            [{"a": 1}, {"a": 2}, {"a": 3}],
            {},
            False,
            id="False",
        ),
        pytest.param(
            True,
            [{"a": 1}, {"a": 2}, {"a": 3}],
            {"a": 1},
            False,
            id="True",
        ),
        pytest.param(
            "override",
            [{"a": 1}, {"a": 2}, {"a": 3}],
            {"a": 1},
            False,
            id="override",
        ),
        pytest.param(
            "drop",
            [{"a": 1}, {"a": 2}, {"a": 3}],
            {},
            False,
            id="drop",
        ),
        pytest.param(
            "drop_conflicts",
            [{"a": 1, "b": 2}, {"b": 1, "c": 3}, {"c": 3, "d": 4}],
            {"a": 1, "c": 3, "d": 4},
            False,
            id="drop_conflicts",
        ),
        pytest.param(
            "no_conflicts",
            [{"a": 1}, {"b": 2}, {"b": 3}],
            None,
            True,
            id="no_conflicts",
        ),
    ),
)
def test_keep_attrs_strategies_dataarray(strategy, attrs, expected, error) -> None:
    a = xr.DataArray(dims="x", data=[0, 1], attrs=attrs[0])
    b = xr.DataArray(dims="x", data=[0, 1], attrs=attrs[1])
    c = xr.DataArray(dims="x", data=[0, 1], attrs=attrs[2])

    if error:
        with pytest.raises(xr.MergeError):
            apply_ufunc(lambda *args: sum(args), a, b, c, keep_attrs=strategy)
    else:
        expected = xr.DataArray(dims="x", data=[0, 3], attrs=expected)
        actual = apply_ufunc(lambda *args: sum(args), a, b, c, keep_attrs=strategy)

        assert_identical(actual, expected)


@pytest.mark.parametrize("variant", ("dim", "coord"))
@pytest.mark.parametrize(
    ["strategy", "attrs", "expected", "error"],
    (
        pytest.param(
            None,
            [{"a": 1}, {"a": 2}, {"a": 3}],
            {},
            False,
            id="default",
        ),
        pytest.param(
            False,
            [{"a": 1}, {"a": 2}, {"a": 3}],
            {},
            False,
            id="False",
        ),
        pytest.param(
            True,
            [{"a": 1}, {"a": 2}, {"a": 3}],
            {"a": 1},
            False,
            id="True",
        ),
        pytest.param(
            "override",
            [{"a": 1}, {"a": 2}, {"a": 3}],
            {"a": 1},
            False,
            id="override",
        ),
        pytest.param(
            "drop",
            [{"a": 1}, {"a": 2}, {"a": 3}],
            {},
            False,
            id="drop",
        ),
        pytest.param(
            "drop_conflicts",
            [{"a": 1, "b": 2}, {"b": 1, "c": 3}, {"c": 3, "d": 4}],
            {"a": 1, "c": 3, "d": 4},
            False,
            id="drop_conflicts",
        ),
        pytest.param(
            "no_conflicts",
            [{"a": 1}, {"b": 2}, {"b": 3}],
            None,
            True,
            id="no_conflicts",
        ),
    ),
)
def test_keep_attrs_strategies_dataarray_variables(
    variant, strategy, attrs, expected, error
):
    compute_attrs = {
        "dim": lambda attrs, default: (attrs, default),
        "coord": lambda attrs, default: (default, attrs),
    }.get(variant)

    dim_attrs, coord_attrs = compute_attrs(attrs, [{}, {}, {}])

    a = xr.DataArray(
        dims="x",
        data=[0, 1],
        coords={"x": ("x", [0, 1], dim_attrs[0]), "u": ("x", [0, 1], coord_attrs[0])},
    )
    b = xr.DataArray(
        dims="x",
        data=[0, 1],
        coords={"x": ("x", [0, 1], dim_attrs[1]), "u": ("x", [0, 1], coord_attrs[1])},
    )
    c = xr.DataArray(
        dims="x",
        data=[0, 1],
        coords={"x": ("x", [0, 1], dim_attrs[2]), "u": ("x", [0, 1], coord_attrs[2])},
    )

    if error:
        with pytest.raises(xr.MergeError):
            apply_ufunc(lambda *args: sum(args), a, b, c, keep_attrs=strategy)
    else:
        dim_attrs, coord_attrs = compute_attrs(expected, {})
        expected = xr.DataArray(
            dims="x",
            data=[0, 3],
            coords={"x": ("x", [0, 1], dim_attrs), "u": ("x", [0, 1], coord_attrs)},
        )
        actual = apply_ufunc(lambda *args: sum(args), a, b, c, keep_attrs=strategy)

        assert_identical(actual, expected)


@pytest.mark.parametrize(
    ["strategy", "attrs", "expected", "error"],
    (
        pytest.param(
            None,
            [{"a": 1}, {"a": 2}, {"a": 3}],
            {},
            False,
            id="default",
        ),
        pytest.param(
            False,
            [{"a": 1}, {"a": 2}, {"a": 3}],
            {},
            False,
            id="False",
        ),
        pytest.param(
            True,
            [{"a": 1}, {"a": 2}, {"a": 3}],
            {"a": 1},
            False,
            id="True",
        ),
        pytest.param(
            "override",
            [{"a": 1}, {"a": 2}, {"a": 3}],
            {"a": 1},
            False,
            id="override",
        ),
        pytest.param(
            "drop",
            [{"a": 1}, {"a": 2}, {"a": 3}],
            {},
            False,
            id="drop",
        ),
        pytest.param(
            "drop_conflicts",
            [{"a": 1, "b": 2}, {"b": 1, "c": 3}, {"c": 3, "d": 4}],
            {"a": 1, "c": 3, "d": 4},
            False,
            id="drop_conflicts",
        ),
        pytest.param(
            "no_conflicts",
            [{"a": 1}, {"b": 2}, {"b": 3}],
            None,
            True,
            id="no_conflicts",
        ),
    ),
)
def test_keep_attrs_strategies_dataset(strategy, attrs, expected, error) -> None:
    a = xr.Dataset({"a": ("x", [0, 1])}, attrs=attrs[0])
    b = xr.Dataset({"a": ("x", [0, 1])}, attrs=attrs[1])
    c = xr.Dataset({"a": ("x", [0, 1])}, attrs=attrs[2])

    if error:
        with pytest.raises(xr.MergeError):
            apply_ufunc(lambda *args: sum(args), a, b, c, keep_attrs=strategy)
    else:
        expected = xr.Dataset({"a": ("x", [0, 3])}, attrs=expected)
        actual = apply_ufunc(lambda *args: sum(args), a, b, c, keep_attrs=strategy)

        assert_identical(actual, expected)


@pytest.mark.parametrize("variant", ("data", "dim", "coord"))
@pytest.mark.parametrize(
    ["strategy", "attrs", "expected", "error"],
    (
        pytest.param(
            None,
            [{"a": 1}, {"a": 2}, {"a": 3}],
            {},
            False,
            id="default",
        ),
        pytest.param(
            False,
            [{"a": 1}, {"a": 2}, {"a": 3}],
            {},
            False,
            id="False",
        ),
        pytest.param(
            True,
            [{"a": 1}, {"a": 2}, {"a": 3}],
            {"a": 1},
            False,
            id="True",
        ),
        pytest.param(
            "override",
            [{"a": 1}, {"a": 2}, {"a": 3}],
            {"a": 1},
            False,
            id="override",
        ),
        pytest.param(
            "drop",
            [{"a": 1}, {"a": 2}, {"a": 3}],
            {},
            False,
            id="drop",
        ),
        pytest.param(
            "drop_conflicts",
            [{"a": 1, "b": 2}, {"b": 1, "c": 3}, {"c": 3, "d": 4}],
            {"a": 1, "c": 3, "d": 4},
            False,
            id="drop_conflicts",
        ),
        pytest.param(
            "no_conflicts",
            [{"a": 1}, {"b": 2}, {"b": 3}],
            None,
            True,
            id="no_conflicts",
        ),
    ),
)
def test_keep_attrs_strategies_dataset_variables(
    variant, strategy, attrs, expected, error
):
    compute_attrs = {
        "data": lambda attrs, default: (attrs, default, default),
        "dim": lambda attrs, default: (default, attrs, default),
        "coord": lambda attrs, default: (default, default, attrs),
    }.get(variant)
    data_attrs, dim_attrs, coord_attrs = compute_attrs(attrs, [{}, {}, {}])

    a = xr.Dataset(
        {"a": ("x", [], data_attrs[0])},
        coords={"x": ("x", [], dim_attrs[0]), "u": ("x", [], coord_attrs[0])},
    )
    b = xr.Dataset(
        {"a": ("x", [], data_attrs[1])},
        coords={"x": ("x", [], dim_attrs[1]), "u": ("x", [], coord_attrs[1])},
    )
    c = xr.Dataset(
        {"a": ("x", [], data_attrs[2])},
        coords={"x": ("x", [], dim_attrs[2]), "u": ("x", [], coord_attrs[2])},
    )

    if error:
        with pytest.raises(xr.MergeError):
            apply_ufunc(lambda *args: sum(args), a, b, c, keep_attrs=strategy)
    else:
        data_attrs, dim_attrs, coord_attrs = compute_attrs(expected, {})
        expected = xr.Dataset(
            {"a": ("x", [], data_attrs)},
            coords={"x": ("x", [], dim_attrs), "u": ("x", [], coord_attrs)},
        )
        actual = apply_ufunc(lambda *args: sum(args), a, b, c, keep_attrs=strategy)

        assert_identical(actual, expected)


def test_dataset_join() -> None:
    ds0 = xr.Dataset({"a": ("x", [1, 2]), "x": [0, 1]})
    ds1 = xr.Dataset({"a": ("x", [99, 3]), "x": [1, 2]})

    # by default, cannot have different labels
    with pytest.raises(ValueError, match=r"cannot align.*join.*exact.*"):
        apply_ufunc(operator.add, ds0, ds1)
    with pytest.raises(TypeError, match=r"must supply"):
        apply_ufunc(operator.add, ds0, ds1, dataset_join="outer")

    def add(a, b, join, dataset_join):
        return apply_ufunc(
            operator.add,
            a,
            b,
            join=join,
            dataset_join=dataset_join,
            dataset_fill_value=np.nan,
        )

    actual = add(ds0, ds1, "outer", "inner")
    expected = xr.Dataset({"a": ("x", [np.nan, 101, np.nan]), "x": [0, 1, 2]})
    assert_identical(actual, expected)

    actual = add(ds0, ds1, "outer", "outer")
    assert_identical(actual, expected)

    with pytest.raises(ValueError, match=r"data variable names"):
        apply_ufunc(operator.add, ds0, xr.Dataset({"b": 1}))

    ds2 = xr.Dataset({"b": ("x", [99, 3]), "x": [1, 2]})
    actual = add(ds0, ds2, "outer", "inner")
    expected = xr.Dataset({"x": [0, 1, 2]})
    assert_identical(actual, expected)

    # we used np.nan as the fill_value in add() above
    actual = add(ds0, ds2, "outer", "outer")
    expected = xr.Dataset(
        {
            "a": ("x", [np.nan, np.nan, np.nan]),
            "b": ("x", [np.nan, np.nan, np.nan]),
            "x": [0, 1, 2],
        }
    )
    assert_identical(actual, expected)


@requires_dask
def test_apply_dask() -> None:
    import dask.array as da

    array = da.ones((2,), chunks=2)
    variable = xr.Variable("x", array)
    coords = xr.DataArray(variable).coords.variables
    data_array = xr.DataArray(variable, dims=["x"], coords=coords)
    dataset = xr.Dataset({"y": variable})

    # encountered dask array, but did not set dask='allowed'
    with pytest.raises(ValueError):
        apply_ufunc(identity, array)
    with pytest.raises(ValueError):
        apply_ufunc(identity, variable)
    with pytest.raises(ValueError):
        apply_ufunc(identity, data_array)
    with pytest.raises(ValueError):
        apply_ufunc(identity, dataset)

    # unknown setting for dask array handling
    with pytest.raises(ValueError):
        apply_ufunc(identity, array, dask="unknown")  # type: ignore[arg-type]

    def dask_safe_identity(x):
        return apply_ufunc(identity, x, dask="allowed")

    assert array is dask_safe_identity(array)

    actual = dask_safe_identity(variable)
    assert isinstance(actual.data, da.Array)
    assert_identical(variable, actual)

    actual = dask_safe_identity(data_array)
    assert isinstance(actual.data, da.Array)
    assert_identical(data_array, actual)

    actual = dask_safe_identity(dataset)
    assert isinstance(actual["y"].data, da.Array)
    assert_identical(dataset, actual)


@requires_dask
def test_apply_dask_parallelized_one_arg() -> None:
    import dask.array as da

    array = da.ones((2, 2), chunks=(1, 1))
    data_array = xr.DataArray(array, dims=("x", "y"))

    def parallel_identity(x):
        return apply_ufunc(identity, x, dask="parallelized", output_dtypes=[x.dtype])

    actual = parallel_identity(data_array)
    assert isinstance(actual.data, da.Array)
    assert actual.data.chunks == array.chunks
    assert_identical(data_array, actual)

    computed = data_array.compute()
    actual = parallel_identity(computed)
    assert_identical(computed, actual)


@requires_dask
def test_apply_dask_parallelized_two_args() -> None:
    import dask.array as da

    array = da.ones((2, 2), chunks=(1, 1), dtype=np.int64)
    data_array = xr.DataArray(array, dims=("x", "y"))
    data_array.name = None

    def parallel_add(x, y):
        return apply_ufunc(
            operator.add, x, y, dask="parallelized", output_dtypes=[np.int64]
        )

    def check(x, y):
        actual = parallel_add(x, y)
        assert isinstance(actual.data, da.Array)
        assert actual.data.chunks == array.chunks
        assert_identical(data_array, actual)

    check(data_array, 0)
    check(0, data_array)
    check(data_array, xr.DataArray(0))
    check(data_array, 0 * data_array)
    check(data_array, 0 * data_array[0])
    check(data_array[:, 0], 0 * data_array[0])
    check(data_array, 0 * data_array.compute())


@requires_dask
def test_apply_dask_parallelized_errors() -> None:
    import dask.array as da

    array = da.ones((2, 2), chunks=(1, 1))
    data_array = xr.DataArray(array, dims=("x", "y"))

    # from apply_array_ufunc
    with pytest.raises(ValueError, match=r"at least one input is an xarray object"):
        apply_ufunc(identity, array, dask="parallelized")

    # formerly from _apply_blockwise, now from apply_variable_ufunc
    with pytest.raises(ValueError, match=r"consists of multiple chunks"):
        apply_ufunc(
            identity,
            data_array,
            dask="parallelized",
            output_dtypes=[float],
            input_core_dims=[("y",)],
            output_core_dims=[("y",)],
        )


# it's currently impossible to silence these warnings from inside dask.array:
# https://github.com/dask/dask/issues/3245
@requires_dask
@pytest.mark.filterwarnings("ignore:Mean of empty slice")
def test_apply_dask_multiple_inputs() -> None:
    import dask.array as da

    def covariance(x, y):
        return (
            (x - x.mean(axis=-1, keepdims=True)) * (y - y.mean(axis=-1, keepdims=True))
        ).mean(axis=-1)

    rs = np.random.default_rng(42)
    array1 = da.from_array(rs.random((4, 4)), chunks=(2, 4))
    array2 = da.from_array(rs.random((4, 4)), chunks=(2, 4))
    data_array_1 = xr.DataArray(array1, dims=("x", "z"))
    data_array_2 = xr.DataArray(array2, dims=("y", "z"))

    expected = apply_ufunc(
        covariance,
        data_array_1.compute(),
        data_array_2.compute(),
        input_core_dims=[["z"], ["z"]],
    )
    allowed = apply_ufunc(
        covariance,
        data_array_1,
        data_array_2,
        input_core_dims=[["z"], ["z"]],
        dask="allowed",
    )
    assert isinstance(allowed.data, da.Array)
    xr.testing.assert_allclose(expected, allowed.compute())

    parallelized = apply_ufunc(
        covariance,
        data_array_1,
        data_array_2,
        input_core_dims=[["z"], ["z"]],
        dask="parallelized",
        output_dtypes=[float],
    )
    assert isinstance(parallelized.data, da.Array)
    xr.testing.assert_allclose(expected, parallelized.compute())


@requires_dask
def test_apply_dask_new_output_dimension() -> None:
    import dask.array as da

    array = da.ones((2, 2), chunks=(1, 1))
    data_array = xr.DataArray(array, dims=("x", "y"))

    def stack_negative(obj):
        def func(x):
            return np.stack([x, -x], axis=-1)

        return apply_ufunc(
            func,
            obj,
            output_core_dims=[["sign"]],
            dask="parallelized",
            output_dtypes=[obj.dtype],
            dask_gufunc_kwargs=dict(output_sizes={"sign": 2}),
        )

    expected = stack_negative(data_array.compute())

    actual = stack_negative(data_array)
    assert actual.dims == ("x", "y", "sign")
    assert actual.shape == (2, 2, 2)
    assert isinstance(actual.data, da.Array)
    assert_identical(expected, actual)


@requires_dask
def test_apply_dask_new_output_sizes() -> None:
    ds = xr.Dataset({"foo": (["lon", "lat"], np.arange(10 * 10).reshape((10, 10)))})
    ds["bar"] = ds["foo"]
    newdims = {"lon_new": 3, "lat_new": 6}

    def extract(obj):
        def func(da):
            return da[1:4, 1:7]

        return apply_ufunc(
            func,
            obj,
            dask="parallelized",
            input_core_dims=[["lon", "lat"]],
            output_core_dims=[["lon_new", "lat_new"]],
            dask_gufunc_kwargs=dict(output_sizes=newdims),
        )

    expected = extract(ds)

    actual = extract(ds.chunk())
    assert actual.sizes == {"lon_new": 3, "lat_new": 6}
    assert_identical(expected.chunk(), actual)


@requires_dask
def test_apply_dask_new_output_sizes_not_supplied_same_dim_names() -> None:
    # test for missing output_sizes kwarg sneaking through
    # see GH discussion 7503

    data = np.random.randn(4, 4, 3, 2)
    da = xr.DataArray(data=data, dims=("x", "y", "i", "j")).chunk(x=1, y=1)

    with pytest.raises(ValueError, match="output_sizes"):
        xr.apply_ufunc(
            np.linalg.pinv,
            da,
            input_core_dims=[["i", "j"]],
            output_core_dims=[["i", "j"]],
            exclude_dims={"i", "j"},
            dask="parallelized",
        )


def pandas_median(x):
    return pd.Series(x).median()


def test_vectorize() -> None:
    data_array = xr.DataArray([[0, 1, 2], [1, 2, 3]], dims=("x", "y"))
    expected = xr.DataArray([1, 2], dims=["x"])
    actual = apply_ufunc(
        pandas_median, data_array, input_core_dims=[["y"]], vectorize=True
    )
    assert_identical(expected, actual)


@requires_dask
def test_vectorize_dask() -> None:
    # run vectorization in dask.array.gufunc by using `dask='parallelized'`
    data_array = xr.DataArray([[0, 1, 2], [1, 2, 3]], dims=("x", "y"))
    expected = xr.DataArray([1, 2], dims=["x"])
    actual = apply_ufunc(
        pandas_median,
        data_array.chunk({"x": 1}),
        input_core_dims=[["y"]],
        vectorize=True,
        dask="parallelized",
        output_dtypes=[float],
    )
    assert_identical(expected, actual)


@requires_dask
def test_vectorize_dask_dtype() -> None:
    # ensure output_dtypes is preserved with vectorize=True
    # GH4015

    # integer
    data_array = xr.DataArray([[0, 1, 2], [1, 2, 3]], dims=("x", "y"))
    expected = xr.DataArray([1, 2], dims=["x"])
    actual = apply_ufunc(
        pandas_median,
        data_array.chunk({"x": 1}),
        input_core_dims=[["y"]],
        vectorize=True,
        dask="parallelized",
        output_dtypes=[int],
    )
    assert_identical(expected, actual)
    assert expected.dtype == actual.dtype

    # complex
    data_array = xr.DataArray([[0 + 0j, 1 + 2j, 2 + 1j]], dims=("x", "y"))
    expected = data_array.copy()
    actual = apply_ufunc(
        identity,
        data_array.chunk({"x": 1}),
        vectorize=True,
        dask="parallelized",
        output_dtypes=[complex],
    )
    assert_identical(expected, actual)
    assert expected.dtype == actual.dtype


@requires_dask
@pytest.mark.parametrize(
    "data_array",
    [
        xr.DataArray([[0, 1, 2], [1, 2, 3]], dims=("x", "y")),
        xr.DataArray([[0 + 0j, 1 + 2j, 2 + 1j]], dims=("x", "y")),
    ],
)
def test_vectorize_dask_dtype_without_output_dtypes(data_array) -> None:
    # ensure output_dtypes is preserved with vectorize=True
    # GH4015

    expected = data_array.copy()
    actual = apply_ufunc(
        identity,
        data_array.chunk({"x": 1}),
        vectorize=True,
        dask="parallelized",
    )

    assert_identical(expected, actual)
    assert expected.dtype == actual.dtype


@requires_dask
def test_vectorize_dask_dtype_meta() -> None:
    data_array = xr.DataArray([[0, 1, 2], [1, 2, 3]], dims=("x", "y"))
    expected = xr.DataArray([1, 2], dims=["x"])

    actual = apply_ufunc(
        pandas_median,
        data_array.chunk({"x": 1}),
        input_core_dims=[["y"]],
        vectorize=True,
        dask="parallelized",
        dask_gufunc_kwargs=dict(meta=np.ndarray((0, 0), dtype=float)),
    )

    assert_identical(expected, actual)
    assert float == actual.dtype


def pandas_median_add(x, y):
    # function which can consume input of unequal length
    return pd.Series(x).median() + pd.Series(y).median()


def test_vectorize_exclude_dims() -> None:
    # GH 3890
    data_array_a = xr.DataArray([[0, 1, 2], [1, 2, 3]], dims=("x", "y"))
    data_array_b = xr.DataArray([[0, 1, 2, 3, 4], [1, 2, 3, 4, 5]], dims=("x", "y"))

    expected = xr.DataArray([3, 5], dims=["x"])
    actual = apply_ufunc(
        pandas_median_add,
        data_array_a,
        data_array_b,
        input_core_dims=[["y"], ["y"]],
        vectorize=True,
        exclude_dims=set("y"),
    )
    assert_identical(expected, actual)


@requires_dask
def test_vectorize_exclude_dims_dask() -> None:
    # GH 3890
    data_array_a = xr.DataArray([[0, 1, 2], [1, 2, 3]], dims=("x", "y"))
    data_array_b = xr.DataArray([[0, 1, 2, 3, 4], [1, 2, 3, 4, 5]], dims=("x", "y"))

    expected = xr.DataArray([3, 5], dims=["x"])
    actual = apply_ufunc(
        pandas_median_add,
        data_array_a.chunk({"x": 1}),
        data_array_b.chunk({"x": 1}),
        input_core_dims=[["y"], ["y"]],
        exclude_dims=set("y"),
        vectorize=True,
        dask="parallelized",
        output_dtypes=[float],
    )
    assert_identical(expected, actual)


def test_corr_only_dataarray() -> None:
    with pytest.raises(TypeError, match="Only xr.DataArray is supported"):
        xr.corr(xr.Dataset(), xr.Dataset())  # type: ignore[type-var]


@pytest.fixture(scope="module")
def arrays():
    da = xr.DataArray(
        np.random.random((3, 21, 4)),
        coords={"time": pd.date_range("2000-01-01", freq="1D", periods=21)},
        dims=("a", "time", "x"),
    )

    return [
        da.isel(time=range(18)),
        da.isel(time=range(2, 20)).rolling(time=3, center=True).mean(),
        xr.DataArray([[1, 2], [1, np.nan]], dims=["x", "time"]),
        xr.DataArray([[1, 2], [np.nan, np.nan]], dims=["x", "time"]),
        xr.DataArray([[1, 2], [2, 1]], dims=["x", "time"]),
    ]


@pytest.fixture(scope="module")
def array_tuples(arrays):
    return [
        (arrays[0], arrays[0]),
        (arrays[0], arrays[1]),
        (arrays[1], arrays[1]),
        (arrays[2], arrays[2]),
        (arrays[2], arrays[3]),
        (arrays[2], arrays[4]),
        (arrays[4], arrays[2]),
        (arrays[3], arrays[3]),
        (arrays[4], arrays[4]),
    ]


@pytest.mark.parametrize("ddof", [0, 1])
@pytest.mark.parametrize("n", [3, 4, 5, 6, 7, 8])
@pytest.mark.parametrize("dim", [None, "x", "time"])
@requires_dask
def test_lazy_corrcov(
    n: int, dim: str | None, ddof: int, array_tuples: tuple[xr.DataArray, xr.DataArray]
) -> None:
    # GH 5284
    from dask import is_dask_collection

    da_a, da_b = array_tuples[n]

    with raise_if_dask_computes():
        cov = xr.cov(da_a.chunk(), da_b.chunk(), dim=dim, ddof=ddof)
        assert is_dask_collection(cov)

        corr = xr.corr(da_a.chunk(), da_b.chunk(), dim=dim)
        assert is_dask_collection(corr)


@pytest.mark.parametrize("ddof", [0, 1])
@pytest.mark.parametrize("n", [0, 1, 2])
@pytest.mark.parametrize("dim", [None, "time"])
def test_cov(
    n: int, dim: str | None, ddof: int, array_tuples: tuple[xr.DataArray, xr.DataArray]
) -> None:
    da_a, da_b = array_tuples[n]

    if dim is not None:

        def np_cov_ind(ts1, ts2, a, x):
            # Ensure the ts are aligned and missing values ignored
            ts1, ts2 = broadcast(ts1, ts2)
            valid_values = ts1.notnull() & ts2.notnull()

            # While dropping isn't ideal here, numpy will return nan
            # if any segment contains a NaN.
            ts1 = ts1.where(valid_values)
            ts2 = ts2.where(valid_values)

            return np.ma.cov(
                np.ma.masked_invalid(ts1.sel(a=a, x=x).data.flatten()),
                np.ma.masked_invalid(ts2.sel(a=a, x=x).data.flatten()),
                ddof=ddof,
            )[0, 1]

        expected = np.zeros((3, 4))
        for a in [0, 1, 2]:
            for x in [0, 1, 2, 3]:
                expected[a, x] = np_cov_ind(da_a, da_b, a=a, x=x)
        actual = xr.cov(da_a, da_b, dim=dim, ddof=ddof)
        assert_allclose(actual, expected)

    else:

        def np_cov(ts1, ts2):
            # Ensure the ts are aligned and missing values ignored
            ts1, ts2 = broadcast(ts1, ts2)
            valid_values = ts1.notnull() & ts2.notnull()

            ts1 = ts1.where(valid_values)
            ts2 = ts2.where(valid_values)

            return np.ma.cov(
                np.ma.masked_invalid(ts1.data.flatten()),
                np.ma.masked_invalid(ts2.data.flatten()),
                ddof=ddof,
            )[0, 1]

        expected = np_cov(da_a, da_b)
        actual = xr.cov(da_a, da_b, dim=dim, ddof=ddof)
        assert_allclose(actual, expected)


@pytest.mark.parametrize("n", [0, 1, 2])
@pytest.mark.parametrize("dim", [None, "time"])
def test_corr(
    n: int, dim: str | None, array_tuples: tuple[xr.DataArray, xr.DataArray]
) -> None:
    da_a, da_b = array_tuples[n]

    if dim is not None:

        def np_corr_ind(ts1, ts2, a, x):
            # Ensure the ts are aligned and missing values ignored
            ts1, ts2 = broadcast(ts1, ts2)
            valid_values = ts1.notnull() & ts2.notnull()

            ts1 = ts1.where(valid_values)
            ts2 = ts2.where(valid_values)

            return np.ma.corrcoef(
                np.ma.masked_invalid(ts1.sel(a=a, x=x).data.flatten()),
                np.ma.masked_invalid(ts2.sel(a=a, x=x).data.flatten()),
            )[0, 1]

        expected = np.zeros((3, 4))
        for a in [0, 1, 2]:
            for x in [0, 1, 2, 3]:
                expected[a, x] = np_corr_ind(da_a, da_b, a=a, x=x)
        actual = xr.corr(da_a, da_b, dim)
        assert_allclose(actual, expected)

    else:

        def np_corr(ts1, ts2):
            # Ensure the ts are aligned and missing values ignored
            ts1, ts2 = broadcast(ts1, ts2)
            valid_values = ts1.notnull() & ts2.notnull()

            ts1 = ts1.where(valid_values)
            ts2 = ts2.where(valid_values)

            return np.ma.corrcoef(
                np.ma.masked_invalid(ts1.data.flatten()),
                np.ma.masked_invalid(ts2.data.flatten()),
            )[0, 1]

        expected = np_corr(da_a, da_b)
        actual = xr.corr(da_a, da_b, dim)
        assert_allclose(actual, expected)


@pytest.mark.parametrize("n", range(9))
@pytest.mark.parametrize("dim", [None, "time", "x"])
def test_covcorr_consistency(
    n: int, dim: str | None, array_tuples: tuple[xr.DataArray, xr.DataArray]
) -> None:
    da_a, da_b = array_tuples[n]
    # Testing that xr.corr and xr.cov are consistent with each other
    # 1. Broadcast the two arrays
    da_a, da_b = broadcast(da_a, da_b)
    # 2. Ignore the nans
    valid_values = da_a.notnull() & da_b.notnull()
    da_a = da_a.where(valid_values)
    da_b = da_b.where(valid_values)

    expected = xr.cov(da_a, da_b, dim=dim, ddof=0) / (
        da_a.std(dim=dim) * da_b.std(dim=dim)
    )
    actual = xr.corr(da_a, da_b, dim=dim)
    assert_allclose(actual, expected)


@requires_dask
@pytest.mark.parametrize("n", range(9))
@pytest.mark.parametrize("dim", [None, "time", "x"])
@pytest.mark.filterwarnings("ignore:invalid value encountered in .*divide")
def test_corr_lazycorr_consistency(
    n: int, dim: str | None, array_tuples: tuple[xr.DataArray, xr.DataArray]
) -> None:
    da_a, da_b = array_tuples[n]
    da_al = da_a.chunk()
    da_bl = da_b.chunk()
    c_abl = xr.corr(da_al, da_bl, dim=dim)
    c_ab = xr.corr(da_a, da_b, dim=dim)
    c_ab_mixed = xr.corr(da_a, da_bl, dim=dim)
    assert_allclose(c_ab, c_abl)
    assert_allclose(c_ab, c_ab_mixed)


@requires_dask
def test_corr_dtype_error():
    da_a = xr.DataArray([[1, 2], [2, 1]], dims=["x", "time"])
    da_b = xr.DataArray([[1, 2], [1, np.nan]], dims=["x", "time"])

    xr.testing.assert_equal(xr.corr(da_a, da_b), xr.corr(da_a.chunk(), da_b.chunk()))
    xr.testing.assert_equal(xr.corr(da_a, da_b), xr.corr(da_a, da_b.chunk()))


@pytest.mark.parametrize("n", range(5))
@pytest.mark.parametrize("dim", [None, "time", "x", ["time", "x"]])
def test_autocov(n: int, dim: str | None, arrays) -> None:
    da = arrays[n]

    # Testing that the autocovariance*(N-1) is ~=~ to the variance matrix
    # 1. Ignore the nans
    valid_values = da.notnull()
    # Because we're using ddof=1, this requires > 1 value in each sample
    da = da.where(valid_values.sum(dim=dim) > 1)
    expected = ((da - da.mean(dim=dim)) ** 2).sum(dim=dim, skipna=True, min_count=1)
    actual = xr.cov(da, da, dim=dim) * (valid_values.sum(dim) - 1)
    assert_allclose(actual, expected)


def test_complex_cov() -> None:
    da = xr.DataArray([1j, -1j])
    actual = xr.cov(da, da)
    assert abs(actual.item()) == 2


@pytest.mark.parametrize("weighted", [True, False])
def test_bilinear_cov_corr(weighted: bool) -> None:
    # Test the bilinear properties of covariance and correlation
    da = xr.DataArray(
        np.random.random((3, 21, 4)),
        coords={"time": pd.date_range("2000-01-01", freq="1D", periods=21)},
        dims=("a", "time", "x"),
    )
    db = xr.DataArray(
        np.random.random((3, 21, 4)),
        coords={"time": pd.date_range("2000-01-01", freq="1D", periods=21)},
        dims=("a", "time", "x"),
    )
    dc = xr.DataArray(
        np.random.random((3, 21, 4)),
        coords={"time": pd.date_range("2000-01-01", freq="1D", periods=21)},
        dims=("a", "time", "x"),
    )
    if weighted:
        weights = xr.DataArray(
            np.abs(np.random.random(4)),
            dims=("x"),
        )
    else:
        weights = None
    k = np.random.random(1)[0]

    # Test covariance properties
    assert_allclose(
        xr.cov(da + k, db, weights=weights), xr.cov(da, db, weights=weights)
    )
    assert_allclose(
        xr.cov(da, db + k, weights=weights), xr.cov(da, db, weights=weights)
    )
    assert_allclose(
        xr.cov(da + dc, db, weights=weights),
        xr.cov(da, db, weights=weights) + xr.cov(dc, db, weights=weights),
    )
    assert_allclose(
        xr.cov(da, db + dc, weights=weights),
        xr.cov(da, db, weights=weights) + xr.cov(da, dc, weights=weights),
    )
    assert_allclose(
        xr.cov(k * da, db, weights=weights), k * xr.cov(da, db, weights=weights)
    )
    assert_allclose(
        xr.cov(da, k * db, weights=weights), k * xr.cov(da, db, weights=weights)
    )

    # Test correlation properties
    assert_allclose(
        xr.corr(da + k, db, weights=weights), xr.corr(da, db, weights=weights)
    )
    assert_allclose(
        xr.corr(da, db + k, weights=weights), xr.corr(da, db, weights=weights)
    )
    assert_allclose(
        xr.corr(k * da, db, weights=weights), xr.corr(da, db, weights=weights)
    )
    assert_allclose(
        xr.corr(da, k * db, weights=weights), xr.corr(da, db, weights=weights)
    )


def test_equally_weighted_cov_corr() -> None:
    # Test that equal weights for all values produces same results as weights=None
    da = xr.DataArray(
        np.random.random((3, 21, 4)),
        coords={"time": pd.date_range("2000-01-01", freq="1D", periods=21)},
        dims=("a", "time", "x"),
    )
    db = xr.DataArray(
        np.random.random((3, 21, 4)),
        coords={"time": pd.date_range("2000-01-01", freq="1D", periods=21)},
        dims=("a", "time", "x"),
    )
    assert_allclose(
        xr.cov(da, db, weights=None), xr.cov(da, db, weights=xr.DataArray(1))
    )
    assert_allclose(
        xr.cov(da, db, weights=None), xr.cov(da, db, weights=xr.DataArray(2))
    )
    assert_allclose(
        xr.corr(da, db, weights=None), xr.corr(da, db, weights=xr.DataArray(1))
    )
    assert_allclose(
        xr.corr(da, db, weights=None), xr.corr(da, db, weights=xr.DataArray(2))
    )


@requires_dask
def test_vectorize_dask_new_output_dims() -> None:
    # regression test for GH3574
    # run vectorization in dask.array.gufunc by using `dask='parallelized'`
    data_array = xr.DataArray([[0, 1, 2], [1, 2, 3]], dims=("x", "y"))
    func = lambda x: x[np.newaxis, ...]
    expected = data_array.expand_dims("z")
    actual = apply_ufunc(
        func,
        data_array.chunk({"x": 1}),
        output_core_dims=[["z"]],
        vectorize=True,
        dask="parallelized",
        output_dtypes=[float],
        dask_gufunc_kwargs=dict(output_sizes={"z": 1}),
    ).transpose(*expected.dims)
    assert_identical(expected, actual)

    with pytest.raises(
        ValueError, match=r"dimension 'z1' in 'output_sizes' must correspond"
    ):
        apply_ufunc(
            func,
            data_array.chunk({"x": 1}),
            output_core_dims=[["z"]],
            vectorize=True,
            dask="parallelized",
            output_dtypes=[float],
            dask_gufunc_kwargs=dict(output_sizes={"z1": 1}),
        )

    with pytest.raises(
        ValueError, match=r"dimension 'z' in 'output_core_dims' needs corresponding"
    ):
        apply_ufunc(
            func,
            data_array.chunk({"x": 1}),
            output_core_dims=[["z"]],
            vectorize=True,
            dask="parallelized",
            output_dtypes=[float],
        )


def test_output_wrong_number() -> None:
    variable = xr.Variable("x", np.arange(10))

    def identity(x):
        return x

    def tuple3x(x):
        return (x, x, x)

    with pytest.raises(
        ValueError,
        match=r"number of outputs.* Received a <class 'numpy.ndarray'> with 10 elements. Expected a tuple of 2 elements:\n\narray\(\[0",
    ):
        apply_ufunc(identity, variable, output_core_dims=[(), ()])

    with pytest.raises(ValueError, match=r"number of outputs"):
        apply_ufunc(tuple3x, variable, output_core_dims=[(), ()])


def test_output_wrong_dims() -> None:
    variable = xr.Variable("x", np.arange(10))

    def add_dim(x):
        return x[..., np.newaxis]

    def remove_dim(x):
        return x[..., 0]

    with pytest.raises(
        ValueError,
        match=r"unexpected number of dimensions.*from:\n\n.*array\(\[\[0",
    ):
        apply_ufunc(add_dim, variable, output_core_dims=[("y", "z")])

    with pytest.raises(ValueError, match=r"unexpected number of dimensions"):
        apply_ufunc(add_dim, variable)

    with pytest.raises(ValueError, match=r"unexpected number of dimensions"):
        apply_ufunc(remove_dim, variable)


def test_output_wrong_dim_size() -> None:
    array = np.arange(10)
    variable = xr.Variable("x", array)
    data_array = xr.DataArray(variable, [("x", -array)])
    dataset = xr.Dataset({"y": variable}, {"x": -array})

    def truncate(array):
        return array[:5]

    def apply_truncate_broadcast_invalid(obj):
        return apply_ufunc(truncate, obj)

    with pytest.raises(ValueError, match=r"size of dimension"):
        apply_truncate_broadcast_invalid(variable)
    with pytest.raises(ValueError, match=r"size of dimension"):
        apply_truncate_broadcast_invalid(data_array)
    with pytest.raises(ValueError, match=r"size of dimension"):
        apply_truncate_broadcast_invalid(dataset)

    def apply_truncate_x_x_invalid(obj):
        return apply_ufunc(
            truncate, obj, input_core_dims=[["x"]], output_core_dims=[["x"]]
        )

    with pytest.raises(ValueError, match=r"size of dimension"):
        apply_truncate_x_x_invalid(variable)
    with pytest.raises(ValueError, match=r"size of dimension"):
        apply_truncate_x_x_invalid(data_array)
    with pytest.raises(ValueError, match=r"size of dimension"):
        apply_truncate_x_x_invalid(dataset)

    def apply_truncate_x_z(obj):
        return apply_ufunc(
            truncate, obj, input_core_dims=[["x"]], output_core_dims=[["z"]]
        )

    assert_identical(xr.Variable("z", array[:5]), apply_truncate_x_z(variable))
    assert_identical(
        xr.DataArray(array[:5], dims=["z"]), apply_truncate_x_z(data_array)
    )
    assert_identical(xr.Dataset({"y": ("z", array[:5])}), apply_truncate_x_z(dataset))

    def apply_truncate_x_x_valid(obj):
        return apply_ufunc(
            truncate,
            obj,
            input_core_dims=[["x"]],
            output_core_dims=[["x"]],
            exclude_dims={"x"},
        )

    assert_identical(xr.Variable("x", array[:5]), apply_truncate_x_x_valid(variable))
    assert_identical(
        xr.DataArray(array[:5], dims=["x"]), apply_truncate_x_x_valid(data_array)
    )
    assert_identical(
        xr.Dataset({"y": ("x", array[:5])}), apply_truncate_x_x_valid(dataset)
    )


@pytest.mark.parametrize("use_dask", [True, False])
def test_dot(use_dask: bool) -> None:
    if use_dask and not has_dask:
        pytest.skip("test for dask.")

    a = np.arange(30 * 4).reshape(30, 4)
    b = np.arange(30 * 4 * 5).reshape(30, 4, 5)
    c = np.arange(5 * 60).reshape(5, 60)
    da_a = xr.DataArray(a, dims=["a", "b"], coords={"a": np.linspace(0, 1, 30)})
    da_b = xr.DataArray(b, dims=["a", "b", "c"], coords={"a": np.linspace(0, 1, 30)})
    da_c = xr.DataArray(c, dims=["c", "e"])
    if use_dask:
        da_a = da_a.chunk({"a": 3})
        da_b = da_b.chunk({"a": 3})
        da_c = da_c.chunk({"c": 3})
    actual = xr.dot(da_a, da_b, dim=["a", "b"])
    assert actual.dims == ("c",)
    assert (actual.data == np.einsum("ij,ijk->k", a, b)).all()
    assert isinstance(actual.variable.data, type(da_a.variable.data))

    actual = xr.dot(da_a, da_b)
    assert actual.dims == ("c",)
    assert (actual.data == np.einsum("ij,ijk->k", a, b)).all()
    assert isinstance(actual.variable.data, type(da_a.variable.data))

    # for only a single array is passed without dims argument, just return
    # as is
    actual = xr.dot(da_a)
    assert_identical(da_a, actual)

    # test for variable
    actual = xr.dot(da_a.variable, da_b.variable)
    assert actual.dims == ("c",)
    assert (actual.data == np.einsum("ij,ijk->k", a, b)).all()
    assert isinstance(actual.data, type(da_a.variable.data))

    if use_dask:
        da_a = da_a.chunk({"a": 3})
        da_b = da_b.chunk({"a": 3})
        actual = xr.dot(da_a, da_b, dim=["b"])
        assert actual.dims == ("a", "c")
        assert (actual.data == np.einsum("ij,ijk->ik", a, b)).all()
        assert isinstance(actual.variable.data, type(da_a.variable.data))

    actual = xr.dot(da_a, da_b, dim=["b"])
    assert actual.dims == ("a", "c")
    assert (actual.data == np.einsum("ij,ijk->ik", a, b)).all()

    actual = xr.dot(da_a, da_b, dim="b")
    assert actual.dims == ("a", "c")
    assert (actual.data == np.einsum("ij,ijk->ik", a, b)).all()

    actual = xr.dot(da_a, da_b, dim="a")
    assert actual.dims == ("b", "c")
    assert (actual.data == np.einsum("ij,ijk->jk", a, b)).all()

    actual = xr.dot(da_a, da_b, dim="c")
    assert actual.dims == ("a", "b")
    assert (actual.data == np.einsum("ij,ijk->ij", a, b)).all()

    actual = xr.dot(da_a, da_b, da_c, dim=["a", "b"])
    assert actual.dims == ("c", "e")
    assert (actual.data == np.einsum("ij,ijk,kl->kl ", a, b, c)).all()

    # should work with tuple
    actual = xr.dot(da_a, da_b, dim=("c",))
    assert actual.dims == ("a", "b")
    assert (actual.data == np.einsum("ij,ijk->ij", a, b)).all()

    # default dims
    actual = xr.dot(da_a, da_b, da_c)
    assert actual.dims == ("e",)
    assert (actual.data == np.einsum("ij,ijk,kl->l ", a, b, c)).all()

    # 1 array summation
    actual = xr.dot(da_a, dim="a")
    assert actual.dims == ("b",)
    assert (actual.data == np.einsum("ij->j ", a)).all()

    # empty dim
    actual = xr.dot(da_a.sel(a=[]), da_a.sel(a=[]), dim="a")
    assert actual.dims == ("b",)
    assert (actual.data == np.zeros(actual.shape)).all()

    # Ellipsis (...) sums over all dimensions
    actual = xr.dot(da_a, da_b, dim=...)
    assert actual.dims == ()
    assert (actual.data == np.einsum("ij,ijk->", a, b)).all()

    actual = xr.dot(da_a, da_b, da_c, dim=...)
    assert actual.dims == ()
    assert (actual.data == np.einsum("ij,ijk,kl-> ", a, b, c)).all()

    actual = xr.dot(da_a, dim=...)
    assert actual.dims == ()
    assert (actual.data == np.einsum("ij-> ", a)).all()

    actual = xr.dot(da_a.sel(a=[]), da_a.sel(a=[]), dim=...)
    assert actual.dims == ()
    assert (actual.data == np.zeros(actual.shape)).all()

    # Invalid cases
    if not use_dask:
        with pytest.raises(TypeError):
            xr.dot(da_a, dim="a", invalid=None)
    with pytest.raises(TypeError):
        xr.dot(da_a.to_dataset(name="da"), dim="a")
    with pytest.raises(TypeError):
        xr.dot(dim="a")

    # einsum parameters
    actual = xr.dot(da_a, da_b, dim=["b"], order="C")
    assert (actual.data == np.einsum("ij,ijk->ik", a, b)).all()
    assert actual.values.flags["C_CONTIGUOUS"]
    assert not actual.values.flags["F_CONTIGUOUS"]
    actual = xr.dot(da_a, da_b, dim=["b"], order="F")
    assert (actual.data == np.einsum("ij,ijk->ik", a, b)).all()
    # dask converts Fortran arrays to C order when merging the final array
    if not use_dask:
        assert not actual.values.flags["C_CONTIGUOUS"]
        assert actual.values.flags["F_CONTIGUOUS"]

    # einsum has a constant string as of the first parameter, which makes
    # it hard to pass to xarray.apply_ufunc.
    # make sure dot() uses functools.partial(einsum, subscripts), which
    # can be pickled, and not a lambda, which can't.
    pickle.loads(pickle.dumps(xr.dot(da_a)))


@pytest.mark.parametrize("use_dask", [True, False])
def test_dot_align_coords(use_dask: bool) -> None:
    # GH 3694

    if use_dask and not has_dask:
        pytest.skip("test for dask.")

    a = np.arange(30 * 4).reshape(30, 4)
    b = np.arange(30 * 4 * 5).reshape(30, 4, 5)

    # use partially overlapping coords
    coords_a = {"a": np.arange(30), "b": np.arange(4)}
    coords_b = {"a": np.arange(5, 35), "b": np.arange(1, 5)}

    da_a = xr.DataArray(a, dims=["a", "b"], coords=coords_a)
    da_b = xr.DataArray(b, dims=["a", "b", "c"], coords=coords_b)

    if use_dask:
        da_a = da_a.chunk({"a": 3})
        da_b = da_b.chunk({"a": 3})

    # join="inner" is the default
    actual = xr.dot(da_a, da_b)
    # `dot` sums over the common dimensions of the arguments
    expected = (da_a * da_b).sum(["a", "b"])
    xr.testing.assert_allclose(expected, actual)

    actual = xr.dot(da_a, da_b, dim=...)
    expected = (da_a * da_b).sum()
    xr.testing.assert_allclose(expected, actual)

    with xr.set_options(arithmetic_join="exact"):
        with pytest.raises(ValueError, match=r"cannot align.*join.*exact.*not equal.*"):
            xr.dot(da_a, da_b)

    # NOTE: dot always uses `join="inner"` because `(a * b).sum()` yields the same for all
    # join method (except "exact")
    with xr.set_options(arithmetic_join="left"):
        actual = xr.dot(da_a, da_b)
        expected = (da_a * da_b).sum(["a", "b"])
        xr.testing.assert_allclose(expected, actual)

    with xr.set_options(arithmetic_join="right"):
        actual = xr.dot(da_a, da_b)
        expected = (da_a * da_b).sum(["a", "b"])
        xr.testing.assert_allclose(expected, actual)

    with xr.set_options(arithmetic_join="outer"):
        actual = xr.dot(da_a, da_b)
        expected = (da_a * da_b).sum(["a", "b"])
        xr.testing.assert_allclose(expected, actual)


def test_where() -> None:
    cond = xr.DataArray([True, False], dims="x")
    actual = xr.where(cond, 1, 0)
    expected = xr.DataArray([1, 0], dims="x")
    assert_identical(expected, actual)


def test_where_attrs() -> None:
    cond = xr.DataArray([True, False], coords={"a": [0, 1]}, attrs={"attr": "cond_da"})
    cond["a"].attrs = {"attr": "cond_coord"}
    input_cond = cond.copy()
    x = xr.DataArray([1, 1], coords={"a": [0, 1]}, attrs={"attr": "x_da"})
    x["a"].attrs = {"attr": "x_coord"}
    y = xr.DataArray([0, 0], coords={"a": [0, 1]}, attrs={"attr": "y_da"})
    y["a"].attrs = {"attr": "y_coord"}

    # 3 DataArrays, takes attrs from x
    actual = xr.where(cond, x, y, keep_attrs=True)
    expected = xr.DataArray([1, 0], coords={"a": [0, 1]}, attrs={"attr": "x_da"})
    expected["a"].attrs = {"attr": "x_coord"}
    assert_identical(expected, actual)
    # Check also that input coordinate attributes weren't modified by reference
    assert x["a"].attrs == {"attr": "x_coord"}
    assert y["a"].attrs == {"attr": "y_coord"}
    assert cond["a"].attrs == {"attr": "cond_coord"}
    assert_identical(cond, input_cond)

    # 3 DataArrays, drop attrs
    actual = xr.where(cond, x, y, keep_attrs=False)
    expected = xr.DataArray([1, 0], coords={"a": [0, 1]})
    assert_identical(expected, actual)
    assert_identical(expected.coords["a"], actual.coords["a"])
    # Check also that input coordinate attributes weren't modified by reference
    assert x["a"].attrs == {"attr": "x_coord"}
    assert y["a"].attrs == {"attr": "y_coord"}
    assert cond["a"].attrs == {"attr": "cond_coord"}
    assert_identical(cond, input_cond)

    # x as a scalar, takes no attrs
    actual = xr.where(cond, 0, y, keep_attrs=True)
    expected = xr.DataArray([0, 0], coords={"a": [0, 1]})
    assert_identical(expected, actual)

    # y as a scalar, takes attrs from x
    actual = xr.where(cond, x, 0, keep_attrs=True)
    expected = xr.DataArray([1, 0], coords={"a": [0, 1]}, attrs={"attr": "x_da"})
    expected["a"].attrs = {"attr": "x_coord"}
    assert_identical(expected, actual)

    # x and y as a scalar, takes no attrs
    actual = xr.where(cond, 1, 0, keep_attrs=True)
    expected = xr.DataArray([1, 0], coords={"a": [0, 1]})
    assert_identical(expected, actual)

    # cond and y as a scalar, takes attrs from x
    actual = xr.where(True, x, y, keep_attrs=True)
    expected = xr.DataArray([1, 1], coords={"a": [0, 1]}, attrs={"attr": "x_da"})
    expected["a"].attrs = {"attr": "x_coord"}
    assert_identical(expected, actual)

    # no xarray objects, handle no attrs
    actual_np = xr.where(True, 0, 1, keep_attrs=True)
    expected_np = np.array(0)
    assert_identical(expected_np, actual_np)

    # DataArray and 2 Datasets, takes attrs from x
    ds_x = xr.Dataset(data_vars={"x": x}, attrs={"attr": "x_ds"})
    ds_y = xr.Dataset(data_vars={"x": y}, attrs={"attr": "y_ds"})
    ds_actual = xr.where(cond, ds_x, ds_y, keep_attrs=True)
    ds_expected = xr.Dataset(
        data_vars={
            "x": xr.DataArray([1, 0], coords={"a": [0, 1]}, attrs={"attr": "x_da"})
        },
        attrs={"attr": "x_ds"},
    )
    ds_expected["a"].attrs = {"attr": "x_coord"}
    assert_identical(ds_expected, ds_actual)

    # 2 DataArrays and 1 Dataset, takes attrs from x
    ds_actual = xr.where(cond, x.rename("x"), ds_y, keep_attrs=True)
    ds_expected = xr.Dataset(
        data_vars={
            "x": xr.DataArray([1, 0], coords={"a": [0, 1]}, attrs={"attr": "x_da"})
        },
    )
    ds_expected["a"].attrs = {"attr": "x_coord"}
    assert_identical(ds_expected, ds_actual)


@pytest.mark.parametrize(
    "use_dask", [pytest.param(False, id="nodask"), pytest.param(True, id="dask")]
)
@pytest.mark.parametrize(
    ["x", "coeffs", "expected"],
    [
        pytest.param(
            xr.DataArray([1, 2, 3], dims="x"),
            xr.DataArray([2, 3, 4], dims="degree", coords={"degree": [0, 1, 2]}),
            xr.DataArray([9, 2 + 6 + 16, 2 + 9 + 36], dims="x"),
            id="simple",
        ),
        pytest.param(
            xr.DataArray([1, 2, 3], dims="x"),
            xr.DataArray(
                [[0, 1], [0, 1]], dims=("y", "degree"), coords={"degree": [0, 1]}
            ),
            xr.DataArray([[1, 1], [2, 2], [3, 3]], dims=("x", "y")),
            id="broadcast-x",
        ),
        pytest.param(
            xr.DataArray([1, 2, 3], dims="x"),
            xr.DataArray(
                [[0, 1], [1, 0], [1, 1]],
                dims=("x", "degree"),
                coords={"degree": [0, 1]},
            ),
            xr.DataArray([1, 1, 1 + 3], dims="x"),
            id="shared-dim",
        ),
        pytest.param(
            xr.DataArray([1, 2, 3], dims="x"),
            xr.DataArray([1, 0, 0], dims="degree", coords={"degree": [2, 1, 0]}),
            xr.DataArray([1, 2**2, 3**2], dims="x"),
            id="reordered-index",
        ),
        pytest.param(
            xr.DataArray([1, 2, 3], dims="x"),
            xr.DataArray([5], dims="degree", coords={"degree": [3]}),
            xr.DataArray([5, 5 * 2**3, 5 * 3**3], dims="x"),
            id="sparse-index",
        ),
        pytest.param(
            xr.DataArray([1, 2, 3], dims="x"),
            xr.Dataset(
                {"a": ("degree", [0, 1]), "b": ("degree", [1, 0])},
                coords={"degree": [0, 1]},
            ),
            xr.Dataset({"a": ("x", [1, 2, 3]), "b": ("x", [1, 1, 1])}),
            id="array-dataset",
        ),
        pytest.param(
            xr.Dataset({"a": ("x", [1, 2, 3]), "b": ("x", [2, 3, 4])}),
            xr.DataArray([1, 1], dims="degree", coords={"degree": [0, 1]}),
            xr.Dataset({"a": ("x", [2, 3, 4]), "b": ("x", [3, 4, 5])}),
            id="dataset-array",
        ),
        pytest.param(
            xr.Dataset({"a": ("x", [1, 2, 3]), "b": ("y", [2, 3, 4])}),
            xr.Dataset(
                {"a": ("degree", [0, 1]), "b": ("degree", [1, 1])},
                coords={"degree": [0, 1]},
            ),
            xr.Dataset({"a": ("x", [1, 2, 3]), "b": ("y", [3, 4, 5])}),
            id="dataset-dataset",
        ),
        pytest.param(
            xr.DataArray(pd.date_range("1970-01-01", freq="s", periods=3), dims="x"),
            xr.DataArray([0, 1], dims="degree", coords={"degree": [0, 1]}),
            xr.DataArray(
                [0, 1e9, 2e9],
                dims="x",
                coords={"x": pd.date_range("1970-01-01", freq="s", periods=3)},
            ),
            id="datetime",
        ),
        pytest.param(
            # Force a non-ns unit for the coordinate, make sure we convert to `ns`
            # for backwards compatibility at the moment. This can be relaxed in the future.
            xr.DataArray(
                pd.date_range("1970-01-01", freq="s", periods=3, unit="s"), dims="x"
            ),
            xr.DataArray([0, 1], dims="degree", coords={"degree": [0, 1]}),
            xr.DataArray(
                [0, 1e9, 2e9],
                dims="x",
                coords={"x": pd.date_range("1970-01-01", freq="s", periods=3)},
            ),
            id="datetime-non-ns",
        ),
        pytest.param(
            xr.DataArray(
                np.array([1000, 2000, 3000], dtype="timedelta64[ns]"), dims="x"
            ),
            xr.DataArray([0, 1], dims="degree", coords={"degree": [0, 1]}),
            xr.DataArray([1000.0, 2000.0, 3000.0], dims="x"),
            id="timedelta",
        ),
        pytest.param(
            xr.DataArray([1, 2, 3], dims="x"),
            xr.DataArray(
                [2, 3, 4],
                dims="degree",
                coords={"degree": np.array([0, 1, 2], dtype=np.int64)},
            ),
            xr.DataArray([9, 2 + 6 + 16, 2 + 9 + 36], dims="x"),
            id="int64-degree",
        ),
        pytest.param(
            xr.DataArray([1, 2, 3], dims="x"),
            xr.DataArray(
                [2, 3, 4],
                dims="degree",
                coords={"degree": np.array([0, 1, 2], dtype=np.int32)},
            ),
            xr.DataArray([9, 2 + 6 + 16, 2 + 9 + 36], dims="x"),
            id="int32-degree",
        ),
        pytest.param(
            xr.DataArray([1, 2, 3], dims="x"),
            xr.DataArray(
                [2, 3, 4],
                dims="degree",
                coords={"degree": np.array([0, 1, 2], dtype=np.uint8)},
            ),
            xr.DataArray([9, 2 + 6 + 16, 2 + 9 + 36], dims="x"),
            id="uint8-degree",
        ),
    ],
)
def test_polyval(
    use_dask: bool,
    x: xr.DataArray | xr.Dataset,
    coeffs: xr.DataArray | xr.Dataset,
    expected: xr.DataArray | xr.Dataset,
) -> None:
    if use_dask:
        if not has_dask:
            pytest.skip("requires dask")
        coeffs = coeffs.chunk({"degree": 2})
        x = x.chunk({"x": 2})

    with raise_if_dask_computes():
        actual = xr.polyval(coord=x, coeffs=coeffs)

    xr.testing.assert_allclose(actual, expected)


@requires_cftime
@pytest.mark.parametrize(
    "use_dask", [pytest.param(False, id="nodask"), pytest.param(True, id="dask")]
)
@pytest.mark.parametrize("date", ["1970-01-01", "0753-04-21"])
def test_polyval_cftime(use_dask: bool, date: str) -> None:
    import cftime

    x = xr.DataArray(
        xr.date_range(date, freq="1s", periods=3, use_cftime=True),
        dims="x",
    )
    coeffs = xr.DataArray([0, 1], dims="degree", coords={"degree": [0, 1]})

    if use_dask:
        if not has_dask:
            pytest.skip("requires dask")
        coeffs = coeffs.chunk({"degree": 2})
        x = x.chunk({"x": 2})

    with raise_if_dask_computes(max_computes=1):
        actual = xr.polyval(coord=x, coeffs=coeffs)

    t0 = xr.date_range(date, periods=1)[0]
    offset = (t0 - cftime.DatetimeGregorian(1970, 1, 1)).total_seconds() * 1e9
    expected = (
        xr.DataArray(
            [0, 1e9, 2e9],
            dims="x",
            coords={"x": xr.date_range(date, freq="1s", periods=3, use_cftime=True)},
        )
        + offset
    )
    xr.testing.assert_allclose(actual, expected)


def test_polyval_degree_dim_checks() -> None:
    x = xr.DataArray([1, 2, 3], dims="x")
    coeffs = xr.DataArray([2, 3, 4], dims="degree", coords={"degree": [0, 1, 2]})
    with pytest.raises(ValueError):
        xr.polyval(x, coeffs.drop_vars("degree"))
    with pytest.raises(ValueError):
        xr.polyval(x, coeffs.assign_coords(degree=coeffs.degree.astype(float)))


@pytest.mark.parametrize(
    "use_dask", [pytest.param(False, id="nodask"), pytest.param(True, id="dask")]
)
@pytest.mark.parametrize(
    "x",
    [
        pytest.param(xr.DataArray([0, 1, 2], dims="x"), id="simple"),
        pytest.param(
            xr.DataArray(pd.date_range("1970-01-01", freq="ns", periods=3), dims="x"),
            id="datetime",
        ),
        # Force a non-ns unit for the coordinate, make sure we convert to `ns` in both polyfit & polval
        # for backwards compatibility at the moment. This can be relaxed in the future.
        pytest.param(
            xr.DataArray(
                pd.date_range("1970-01-01", freq="s", unit="s", periods=3), dims="x"
            ),
            id="datetime-non-ns",
        ),
        pytest.param(
            xr.DataArray(np.array([0, 1, 2], dtype="timedelta64[ns]"), dims="x"),
            id="timedelta",
        ),
    ],
)
@pytest.mark.parametrize(
    "y",
    [
        pytest.param(xr.DataArray([1, 6, 17], dims="x"), id="1D"),
        pytest.param(
            xr.DataArray([[1, 6, 17], [34, 57, 86]], dims=("y", "x")), id="2D"
        ),
    ],
)
def test_polyfit_polyval_integration(
    use_dask: bool, x: xr.DataArray, y: xr.DataArray
) -> None:
    y.coords["x"] = x
    if use_dask:
        if not has_dask:
            pytest.skip("requires dask")
        y = y.chunk({"x": 2})

    fit = y.polyfit(dim="x", deg=2)
    evaluated = xr.polyval(y.x, fit.polyfit_coefficients)
    expected = y.transpose(*evaluated.dims)
    xr.testing.assert_allclose(evaluated.variable, expected.variable)


@pytest.mark.parametrize("use_dask", [False, True])
@pytest.mark.parametrize(
    "a, b, ae, be, dim, axis",
    [
        [
            xr.DataArray([1, 2, 3]),
            xr.DataArray([4, 5, 6]),
            np.array([1, 2, 3]),
            np.array([4, 5, 6]),
            "dim_0",
            -1,
        ],
        [
            xr.DataArray([1, 2]),
            xr.DataArray([4, 5, 6]),
            np.array([1, 2, 0]),
            np.array([4, 5, 6]),
            "dim_0",
            -1,
        ],
        [
            xr.Variable(dims=["dim_0"], data=[1, 2, 3]),
            xr.Variable(dims=["dim_0"], data=[4, 5, 6]),
            np.array([1, 2, 3]),
            np.array([4, 5, 6]),
            "dim_0",
            -1,
        ],
        [
            xr.Variable(dims=["dim_0"], data=[1, 2]),
            xr.Variable(dims=["dim_0"], data=[4, 5, 6]),
            np.array([1, 2, 0]),
            np.array([4, 5, 6]),
            "dim_0",
            -1,
        ],
        [  # Test dim in the middle:
            xr.DataArray(
                np.arange(0, 5 * 3 * 4).reshape((5, 3, 4)),
                dims=["time", "cartesian", "var"],
                coords=dict(
                    time=(["time"], np.arange(0, 5)),
                    cartesian=(["cartesian"], ["x", "y", "z"]),
                    var=(["var"], [1, 1.5, 2, 2.5]),
                ),
            ),
            xr.DataArray(
                np.arange(0, 5 * 3 * 4).reshape((5, 3, 4)) + 1,
                dims=["time", "cartesian", "var"],
                coords=dict(
                    time=(["time"], np.arange(0, 5)),
                    cartesian=(["cartesian"], ["x", "y", "z"]),
                    var=(["var"], [1, 1.5, 2, 2.5]),
                ),
            ),
            np.arange(0, 5 * 3 * 4).reshape((5, 3, 4)),
            np.arange(0, 5 * 3 * 4).reshape((5, 3, 4)) + 1,
            "cartesian",
            1,
        ],
        # Test 1 sized arrays with coords:
        pytest.param(
            xr.DataArray(
                np.array([1]),
                dims=["cartesian"],
                coords=dict(cartesian=(["cartesian"], ["z"])),
            ),
            xr.DataArray(
                np.array([4, 5, 6]),
                dims=["cartesian"],
                coords=dict(cartesian=(["cartesian"], ["x", "y", "z"])),
            ),
            np.array([0, 0, 1]),
            np.array([4, 5, 6]),
            "cartesian",
            -1,
            marks=(pytest.mark.xfail(),),
        ),
        # Test filling in between with coords:
        pytest.param(
            xr.DataArray(
                [1, 2],
                dims=["cartesian"],
                coords=dict(cartesian=(["cartesian"], ["x", "z"])),
            ),
            xr.DataArray(
                [4, 5, 6],
                dims=["cartesian"],
                coords=dict(cartesian=(["cartesian"], ["x", "y", "z"])),
            ),
            np.array([1, 0, 2]),
            np.array([4, 5, 6]),
            "cartesian",
            -1,
            marks=(pytest.mark.xfail(),),
        ),
    ],
)
def test_cross(a, b, ae, be, dim: str, axis: int, use_dask: bool) -> None:
    expected = np.cross(ae, be, axis=axis)

    if use_dask:
        if not has_dask:
            pytest.skip("test for dask.")
        a = a.chunk()
        b = b.chunk()

    actual = xr.cross(a, b, dim=dim)
    xr.testing.assert_duckarray_allclose(expected, actual)


@pytest.mark.parametrize("compute_backend", ["numbagg"], indirect=True)
def test_complex_number_reduce(compute_backend):
    da = xr.DataArray(np.ones((2,), dtype=np.complex64), dims=["x"])
    # Check that xarray doesn't call into numbagg, which doesn't compile for complex
    # numbers at the moment (but will when numba supports dynamic compilation)
    da.min()


def test_fix() -> None:
    val = 3.0
    val_fixed = np.fix(val)

    da = xr.DataArray([val])
    expected = xr.DataArray([val_fixed])

    actual = np.fix(da)
    assert_identical(expected, actual)
