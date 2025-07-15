from __future__ import annotations

from collections.abc import Hashable
from types import EllipsisType

import numpy as np
import pandas as pd
import pytest

from xarray.core import duck_array_ops, utils
from xarray.core.utils import (
    attempt_import,
    either_dict_or_kwargs,
    infix_dims,
    iterate_nested,
)
from xarray.tests import assert_array_equal, requires_dask


class TestAlias:
    def test(self):
        def new_method():
            pass

        old_method = utils.alias(new_method, "old_method")
        assert "deprecated" in old_method.__doc__
        with pytest.warns(Warning, match="deprecated"):
            old_method()


@pytest.mark.parametrize(
    ["a", "b", "expected"],
    [
        [np.array(["a"]), np.array(["b"]), np.array(["a", "b"])],
        [np.array([1], dtype="int64"), np.array([2], dtype="int64"), pd.Index([1, 2])],
    ],
)
def test_maybe_coerce_to_str(a, b, expected):
    index = pd.Index(a).append(pd.Index(b))

    actual = utils.maybe_coerce_to_str(index, [a, b])

    assert_array_equal(expected, actual)
    assert expected.dtype == actual.dtype


def test_maybe_coerce_to_str_minimal_str_dtype():
    a = np.array(["a", "a_long_string"])
    index = pd.Index(["a"])

    actual = utils.maybe_coerce_to_str(index, [a])
    expected = np.array("a")

    assert_array_equal(expected, actual)
    assert expected.dtype == actual.dtype


class TestArrayEquiv:
    def test_0d(self):
        # verify our work around for pd.isnull not working for 0-dimensional
        # object arrays
        assert duck_array_ops.array_equiv(0, np.array(0, dtype=object))
        assert duck_array_ops.array_equiv(np.nan, np.array(np.nan, dtype=object))
        assert not duck_array_ops.array_equiv(0, np.array(1, dtype=object))


class TestDictionaries:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.x = {"a": "A", "b": "B"}
        self.y = {"c": "C", "b": "B"}
        self.z = {"a": "Z"}

    def test_equivalent(self):
        assert utils.equivalent(0, 0)
        assert utils.equivalent(np.nan, np.nan)
        assert utils.equivalent(0, np.array(0.0))
        assert utils.equivalent([0], np.array([0]))
        assert utils.equivalent(np.array([0]), [0])
        assert utils.equivalent(np.arange(3), 1.0 * np.arange(3))
        assert not utils.equivalent(0, np.zeros(3))

    def test_safe(self):
        # should not raise exception:
        utils.update_safety_check(self.x, self.y)

    def test_unsafe(self):
        with pytest.raises(ValueError):
            utils.update_safety_check(self.x, self.z)

    def test_compat_dict_intersection(self):
        assert {"b": "B"} == utils.compat_dict_intersection(self.x, self.y)
        assert {} == utils.compat_dict_intersection(self.x, self.z)

    def test_compat_dict_union(self):
        assert {"a": "A", "b": "B", "c": "C"} == utils.compat_dict_union(self.x, self.y)
        with pytest.raises(
            ValueError,
            match=r"unsafe to merge dictionaries without "
            "overriding values; conflicting key",
        ):
            utils.compat_dict_union(self.x, self.z)

    def test_dict_equiv(self):
        x = {}
        x["a"] = 3
        x["b"] = np.array([1, 2, 3])
        y = {}
        y["b"] = np.array([1.0, 2.0, 3.0])
        y["a"] = 3
        assert utils.dict_equiv(x, y)  # two nparrays are equal
        y["b"] = [1, 2, 3]  # np.array not the same as a list
        assert utils.dict_equiv(x, y)  # nparray == list
        x["b"] = [1.0, 2.0, 3.0]
        assert utils.dict_equiv(x, y)  # list vs. list
        x["c"] = None
        assert not utils.dict_equiv(x, y)  # new key in x
        x["c"] = np.nan
        y["c"] = np.nan
        assert utils.dict_equiv(x, y)  # as intended, nan is nan
        x["c"] = np.inf
        y["c"] = np.inf
        assert utils.dict_equiv(x, y)  # inf == inf
        y = dict(y)
        assert utils.dict_equiv(x, y)  # different dictionary types are fine
        y["b"] = 3 * np.arange(3)
        assert not utils.dict_equiv(x, y)  # not equal when arrays differ

    def test_frozen(self):
        x = utils.Frozen(self.x)
        with pytest.raises(TypeError):
            x["foo"] = "bar"
        with pytest.raises(TypeError):
            del x["a"]
        with pytest.raises(AttributeError):
            x.update(self.y)
        assert x.mapping == self.x
        assert repr(x) in (
            "Frozen({'a': 'A', 'b': 'B'})",
            "Frozen({'b': 'B', 'a': 'A'})",
        )

    def test_filtered(self):
        x = utils.FilteredMapping(keys={"a"}, mapping={"a": 1, "b": 2})
        assert "a" in x
        assert "b" not in x
        assert x["a"] == 1
        assert list(x) == ["a"]
        assert len(x) == 1
        assert repr(x) == "FilteredMapping(keys={'a'}, mapping={'a': 1, 'b': 2})"
        assert dict(x) == {"a": 1}


def test_repr_object():
    obj = utils.ReprObject("foo")
    assert repr(obj) == "foo"
    assert isinstance(obj, Hashable)
    assert not isinstance(obj, str)


def test_repr_object_magic_methods():
    o1 = utils.ReprObject("foo")
    o2 = utils.ReprObject("foo")
    o3 = utils.ReprObject("bar")
    o4 = "foo"
    assert o1 == o2
    assert o1 != o3
    assert o1 != o4
    assert hash(o1) == hash(o2)
    assert hash(o1) != hash(o3)
    assert hash(o1) != hash(o4)


def test_is_remote_uri():
    assert utils.is_remote_uri("http://example.com")
    assert utils.is_remote_uri("https://example.com")
    assert not utils.is_remote_uri(" http://example.com")
    assert not utils.is_remote_uri("example.nc")


class Test_is_uniform_and_sorted:
    def test_sorted_uniform(self):
        assert utils.is_uniform_spaced(np.arange(5))

    def test_sorted_not_uniform(self):
        assert not utils.is_uniform_spaced([-2, 1, 89])

    def test_not_sorted_uniform(self):
        assert not utils.is_uniform_spaced([1, -1, 3])

    def test_not_sorted_not_uniform(self):
        assert not utils.is_uniform_spaced([4, 1, 89])

    def test_two_numbers(self):
        assert utils.is_uniform_spaced([0, 1.7])

    def test_relative_tolerance(self):
        assert utils.is_uniform_spaced([0, 0.97, 2], rtol=0.1)


class Test_hashable:
    def test_hashable(self):
        for v in [False, 1, (2,), (3, 4), "four"]:
            assert utils.hashable(v)
        for v in [[5, 6], ["seven", "8"], {9: "ten"}]:
            assert not utils.hashable(v)


@requires_dask
def test_dask_array_is_scalar():
    # regression test for GH1684
    import dask.array as da

    y = da.arange(8, chunks=4)
    assert not utils.is_scalar(y)


def test_hidden_key_dict():
    hidden_key = "_hidden_key"
    data = {"a": 1, "b": 2, hidden_key: 3}
    data_expected = {"a": 1, "b": 2}
    hkd = utils.HiddenKeyDict(data, [hidden_key])
    assert len(hkd) == 2
    assert hidden_key not in hkd
    for k, v in data_expected.items():
        assert hkd[k] == v
    with pytest.raises(KeyError):
        hkd[hidden_key]
    with pytest.raises(KeyError):
        del hkd[hidden_key]


def test_either_dict_or_kwargs():
    result = either_dict_or_kwargs(dict(a=1), None, "foo")
    expected = dict(a=1)
    assert result == expected

    result = either_dict_or_kwargs(None, dict(a=1), "foo")
    expected = dict(a=1)
    assert result == expected

    with pytest.raises(ValueError, match=r"foo"):
        result = either_dict_or_kwargs(dict(a=1), dict(a=1), "foo")


@pytest.mark.parametrize(
    ["supplied", "all_", "expected"],
    [
        (list("abc"), list("abc"), list("abc")),
        (["a", ..., "c"], list("abc"), list("abc")),
        (["a", ...], list("abc"), list("abc")),
        (["c", ...], list("abc"), list("cab")),
        ([..., "b"], list("abc"), list("acb")),
        ([...], list("abc"), list("abc")),
    ],
)
def test_infix_dims(supplied, all_, expected):
    result = list(infix_dims(supplied, all_))
    assert result == expected


@pytest.mark.parametrize(
    ["supplied", "all_"], [([..., ...], list("abc")), ([...], list("aac"))]
)
def test_infix_dims_errors(supplied, all_):
    with pytest.raises(ValueError):
        list(infix_dims(supplied, all_))


@pytest.mark.parametrize(
    ["dim", "expected"],
    [
        pytest.param("a", ("a",), id="str"),
        pytest.param(["a", "b"], ("a", "b"), id="list_of_str"),
        pytest.param(["a", 1], ("a", 1), id="list_mixed"),
        pytest.param(["a", ...], ("a", ...), id="list_with_ellipsis"),
        pytest.param(("a", "b"), ("a", "b"), id="tuple_of_str"),
        pytest.param(["a", ("b", "c")], ("a", ("b", "c")), id="list_with_tuple"),
        pytest.param((("b", "c"),), (("b", "c"),), id="tuple_of_tuple"),
        pytest.param({"a", 1}, tuple({"a", 1}), id="non_sequence_collection"),
        pytest.param((), (), id="empty_tuple"),
        pytest.param(set(), (), id="empty_collection"),
        pytest.param(None, None, id="None"),
        pytest.param(..., ..., id="ellipsis"),
    ],
)
def test_parse_dims_as_tuple(dim, expected) -> None:
    all_dims = ("a", "b", 1, ("b", "c"))  # selection of different Hashables
    actual = utils.parse_dims_as_tuple(dim, all_dims, replace_none=False)
    assert actual == expected


def test_parse_dims_set() -> None:
    all_dims = ("a", "b", 1, ("b", "c"))  # selection of different Hashables
    dim = {"a", 1}
    actual = utils.parse_dims_as_tuple(dim, all_dims)
    assert set(actual) == dim


@pytest.mark.parametrize(
    "dim", [pytest.param(None, id="None"), pytest.param(..., id="ellipsis")]
)
def test_parse_dims_replace_none(dim: EllipsisType | None) -> None:
    all_dims = ("a", "b", 1, ("b", "c"))  # selection of different Hashables
    actual = utils.parse_dims_as_tuple(dim, all_dims, replace_none=True)
    assert actual == all_dims


@pytest.mark.parametrize(
    "dim",
    [
        pytest.param("x", id="str_missing"),
        pytest.param(["a", "x"], id="list_missing_one"),
        pytest.param(["x", 2], id="list_missing_all"),
    ],
)
def test_parse_dims_raises(dim) -> None:
    all_dims = ("a", "b", 1, ("b", "c"))  # selection of different Hashables
    with pytest.raises(ValueError, match="'x'"):
        utils.parse_dims_as_tuple(dim, all_dims, check_exists=True)


@pytest.mark.parametrize(
    ["dim", "expected"],
    [
        pytest.param("a", ("a",), id="str"),
        pytest.param(["a", "b"], ("a", "b"), id="list"),
        pytest.param([...], ("a", "b", "c"), id="list_only_ellipsis"),
        pytest.param(["a", ...], ("a", "b", "c"), id="list_with_ellipsis"),
        pytest.param(["a", ..., "b"], ("a", "c", "b"), id="list_with_middle_ellipsis"),
    ],
)
def test_parse_ordered_dims(dim, expected) -> None:
    all_dims = ("a", "b", "c")
    actual = utils.parse_ordered_dims(dim, all_dims)
    assert actual == expected


def test_parse_ordered_dims_raises() -> None:
    all_dims = ("a", "b", "c")

    with pytest.raises(ValueError, match="'x' do not exist"):
        utils.parse_ordered_dims("x", all_dims, check_exists=True)

    with pytest.raises(ValueError, match="repeated dims"):
        utils.parse_ordered_dims(["a", ...], all_dims + ("a",))

    with pytest.raises(ValueError, match="More than one ellipsis"):
        utils.parse_ordered_dims(["a", ..., "b", ...], all_dims)


@pytest.mark.parametrize(
    "nested_list, expected",
    [
        ([], []),
        ([1], [1]),
        ([1, 2, 3], [1, 2, 3]),
        ([[1]], [1]),
        ([[1, 2], [3, 4]], [1, 2, 3, 4]),
        ([[[1, 2, 3], [4]], [5, 6]], [1, 2, 3, 4, 5, 6]),
    ],
)
def test_iterate_nested(nested_list, expected):
    assert list(iterate_nested(nested_list)) == expected


def test_find_stack_level():
    assert utils.find_stack_level() == 1
    assert utils.find_stack_level(test_mode=True) == 2

    def f():
        return utils.find_stack_level(test_mode=True)

    assert f() == 3


def test_attempt_import() -> None:
    """Test optional dependency handling."""
    np = attempt_import("numpy")
    assert np.__name__ == "numpy"

    with pytest.raises(ImportError, match="The foo package is required"):
        attempt_import(module="foo")
    with pytest.raises(ImportError, match="The foo package is required"):
        attempt_import(module="foo.bar")
