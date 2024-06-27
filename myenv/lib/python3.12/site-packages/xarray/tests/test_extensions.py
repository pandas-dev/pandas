from __future__ import annotations

import pickle

import pytest

import xarray as xr

# TODO: Remove imports in favour of xr.DataTree etc, once part of public API
from xarray.core.datatree import DataTree
from xarray.core.extensions import register_datatree_accessor
from xarray.tests import assert_identical


@register_datatree_accessor("example_accessor")
@xr.register_dataset_accessor("example_accessor")
@xr.register_dataarray_accessor("example_accessor")
class ExampleAccessor:
    """For the pickling tests below."""

    def __init__(self, xarray_obj):
        self.obj = xarray_obj


class TestAccessor:
    def test_register(self) -> None:
        @register_datatree_accessor("demo")
        @xr.register_dataset_accessor("demo")
        @xr.register_dataarray_accessor("demo")
        class DemoAccessor:
            """Demo accessor."""

            def __init__(self, xarray_obj):
                self._obj = xarray_obj

            @property
            def foo(self):
                return "bar"

        dt: DataTree = DataTree()
        assert dt.demo.foo == "bar"

        ds = xr.Dataset()
        assert ds.demo.foo == "bar"

        da = xr.DataArray(0)
        assert da.demo.foo == "bar"
        # accessor is cached
        assert ds.demo is ds.demo

        # check descriptor
        assert ds.demo.__doc__ == "Demo accessor."
        # TODO: typing doesn't seem to work with accessors
        assert xr.Dataset.demo.__doc__ == "Demo accessor."  # type: ignore
        assert isinstance(ds.demo, DemoAccessor)
        assert xr.Dataset.demo is DemoAccessor  # type: ignore

        # ensure we can remove it
        del xr.Dataset.demo  # type: ignore
        assert not hasattr(xr.Dataset, "demo")

        with pytest.warns(Warning, match="overriding a preexisting attribute"):

            @xr.register_dataarray_accessor("demo")
            class Foo:
                pass

        # it didn't get registered again
        assert not hasattr(xr.Dataset, "demo")

    def test_pickle_dataset(self) -> None:
        ds = xr.Dataset()
        ds_restored = pickle.loads(pickle.dumps(ds))
        assert_identical(ds, ds_restored)

        # state save on the accessor is restored
        assert ds.example_accessor is ds.example_accessor
        ds.example_accessor.value = "foo"
        ds_restored = pickle.loads(pickle.dumps(ds))
        assert_identical(ds, ds_restored)
        assert ds_restored.example_accessor.value == "foo"

    def test_pickle_dataarray(self) -> None:
        array = xr.Dataset()
        assert array.example_accessor is array.example_accessor
        array_restored = pickle.loads(pickle.dumps(array))
        assert_identical(array, array_restored)

    def test_broken_accessor(self) -> None:
        # regression test for GH933

        @xr.register_dataset_accessor("stupid_accessor")
        class BrokenAccessor:
            def __init__(self, xarray_obj):
                raise AttributeError("broken")

        with pytest.raises(RuntimeError, match=r"error initializing"):
            xr.Dataset().stupid_accessor
