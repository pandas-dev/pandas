from __future__ import annotations

import pytest

import xarray
from xarray import concat, merge
from xarray.backends.file_manager import FILE_CACHE
from xarray.core.options import OPTIONS, _get_keep_attrs
from xarray.tests.test_dataset import create_test_data


def test_invalid_option_raises() -> None:
    with pytest.raises(ValueError):
        xarray.set_options(not_a_valid_options=True)


def test_display_width() -> None:
    with pytest.raises(ValueError):
        xarray.set_options(display_width=0)
    with pytest.raises(ValueError):
        xarray.set_options(display_width=-10)
    with pytest.raises(ValueError):
        xarray.set_options(display_width=3.5)


def test_arithmetic_join() -> None:
    with pytest.raises(ValueError):
        xarray.set_options(arithmetic_join="invalid")
    with xarray.set_options(arithmetic_join="exact"):
        assert OPTIONS["arithmetic_join"] == "exact"


def test_enable_cftimeindex() -> None:
    with pytest.raises(ValueError):
        xarray.set_options(enable_cftimeindex=None)
    with pytest.warns(FutureWarning, match="no-op"):
        with xarray.set_options(enable_cftimeindex=True):
            assert OPTIONS["enable_cftimeindex"]


def test_file_cache_maxsize() -> None:
    with pytest.raises(ValueError):
        xarray.set_options(file_cache_maxsize=0)
    original_size = FILE_CACHE.maxsize
    with xarray.set_options(file_cache_maxsize=123):
        assert FILE_CACHE.maxsize == 123
    assert FILE_CACHE.maxsize == original_size


def test_keep_attrs() -> None:
    with pytest.raises(ValueError):
        xarray.set_options(keep_attrs="invalid_str")
    with xarray.set_options(keep_attrs=True):
        assert OPTIONS["keep_attrs"]
    with xarray.set_options(keep_attrs=False):
        assert not OPTIONS["keep_attrs"]
    with xarray.set_options(keep_attrs="default"):
        assert _get_keep_attrs(default=True)
        assert not _get_keep_attrs(default=False)


def test_nested_options() -> None:
    original = OPTIONS["display_width"]
    with xarray.set_options(display_width=1):
        assert OPTIONS["display_width"] == 1
        with xarray.set_options(display_width=2):
            assert OPTIONS["display_width"] == 2
        assert OPTIONS["display_width"] == 1
    assert OPTIONS["display_width"] == original


def test_display_style() -> None:
    original = "html"
    assert OPTIONS["display_style"] == original
    with pytest.raises(ValueError):
        xarray.set_options(display_style="invalid_str")
    with xarray.set_options(display_style="text"):
        assert OPTIONS["display_style"] == "text"
    assert OPTIONS["display_style"] == original


def create_test_dataset_attrs(seed=0):
    ds = create_test_data(seed)
    ds.attrs = {"attr1": 5, "attr2": "history", "attr3": {"nested": "more_info"}}
    return ds


def create_test_dataarray_attrs(seed=0, var="var1"):
    da = create_test_data(seed)[var]
    da.attrs = {"attr1": 5, "attr2": "history", "attr3": {"nested": "more_info"}}
    return da


class TestAttrRetention:
    def test_dataset_attr_retention(self) -> None:
        # Use .mean() for all tests: a typical reduction operation
        ds = create_test_dataset_attrs()
        original_attrs = ds.attrs

        # Test default behaviour
        result = ds.mean()
        assert result.attrs == {}
        with xarray.set_options(keep_attrs="default"):
            result = ds.mean()
            assert result.attrs == {}

        with xarray.set_options(keep_attrs=True):
            result = ds.mean()
            assert result.attrs == original_attrs

        with xarray.set_options(keep_attrs=False):
            result = ds.mean()
            assert result.attrs == {}

    def test_dataarray_attr_retention(self) -> None:
        # Use .mean() for all tests: a typical reduction operation
        da = create_test_dataarray_attrs()
        original_attrs = da.attrs

        # Test default behaviour
        result = da.mean()
        assert result.attrs == {}
        with xarray.set_options(keep_attrs="default"):
            result = da.mean()
            assert result.attrs == {}

        with xarray.set_options(keep_attrs=True):
            result = da.mean()
            assert result.attrs == original_attrs

        with xarray.set_options(keep_attrs=False):
            result = da.mean()
            assert result.attrs == {}

    def test_groupby_attr_retention(self) -> None:
        da = xarray.DataArray([1, 2, 3], [("x", [1, 1, 2])])
        da.attrs = {"attr1": 5, "attr2": "history", "attr3": {"nested": "more_info"}}
        original_attrs = da.attrs

        # Test default behaviour
        result = da.groupby("x").sum(keep_attrs=True)
        assert result.attrs == original_attrs
        with xarray.set_options(keep_attrs="default"):
            result = da.groupby("x").sum(keep_attrs=True)
            assert result.attrs == original_attrs

        with xarray.set_options(keep_attrs=True):
            result1 = da.groupby("x")
            result = result1.sum()
            assert result.attrs == original_attrs

        with xarray.set_options(keep_attrs=False):
            result = da.groupby("x").sum()
            assert result.attrs == {}

    def test_concat_attr_retention(self) -> None:
        ds1 = create_test_dataset_attrs()
        ds2 = create_test_dataset_attrs()
        ds2.attrs = {"wrong": "attributes"}
        original_attrs = ds1.attrs

        # Test default behaviour of keeping the attrs of the first
        # dataset in the supplied list
        # global keep_attrs option current doesn't affect concat
        result = concat([ds1, ds2], dim="dim1")
        assert result.attrs == original_attrs

    def test_merge_attr_retention(self) -> None:
        da1 = create_test_dataarray_attrs(var="var1")
        da2 = create_test_dataarray_attrs(var="var2")
        da2.attrs = {"wrong": "attributes"}
        original_attrs = da1.attrs

        # merge currently discards attrs, and the global keep_attrs
        # option doesn't affect this
        result = merge([da1, da2])
        assert result.attrs == original_attrs

    def test_display_style_text(self) -> None:
        ds = create_test_dataset_attrs()
        with xarray.set_options(display_style="text"):
            text = ds._repr_html_()
            assert text.startswith("<pre>")
            assert "&#x27;nested&#x27;" in text
            assert "&lt;xarray.Dataset&gt;" in text

    def test_display_style_html(self) -> None:
        ds = create_test_dataset_attrs()
        with xarray.set_options(display_style="html"):
            html = ds._repr_html_()
            assert html.startswith("<div>")
            assert "&#x27;nested&#x27;" in html

    def test_display_dataarray_style_text(self) -> None:
        da = create_test_dataarray_attrs()
        with xarray.set_options(display_style="text"):
            text = da._repr_html_()
            assert text.startswith("<pre>")
            assert "&lt;xarray.DataArray &#x27;var1&#x27;" in text

    def test_display_dataarray_style_html(self) -> None:
        da = create_test_dataarray_attrs()
        with xarray.set_options(display_style="html"):
            html = da._repr_html_()
            assert html.startswith("<div>")
            assert "#x27;nested&#x27;" in html


@pytest.mark.parametrize(
    "set_value",
    [("left"), ("exact")],
)
def test_get_options_retention(set_value):
    """Test to check if get_options will return changes made by set_options"""
    with xarray.set_options(arithmetic_join=set_value):
        get_options = xarray.get_options()
        assert get_options["arithmetic_join"] == set_value
