from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import xarray as xr
from xarray.core import formatting_html as fh
from xarray.core.coordinates import Coordinates


@pytest.fixture
def dataarray() -> xr.DataArray:
    return xr.DataArray(np.random.default_rng(0).random((4, 6)))


@pytest.fixture
def dask_dataarray(dataarray: xr.DataArray) -> xr.DataArray:
    pytest.importorskip("dask")
    return dataarray.chunk()


@pytest.fixture
def multiindex() -> xr.Dataset:
    midx = pd.MultiIndex.from_product(
        [["a", "b"], [1, 2]], names=("level_1", "level_2")
    )
    midx_coords = Coordinates.from_pandas_multiindex(midx, "x")
    return xr.Dataset({}, midx_coords)


@pytest.fixture
def dataset() -> xr.Dataset:
    times = pd.date_range("2000-01-01", "2001-12-31", name="time")
    annual_cycle = np.sin(2 * np.pi * (times.dayofyear.values / 365.25 - 0.28))

    base = 10 + 15 * annual_cycle.reshape(-1, 1)
    tmin_values = base + 3 * np.random.randn(annual_cycle.size, 3)
    tmax_values = base + 10 + 3 * np.random.randn(annual_cycle.size, 3)

    return xr.Dataset(
        {
            "tmin": (("time", "location"), tmin_values),
            "tmax": (("time", "location"), tmax_values),
        },
        {"time": times, "location": ["<IA>", "IN", "IL"]},
        attrs={"description": "Test data."},
    )


def test_short_data_repr_html(dataarray: xr.DataArray) -> None:
    data_repr = fh.short_data_repr_html(dataarray)
    assert data_repr.startswith("<pre>array")


def test_short_data_repr_html_non_str_keys(dataset: xr.Dataset) -> None:
    ds = dataset.assign({2: lambda x: x["tmin"]})
    fh.dataset_repr(ds)


def test_short_data_repr_html_dask(dask_dataarray: xr.DataArray) -> None:
    assert hasattr(dask_dataarray.data, "_repr_html_")
    data_repr = fh.short_data_repr_html(dask_dataarray)
    assert data_repr == dask_dataarray.data._repr_html_()


def test_format_dims_no_dims() -> None:
    dims: dict = {}
    dims_with_index: list = []
    formatted = fh.format_dims(dims, dims_with_index)
    assert formatted == ""


def test_format_dims_unsafe_dim_name() -> None:
    dims = {"<x>": 3, "y": 2}
    dims_with_index: list = []
    formatted = fh.format_dims(dims, dims_with_index)
    assert "&lt;x&gt;" in formatted


def test_format_dims_non_index() -> None:
    dims, dims_with_index = {"x": 3, "y": 2}, ["time"]
    formatted = fh.format_dims(dims, dims_with_index)
    assert "class='xr-has-index'" not in formatted


def test_format_dims_index() -> None:
    dims, dims_with_index = {"x": 3, "y": 2}, ["x"]
    formatted = fh.format_dims(dims, dims_with_index)
    assert "class='xr-has-index'" in formatted


def test_summarize_attrs_with_unsafe_attr_name_and_value() -> None:
    attrs = {"<x>": 3, "y": "<pd.DataFrame>"}
    formatted = fh.summarize_attrs(attrs)
    assert "<dt><span>&lt;x&gt; :</span></dt>" in formatted
    assert "<dt><span>y :</span></dt>" in formatted
    assert "<dd>3</dd>" in formatted
    assert "<dd>&lt;pd.DataFrame&gt;</dd>" in formatted


def test_repr_of_dataarray(dataarray: xr.DataArray) -> None:
    formatted = fh.array_repr(dataarray)
    assert "dim_0" in formatted
    # has an expanded data section
    assert formatted.count("class='xr-array-in' type='checkbox' checked>") == 1
    # coords, indexes and attrs don't have an items so they'll be be disabled and collapsed
    assert (
        formatted.count("class='xr-section-summary-in' type='checkbox' disabled >") == 3
    )

    with xr.set_options(display_expand_data=False):
        formatted = fh.array_repr(dataarray)
        assert "dim_0" in formatted
        # has a collapsed data section
        assert formatted.count("class='xr-array-in' type='checkbox' checked>") == 0
        # coords, indexes and attrs don't have an items so they'll be be disabled and collapsed
        assert (
            formatted.count("class='xr-section-summary-in' type='checkbox' disabled >")
            == 3
        )


def test_repr_of_multiindex(multiindex: xr.Dataset) -> None:
    formatted = fh.dataset_repr(multiindex)
    assert "(x)" in formatted


def test_repr_of_dataset(dataset: xr.Dataset) -> None:
    formatted = fh.dataset_repr(dataset)
    # coords, attrs, and data_vars are expanded
    assert (
        formatted.count("class='xr-section-summary-in' type='checkbox'  checked>") == 3
    )
    # indexes is collapsed
    assert formatted.count("class='xr-section-summary-in' type='checkbox'  >") == 1
    assert "&lt;U4" in formatted or "&gt;U4" in formatted
    assert "&lt;IA&gt;" in formatted

    with xr.set_options(
        display_expand_coords=False,
        display_expand_data_vars=False,
        display_expand_attrs=False,
        display_expand_indexes=True,
    ):
        formatted = fh.dataset_repr(dataset)
        # coords, attrs, and data_vars are collapsed, indexes is expanded
        assert (
            formatted.count("class='xr-section-summary-in' type='checkbox'  checked>")
            == 1
        )
        assert "&lt;U4" in formatted or "&gt;U4" in formatted
        assert "&lt;IA&gt;" in formatted


def test_repr_text_fallback(dataset: xr.Dataset) -> None:
    formatted = fh.dataset_repr(dataset)

    # Just test that the "pre" block used for fallback to plain text is present.
    assert "<pre class='xr-text-repr-fallback'>" in formatted


def test_variable_repr_html() -> None:
    v = xr.Variable(["time", "x"], [[1, 2, 3], [4, 5, 6]], {"foo": "bar"})
    assert hasattr(v, "_repr_html_")
    with xr.set_options(display_style="html"):
        html = v._repr_html_().strip()
    # We don't do a complete string identity since
    # html output is probably subject to change, is long and... reasons.
    # Just test that something reasonable was produced.
    assert html.startswith("<div") and html.endswith("</div>")
    assert "xarray.Variable" in html


def test_repr_of_nonstr_dataset(dataset: xr.Dataset) -> None:
    ds = dataset.copy()
    ds.attrs[1] = "Test value"
    ds[2] = ds["tmin"]
    formatted = fh.dataset_repr(ds)
    assert "<dt><span>1 :</span></dt><dd>Test value</dd>" in formatted
    assert "<div class='xr-var-name'><span>2</span>" in formatted


def test_repr_of_nonstr_dataarray(dataarray: xr.DataArray) -> None:
    da = dataarray.rename(dim_0=15)
    da.attrs[1] = "value"
    formatted = fh.array_repr(da)
    assert "<dt><span>1 :</span></dt><dd>value</dd>" in formatted
    assert "<li><span>15</span>: 4</li>" in formatted


def test_nonstr_variable_repr_html() -> None:
    v = xr.Variable(["time", 10], [[1, 2, 3], [4, 5, 6]], {22: "bar"})
    assert hasattr(v, "_repr_html_")
    with xr.set_options(display_style="html"):
        html = v._repr_html_().strip()
    assert "<dt><span>22 :</span></dt><dd>bar</dd>" in html
    assert "<li><span>10</span>: 3</li></ul>" in html


@pytest.fixture(scope="module", params=["some html", "some other html"])
def repr(request):
    return request.param


class Test_summarize_datatree_children:
    """
    Unit tests for summarize_datatree_children.
    """

    func = staticmethod(fh.summarize_datatree_children)

    @pytest.fixture(scope="class")
    def childfree_tree_factory(self):
        """
        Fixture for a child-free DataTree factory.
        """
        from random import randint

        def _childfree_tree_factory():
            return xr.DataTree(
                dataset=xr.Dataset({"z": ("y", [randint(1, 100) for _ in range(3)])})
            )

        return _childfree_tree_factory

    @pytest.fixture(scope="class")
    def childfree_tree(self, childfree_tree_factory):
        """
        Fixture for a child-free DataTree.
        """
        return childfree_tree_factory()

    @pytest.fixture
    def mock_datatree_node_repr(self, monkeypatch):
        """
        Apply mocking for datatree_node_repr.
        """

        def mock(group_title, dt):
            """
            Mock with a simple result
            """
            return group_title + " " + str(id(dt))

        monkeypatch.setattr(fh, "datatree_node_repr", mock)

    @pytest.fixture
    def mock_wrap_datatree_repr(self, monkeypatch):
        """
        Apply mocking for _wrap_datatree_repr.
        """

        def mock(r, *, end, **kwargs):
            """
            Mock by appending "end" or "not end".
            """
            return r + " " + ("end" if end else "not end") + "//"

        monkeypatch.setattr(fh, "_wrap_datatree_repr", mock)

    def test_empty_mapping(self):
        """
        Test with an empty mapping of children.
        """
        children: dict[str, xr.DataTree] = {}
        assert self.func(children) == (
            "<div style='display: inline-grid; grid-template-columns: 100%; grid-column: 1 / -1'>"
            "</div>"
        )

    def test_one_child(
        self, childfree_tree, mock_wrap_datatree_repr, mock_datatree_node_repr
    ):
        """
        Test with one child.

        Uses a mock of _wrap_datatree_repr and _datatree_node_repr to essentially mock
        the inline lambda function "lines_callback".
        """
        # Create mapping of children
        children = {"a": childfree_tree}

        # Expect first line to be produced from the first child, and
        # wrapped as the last child
        first_line = f"a {id(children['a'])} end//"

        assert self.func(children) == (
            "<div style='display: inline-grid; grid-template-columns: 100%; grid-column: 1 / -1'>"
            f"{first_line}"
            "</div>"
        )

    def test_two_children(
        self, childfree_tree_factory, mock_wrap_datatree_repr, mock_datatree_node_repr
    ):
        """
        Test with two level deep children.

        Uses a mock of _wrap_datatree_repr and datatree_node_repr to essentially mock
        the inline lambda function "lines_callback".
        """

        # Create mapping of children
        children = {"a": childfree_tree_factory(), "b": childfree_tree_factory()}

        # Expect first line to be produced from the first child, and
        # wrapped as _not_ the last child
        first_line = f"a {id(children['a'])} not end//"

        # Expect second line to be produced from the second child, and
        # wrapped as the last child
        second_line = f"b {id(children['b'])} end//"

        assert self.func(children) == (
            "<div style='display: inline-grid; grid-template-columns: 100%; grid-column: 1 / -1'>"
            f"{first_line}"
            f"{second_line}"
            "</div>"
        )


class TestDataTreeTruncatesNodes:
    def test_many_nodes(self) -> None:
        # construct a datatree with 500 nodes
        number_of_files = 20
        number_of_groups = 25
        tree_dict = {}
        for f in range(number_of_files):
            for g in range(number_of_groups):
                tree_dict[f"file_{f}/group_{g}"] = xr.Dataset({"g": f * g})

        tree = xr.DataTree.from_dict(tree_dict)
        with xr.set_options(display_style="html"):
            result = tree._repr_html_()

        assert "6/20" in result
        for i in range(number_of_files):
            if i < 3 or i >= (number_of_files - 3):
                assert f"file_{i}</div>" in result
            else:
                assert f"file_{i}</div>" not in result

        assert "6/25" in result
        for i in range(number_of_groups):
            if i < 3 or i >= (number_of_groups - 3):
                assert f"group_{i}</div>" in result
            else:
                assert f"group_{i}</div>" not in result

        with xr.set_options(display_style="html", display_max_children=3):
            result = tree._repr_html_()

        assert "3/20" in result
        for i in range(number_of_files):
            if i < 2 or i >= (number_of_files - 1):
                assert f"file_{i}</div>" in result
            else:
                assert f"file_{i}</div>" not in result

        assert "3/25" in result
        for i in range(number_of_groups):
            if i < 2 or i >= (number_of_groups - 1):
                assert f"group_{i}</div>" in result
            else:
                assert f"group_{i}</div>" not in result


class TestDataTreeInheritance:
    def test_inherited_section_present(self) -> None:
        dt = xr.DataTree.from_dict(
            {
                "/": None,
                "a": None,
            }
        )
        with xr.set_options(display_style="html"):
            html = dt._repr_html_().strip()
        # checks that the section appears somewhere
        assert "Inherited coordinates" in html

        # TODO how can we assert that the Inherited coordinates section does not appear in the child group?
        # with xr.set_options(display_style="html"):
        #     child_html = dt["a"]._repr_html_().strip()
        # assert "Inherited coordinates" not in child_html


class Test__wrap_datatree_repr:
    """
    Unit tests for _wrap_datatree_repr.
    """

    func = staticmethod(fh._wrap_datatree_repr)

    def test_end(self, repr):
        """
        Test with end=True.
        """
        r = self.func(repr, end=True)
        assert r == (
            "<div style='display: inline-grid; grid-template-columns: 0px 20px auto; width: 100%;'>"
            "<div style='"
            "grid-column-start: 1;"
            "border-right: 0.2em solid;"
            "border-color: var(--xr-border-color);"
            "height: 1.2em;"
            "width: 0px;"
            "'>"
            "</div>"
            "<div style='"
            "grid-column-start: 2;"
            "grid-row-start: 1;"
            "height: 1em;"
            "width: 20px;"
            "border-bottom: 0.2em solid;"
            "border-color: var(--xr-border-color);"
            "'>"
            "</div>"
            "<div style='"
            "grid-column-start: 3;"
            "'>"
            f"{repr}"
            "</div>"
            "</div>"
        )

    def test_not_end(self, repr):
        """
        Test with end=False.
        """
        r = self.func(repr, end=False)
        assert r == (
            "<div style='display: inline-grid; grid-template-columns: 0px 20px auto; width: 100%;'>"
            "<div style='"
            "grid-column-start: 1;"
            "border-right: 0.2em solid;"
            "border-color: var(--xr-border-color);"
            "height: 100%;"
            "width: 0px;"
            "'>"
            "</div>"
            "<div style='"
            "grid-column-start: 2;"
            "grid-row-start: 1;"
            "height: 1em;"
            "width: 20px;"
            "border-bottom: 0.2em solid;"
            "border-color: var(--xr-border-color);"
            "'>"
            "</div>"
            "<div style='"
            "grid-column-start: 3;"
            "'>"
            f"{repr}"
            "</div>"
            "</div>"
        )
