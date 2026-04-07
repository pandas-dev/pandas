import re

import pytest

from pandas import (
    ArrowDtype,
    DataFrame,
    Index,
    Series,
)
import pandas._testing as tm

pa = pytest.importorskip("pyarrow")
pc = pytest.importorskip("pyarrow.compute")


def test_struct_accessor_dtypes():
    ser = Series(
        [],
        dtype=ArrowDtype(
            pa.struct(
                [
                    ("int_col", pa.int64()),
                    ("string_col", pa.string()),
                    (
                        "struct_col",
                        pa.struct(
                            [
                                ("int_col", pa.int64()),
                                ("float_col", pa.float64()),
                            ]
                        ),
                    ),
                ]
            )
        ),
    )
    actual = ser.struct.dtypes
    expected = Series(
        [
            ArrowDtype(pa.int64()),
            ArrowDtype(pa.string()),
            ArrowDtype(
                pa.struct(
                    [
                        ("int_col", pa.int64()),
                        ("float_col", pa.float64()),
                    ]
                )
            ),
        ],
        index=Index(["int_col", "string_col", "struct_col"]),
    )
    tm.assert_series_equal(actual, expected)


def test_struct_accessor_field():
    index = Index([-100, 42, 123])
    ser = Series(
        [
            {"rice": 1.0, "maize": -1, "wheat": "a"},
            {"rice": 2.0, "maize": 0, "wheat": "b"},
            {"rice": 3.0, "maize": 1, "wheat": "c"},
        ],
        dtype=ArrowDtype(
            pa.struct(
                [
                    ("rice", pa.float64()),
                    ("maize", pa.int64()),
                    ("wheat", pa.string()),
                ]
            )
        ),
        index=index,
    )
    by_name = ser.struct.field("maize")
    by_name_expected = Series(
        [-1, 0, 1],
        dtype=ArrowDtype(pa.int64()),
        index=index,
        name="maize",
    )
    tm.assert_series_equal(by_name, by_name_expected)

    by_index = ser.struct.field(2)
    by_index_expected = Series(
        ["a", "b", "c"],
        dtype=ArrowDtype(pa.string()),
        index=index,
        name="wheat",
    )
    tm.assert_series_equal(by_index, by_index_expected)


def test_struct_accessor_field_with_invalid_name_or_index():
    ser = Series([], dtype=ArrowDtype(pa.struct([("field", pa.int64())])))

    with pytest.raises(ValueError, match="name_or_index must be an int, str,"):
        ser.struct.field(1.1)


def test_struct_accessor_explode():
    index = Index([-100, 42, 123])
    ser = Series(
        [
            {"painted": 1, "snapping": {"sea": "green"}},
            {"painted": 2, "snapping": {"sea": "leatherback"}},
            {"painted": 3, "snapping": {"sea": "hawksbill"}},
        ],
        dtype=ArrowDtype(
            pa.struct(
                [
                    ("painted", pa.int64()),
                    ("snapping", pa.struct([("sea", pa.string())])),
                ]
            )
        ),
        index=index,
    )
    actual = ser.struct.explode()
    expected = DataFrame(
        {
            "painted": Series([1, 2, 3], index=index, dtype=ArrowDtype(pa.int64())),
            "snapping": Series(
                [{"sea": "green"}, {"sea": "leatherback"}, {"sea": "hawksbill"}],
                index=index,
                dtype=ArrowDtype(pa.struct([("sea", pa.string())])),
            ),
        },
    )
    tm.assert_frame_equal(actual, expected)


def test_struct_accessor_explode_recursive():
    # Test for issue #64915 - recursive flattening of nested structs
    index = Index([0, 1, 2])
    ser = Series(
        [
            {"version": 1, "project": "pandas", "bar": {"a": [1, 2, 3], "b": 10}},
            {"version": 2, "project": "pandas", "bar": {"a": [1, 2], "b": 20}},
            {"version": 1, "project": "numpy", "bar": {"a": [1], "b": 30}},
        ],
        dtype=ArrowDtype(
            pa.struct(
                [
                    ("version", pa.int64()),
                    ("project", pa.string()),
                    (
                        "bar",
                        pa.struct(
                            [
                                ("a", pa.list_(pa.int64())),
                                ("b", pa.int64()),
                            ]
                        ),
                    ),
                ]
            )
        ),
        index=index,
    )

    # Test recursive=True flattens nested structs
    actual = ser.struct.explode(recursive=True)
    expected = DataFrame(
        {
            "version": Series([1, 2, 1], index=index, dtype=ArrowDtype(pa.int64())),
            "project": Series(
                ["pandas", "pandas", "numpy"], index=index, dtype=ArrowDtype(pa.string())
            ),
            "bar_a": Series(
                [[1, 2, 3], [1, 2], [1]], index=index, dtype=ArrowDtype(pa.list_(pa.int64()))
            ),
            "bar_b": Series([10, 20, 30], index=index, dtype=ArrowDtype(pa.int64())),
        },
    )
    tm.assert_frame_equal(actual, expected)


def test_struct_accessor_explode_recursive_separator():
    # Test custom separator
    index = Index([0, 1])
    ser = Series(
        [
            {"bar": {"a": 1, "b": "x"}},
            {"bar": {"a": 2, "b": "y"}},
        ],
        dtype=ArrowDtype(
            pa.struct(
                [
                    (
                        "bar",
                        pa.struct(
                            [
                                ("a", pa.int64()),
                                ("b", pa.string()),
                            ]
                        ),
                    ),
                ]
            )
        ),
        index=index,
    )

    # Default separator "_"
    actual_default = ser.struct.explode(recursive=True)
    assert list(actual_default.columns) == ["bar_a", "bar_b"]

    # Custom separator "."
    actual_dot = ser.struct.explode(recursive=True, separator=".")
    assert list(actual_dot.columns) == ["bar.a", "bar.b"]

    # Custom separator "/"
    actual_slash = ser.struct.explode(recursive=True, separator="/")
    assert list(actual_slash.columns) == ["bar/a", "bar/b"]


def test_struct_accessor_explode_recursive_deeply_nested():
    # Test deeply nested structs (3 levels)
    index = Index([0, 1])
    ser = Series(
        [
            {"a": {"b": {"c": 1, "d": "x"}}},
            {"a": {"b": {"c": 2, "d": "y"}}},
        ],
        dtype=ArrowDtype(
            pa.struct(
                [
                    (
                        "a",
                        pa.struct(
                            [
                                (
                                    "b",
                                    pa.struct(
                                        [
                                            ("c", pa.int64()),
                                            ("d", pa.string()),
                                        ]
                                    ),
                                ),
                            ]
                        ),
                    ),
                ]
            )
        ),
        index=index,
    )

    actual = ser.struct.explode(recursive=True)
    expected = DataFrame(
        {
            "a_b_c": Series([1, 2], index=index, dtype=ArrowDtype(pa.int64())),
            "a_b_d": Series(["x", "y"], index=index, dtype=ArrowDtype(pa.string())),
        },
    )
    tm.assert_frame_equal(actual, expected)


def test_struct_accessor_explode_recursive_empty():
    # Test recursive explode on empty series
    ser = Series(
        [],
        dtype=ArrowDtype(
            pa.struct(
                [
                    ("a", pa.int64()),
                    (
                        "b",
                        pa.struct([("c", pa.string())]),
                    ),
                ]
            )
        ),
    )

    actual = ser.struct.explode(recursive=True)
    expected = DataFrame(
        {
            "a": Series([], dtype=ArrowDtype(pa.int64())),
            "b_c": Series([], dtype=ArrowDtype(pa.string())),
        },
        index=Index([], dtype="object"),
    )
    tm.assert_frame_equal(actual, expected)


def test_struct_accessor_explode_recursive_bytes_field():
    # Test that bytes field names are properly handled (converted to str)
    index = Index([0, 1])
    # Arrow supports bytes as field names
    arrow_type = pa.struct(
        [
            (b"bytes_field", pa.int64()),  # bytes field name
            (
                "nested",
                pa.struct(
                    [
                        (b"bytes_nested", pa.string()),  # bytes in nested struct
                    ]
                ),
            ),
        ]
    )
    ser = Series(
        [
            {"bytes_field": 1, "nested": {b"bytes_nested": "x"}},
            {"bytes_field": 2, "nested": {b"bytes_nested": "y"}},
        ],
        dtype=ArrowDtype(arrow_type),
        index=index,
    )

    actual = ser.struct.explode(recursive=True)
    # bytes field names should be decoded to str
    assert list(actual.columns) == ["bytes_field", "nested_bytes_nested"]
    expected = DataFrame(
        {
            "bytes_field": Series([1, 2], index=index, dtype=ArrowDtype(pa.int64())),
            "nested_bytes_nested": Series(
                ["x", "y"], index=index, dtype=ArrowDtype(pa.string())
            ),
        },
    )
    tm.assert_frame_equal(actual, expected)


@pytest.mark.parametrize(
    "invalid",
    [
        pytest.param(Series([1, 2, 3], dtype="int64"), id="int64"),
        pytest.param(
            Series(["a", "b", "c"], dtype="string[pyarrow]"), id="string-pyarrow"
        ),
    ],
)
def test_struct_accessor_api_for_invalid(invalid):
    with pytest.raises(
        AttributeError,
        match=re.escape(
            "Can only use the '.struct' accessor with 'struct[pyarrow]' dtype, "
            f"not {invalid.dtype}."
        ),
    ):
        invalid.struct


@pytest.mark.parametrize(
    ["indices", "name"],
    [
        (0, "int_col"),
        ([1, 2], "str_col"),
        (pc.field("int_col"), "int_col"),
        ("int_col", "int_col"),
        (b"string_col", b"string_col"),
        ([b"string_col"], "string_col"),
    ],
)
def test_struct_accessor_field_expanded(indices, name):
    arrow_type = pa.struct(
        [
            ("int_col", pa.int64()),
            (
                "struct_col",
                pa.struct(
                    [
                        ("int_col", pa.int64()),
                        ("float_col", pa.float64()),
                        ("str_col", pa.string()),
                    ]
                ),
            ),
            (b"string_col", pa.string()),
        ]
    )

    data = pa.array([], type=arrow_type)
    ser = Series(data, dtype=ArrowDtype(arrow_type))
    expected = pc.struct_field(data, indices)
    result = ser.struct.field(indices)
    tm.assert_equal(result.array._pa_array.combine_chunks(), expected)
    assert result.name == name
