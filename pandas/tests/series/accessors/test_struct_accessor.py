import re

import pytest

import pandas as pd
import pandas._testing as tm

pa = pytest.importorskip("pyarrow")
pc = pytest.importorskip("pyarrow.compute")


def test_struct_accessor_dtypes():
    ser = pd.Series(
        [],
        dtype=pd.ArrowDtype(
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
    expected = pd.Series(
        [
            pd.ArrowDtype(pa.int64()),
            pd.ArrowDtype(pa.string()),
            pd.ArrowDtype(
                pa.struct(
                    [
                        ("int_col", pa.int64()),
                        ("float_col", pa.float64()),
                    ]
                )
            ),
        ],
        index=pd.Index(["int_col", "string_col", "struct_col"]),
    )
    tm.assert_series_equal(actual, expected)


def test_struct_accessor_field():
    index = pd.Index([-100, 42, 123])
    ser = pd.Series(
        [
            {"rice": 1.0, "maize": -1, "wheat": "a"},
            {"rice": 2.0, "maize": 0, "wheat": "b"},
            {"rice": 3.0, "maize": 1, "wheat": "c"},
        ],
        dtype=pd.ArrowDtype(
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
    by_name_expected = pd.Series(
        [-1, 0, 1],
        dtype=pd.ArrowDtype(pa.int64()),
        index=index,
        name="maize",
    )
    tm.assert_series_equal(by_name, by_name_expected)

    by_index = ser.struct.field(2)
    by_index_expected = pd.Series(
        ["a", "b", "c"],
        dtype=pd.ArrowDtype(pa.string()),
        index=index,
        name="wheat",
    )
    tm.assert_series_equal(by_index, by_index_expected)


def test_struct_accessor_field_with_invalid_name_or_index():
    ser = pd.Series([], dtype=pd.ArrowDtype(pa.struct([("field", pa.int64())])))

    with pytest.raises(ValueError, match="name_or_index must be an int, str,"):
        ser.struct.field(1.1)


def test_struct_accessor_explode():
    index = pd.Index([-100, 42, 123])
    ser = pd.Series(
        [
            {"painted": 1, "snapping": {"sea": "green"}},
            {"painted": 2, "snapping": {"sea": "leatherback"}},
            {"painted": 3, "snapping": {"sea": "hawksbill"}},
        ],
        dtype=pd.ArrowDtype(
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
    expected = pd.DataFrame(
        {
            "painted": pd.Series(
                [1, 2, 3], index=index, dtype=pd.ArrowDtype(pa.int64())
            ),
            "snapping": pd.Series(
                [{"sea": "green"}, {"sea": "leatherback"}, {"sea": "hawksbill"}],
                index=index,
                dtype=pd.ArrowDtype(pa.struct([("sea", pa.string())])),
            ),
        },
    )
    tm.assert_frame_equal(actual, expected)


@pytest.mark.parametrize(
    "invalid",
    [
        pytest.param(pd.Series([1, 2, 3], dtype="int64"), id="int64"),
        pytest.param(
            pd.Series(["a", "b", "c"], dtype="string[pyarrow]"), id="string-pyarrow"
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
    ser = pd.Series(data, dtype=pd.ArrowDtype(arrow_type))
    expected = pc.struct_field(data, indices)
    result = ser.struct.field(indices)
    tm.assert_equal(result.array._pa_array.combine_chunks(), expected)
    assert result.name == name
