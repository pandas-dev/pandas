from __future__ import annotations

from contextlib import suppress

import numpy as np
import pytest

from xarray import Variable
from xarray.coding import strings
from xarray.core import indexing
from xarray.tests import (
    IndexerMaker,
    assert_array_equal,
    assert_identical,
    requires_dask,
)

with suppress(ImportError):
    import dask.array as da


def test_vlen_dtype() -> None:
    dtype = strings.create_vlen_dtype(str)
    assert dtype.metadata["element_type"] is str
    assert strings.is_unicode_dtype(dtype)
    assert not strings.is_bytes_dtype(dtype)
    assert strings.check_vlen_dtype(dtype) is str

    dtype = strings.create_vlen_dtype(bytes)
    assert dtype.metadata["element_type"] is bytes
    assert not strings.is_unicode_dtype(dtype)
    assert strings.is_bytes_dtype(dtype)
    assert strings.check_vlen_dtype(dtype) is bytes

    # check h5py variant ("vlen")
    dtype = np.dtype("O", metadata={"vlen": str})  # type: ignore[call-overload,unused-ignore]
    assert strings.check_vlen_dtype(dtype) is str

    assert strings.check_vlen_dtype(np.dtype(object)) is None


@pytest.mark.parametrize("numpy_str_type", (np.str_, np.bytes_))
def test_numpy_subclass_handling(numpy_str_type) -> None:
    with pytest.raises(TypeError, match="unsupported type for vlen_dtype"):
        strings.create_vlen_dtype(numpy_str_type)


def test_EncodedStringCoder_decode() -> None:
    coder = strings.EncodedStringCoder()

    raw_data = np.array([b"abc", "ß∂µ∆".encode()])
    raw = Variable(("x",), raw_data, {"_Encoding": "utf-8"})
    actual = coder.decode(raw)

    expected = Variable(("x",), np.array(["abc", "ß∂µ∆"], dtype=object))
    assert_identical(actual, expected)

    assert_identical(coder.decode(actual[0]), expected[0])


@requires_dask
def test_EncodedStringCoder_decode_dask() -> None:
    coder = strings.EncodedStringCoder()

    raw_data = np.array([b"abc", "ß∂µ∆".encode()])
    raw = Variable(("x",), raw_data, {"_Encoding": "utf-8"}).chunk()
    actual = coder.decode(raw)
    assert isinstance(actual.data, da.Array)

    expected = Variable(("x",), np.array(["abc", "ß∂µ∆"], dtype=object))
    assert_identical(actual, expected)

    actual_indexed = coder.decode(actual[0])
    assert isinstance(actual_indexed.data, da.Array)
    assert_identical(actual_indexed, expected[0])


def test_EncodedStringCoder_encode() -> None:
    dtype = strings.create_vlen_dtype(str)
    raw_data = np.array(["abc", "ß∂µ∆"], dtype=dtype)
    expected_data = np.array([r.encode("utf-8") for r in raw_data], dtype=object)

    coder = strings.EncodedStringCoder(allows_unicode=True)
    raw = Variable(("x",), raw_data, encoding={"dtype": "S1"})
    actual = coder.encode(raw)
    expected = Variable(("x",), expected_data, attrs={"_Encoding": "utf-8"})
    assert_identical(actual, expected)

    raw = Variable(("x",), raw_data)
    assert_identical(coder.encode(raw), raw)

    coder = strings.EncodedStringCoder(allows_unicode=False)
    assert_identical(coder.encode(raw), expected)


@pytest.mark.parametrize(
    "original",
    [
        Variable(("x",), [b"ab", b"cdef"]),
        Variable((), b"ab"),
        Variable(("x",), [b"a", b"b"]),
        Variable((), b"a"),
    ],
)
def test_CharacterArrayCoder_roundtrip(original) -> None:
    coder = strings.CharacterArrayCoder()
    roundtripped = coder.decode(coder.encode(original))
    assert_identical(original, roundtripped)


@pytest.mark.parametrize(
    "data",
    [
        np.array([b"a", b"bc"]),
        np.array([b"a", b"bc"], dtype=strings.create_vlen_dtype(bytes)),
    ],
)
def test_CharacterArrayCoder_encode(data) -> None:
    coder = strings.CharacterArrayCoder()
    raw = Variable(("x",), data)
    actual = coder.encode(raw)
    expected = Variable(("x", "string2"), np.array([[b"a", b""], [b"b", b"c"]]))
    assert_identical(actual, expected)


@pytest.mark.parametrize(
    ["original", "expected_char_dim_name"],
    [
        (Variable(("x",), [b"ab", b"cdef"]), "string4"),
        (Variable(("x",), [b"ab", b"cdef"], encoding={"char_dim_name": "foo"}), "foo"),
    ],
)
def test_CharacterArrayCoder_char_dim_name(original, expected_char_dim_name) -> None:
    coder = strings.CharacterArrayCoder()
    encoded = coder.encode(original)
    roundtripped = coder.decode(encoded)
    assert encoded.dims[-1] == expected_char_dim_name
    assert roundtripped.encoding["char_dim_name"] == expected_char_dim_name
    assert roundtripped.dims[-1] == original.dims[-1]


@pytest.mark.parametrize(
    [
        "original",
        "expected_char_dim_name",
        "expected_char_dim_length",
        "warning_message",
    ],
    [
        (
            Variable(("x",), [b"ab", b"cde"], encoding={"char_dim_name": "foo4"}),
            "foo3",
            3,
            "String dimension naming mismatch",
        ),
        (
            Variable(
                ("x",),
                [b"ab", b"cde"],
                encoding={"original_shape": (2, 4), "char_dim_name": "foo"},
            ),
            "foo3",
            3,
            "String dimension length mismatch",
        ),
    ],
)
def test_CharacterArrayCoder_dim_mismatch_warnings(
    original, expected_char_dim_name, expected_char_dim_length, warning_message
) -> None:
    coder = strings.CharacterArrayCoder()
    with pytest.warns(UserWarning, match=warning_message):
        encoded = coder.encode(original)
    roundtripped = coder.decode(encoded)
    assert encoded.dims[-1] == expected_char_dim_name
    assert encoded.sizes[expected_char_dim_name] == expected_char_dim_length
    assert roundtripped.encoding["char_dim_name"] == expected_char_dim_name
    assert roundtripped.dims[-1] == original.dims[-1]


def test_StackedBytesArray() -> None:
    array = np.array([[b"a", b"b", b"c"], [b"d", b"e", b"f"]], dtype="S")
    actual = strings.StackedBytesArray(array)
    expected = np.array([b"abc", b"def"], dtype="S")
    assert actual.dtype == expected.dtype
    assert actual.shape == expected.shape
    assert actual.size == expected.size
    assert actual.ndim == expected.ndim
    assert len(actual) == len(expected)
    assert_array_equal(expected, actual)

    B = IndexerMaker(indexing.BasicIndexer)
    assert_array_equal(expected[:1], actual[B[:1]])
    with pytest.raises(IndexError):
        actual[B[:, :2]]


def test_StackedBytesArray_scalar() -> None:
    array = np.array([b"a", b"b", b"c"], dtype="S")
    actual = strings.StackedBytesArray(array)

    expected = np.array(b"abc")
    assert actual.dtype == expected.dtype
    assert actual.shape == expected.shape
    assert actual.size == expected.size
    assert actual.ndim == expected.ndim
    with pytest.raises(TypeError):
        len(actual)
    np.testing.assert_array_equal(expected, actual)

    B = IndexerMaker(indexing.BasicIndexer)
    with pytest.raises(IndexError):
        actual[B[:2]]


def test_StackedBytesArray_vectorized_indexing() -> None:
    array = np.array([[b"a", b"b", b"c"], [b"d", b"e", b"f"]], dtype="S")
    stacked = strings.StackedBytesArray(array)
    expected = np.array([[b"abc", b"def"], [b"def", b"abc"]])

    V = IndexerMaker(indexing.VectorizedIndexer)
    indexer = V[np.array([[0, 1], [1, 0]])]
    actual = stacked.vindex[indexer]
    assert_array_equal(actual, expected)


def test_char_to_bytes() -> None:
    array = np.array([[b"a", b"b", b"c"], [b"d", b"e", b"f"]])
    expected = np.array([b"abc", b"def"])
    actual = strings.char_to_bytes(array)
    assert_array_equal(actual, expected)

    expected = np.array([b"ad", b"be", b"cf"])
    actual = strings.char_to_bytes(array.T)  # non-contiguous
    assert_array_equal(actual, expected)


def test_char_to_bytes_ndim_zero() -> None:
    expected = np.array(b"a")
    actual = strings.char_to_bytes(expected)
    assert_array_equal(actual, expected)


def test_char_to_bytes_size_zero() -> None:
    array = np.zeros((3, 0), dtype="S1")
    expected = np.array([b"", b"", b""])
    actual = strings.char_to_bytes(array)
    assert_array_equal(actual, expected)


@requires_dask
def test_char_to_bytes_dask() -> None:
    numpy_array = np.array([[b"a", b"b", b"c"], [b"d", b"e", b"f"]])
    array = da.from_array(numpy_array, ((2,), (3,)))
    expected = np.array([b"abc", b"def"])
    actual = strings.char_to_bytes(array)
    assert isinstance(actual, da.Array)
    assert actual.chunks == ((2,),)
    assert actual.dtype == "S3"
    assert_array_equal(np.array(actual), expected)

    with pytest.raises(ValueError, match=r"stacked dask character array"):
        strings.char_to_bytes(array.rechunk(1))


def test_bytes_to_char() -> None:
    array = np.array([[b"ab", b"cd"], [b"ef", b"gh"]])
    expected = np.array([[[b"a", b"b"], [b"c", b"d"]], [[b"e", b"f"], [b"g", b"h"]]])
    actual = strings.bytes_to_char(array)
    assert_array_equal(actual, expected)

    expected = np.array([[[b"a", b"b"], [b"e", b"f"]], [[b"c", b"d"], [b"g", b"h"]]])
    actual = strings.bytes_to_char(array.T)  # non-contiguous
    assert_array_equal(actual, expected)


@requires_dask
def test_bytes_to_char_dask() -> None:
    numpy_array = np.array([b"ab", b"cd"])
    array = da.from_array(numpy_array, ((1, 1),))
    expected = np.array([[b"a", b"b"], [b"c", b"d"]])
    actual = strings.bytes_to_char(array)
    assert isinstance(actual, da.Array)
    assert actual.chunks == ((1, 1), ((2,)))
    assert actual.dtype == "S1"
    assert_array_equal(np.array(actual), expected)
