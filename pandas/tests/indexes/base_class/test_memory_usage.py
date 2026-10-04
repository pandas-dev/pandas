from decimal import Decimal
import sys

import numpy as np
import pytest

from pandas.compat import PYPY

import pandas as pd
import pandas._testing as tm


@pytest.fixture(
    params=["str", "string[pyarrow]", "arrow_string", "large_string[pyarrow]"]
)
def arrow_string_dtype(request):
    pa = pytest.importorskip("pyarrow")
    if request.param == "str":
        return pd.StringDtype(storage="pyarrow", na_value=np.nan)
    if request.param == "arrow_string":
        return pd.ArrowDtype(pa.string())
    return request.param


@pytest.mark.parametrize("deep", [False, True])
@pytest.mark.parametrize("state", ["empty_mapping", "populated", "cleared"])
@pytest.mark.parametrize("values", [[], ["first label", "second label", None]])
def test_memory_usage_arrow_engine_copy(arrow_string_dtype, deep, state, values):
    # GH#66593: the object copy exists before the hash table is populated.
    index = pd.Index(values, dtype=arrow_string_dtype)
    before = index.memory_usage(deep=deep)
    assert "_engine" not in index._cache

    engine = index._engine
    if state != "empty_mapping" and len(index):
        index.get_loc(index[0])
    if state == "cleared":
        engine.clear_mapping()

    expected = before + engine.sizeof(deep=deep) + engine.values.nbytes
    if deep and not PYPY:
        expected += sum(sys.getsizeof(value) for value in engine.values)
    assert index.memory_usage(deep=deep) == expected
    assert index._engine is engine


@pytest.mark.parametrize("dtype", ["object", "string[python]"])
@pytest.mark.parametrize("deep", [False, True])
@pytest.mark.parametrize("populate", [False, True])
def test_memory_usage_engine_shared_values(dtype, deep, populate):
    # GH#66593: the original object array and its values are already counted.
    index = pd.Index(["first label", "second label", None], dtype=dtype)
    before = index.memory_usage(deep=deep)
    assert "_engine" not in index._cache

    engine = index._engine
    values = index._values if dtype == "object" else index._values._ndarray
    assert engine.values is values
    if populate:
        index.get_loc(index[0])

    expected = before + engine.sizeof(deep=deep)
    assert index.memory_usage(deep=deep) == expected


@pytest.mark.parametrize("deep", [False, True])
def test_memory_usage_frame_arrow_index(arrow_string_dtype, deep):
    # GH#66593: DataFrame.memory_usage must include the index's retained copy.
    index = pd.Index(["first label", "second label"], dtype=arrow_string_dtype)
    frame = pd.DataFrame({"value": [1, 2]}, index=index)
    before = frame.memory_usage(index=False, deep=deep).sum()
    before += index.memory_usage(deep=deep)
    frame.loc["first label"]

    engine = index._engine
    expected = before + engine.sizeof(deep=deep) + engine.values.nbytes
    if deep and not PYPY:
        expected += sum(sys.getsizeof(value) for value in engine.values)
    assert frame.memory_usage(deep=deep).sum() == expected


@pytest.mark.parametrize("deep", [False, True])
@pytest.mark.parametrize("kind", ["categorical", "multi"])
def test_memory_usage_nested_arrow_index(arrow_string_dtype, deep, kind):
    # GH#66593: category and level indexes retain their own lookup engines.
    index = pd.Index(["first label", "second label"], dtype=arrow_string_dtype)
    if kind == "categorical":
        obj = pd.CategoricalIndex(
            pd.Categorical.from_codes([0, 1, 0], categories=index)
        )
        retained = obj.categories
    else:
        obj = pd.MultiIndex(levels=[index], codes=[[0, 1, 0]])
        retained = obj.levels[0]
    retained._cache.clear()
    before = obj.memory_usage(deep=deep)
    assert "_engine" not in retained._cache

    retained.get_loc("first label")
    engine = retained._engine
    expected = before + engine.sizeof(deep=deep) + engine.values.nbytes
    if deep and not PYPY:
        expected += sum(sys.getsizeof(value) for value in engine.values)
    assert obj.memory_usage(deep=deep) == expected


@pytest.mark.parametrize("deep", [False, True])
@pytest.mark.parametrize(
    "dtype, values",
    [
        ("bool[pyarrow]", [True, False, None]),
        ("binary[pyarrow]", [b"first label", b"second label", None]),
        ("decimal", [Decimal("1.25"), Decimal("2.50"), None]),
    ],
)
def test_memory_usage_arrow_object_engine(dtype, values, deep):
    # GH#66593: other Arrow types also retain an object copy for lookups.
    pa = pytest.importorskip("pyarrow")
    if dtype == "decimal":
        dtype = pd.ArrowDtype(pa.decimal128(5, 2))
    index = pd.Index(values, dtype=dtype)
    before = index.memory_usage(deep=deep)
    index.get_loc(index[0])

    engine = index._engine
    assert engine.values.dtype == object
    expected = before + engine.sizeof(deep=deep) + engine.values.nbytes
    if deep and not PYPY:
        expected += sum(sys.getsizeof(value) for value in engine.values)
    assert index.memory_usage(deep=deep) == expected


@pytest.mark.parametrize("deep", [False, True])
@pytest.mark.parametrize(
    "dtype, values",
    [
        ("int64[pyarrow]", [1, 2, None]),
        ("float64[pyarrow]", [1.0, 2.0, None]),
        ("timestamp[ns][pyarrow]", [pd.Timestamp("2020-01-01"), None]),
        ("duration[ns][pyarrow]", [pd.Timedelta(1), None]),
    ],
)
def test_memory_usage_arrow_non_object_engine(dtype, values, deep):
    # GH#66593: the new accounting only applies to object-dtype engine arrays.
    pytest.importorskip("pyarrow")
    index = pd.Index(values, dtype=dtype)
    before = index.memory_usage(deep=deep)
    index.get_loc(index[0])

    engine = index._engine
    assert engine.values.dtype != object
    assert index.memory_usage(deep=deep) == before + engine.sizeof(deep=deep)


@pytest.mark.parametrize("deep", [False, True])
def test_memory_usage_arrow_engine_repeated_calls(arrow_string_dtype, deep):
    # GH#66593: repeated measurements and mapping resets must not accumulate bytes.
    index = pd.Index(["first label", "second label"], dtype=arrow_string_dtype)
    before = index.memory_usage(deep=deep)
    engine = index._engine
    values = engine.values
    retained = values.nbytes
    if deep and not PYPY:
        retained += sum(sys.getsizeof(value) for value in values)

    for _ in range(3):
        engine.clear_mapping()
        assert engine.values is values
        assert not engine.is_mapping_populated
        assert index.memory_usage(deep=deep) == before + retained
        assert index.memory_usage(deep=deep) == before + retained
        assert not engine.is_mapping_populated

        index.get_loc("first label")
        assert index._engine is engine
        assert engine.values is values
        expected = before + retained + engine.sizeof(deep=deep)
        assert index.memory_usage(deep=deep) == expected
        assert index.memory_usage(deep=deep) == expected


@pytest.mark.parametrize("deep", [False, True])
@pytest.mark.parametrize("operation", ["view", "shallow_copy", "rename"])
def test_memory_usage_arrow_shared_engine(arrow_string_dtype, deep, operation):
    # GH#66593: views report their retained engine once per index, not globally.
    index = pd.Index(["first label", "second label"], dtype=arrow_string_dtype)
    index.get_loc("first label")
    expected = index.memory_usage(deep=deep)

    if operation == "view":
        result = index.view()
    elif operation == "shallow_copy":
        result = index.copy(deep=False)
    else:
        result = index.rename("renamed")

    assert result is not index
    assert result._engine is index._engine
    assert result.memory_usage(deep=deep) == expected
    assert result.memory_usage(deep=deep) == expected
    assert index.memory_usage(deep=deep) == expected


@pytest.mark.parametrize("deep", [False, True])
@pytest.mark.parametrize("operation", ["take", "slice", "insert", "deep_copy"])
def test_memory_usage_arrow_replaced_values(arrow_string_dtype, deep, operation):
    # GH#66593: new values must not inherit an engine for the old labels.
    index = pd.Index(
        ["first label", "second label", "third label"], dtype=arrow_string_dtype
    )
    index.get_loc("first label")
    original_engine = index._engine
    original_usage = index.memory_usage(deep=deep)

    if operation == "take":
        result = index.take([2, 0])
    elif operation == "slice":
        result = index[1:]
    elif operation == "insert":
        result = index.insert(1, "inserted label")
    else:
        result = index.copy(deep=True)

    # Slicing propagates cached properties into a new engine immediately.
    assert ("_engine" in result._cache) == (operation == "slice")
    before = result._memory_usage(deep=deep)
    result.get_loc(result[0])
    engine = result._engine
    assert engine is not original_engine
    tm.assert_numpy_array_equal(engine.values, result.to_numpy(dtype=object))
    expected = before + engine.sizeof(deep=deep) + engine.values.nbytes
    if deep and not PYPY:
        expected += sum(sys.getsizeof(value) for value in engine.values)
    assert result.memory_usage(deep=deep) == expected
    assert result.memory_usage(deep=deep) == expected
    assert index._engine is original_engine
    assert index.memory_usage(deep=deep) == original_usage
