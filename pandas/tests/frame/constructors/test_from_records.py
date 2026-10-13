from collections.abc import Iterator
from datetime import (
    UTC,
    datetime,
)
from decimal import Decimal

import numpy as np
import pytest

from pandas._config import using_string_dtype

from pandas.compat import is_platform_little_endian
from pandas.errors import Pandas4Warning

import pandas as pd
import pandas._testing as tm


class TestFromRecords:
    def test_from_records_dt64tz_frame(self):
        # GH#51697
        df = pd.DataFrame({"a": [1, 2, 3]})
        with pytest.raises(TypeError, match="not supported"):
            pd.DataFrame.from_records(df)

    def test_from_records_with_datetimes(self):
        # this may fail on certain platforms because of a numpy issue
        # related GH#6140
        if not is_platform_little_endian():
            pytest.skip("known failure of test on non-little endian")

        # construction with a null in a recarray
        # GH#6140
        expected = pd.DataFrame({"EXPIRY": [datetime(2005, 3, 1, 0, 0), None]})

        arrdata = [np.array([datetime(2005, 3, 1, 0, 0), None])]
        dtypes = [("EXPIRY", "<M8[us]")]

        recarray = np.rec.fromarrays(arrdata, dtype=dtypes)

        result = pd.DataFrame.from_records(recarray)
        tm.assert_frame_equal(result, expected)

        # coercion should work too
        arrdata = [np.array([datetime(2005, 3, 1, 0, 0), None])]
        dtypes = [("EXPIRY", "<M8[m]")]
        recarray = np.rec.fromarrays(arrdata, dtype=dtypes)
        result = pd.DataFrame.from_records(recarray)
        # we get the closest supported unit, "s"
        expected["EXPIRY"] = expected["EXPIRY"].astype("M8[s]")
        tm.assert_frame_equal(result, expected)

    @pytest.mark.xfail(using_string_dtype(), reason="dtype checking logic doesn't work")
    def test_from_records_sequencelike(self):
        df = pd.DataFrame(
            {
                "A": np.array(
                    np.random.default_rng(2).standard_normal(6), dtype=np.float64
                ),
                "A1": np.array(
                    np.random.default_rng(2).standard_normal(6), dtype=np.float64
                ),
                "B": np.array(np.arange(6), dtype=np.int64),
                "C": ["foo"] * 6,
                "D": np.array([True, False] * 3, dtype=bool),
                "E": np.array(
                    np.random.default_rng(2).standard_normal(6), dtype=np.float32
                ),
                "E1": np.array(
                    np.random.default_rng(2).standard_normal(6), dtype=np.float32
                ),
                "F": np.array(np.arange(6), dtype=np.int32),
            }
        )

        # this is actually tricky to create the recordlike arrays and
        # have the dtypes be intact
        blocks = df._to_dict_of_blocks()
        tuples = []
        columns = []
        dtypes = []
        for dtype, b in blocks.items():
            columns.extend(b.columns)
            dtypes.extend([(c, np.dtype(dtype).descr[0][1]) for c in b.columns])
        for i in range(len(df.index)):
            tup = []
            for _, b in blocks.items():
                tup.extend(b.iloc[i].values)
            tuples.append(tuple(tup))

        recarray = np.array(tuples, dtype=dtypes).view(np.rec.recarray)
        recarray2 = df.to_records()
        lists = [list(x) for x in tuples]

        # tuples (lose the dtype info)
        result = pd.DataFrame.from_records(tuples, columns=columns).reindex(
            columns=df.columns
        )

        # created recarray and with to_records recarray (have dtype info)
        result2 = pd.DataFrame.from_records(recarray, columns=columns).reindex(
            columns=df.columns
        )
        result3 = pd.DataFrame.from_records(recarray2, columns=columns).reindex(
            columns=df.columns
        )

        # list of tuples (no dtype info)
        result4 = pd.DataFrame.from_records(lists, columns=columns).reindex(
            columns=df.columns
        )

        tm.assert_frame_equal(result, df, check_dtype=False)
        tm.assert_frame_equal(result2, df)
        tm.assert_frame_equal(result3, df)
        tm.assert_frame_equal(result4, df, check_dtype=False)

        # tuples is in the order of the columns
        result = pd.DataFrame.from_records(tuples)
        tm.assert_index_equal(result.columns, pd.RangeIndex(8))

        # test exclude parameter & we are casting the results here (as we don't
        # have dtype info to recover)
        columns_to_test = [columns.index("C"), columns.index("E1")]

        exclude = list(set(range(8)) - set(columns_to_test))
        result = pd.DataFrame.from_records(tuples, exclude=exclude)
        result.columns = [columns[i] for i in sorted(columns_to_test)]
        tm.assert_series_equal(result["C"], df["C"])
        tm.assert_series_equal(result["E1"], df["E1"])

    def test_from_records_sequencelike_empty(self):
        # empty case
        result = pd.DataFrame.from_records([], columns=["foo", "bar", "baz"])
        assert len(result) == 0
        tm.assert_index_equal(result.columns, pd.Index(["foo", "bar", "baz"]))

        result = pd.DataFrame.from_records([])
        assert len(result) == 0
        assert len(result.columns) == 0

    def test_from_records_sequencelike_empty_index(self):
        result = pd.DataFrame.from_records([], index=[])
        assert len(result) == 0
        assert len(result.columns) == 0
        assert len(result.index) == 0

    def test_from_records_dictlike(self):
        # test the dict methods
        df = pd.DataFrame(
            {
                "A": np.array(
                    np.random.default_rng(2).standard_normal(6), dtype=np.float64
                ),
                "A1": np.array(
                    np.random.default_rng(2).standard_normal(6), dtype=np.float64
                ),
                "B": np.array(np.arange(6), dtype=np.int64),
                "C": ["foo"] * 6,
                "D": np.array([True, False] * 3, dtype=bool),
                "E": np.array(
                    np.random.default_rng(2).standard_normal(6), dtype=np.float32
                ),
                "E1": np.array(
                    np.random.default_rng(2).standard_normal(6), dtype=np.float32
                ),
                "F": np.array(np.arange(6), dtype=np.int32),
            }
        )

        # columns is in a different order here than the actual items iterated
        # from the dict
        blocks = df._to_dict_of_blocks()
        columns = []
        for b in blocks.values():
            columns.extend(b.columns)

        asdict = dict(df.items())
        asdict2 = {x: y.values for x, y in df.items()}

        # dict of series & dict of ndarrays (have dtype info)
        msg = "Passing a dict to DataFrame.from_records is deprecated"
        results = []
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            results.append(
                pd.DataFrame.from_records(asdict).reindex(columns=df.columns)
            )
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            results.append(
                pd.DataFrame.from_records(asdict, columns=columns).reindex(
                    columns=df.columns
                )
            )
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            results.append(
                pd.DataFrame.from_records(asdict2, columns=columns).reindex(
                    columns=df.columns
                )
            )

        for r in results:
            tm.assert_frame_equal(r, df)

    def test_from_records_non_tuple(self):
        class Record:
            def __init__(self, *args) -> None:
                self.args = args

            def __getitem__(self, i):
                return self.args[i]

            def __iter__(self) -> Iterator:
                return iter(self.args)

        recs = [Record(1, 2, 3), Record(4, 5, 6), Record(7, 8, 9)]
        tups = [tuple(rec) for rec in recs]

        result = pd.DataFrame.from_records(recs)
        expected = pd.DataFrame.from_records(tups)
        tm.assert_frame_equal(result, expected)

    def test_from_records_len0_with_columns(self):
        # GH#2633
        result = pd.DataFrame.from_records([], index="foo", columns=["foo", "bar"])
        expected = pd.Index(["bar"])

        assert len(result) == 0
        assert result.index.name == "foo"
        tm.assert_index_equal(result.columns, expected)

    def test_from_records_series_list_dict(self):
        # GH#27358
        expected = pd.DataFrame([[{"a": 1, "b": 2}, {"a": 3, "b": 4}]]).T
        data = pd.Series([[{"a": 1, "b": 2}], [{"a": 3, "b": 4}]])
        result = pd.DataFrame.from_records(data)
        tm.assert_frame_equal(result, expected)

    def test_from_records_series_categorical_index(self):
        # GH#32805
        index = pd.CategoricalIndex(
            [pd.Interval(-20, -10), pd.Interval(-10, 0), pd.Interval(0, 10)]
        )
        series_of_dicts = pd.Series([{"a": 1}, {"a": 2}, {"b": 3}], index=index)
        frame = pd.DataFrame.from_records(series_of_dicts, index=index)
        expected = pd.DataFrame(
            {"a": [1, 2, np.nan], "b": [np.nan, np.nan, 3]}, index=index
        )
        tm.assert_frame_equal(frame, expected)

    def test_frame_from_records_utc(self):
        rec = {"datum": 1.5, "begin_time": datetime(2006, 4, 27, tzinfo=UTC)}

        # it works
        pd.DataFrame.from_records([rec], index="begin_time")

    def test_from_records_to_records(self):
        # from numpy documentation
        arr = np.zeros((2,), dtype=("i4,f4,S10"))
        arr[:] = [(1, 2.0, "Hello"), (2, 3.0, "World")]

        pd.DataFrame.from_records(arr)

        index = pd.Index(np.arange(len(arr))[::-1])
        indexed_frame = pd.DataFrame.from_records(arr, index=index)
        tm.assert_index_equal(indexed_frame.index, index)

        # without names, it should go to last ditch
        arr2 = np.zeros((2, 3))
        tm.assert_frame_equal(pd.DataFrame.from_records(arr2), pd.DataFrame(arr2))

        # wrong length
        msg = "|".join(
            [
                r"Length of values \(2\) does not match length of index \(1\)",
            ]
        )
        with pytest.raises(ValueError, match=msg):
            pd.DataFrame.from_records(arr, index=index[:-1])

        indexed_frame = pd.DataFrame.from_records(arr, index="f1")

        # what to do?
        records = indexed_frame.to_records()
        assert len(records.dtype.names) == 3

        records = indexed_frame.to_records(index=False)
        assert len(records.dtype.names) == 2
        assert "index" not in records.dtype.names

    def test_from_records_nones(self):
        tuples = [(1, 2, None, 3), (1, 2, None, 3), (None, 2, 5, 3)]

        df = pd.DataFrame.from_records(tuples, columns=["a", "b", "c", "d"])
        assert np.isnan(df["c"][0])

    def test_from_records_iterator(self):
        arr = np.array(
            [(1.0, 1.0, 2, 2), (3.0, 3.0, 4, 4), (5.0, 5.0, 6, 6), (7.0, 7.0, 8, 8)],
            dtype=[
                ("x", np.float64),
                ("u", np.float32),
                ("y", np.int64),
                ("z", np.int32),
            ],
        )
        df = pd.DataFrame.from_records(iter(arr), nrows=2)
        xp = pd.DataFrame(
            {
                "x": np.array([1.0, 3.0], dtype=np.float64),
                "u": np.array([1.0, 3.0], dtype=np.float32),
                "y": np.array([2, 4], dtype=np.int64),
                "z": np.array([2, 4], dtype=np.int32),
            }
        )
        tm.assert_frame_equal(df.reindex_like(xp), xp)

        # no dtypes specified here, so just compare with the default
        arr = [(1.0, 2), (3.0, 4), (5.0, 6), (7.0, 8)]
        df = pd.DataFrame.from_records(iter(arr), columns=["x", "y"], nrows=2)
        tm.assert_frame_equal(df, xp.reindex(columns=["x", "y"]), check_dtype=False)

    def test_from_records_tuples_generator(self):
        def tuple_generator(length):
            for i in range(length):
                letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                yield (i, letters[i % len(letters)], i / length)

        columns_names = ["Integer", "String", "Float"]
        columns = [
            [i[j] for i in tuple_generator(10)] for j in range(len(columns_names))
        ]
        data = {"Integer": columns[0], "String": columns[1], "Float": columns[2]}
        expected = pd.DataFrame(data, columns=columns_names)

        generator = tuple_generator(10)
        result = pd.DataFrame.from_records(generator, columns=columns_names)
        tm.assert_frame_equal(result, expected)

    def test_from_records_lists_generator(self):
        def list_generator(length):
            for i in range(length):
                letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                yield [i, letters[i % len(letters)], i / length]

        columns_names = ["Integer", "String", "Float"]
        columns = [
            [i[j] for i in list_generator(10)] for j in range(len(columns_names))
        ]
        data = {"Integer": columns[0], "String": columns[1], "Float": columns[2]}
        expected = pd.DataFrame(data, columns=columns_names)

        generator = list_generator(10)
        result = pd.DataFrame.from_records(generator, columns=columns_names)
        tm.assert_frame_equal(result, expected)

    def test_from_records_columns_not_modified(self):
        tuples = [(1, 2, 3), (1, 2, 3), (2, 5, 3)]

        columns = ["a", "b", "c"]
        original_columns = list(columns)

        pd.DataFrame.from_records(tuples, columns=columns, index="a")

        assert columns == original_columns

    def test_from_records_decimal(self):
        tuples = [(Decimal("1.5"),), (Decimal("2.5"),), (None,)]

        df = pd.DataFrame.from_records(tuples, columns=["a"])
        assert df["a"].dtype == object

        df = pd.DataFrame.from_records(tuples, columns=["a"], coerce_float=True)
        assert df["a"].dtype == np.float64
        assert np.isnan(df["a"].values[-1])

    def test_from_records_duplicates(self):
        result = pd.DataFrame.from_records(
            [(1, 2, 3), (4, 5, 6)], columns=["a", "b", "a"]
        )

        expected = pd.DataFrame([(1, 2, 3), (4, 5, 6)], columns=["a", "b", "a"])

        tm.assert_frame_equal(result, expected)

    def test_from_records_set_index_name(self):
        def create_dict(order_id):
            return {
                "order_id": order_id,
                "quantity": np.random.default_rng(2).integers(1, 10),
                "price": np.random.default_rng(2).integers(1, 10),
            }

        documents = [create_dict(i) for i in range(10)]
        # demo missing data
        documents.append({"order_id": 10, "quantity": 5})

        result = pd.DataFrame.from_records(documents, index="order_id")
        assert result.index.name == "order_id"

        # MultiIndex
        result = pd.DataFrame.from_records(documents, index=["order_id", "quantity"])
        assert result.index.names == ("order_id", "quantity")

    def test_from_records_misc_brokenness(self):
        # GH#2179
        msg = "Passing a dict to DataFrame.from_records is deprecated"

        data = {1: ["foo"], 2: ["bar"]}

        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            result = pd.DataFrame.from_records(data, columns=["a", "b"])
        exp = pd.DataFrame(data, columns=["a", "b"])
        tm.assert_frame_equal(result, exp)

        # overlap in index/index_names

        data = {"a": [1, 2, 3], "b": [4, 5, 6]}

        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            result = pd.DataFrame.from_records(data, index=["a", "b", "c"])
        exp = pd.DataFrame(data, index=["a", "b", "c"])
        tm.assert_frame_equal(result, exp)

    def test_from_records_misc_brokenness2(self):
        # GH#2623
        rows = []
        rows.append([datetime(2010, 1, 1), 1])
        rows.append([datetime(2010, 1, 2), "hi"])  # test col upconverts to obj
        result = pd.DataFrame.from_records(rows, columns=["date", "test"])
        expected = pd.DataFrame(
            {"date": [row[0] for row in rows], "test": [row[1] for row in rows]}
        )
        tm.assert_frame_equal(result, expected)
        assert result.dtypes["test"] == np.dtype(object)

    def test_from_records_misc_brokenness3(self):
        rows = []
        rows.append([datetime(2010, 1, 1), 1])
        rows.append([datetime(2010, 1, 2), 1])
        result = pd.DataFrame.from_records(rows, columns=["date", "test"])
        expected = pd.DataFrame(
            {"date": [row[0] for row in rows], "test": [row[1] for row in rows]}
        )
        tm.assert_frame_equal(result, expected)

    def test_from_records_empty(self):
        # GH#3562
        result = pd.DataFrame.from_records([], columns=["a", "b", "c"])
        expected = pd.DataFrame(columns=["a", "b", "c"])
        tm.assert_frame_equal(result, expected)

        result = pd.DataFrame.from_records([], columns=["a", "b", "b"])
        expected = pd.DataFrame(columns=["a", "b", "b"])
        tm.assert_frame_equal(result, expected)

    def test_from_records_empty_with_nonempty_fields_gh3682(self):
        a = np.array([(1, 2)], dtype=[("id", np.int64), ("value", np.int64)])
        df = pd.DataFrame.from_records(a, index="id")

        ex_index = pd.Index([1], name="id")
        expected = pd.DataFrame({"value": [2]}, index=ex_index, columns=["value"])
        tm.assert_frame_equal(df, expected)

        b = a[:0]
        df2 = pd.DataFrame.from_records(b, index="id")
        tm.assert_frame_equal(df2, df.iloc[:0])

    @pytest.mark.parametrize("dtype", [np.int64, np.float64, np.bool_])
    def test_from_records_ndarray_2d_dtype_preserved(self, dtype):
        # GH#22025 - 2D ndarray should slice columns directly
        arr = np.array([[1, 2, 3], [4, 5, 6]], dtype=dtype)
        result = pd.DataFrame.from_records(arr)
        expected = pd.DataFrame(arr)
        tm.assert_frame_equal(result, expected)

    def test_from_records_ndarray_2d_with_columns(self):
        # GH#22025
        arr = np.array([[1, 2, 3], [4, 5, 6]])
        result = pd.DataFrame.from_records(arr, columns=["a", "b", "c"])
        expected = pd.DataFrame(arr, columns=["a", "b", "c"])
        tm.assert_frame_equal(result, expected)

    def test_from_records_ndarray_2d_with_index(self):
        # GH#22025
        arr = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        result = pd.DataFrame.from_records(arr, columns=["a", "b", "c"], index="a")
        expected = pd.DataFrame(
            {"b": arr[:, 1], "c": arr[:, 2]},
            index=pd.Index(arr[:, 0], name="a"),
        )
        tm.assert_frame_equal(result, expected)

    def test_from_records_ndarray_2d_with_exclude(self):
        # GH#22025
        arr = np.array([[1, 2, 3], [4, 5, 6]])
        result = pd.DataFrame.from_records(arr, columns=["a", "b", "c"], exclude=["b"])
        expected = pd.DataFrame(arr, columns=["a", "b", "c"]).drop(columns=["b"])
        tm.assert_frame_equal(result, expected)

    def test_from_records_ndarray_2d_empty(self):
        # GH#22025
        arr = np.empty((0, 3), dtype=np.int64)
        result = pd.DataFrame.from_records(arr, columns=["a", "b", "c"])
        expected = pd.DataFrame(arr, columns=["a", "b", "c"])
        tm.assert_frame_equal(result, expected)

    def test_from_records_ndarray_2d_single_row(self):
        # GH#22025
        arr = np.array([[10, 20, 30]])
        result = pd.DataFrame.from_records(arr, columns=["x", "y", "z"])
        expected = pd.DataFrame(arr, columns=["x", "y", "z"])
        tm.assert_frame_equal(result, expected)

    def test_from_records_ndarray_2d_single_column(self):
        # GH#22025
        arr = np.array([[1], [2], [3]])
        result = pd.DataFrame.from_records(arr, columns=["a"])
        expected = pd.DataFrame(arr, columns=["a"])
        tm.assert_frame_equal(result, expected)

    def test_from_records_empty2(self):
        # GH#42456
        dtype = [("prop", int)]
        shape = (0, len(dtype))
        arr = np.empty(shape, dtype=dtype)

        result = pd.DataFrame.from_records(arr)
        expected = pd.DataFrame({"prop": np.array([], dtype=int)})
        tm.assert_frame_equal(result, expected)

        alt = pd.DataFrame(arr)
        tm.assert_frame_equal(alt, expected)

    def test_from_records_structured_array(self):
        # GH 59717
        data = np.array(
            [
                ("John", 25, "New York", 50000),
                ("Jane", 30, "San Francisco", 75000),
                ("Bob", 35, "Chicago", 65000),
                ("Alice", 28, "Los Angeles", 60000),
            ],
            dtype=[("name", "U10"), ("age", "i4"), ("city", "U15"), ("salary", "i4")],
        )

        actual_result = pd.DataFrame.from_records(
            data, columns=["name", "salary", "city"]
        )

        modified_data = {
            "name": ["John", "Jane", "Bob", "Alice"],
            "salary": np.array([50000, 75000, 65000, 60000], dtype="int32"),
            "city": ["New York", "San Francisco", "Chicago", "Los Angeles"],
        }
        expected_result = pd.DataFrame(modified_data)

        tm.assert_frame_equal(actual_result, expected_result)

    def test_from_records_empty_iterator_with_preserve_columns(self):
        # GH#61140
        rows = []
        result = pd.DataFrame.from_records(
            iter(rows), index=[0, 1], columns=["col_1", "Col_2"], nrows=0
        )
        expected = pd.DataFrame([], index=[0, 1], columns=["col_1", "Col_2"])
        tm.assert_frame_equal(result, expected)

    def test_from_records_empty_iterator_applies_exclude(self):
        # GH#63774
        rows = []
        result = pd.DataFrame.from_records(
            iter(rows), columns=["a", "b", "c"], exclude=["b"], nrows=0
        )
        expected = pd.DataFrame([], columns=["a", "c"])
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("missing_value", [None, np.nan, pd.NA])
    def test_from_records_missing_value_key(self, missing_value, using_infer_string):
        # https://github.com/pandas-dev/pandas/issues/63889
        # preserve values when None key is converted to NaN column name
        dict_data = [
            {"colA": 1, missing_value: 2},
            {"colA": 3, missing_value: 4},
        ]
        result = pd.DataFrame.from_records(dict_data)
        expected = pd.DataFrame(
            [[1, 2], [3, 4]],
            columns=["colA", np.nan if using_infer_string else missing_value],
        )
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("missing_value", [None, np.nan, pd.NA])
    def test_from_records_missing_value_key_only(self, missing_value):
        dict_data = [
            {missing_value: 1},
            {missing_value: 2},
        ]
        result = pd.DataFrame.from_records(dict_data)
        expected = pd.DataFrame([[1], [2]], columns=pd.Index([missing_value]))
        tm.assert_frame_equal(result, expected)
