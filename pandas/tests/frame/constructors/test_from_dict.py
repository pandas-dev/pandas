from collections import OrderedDict

import numpy as np
import pytest

from pandas.errors import Pandas4Warning

import pandas as pd
import pandas._testing as tm


class TestFromDict:
    # Note: these tests are specific to the from_dict method, not for
    #  passing dictionaries to DataFrame.__init__

    def test_constructor_list_of_odicts(self):
        data = [
            OrderedDict([["a", 1.5], ["b", 3], ["c", 4], ["d", 6]]),
            OrderedDict([["a", 1.5], ["b", 3], ["d", 6]]),
            OrderedDict([["a", 1.5], ["d", 6]]),
            OrderedDict(),
            OrderedDict([["a", 1.5], ["b", 3], ["c", 4]]),
            OrderedDict([["b", 3], ["c", 4], ["d", 6]]),
        ]

        result = pd.DataFrame(data)
        expected = pd.DataFrame.from_dict(
            dict(zip(range(len(data)), data, strict=True)), orient="index"
        )
        tm.assert_frame_equal(result, expected.reindex(result.index))

    def test_constructor_single_row(self):
        data = [OrderedDict([["a", 1.5], ["b", 3], ["c", 4], ["d", 6]])]

        result = pd.DataFrame(data)
        expected = pd.DataFrame.from_dict(
            dict(zip([0], data, strict=True)), orient="index"
        ).reindex(result.index)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("blank", [pd.Series(), {}])
    def test_from_nested_dict_with_empty_entry(self, blank):
        # GH#62775 an empty Series/dict should give an all-NA row, not be dropped
        data = {"good": pd.Series({"a": 1, "b": 2}), "blank": blank}

        result = pd.DataFrame.from_dict(data, orient="index")
        expected = pd.DataFrame(
            {
                "a": {"good": 1.0, "blank": np.nan},
                "b": {"good": 2.0, "blank": np.nan},
            }
        )
        tm.assert_frame_equal(result, expected)

    def test_constructor_list_of_series(self):
        data = [
            OrderedDict([["a", 1.5], ["b", 3.0], ["c", 4.0]]),
            OrderedDict([["a", 1.5], ["b", 3.0], ["c", 6.0]]),
        ]
        sdict = OrderedDict(zip(["x", "y"], data, strict=True))
        idx = pd.Index(["a", "b", "c"])

        # all named
        data2 = [
            pd.Series([1.5, 3, 4], idx, dtype="O", name="x"),
            pd.Series([1.5, 3, 6], idx, name="y"),
        ]
        result = pd.DataFrame(data2)
        expected = pd.DataFrame.from_dict(sdict, orient="index")
        tm.assert_frame_equal(result, expected)

        # some unnamed
        data2 = [
            pd.Series([1.5, 3, 4], idx, dtype="O", name="x"),
            pd.Series([1.5, 3, 6], idx),
        ]
        result = pd.DataFrame(data2)

        sdict = OrderedDict(zip(["x", "Unnamed 0"], data, strict=True))
        expected = pd.DataFrame.from_dict(sdict, orient="index")
        tm.assert_frame_equal(result, expected)

        # none named
        data = [
            OrderedDict([["a", 1.5], ["b", 3], ["c", 4], ["d", 6]]),
            OrderedDict([["a", 1.5], ["b", 3], ["d", 6]]),
            OrderedDict([["a", 1.5], ["d", 6]]),
            OrderedDict(),
            OrderedDict([["a", 1.5], ["b", 3], ["c", 4]]),
            OrderedDict([["b", 3], ["c", 4], ["d", 6]]),
        ]
        data = [pd.Series(d) for d in data]

        result = pd.DataFrame(data)
        sdict = OrderedDict(zip(range(len(data)), data, strict=True))
        expected = pd.DataFrame.from_dict(sdict, orient="index")
        tm.assert_frame_equal(result, expected.reindex(result.index))

        result2 = pd.DataFrame(data, index=np.arange(6, dtype=np.int64))
        tm.assert_frame_equal(result, result2)

        result = pd.DataFrame([pd.Series(dtype=object)])
        expected = pd.DataFrame(index=[0])
        tm.assert_frame_equal(result, expected)

        data = [
            OrderedDict([["a", 1.5], ["b", 3.0], ["c", 4.0]]),
            OrderedDict([["a", 1.5], ["b", 3.0], ["c", 6.0]]),
        ]
        sdict = OrderedDict(zip(range(len(data)), data, strict=True))

        idx = pd.Index(["a", "b", "c"])
        data2 = [pd.Series([1.5, 3, 4], idx, dtype="O"), pd.Series([1.5, 3, 6], idx)]
        result = pd.DataFrame(data2)
        expected = pd.DataFrame.from_dict(sdict, orient="index")
        tm.assert_frame_equal(result, expected)

    def test_constructor_orient(self, float_string_frame):
        data_dict = float_string_frame.T._series
        recons = pd.DataFrame.from_dict(data_dict, orient="index")
        expected = float_string_frame.reindex(index=recons.index)
        tm.assert_frame_equal(recons, expected)

        # dict of sequence
        a = {"hi": [32, 3, 3], "there": [3, 5, 3]}
        rs = pd.DataFrame.from_dict(a, orient="index")
        xp = pd.DataFrame.from_dict(a).T.reindex(list(a.keys()))
        tm.assert_frame_equal(rs, xp)

    def test_constructor_from_ordered_dict(self):
        # GH#8425
        a = OrderedDict(
            [
                ("one", OrderedDict([("col_a", "foo1"), ("col_b", "bar1")])),
                ("two", OrderedDict([("col_a", "foo2"), ("col_b", "bar2")])),
                ("three", OrderedDict([("col_a", "foo3"), ("col_b", "bar3")])),
            ]
        )
        expected = pd.DataFrame.from_dict(a, orient="columns").T
        result = pd.DataFrame.from_dict(a, orient="index")
        tm.assert_frame_equal(result, expected)

    def test_from_dict_columns_parameter(self):
        # GH#18529
        # Test new columns parameter for from_dict that was added to make
        # from_items(..., orient='index', columns=[...]) easier to replicate
        result = pd.DataFrame.from_dict(
            OrderedDict([("A", [1, 2]), ("B", [4, 5])]),
            orient="index",
            columns=["one", "two"],
        )
        expected = pd.DataFrame(
            [[1, 2], [4, 5]], index=["A", "B"], columns=["one", "two"]
        )
        tm.assert_frame_equal(result, expected)

        msg = "cannot use columns parameter with orient='columns'"
        with pytest.raises(ValueError, match=msg):
            pd.DataFrame.from_dict(
                {"A": [1, 2], "B": [4, 5]},
                orient="columns",
                columns=["one", "two"],
            )
        with pytest.raises(ValueError, match=msg):
            pd.DataFrame.from_dict({"A": [1, 2], "B": [4, 5]}, columns=["one", "two"])

    @pytest.mark.parametrize(
        "data_dict, orient, expected",
        [
            ({}, "index", pd.RangeIndex(0)),
            (
                [{("a",): 1}, {("a",): 2}],
                "columns",
                pd.Index([("a",)], tupleize_cols=False),
            ),
            (
                [OrderedDict([(("a",), 1), (("b",), 2)])],
                "columns",
                pd.Index([("a",), ("b",)], tupleize_cols=False),
            ),
            ([{("a", "b"): 1}], "columns", pd.Index([("a", "b")], tupleize_cols=False)),
        ],
    )
    def test_constructor_from_dict_tuples(self, data_dict, orient, expected):
        # GH#16769
        warn = Pandas4Warning if isinstance(data_dict, list) else None
        with tm.assert_produces_warning(warn, match="from_dict is deprecated"):
            df = pd.DataFrame.from_dict(data_dict, orient)
        result = df.columns
        tm.assert_index_equal(result, expected)

    def test_frame_dict_constructor_empty_series(self):
        s1 = pd.Series(
            [1, 2, 3, 4],
            index=pd.MultiIndex.from_tuples([(1, 2), (1, 3), (2, 2), (2, 4)]),
        )
        s2 = pd.Series(
            [1, 2, 3, 4],
            index=pd.MultiIndex.from_tuples([(1, 2), (1, 3), (3, 2), (3, 4)]),
        )
        s3 = pd.Series(dtype=object)

        # it works!
        pd.DataFrame({"foo": s1, "bar": s2, "baz": s3})
        pd.DataFrame.from_dict({"foo": s1, "baz": s3, "bar": s2})

    def test_from_dict_scalars_requires_index(self):
        # GH#25515 the message must be actionable: from_dict has no index
        #  parameter, so it points to orient="index" / the DataFrame constructor
        msg = "If using all scalar values, pass orient='index'"
        with pytest.raises(ValueError, match=msg):
            pd.DataFrame.from_dict({"a": 0.7})

        with pytest.raises(ValueError, match=msg):
            pd.DataFrame.from_dict({"b": 8, "a": 6})

    def test_from_dict_orient_invalid(self):
        msg = (
            "Expected 'index', 'columns' or 'tight' for orient parameter. "
            "Got 'abc' instead"
        )
        with pytest.raises(ValueError, match=msg):
            pd.DataFrame.from_dict({"foo": 1, "baz": 3, "bar": 2}, orient="abc")

    def test_from_dict_list_of_dicts_deprecated(self):
        # GH#58862
        data = [
            {"key1": "value1", "key2": 42},
            {"key1": "value2", "key2": 123},
        ]
        msg = "Passing a list to DataFrame.from_dict is deprecated"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            result = pd.DataFrame.from_dict(data)
        expected = pd.DataFrame(data)
        tm.assert_frame_equal(result, expected)

    def test_from_dict_order_with_single_column(self):
        data = {
            "alpha": {
                "value2": 123,
                "value1": 532,
                "animal": 222,
                "plant": False,
                "name": "test",
            }
        }
        result = pd.DataFrame.from_dict(
            data,
            orient="columns",
        )
        expected = pd.DataFrame(
            [[123], [532], [222], [False], ["test"]],
            index=["value2", "value1", "animal", "plant", "name"],
            columns=["alpha"],
        )
        tm.assert_frame_equal(result, expected)
