from __future__ import annotations

import pytest

import pandas as pd
from pandas import api
import pandas._testing as tm


class Base:
    def check(self, namespace, expected, ignored=None):
        # see which names are in the namespace, minus optional
        # ignored ones
        # compare vs the expected

        result = sorted(
            f for f in dir(namespace) if not f.startswith("__") and f != "annotations"
        )
        if ignored is not None:
            result = sorted(set(result) - set(ignored))

        expected = sorted(expected)
        tm.assert_almost_equal(result, expected)


class TestPDApi(Base):
    # these are optionally imported based on testing
    # & need to be ignored
    ignored = ["tests", "locale", "conftest"]

    # top-level sub-packages
    public_lib = [
        "api",
        "arrays",
        "options",
        "test",
        "testing",
        "errors",
        "plotting",
        "io",
        "tseries",
    ]
    private_lib = ["compat", "core", "pandas", "util"]

    # misc
    misc = ["IndexSlice", "NaT", "NA"]

    # top-level classes
    classes = [
        "ArrowDtype",
        "Categorical",
        "CategoricalIndex",
        "DataFrame",
        "DateOffset",
        "DatetimeIndex",
        "ExcelFile",
        "ExcelWriter",
        "Flags",
        "Grouper",
        "HDFStore",
        "Index",
        "MultiIndex",
        "Period",
        "PeriodIndex",
        "RangeIndex",
        "Series",
        "SparseDtype",
        "StringDtype",
        "Timedelta",
        "TimedeltaIndex",
        "Timestamp",
        "Interval",
        "IntervalIndex",
        "CategoricalDtype",
        "PeriodDtype",
        "IntervalDtype",
        "DatetimeTZDtype",
        "BooleanDtype",
        "Int8Dtype",
        "Int16Dtype",
        "Int32Dtype",
        "Int64Dtype",
        "UInt8Dtype",
        "UInt16Dtype",
        "UInt32Dtype",
        "UInt64Dtype",
        "Float32Dtype",
        "Float64Dtype",
        "NamedAgg",
    ]

    # these are already deprecated; awaiting removal
    deprecated_classes: list[str] = []

    # external modules exposed in pandas namespace
    modules: list[str] = []

    # top-level functions
    funcs = [
        "array",
        "bdate_range",
        "concat",
        "crosstab",
        "cut",
        "date_range",
        "interval_range",
        "eval",
        "factorize",
        "get_dummies",
        "from_dummies",
        "infer_freq",
        "isna",
        "isnull",
        "lreshape",
        "melt",
        "notna",
        "notnull",
        "offsets",
        "merge",
        "merge_ordered",
        "merge_asof",
        "period_range",
        "pivot",
        "pivot_table",
        "qcut",
        "show_versions",
        "timedelta_range",
        "unique",
        "value_counts",
        "wide_to_long",
    ]

    # top-level option funcs
    funcs_option = [
        "reset_option",
        "describe_option",
        "get_option",
        "option_context",
        "set_option",
        "set_eng_float_format",
    ]

    # top-level read_* funcs
    funcs_read = [
        "read_clipboard",
        "read_csv",
        "read_excel",
        "read_fwf",
        "read_gbq",
        "read_hdf",
        "read_html",
        "read_xml",
        "read_json",
        "read_pickle",
        "read_sas",
        "read_sql",
        "read_sql_query",
        "read_sql_table",
        "read_stata",
        "read_table",
        "read_feather",
        "read_parquet",
        "read_orc",
        "read_spss",
    ]

    # top-level json funcs
    funcs_json = ["json_normalize"]

    # top-level to_* funcs
    funcs_to = ["to_datetime", "to_numeric", "to_pickle", "to_timedelta"]

    # top-level to deprecate in the future
    deprecated_funcs_in_future: list[str] = []

    # these are already deprecated; awaiting removal
    deprecated_funcs: list[str] = []

    # private modules in pandas namespace
    private_modules = [
        "_config",
        "_libs",
        "_is_numpy_dev",
        "_testing",
        "_typing",
        "_version",
    ]

    def test_api(self):
        checkthese = (
            self.public_lib
            + self.private_lib
            + self.misc
            + self.modules
            + self.classes
            + self.funcs
            + self.funcs_option
            + self.funcs_read
            + self.funcs_json
            + self.funcs_to
            + self.private_modules
        )
        self.check(namespace=pd, expected=checkthese, ignored=self.ignored)

    def test_api_all(self):
        expected = set(
            self.public_lib
            + self.misc
            + self.modules
            + self.classes
            + self.funcs
            + self.funcs_option
            + self.funcs_read
            + self.funcs_json
            + self.funcs_to
        ) - set(self.deprecated_classes)
        actual = set(pd.__all__)

        extraneous = actual - expected
        assert not extraneous

        missing = expected - actual
        assert not missing

    def test_depr(self):
        deprecated_list = (
            self.deprecated_classes
            + self.deprecated_funcs
            + self.deprecated_funcs_in_future
        )
        for depr in deprecated_list:
            with tm.assert_produces_warning(FutureWarning):
                _ = getattr(pd, depr)


class TestApi(Base):
    allowed = ["types", "extensions", "indexers", "interchange"]

    def test_api(self):
        self.check(api, self.allowed)


class TestTesting(Base):
    funcs = [
        "assert_frame_equal",
        "assert_series_equal",
        "assert_index_equal",
        "assert_extension_array_equal",
    ]

    def test_testing(self):
        from pandas import testing

        self.check(testing, self.funcs)

    def test_util_in_top_level(self):
        with pytest.raises(AttributeError, match="foo"):
            pd.util.foo
