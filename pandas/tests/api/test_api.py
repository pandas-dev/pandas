import pandas as pd
from pandas import api, compat
from pandas.util import testing as tm


class Base:
    def check(self, namespace, expected, ignored=None):
        # see which names are in the namespace, minus optional
        # ignored ones
        # compare vs the expected

        result = sorted(f for f in dir(namespace) if not f.startswith("__"))
        if ignored is not None:
            result = sorted(list(set(result) - set(ignored)))

        expected = sorted(expected)
        tm.assert_almost_equal(result, expected)


class TestPDApi(Base):

    # these are optionally imported based on testing
    # & need to be ignored
    ignored = ["tests", "locale", "conftest"]

    # top-level sub-packages
    lib = [
        "api",
        "arrays",
        "compat",
        "core",
        "errors",
        "pandas",
        "plotting",
        "test",
        "testing",
        "tseries",
        "util",
        "options",
        "io",
    ]

    # these are already deprecated; awaiting removal
    deprecated_modules = []

    # misc
    misc = ["IndexSlice", "NaT"]

    # top-level classes
    classes = [
        "Categorical",
        "CategoricalIndex",
        "DataFrame",
        "DateOffset",
        "DatetimeIndex",
        "ExcelFile",
        "ExcelWriter",
        "Float64Index",
        "Grouper",
        "HDFStore",
        "Index",
        "Int64Index",
        "MultiIndex",
        "Period",
        "PeriodIndex",
        "RangeIndex",
        "UInt64Index",
        "Series",
        "SparseArray",
        "SparseDataFrame",
        "SparseDtype",
        "SparseSeries",
        "Timedelta",
        "TimedeltaIndex",
        "Timestamp",
        "Interval",
        "IntervalIndex",
        "CategoricalDtype",
        "PeriodDtype",
        "IntervalDtype",
        "DatetimeTZDtype",
        "Int8Dtype",
        "Int16Dtype",
        "Int32Dtype",
        "Int64Dtype",
        "UInt8Dtype",
        "UInt16Dtype",
        "UInt32Dtype",
        "UInt64Dtype",
        "NamedAgg",
    ]
    if not compat.PY37:
        classes.append("Panel")

    # these are already deprecated; awaiting removal
    deprecated_classes = []

    # these should be deprecated in the future
    deprecated_classes_in_future = []

    # external modules exposed in pandas namespace
    modules = ["np", "datetime"]

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
        "read_json",
        "read_msgpack",
        "read_pickle",
        "read_sas",
        "read_sql",
        "read_sql_query",
        "read_sql_table",
        "read_stata",
        "read_table",
        "read_feather",
        "read_parquet",
        "read_spss",
    ]

    # top-level to_* funcs
    funcs_to = ["to_datetime", "to_msgpack", "to_numeric", "to_pickle", "to_timedelta"]

    # top-level to deprecate in the future
    deprecated_funcs_in_future = []

    # these are already deprecated; awaiting removal
    deprecated_funcs = []

    # private modules in pandas namespace
    private_modules = [
        "_config",
        "_hashtable",
        "_lib",
        "_libs",
        "_np_version_under1p14",
        "_np_version_under1p15",
        "_np_version_under1p16",
        "_np_version_under1p17",
        "_tslib",
        "_typing",
        "_version",
    ]

    def test_api(self):

        self.check(
            pd,
            self.lib
            + self.misc
            + self.modules
            + self.deprecated_modules
            + self.classes
            + self.deprecated_classes
            + self.deprecated_classes_in_future
            + self.funcs
            + self.funcs_option
            + self.funcs_read
            + self.funcs_to
            + self.deprecated_funcs_in_future
            + self.deprecated_funcs
            + self.private_modules,
            self.ignored,
        )


class TestApi(Base):

    allowed = ["types", "extensions"]

    def test_api(self):

        self.check(api, self.allowed)


class TestTesting(Base):

    funcs = ["assert_frame_equal", "assert_series_equal", "assert_index_equal"]

    def test_testing(self):

        from pandas import testing

        self.check(testing, self.funcs)
