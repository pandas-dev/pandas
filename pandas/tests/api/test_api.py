from __future__ import annotations

import importlib
import inspect
import pathlib
import pkgutil

import pytest

import pandas as pd
from pandas import api
import pandas._testing as tm
from pandas.api import (
    executors as api_executors,
    extensions as api_extensions,
    indexers as api_indexers,
    interchange as api_interchange,
    types as api_types,
    typing as api_typing,
)
from pandas.api.typing import aliases as api_aliases


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
    ignored = ["tests", "locale", "conftest", "_version_meson"]

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
    private_lib = ["compat", "core", "pandas", "util", "_built_with_meson"]

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
        "col",
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
        "read_iceberg",
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
        "_pandas_datetime_CAPI",
        "_pandas_parser_CAPI",
        "_testing",
        "_typing",
    ]
    if not pd._built_with_meson:
        private_modules.append("_version")

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
    allowed_api_dirs = [
        "executors",
        "types",
        "extensions",
        "indexers",
        "interchange",
        "typing",
        "internals",
    ]
    allowed_typing = [
        "DataFrameGroupBy",
        "DatetimeIndexResamplerGroupby",
        "Expanding",
        "ExpandingGroupby",
        "ExponentialMovingWindow",
        "ExponentialMovingWindowGroupby",
        "Expression",
        "FrozenList",
        "JsonReader",
        "NaTType",
        "NAType",
        "NoDefault",
        "PeriodIndexResamplerGroupby",
        "Resampler",
        "Rolling",
        "RollingGroupby",
        "SeriesGroupBy",
        "StataReader",
        "SASReader",
        "TimedeltaIndexResamplerGroupby",
        "TimeGrouper",
        "Window",
        "aliases",
    ]
    allowed_api_types = [
        "is_any_real_numeric_dtype",
        "is_array_like",
        "is_bool",
        "is_bool_dtype",
        "is_categorical_dtype",
        "is_complex",
        "is_complex_dtype",
        "is_datetime64_any_dtype",
        "is_datetime64_dtype",
        "is_datetime64_ns_dtype",
        "is_datetime64tz_dtype",
        "is_dict_like",
        "is_dtype_equal",
        "is_extension_array_dtype",
        "is_file_like",
        "is_float",
        "is_float_dtype",
        "is_hashable",
        "is_int64_dtype",
        "is_integer",
        "is_integer_dtype",
        "is_interval_dtype",
        "is_iterator",
        "is_list_like",
        "is_named_tuple",
        "is_number",
        "is_numeric_dtype",
        "is_object_dtype",
        "is_period_dtype",
        "is_re",
        "is_re_compilable",
        "is_scalar",
        "is_signed_integer_dtype",
        "is_sparse",
        "is_string_dtype",
        "is_timedelta64_dtype",
        "is_timedelta64_ns_dtype",
        "is_unsigned_integer_dtype",
        "pandas_dtype",
        "infer_dtype",
        "union_categoricals",
        "CategoricalDtype",
        "DatetimeTZDtype",
        "IntervalDtype",
        "PeriodDtype",
    ]
    allowed_api_interchange = ["from_dataframe", "DataFrame"]
    allowed_api_indexers = [
        "check_array_indexer",
        "BaseIndexer",
        "FixedForwardWindowIndexer",
        "VariableOffsetWindowIndexer",
    ]
    allowed_api_extensions = [
        "no_default",
        "ExtensionDtype",
        "register_extension_dtype",
        "register_dataframe_accessor",
        "register_index_accessor",
        "register_series_accessor",
        "take",
        "ExtensionArray",
        "ExtensionScalarOpsMixin",
    ]
    allowed_api_executors = ["BaseExecutionEngine"]
    allowed_api_aliases = [
        "AggFuncType",
        "AlignJoin",
        "AnyAll",
        "AnyArrayLike",
        "ArrayLike",
        "AstypeArg",
        "Axes",
        "Axis",
        "CSVEngine",
        "ColspaceArgType",
        "CompressionOptions",
        "CorrelationMethod",
        "DropKeep",
        "Dtype",
        "DtypeArg",
        "DtypeBackend",
        "DtypeObj",
        "ExcelWriterIfSheetExists",
        "ExcelWriterMergeCells",
        "FilePath",
        "FillnaOptions",
        "FloatFormatType",
        "FormattersType",
        "FromDictOrient",
        "HTMLFlavors",
        "IgnoreRaise",
        "IndexLabel",
        "InterpolateOptions",
        "JSONEngine",
        "JSONSerializable",
        "JoinHow",
        "JoinValidate",
        "MergeHow",
        "MergeValidate",
        "NaPosition",
        "NsmallestNlargestKeep",
        "OpenFileErrors",
        "Ordered",
        "ParquetCompressionOptions",
        "QuantileInterpolation",
        "ReadBuffer",
        "ReadCsvBuffer",
        "ReadPickleBuffer",
        "ReindexMethod",
        "Scalar",
        "SequenceNotStr",
        "SliceType",
        "SortKind",
        "StorageOptions",
        "Suffixes",
        "TakeIndexer",
        "TimeAmbiguous",
        "TimeGrouperOrigin",
        "TimeNonexistent",
        "TimeUnit",
        "TimedeltaConvertibleTypes",
        "TimestampConvertibleTypes",
        "ToStataByteorder",
        "ToTimestampHow",
        "UpdateJoin",
        "UsecolsArgType",
        "WindowingRankType",
        "WriteBuffer",
        "WriteExcelBuffer",
        "XMLParsers",
    ]

    def test_api(self):
        self.check(api, self.allowed_api_dirs)

    def test_api_typing(self):
        self.check(api_typing, self.allowed_typing)

    def test_api_types(self):
        self.check(api_types, self.allowed_api_types)

    def test_api_interchange(self):
        self.check(api_interchange, self.allowed_api_interchange)

    def test_api_indexers(self):
        self.check(api_indexers, self.allowed_api_indexers)

    def test_api_extensions(self):
        self.check(api_extensions, self.allowed_api_extensions)

    def test_api_executors(self):
        self.check(api_executors, self.allowed_api_executors)

    def test_api_typing_aliases(self):
        self.check(api_aliases, self.allowed_api_aliases)


class TestErrors(Base):
    def test_errors(self):
        ignored = ["_CurrentDeprecationWarning", "abc", "ctypes", "cow"]
        self.check(pd.errors, pd.errors.__all__, ignored=ignored)


class TestUtil(Base):
    def test_util(self):
        self.check(
            pd.util,
            ["hash_array", "hash_pandas_object"],
            ignored=[
                "_decorators",
                "_test_decorators",
                "_exceptions",
                "_validators",
                "capitalize_first_letter",
                "version",
                "_print_versions",
                "_tester",
            ],
        )


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


def get_pandas_objects(
    module_name: str, recurse: bool
) -> list[tuple[str, str, object]]:
    """
    Get all pandas objects within a module.

    An object is determined to be part of pandas if it has a string
    __module__ attribute that starts with ``"pandas"``.

    Parameters
    ----------
    module_name : str
        Name of the module to search.
    recurse : bool
        Whether to search submodules.

    Returns
    -------
        List of all objects that are determined to be a part of pandas.
    """
    module = importlib.import_module(module_name)
    objs = []

    for name, obj in inspect.getmembers(module):
        module_dunder = getattr(obj, "__module__", None)
        if isinstance(module_dunder, str) and module_dunder.startswith("pandas"):
            objs.append((module_name, name, obj))

    if not recurse:
        return objs

    # __file__ can, but shouldn't, be None
    assert isinstance(module.__file__, str)
    paths = [pathlib.Path(module.__file__).parent]
    for module_info in pkgutil.walk_packages(paths):
        name = module_info.name
        if name.startswith("_") or name == "internals":
            continue
        objs.extend(
            get_pandas_objects(f"{module.__name__}.{name}", recurse=module_info.ispkg)
        )
    return objs


@pytest.mark.slow
@pytest.mark.parametrize(
    "module_name",
    [
        "pandas",
        "pandas.api",
        "pandas.arrays",
        "pandas.errors",
        pytest.param("pandas.io", marks=pytest.mark.xfail(reason="Private imports")),
        "pandas.plotting",
        "pandas.testing",
    ],
)
def test_attributes_module(module_name):
    """
    Ensures that all public objects have their __module__ set to the public import path.
    """
    recurse = module_name not in ["pandas", "pandas.testing"]
    objs = get_pandas_objects(module_name, recurse=recurse)
    failures = [
        (module_name, name, type(obj), obj.__module__)
        for module_name, name, obj in objs
        if not (
            obj.__module__ == module_name
            # Explicit exceptions
            or ("Dtype" in name and obj.__module__ == "pandas")
            or (name == "Categorical" and obj.__module__ == "pandas")
        )
    ]
    assert len(failures) == 0, "\n".join(str(e) for e in failures)

    # Check that all objects can indeed be imported from their __module__
    failures = []
    for module_name, name, obj in objs:
        module = importlib.import_module(obj.__module__)
        try:
            getattr(module, name)
        except Exception:
            failures.append((module_name, name, type(obj), obj.__module__))
    assert len(failures) == 0, "\n".join(str(e) for e in failures)
