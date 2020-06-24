""" define generic base classes for pandas objects """


# define abstract base classes to enable isinstance type checking on our
# objects
def create_pandas_abc_type(name, attr, comp):

    # https://github.com/python/mypy/issues/1006
    # error: 'classmethod' used with a non-method
    @classmethod  # type: ignore
    def _check(cls, inst) -> bool:
        return getattr(inst, attr, "_typ") in comp

    dct = dict(__instancecheck__=_check, __subclasscheck__=_check)
    meta = type("ABCBase", (type,), dct)
    return meta(name, tuple(), dct)


ABCIndex = create_pandas_abc_type("ABCIndex", "_typ", ("index",))
ABCInt64Index = create_pandas_abc_type("ABCInt64Index", "_typ", ("int64index",))
ABCUInt64Index = create_pandas_abc_type("ABCUInt64Index", "_typ", ("uint64index",))
ABCRangeIndex = create_pandas_abc_type("ABCRangeIndex", "_typ", ("rangeindex",))
ABCFloat64Index = create_pandas_abc_type("ABCFloat64Index", "_typ", ("float64index",))
ABCMultiIndex = create_pandas_abc_type("ABCMultiIndex", "_typ", ("multiindex",))
ABCDatetimeIndex = create_pandas_abc_type(
    "ABCDatetimeIndex", "_typ", ("datetimeindex",)
)
ABCTimedeltaIndex = create_pandas_abc_type(
    "ABCTimedeltaIndex", "_typ", ("timedeltaindex",)
)
ABCPeriodIndex = create_pandas_abc_type("ABCPeriodIndex", "_typ", ("periodindex",))
ABCCategoricalIndex = create_pandas_abc_type(
    "ABCCategoricalIndex", "_typ", ("categoricalindex",)
)
ABCIntervalIndex = create_pandas_abc_type(
    "ABCIntervalIndex", "_typ", ("intervalindex",)
)
ABCIndexClass = create_pandas_abc_type(
    "ABCIndexClass",
    "_typ",
    {
        "index",
        "int64index",
        "rangeindex",
        "float64index",
        "uint64index",
        "multiindex",
        "datetimeindex",
        "timedeltaindex",
        "periodindex",
        "categoricalindex",
        "intervalindex",
    },
)

ABCSeries = create_pandas_abc_type("ABCSeries", "_typ", ("series",))
ABCDataFrame = create_pandas_abc_type("ABCDataFrame", "_typ", ("dataframe",))

ABCCategorical = create_pandas_abc_type("ABCCategorical", "_typ", ("categorical"))
ABCDatetimeArray = create_pandas_abc_type("ABCDatetimeArray", "_typ", ("datetimearray"))
ABCTimedeltaArray = create_pandas_abc_type(
    "ABCTimedeltaArray", "_typ", ("timedeltaarray")
)
ABCPeriodArray = create_pandas_abc_type("ABCPeriodArray", "_typ", ("periodarray",))
ABCExtensionArray = create_pandas_abc_type(
    "ABCExtensionArray",
    "_typ",
    # Note: IntervalArray and SparseArray are included bc they have _typ="extension"
    {"extension", "categorical", "periodarray", "datetimearray", "timedeltaarray"},
)
ABCPandasArray = create_pandas_abc_type("ABCPandasArray", "_typ", ("npy_extension",))
