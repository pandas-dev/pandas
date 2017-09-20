#!/usr/bin/env python
# -*- coding: utf-8 -*-
# cython: profile=False
# cython: linetrace=False
# distutils: define_macros=CYTHON_TRACE=0
# distutils: define_macros=CYTHON_TRACE_NOGIL=0
"""define generic base classes for pandas objects"""


# define abstract base classes to enable isinstance type checking on our
# objects
def create_pandas_abc_type(name, attr, comp):
    comp = set(comp)

    if len(comp) == 1:
        comp = list(comp)[0]
        @classmethod
        def _check(cls, object inst):
            return getattr(inst, attr, "_typ") == comp

    else:
        @classmethod
        def _check(cls, object inst):
            return getattr(inst, attr, "_typ") in comp

    dct = dict(__instancecheck__=_check, __subclasscheck__=_check)
    meta = type("ABCBase", (type, ), dct)
    return meta(name, tuple(), dct)


ABCIndex = create_pandas_abc_type("ABCIndex", "_typ", ("index",))
ABCInt64Index = create_pandas_abc_type("ABCInt64Index", "_typ",
                                       ("int64index",))
ABCUInt64Index = create_pandas_abc_type("ABCUInt64Index", "_typ",
                                        ("uint64index",))
ABCRangeIndex = create_pandas_abc_type("ABCRangeIndex", "_typ",
                                       ("rangeindex",))
ABCFloat64Index = create_pandas_abc_type("ABCFloat64Index", "_typ",
                                         ("float64index",))
ABCMultiIndex = create_pandas_abc_type("ABCMultiIndex", "_typ",
                                       ("multiindex",))
ABCDatetimeIndex = create_pandas_abc_type("ABCDatetimeIndex", "_typ",
                                          ("datetimeindex",))
ABCTimedeltaIndex = create_pandas_abc_type("ABCTimedeltaIndex", "_typ",
                                           ("timedeltaindex",))
ABCPeriodIndex = create_pandas_abc_type("ABCPeriodIndex", "_typ",
                                        ("periodindex",))
ABCCategoricalIndex = create_pandas_abc_type("ABCCategoricalIndex", "_typ",
                                             ("categoricalindex",))
ABCIntervalIndex = create_pandas_abc_type("ABCIntervalIndex", "_typ",
                                          ("intervalindex",))
ABCIndexClass = create_pandas_abc_type("ABCIndexClass", "_typ",
                                       ("index", "int64index", "rangeindex",
                                        "float64index", "uint64index",
                                        "multiindex", "datetimeindex",
                                        "timedeltaindex", "periodindex",
                                        "categoricalindex", "intervalindex"))

ABCSeries = create_pandas_abc_type("ABCSeries", "_typ", ("series",))
ABCDataFrame = create_pandas_abc_type("ABCDataFrame", "_typ", ("dataframe",))
ABCPanel = create_pandas_abc_type("ABCPanel", "_typ", ("panel", "panel4d"))

ABCSparseSeries = create_pandas_abc_type("ABCSparseSeries", "_subtyp",
                                         ("sparse_series",
                                          "sparse_time_series"))
ABCSparseArray = create_pandas_abc_type("ABCSparseArray", "_subtyp",
                                        ("sparse_array", "sparse_series"))

ABCCategorical = create_pandas_abc_type("ABCCategorical", "_typ",
                                        ("categorical",))
ABCPeriod = create_pandas_abc_type("ABCPeriod", "_typ", ("period",))
ABCDateOffset = create_pandas_abc_type("ABCDateOffset", "_typ",
                                       ("dateoffset",))
ABCTick = create_pandas_abc_type("ABCTick", "_typ", ("tick",))

ABCBlockManager = create_pandas_abc_type("ABCBlockManager", "_typ",
                                         ("blockmanager",))
ABCNDFrame = create_pandas_abc_type("ABCNDFrame", "_typ", ("ndframe",))
ABCGroupBy = create_pandas_abc_type("ABCGroupy", "_typ", ("groupby",))
ABCGrouper = create_pandas_abc_type("ABCGrouper", "_typ", ("grouper",))
