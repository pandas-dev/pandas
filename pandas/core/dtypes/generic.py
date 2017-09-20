""" define generic base classes for pandas objects """

from pandas._libs.dtypes.generic import (  # noqa
    create_pandas_abc_type,
    ABCIndex,
    ABCInt64Index, ABCUInt64Index,
    ABCRangeIndex,
    ABCFloat64Index,
    ABCMultiIndex,
    ABCDatetimeIndex, ABCTimedeltaIndex, ABCPeriodIndex,
    ABCCategoricalIndex,
    ABCIntervalIndex,
    ABCIndexClass,
    ABCSeries, ABCDataFrame, ABCPanel,
    ABCSparseSeries, ABCSparseArray,
    ABCCategorical,
    ABCPeriod,
    ABCDateOffset)


class _ABCGeneric(type):

    def __instancecheck__(cls, inst):
        return hasattr(inst, "_data")


ABCGeneric = _ABCGeneric("ABCGeneric", tuple(), {})
