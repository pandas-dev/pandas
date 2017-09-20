""" define generic base classes for pandas objects """

from pandas._libs.dtypes.generic import (
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


# define abstract base classes to enable isinstance type checking on our
# objects
def create_pandas_abc_type(name, attr, comp):
    @classmethod
    def _check(cls, inst):
        return getattr(inst, attr, '_typ') in comp

    dct = dict(__instancecheck__=_check, __subclasscheck__=_check)
    meta = type("ABCBase", (type, ), dct)
    return meta(name, tuple(), dct)


class _ABCGeneric(type):

    def __instancecheck__(cls, inst):
        return hasattr(inst, "_data")


ABCGeneric = _ABCGeneric("ABCGeneric", tuple(), {})
