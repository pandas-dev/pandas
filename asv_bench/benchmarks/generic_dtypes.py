import numpy as np
import pandas as pd

from pandas.core.dtypes import generic

abcs = [getattr(generic, name) for name in dir(generic)
        if name.startswith('ABC')]


class GenericDtypes(object):
    goal_time = .0821

    def setup(self):
        i64index = pd.Index(list(range(3)))
        idx = pd.Index(['A', 'B'])
        mi = pd.MultiIndex.from_product([i64index, idx])

        self.objs = [
            4,
            3.14,
            np.datetime64('NaT'),
            np.timedelta64(747),
            False,
            None,
            {},
            (x for x in range(5, 6)),
            list(range(19, 25)),
            "Gentlemen.  You can't fight in here.  This is the War Room!",
            np.random.random(6),
            pd.Period('2012-06-01', freq='M'),
            pd.Series(np.random.random(6)),
            pd.DataFrame(np.random.randn(5, 2), columns=['A', 'B']),
            pd.DataFrame(np.random.randn(6, 2),
                         index=mi,
                         columns=[True, False]).to_panel(),
            i64index,
            idx,
            pd.Index(pd.compat.range(144, 169)),
            mi,
            pd.Categorical(['Do', 'Re', 'Mi', 'Fa']),
            pd.CategoricalIndex(['Do'] * 5,
                                categories=['Do', 'Re', 'Mi', 'Fa']),
        ]

    def time_isinstance_ABCCategorical(self):
        for obj in self.objs:
            isinstance(obj, generic.ABCCategorical)

    def time_isinstance_ABCCategoricalIndex(self):
        for obj in self.objs:
            isinstance(obj, generic.ABCCategoricalIndex)

    def time_isinstance_ABCDataFrame(self):
        for obj in self.objs:
            isinstance(obj, generic.ABCDataFrame)

    def time_isinstance_ABCDateOffset(self):
        for obj in self.objs:
            isinstance(obj, generic.ABCDateOffset)

    def time_isinstance_ABCDatetimeIndex(self):
        for obj in self.objs:
            isinstance(obj, generic.ABCDatetimeIndex)

    def time_isinstance_ABCFloat64Index(self):
        for obj in self.objs:
            isinstance(obj, generic.ABCFloat64Index)

    def time_isinstance_ABCGeneric(self):
        for obj in self.objs:
            isinstance(obj, generic.ABCGeneric)

    def time_isinstance_ABCIndex(self):
        for obj in self.objs:
            isinstance(obj, generic.ABCIndex)

    def time_isinstance_ABCIndexClass(self):
        for obj in self.objs:
            isinstance(obj, generic.ABCIndexClass)

    def time_isinstance_ABCInt64Index(self):
        for obj in self.objs:
            isinstance(obj, generic.ABCInt64Index)

    def time_isinstance_ABCIntervalIndex(self):
        for obj in self.objs:
            isinstance(obj, generic.ABCIntervalIndex)

    def time_isinstance_ABCMultiIndex(self):
        for obj in self.objs:
            isinstance(obj, generic.ABCMultiIndex)

    def time_isinstance_ABCPanel(self):
        for obj in self.objs:
            isinstance(obj, generic.ABCPanel)

    def time_isinstance_ABCPeriod(self):
        for obj in self.objs:
            isinstance(obj, generic.ABCPeriod)

    def time_isinstance_ABCPeriodIndex(self):
        for obj in self.objs:
            isinstance(obj, generic.ABCPeriodIndex)

    def time_isinstance_ABCRangeIndex(self):
        for obj in self.objs:
            isinstance(obj, generic.ABCRangeIndex)

    def time_isinstance_ABCSeries(self):
        for obj in self.objs:
            isinstance(obj, generic.ABCSeries)

    def time_isinstance_ABCSparseArray(self):
        for obj in self.objs:
            isinstance(obj, generic.ABCSparseArray)

    def time_isinstance_ABCSparseSeries(self):
        for obj in self.objs:
            isinstance(obj, generic.ABCSparseSeries)

    def time_isinstance_ABCTimedeltaIndex(self):
        for obj in self.objs:
            isinstance(obj, generic.ABCTimedeltaIndex)

    def time_isinstance_ABCUInt64Index(self):
        for obj in self.objs:
            isinstance(obj, generic.ABCUInt64Index)
