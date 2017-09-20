import numpy as np
import pandas as pd

from pandas.core.dtypes import generic

abcs = [getattr(generic, name) for name in dir(generic)
        if name.startswith('ABC')]


class GenericDtypes(object):
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

    def time_isinstance(self):
        for obj in self.objs:
            for cls in abcs:
                isinstance(obj, cls)
