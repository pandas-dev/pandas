# -*- coding: utf-8 -*-

import numpy as np

import pandas as pd
from pandas import Categorical, DataFrame


class TestCategorical(object):

    def setup_method(self, method):
        self.factor = Categorical(['a', 'b', 'b', 'a', 'a', 'c', 'c', 'c'],
                                  ordered=True)


class TestCategoricalBlock(object):

    def setup_method(self, method):
        self.factor = Categorical(['a', 'b', 'b', 'a', 'a', 'c', 'c', 'c'])

        df = DataFrame({'value': np.random.randint(0, 10000, 100)})
        labels = ["{0} - {1}".format(i, i + 499) for i in range(0, 10000, 500)]
        cat_labels = Categorical(labels, labels)

        df = df.sort_values(by=['value'], ascending=True)
        df['value_group'] = pd.cut(df.value, range(0, 10500, 500),
                                   right=False, labels=cat_labels)
        self.cat = df
