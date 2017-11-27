# -*- coding: utf-8 -*-


import pandas.util.testing as tm
from pandas import CategoricalIndex
from pandas.tests.categorical.common import TestCategoricalBlock


class TestCategoricalBlockGroubpy(TestCategoricalBlock):

    def test_groupby_sort(self):

        # http://stackoverflow.com/questions/23814368/sorting-pandas-categorical-labels-after-groupby
        # This should result in a properly sorted Series so that the plot
        # has a sorted x axis
        # self.cat.groupby(['value_group'])['value_group'].count().plot(kind='bar')

        res = self.cat.groupby(['value_group'])['value_group'].count()
        exp = res[sorted(res.index, key=lambda x: float(x.split()[0]))]
        exp.index = CategoricalIndex(exp.index, name=exp.index.name)
        tm.assert_series_equal(res, exp)
