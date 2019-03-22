import numpy as np

from pandas._libs.tslibs.parsing import _does_string_look_like_datetime

from pandas.io.parsers import _concat_date_cols


class DoesStringLookLikeDatetime(object):

    params = (['2Q2005', '0.0', '10000'],)
    param_names = ['value']

    def setup(self, value):
        self.objects = [value] * 1000000

    def time_check_datetimes(self, value):
        for obj in self.objects:
            try:
                _does_string_look_like_datetime(obj)
            except ValueError:
                pass


class ConcatDateCols(object):

    params = ([1234567890, 'AAAA'], [1, 2], [np.array, list])
    param_names = ['value', 'dim', 'container']

    def setup(self, value, dim, container):
        count_elem = 10000
        if dim == 1:
            self.object = (container([value] * count_elem),)
        if dim == 2:
            self.object = (container([value] * count_elem),
                           container([value] * count_elem))

    def time_check_concat(self, value, dim, container):
        _concat_date_cols(self.object)
