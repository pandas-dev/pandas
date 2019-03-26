import numpy as np

from pandas._libs.tslib import array_to_datetime

class ParseDateString(object):
    params = (['mY', 'mdY', 'mQY', 'hm'],)
    params_name = ['value']
    objects = {
        'mY':  ['01-2019', '1-2019'],
        'mdY': ['12/02/2010'],
        'mQY': ['1Q09', '1Q2000', '09Q1', '2000Q1'],
        'hm': ['21:34']
    }

    def setup(self, value):
        count_elem = 100000
        self.data = np.array(self.objects[value] * count_elem, dtype=np.object)

    def time_parse_datestring(self, value):
        array_to_datetime(self.data)
