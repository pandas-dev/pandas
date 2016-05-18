# -*- coding: utf-8 -*-
""""""

from __future__ import (absolute_import, division, print_function)

import pandas.util.testing as tm
import pandas as pd
import json


class TestJSON(tm.TestCase):
    testjson = u'''
    [{"Ünicøde":0,"sub":{"A":1, "B":2}},
     {"Ünicøde":1,"sub":{"A":3, "B":4}}]
     '''.encode('utf8')

    testdata = {
        u'sub.A': [1, 3],
        u'sub.B': [2, 4],
        u'Ünicøde': [0, 1]
    }
    testdf = pd.DataFrame(testdata)

    def test_json_normalize(self):
        df = pd.io.json.json_normalize(json.loads(self.testjson))
        tm.assert_frame_equal(df, self.testdf)
