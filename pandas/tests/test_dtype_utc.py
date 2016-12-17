#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 17 09:34:47 2016

@author: rodolfoxps
"""
import nose
import pandas as pd
import datetime, pytz
from pandas.util.testing import assert_frame_equal

def test_dtype_utc(self):
    
    data=pd.Series( [pd.NaT, pd.NaT, datetime.datetime(2016, 12, 12, 22, 24, 6, 100001, tzinfo=pytz.utc) ] )
  
    filled=data.fillna(method='bfill')
    
    expected=pd.Series([datetime.datetime(2016, 12, 12, 22, 24, 6, 100001, tzinfo=pytz.utc) ,
                      datetime.datetime(2016, 12, 12, 22, 24, 6, 100001, tzinfo=pytz.utc) ,
                      datetime.datetime(2016, 12, 12, 22, 24, 6, 100001, tzinfo=pytz.utc) ])

    assert_frame_equal(filled, expected)
    
if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)