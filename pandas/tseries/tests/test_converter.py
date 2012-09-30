from datetime import datetime, time, timedelta
import sys
import os
import unittest

import nose

import numpy as np

from pandas.tseries import converter

def test_timtetonum_accepts_unicode():
    assert(converter.time2num("00:01")==converter.time2num(u"00:01"))

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)
