import nose
import unittest

import numpy as np

from pandas import DataFrame, Series
import pandas.util.testing as tm

from pandas.tools.tile import cut


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)


