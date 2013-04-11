# pylint: disable-msg=W0612,E1101

import unittest
import nose

from pandas.util.testing import assert_series_equal

class TestDgrep(unittest.TestCase):
    def test_dgrep(self):
        import pandas as pd
        from pandas import Series as Series
        from pandas.util.testing import makeCustomDataframe as mkdf

        import re
        pd.options.sandbox.dgrep=True # turn it on
        df=mkdf(30,4,r_idx_nlevels=3)
        df.index=range(30)
        df.iloc[5,0] = "supercool"
        df.iloc[6,0] = "supercool"
        df.iloc[29,0] = "supercool"
        df.iloc[15,1] = "supercool"
        df.iloc[17,2] = "supercool"
        # accepts colname and regex string
        rs = df.dgrep(".cool$","C_l0_g0")
        assert_series_equal(rs.C_l0_g0,Series(["supercool"]*3,index=[5,6,29]))
        # accepts lists of cols, can include a series such as df.series_name
        # (convenient for tab completion on columns)
        rs = df.dgrep(".cool$",['C_l0_g0','C_l0_g1'])
        xp = Series(["supercool","supercool","R15C0","supercool"],index=[5,6,15,29])
        assert_series_equal(rs.C_l0_g0,xp)
        self.assertEqual(rs.iloc[2,1],"supercool")

        # accepts a single named series
        rs = df.dgrep(".cool$",'C_l0_g1')
        xp = Series(["supercool"],index=[15])
        assert_series_equal(rs.C_l0_g1,xp)


        # specifying C=2 (or A/B=) does a grep context , providing
        # context lines around the hit
        # NB overlapping context lines do not cause line duplication (*)
        rs = df.dgrep(".cool$",["C_l0_g0"],C=2)
        xp = Series(['R4C0', 'supercool', 'supercool', 'R28C0', 'supercool'],index=[4,5,6,28,29])
        assert_series_equal(rs.C_l0_g0,xp)

        # also accepts lambda
        # NB, last match is at end, so only previous line of context displayed
        rs=df.dgrep(lambda x: bool(re.search(".cool$",x)),["C_l0_g0"],C=3)
        xp = Series(['R4C0', 'supercool', 'supercool', 'R7C0', 'R28C0', 'supercool'],index=[4,5,6,7,28,29])
        assert_series_equal(xp,rs.C_l0_g0)
        # split=True returns a series of (index_label_matched, dataframe)
        # pairs, similar to groupby
        # NB some lines appear in more then one group in this case (*)
        rs = df.dgrep(".cool$",["C_l0_g0"],split=True,C=3)
        self.assertEqual(len(rs),3)
        xp = Series(['R4C0', 'supercool', 'supercool'],index=[4,5,6])
        assert_series_equal(xp,rs[0][1].C_l0_g0)

        # works on series too
        s = df.C_l0_g0.dgrep(".cool$",C=3)
        xp = Series(['R4C0', 'supercool', 'supercool', 'R7C0', 'R28C0', 'supercool'],index=[4,5,6,7,28,29])
        assert_series_equal(xp,s)

        # can also get the values "applied" onto the function
        df.dgrep(lambda c1,c2: "cool" in c1 or "cool" in c2,df.columns[:2])

        # which also works with *args
        df.dgrep(lambda *args: "supercool" in args,df.columns[:3])
