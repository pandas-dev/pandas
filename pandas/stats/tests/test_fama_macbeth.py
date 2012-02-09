from pandas import DataFrame, Panel
from pandas.stats.api import fama_macbeth
from common import assert_almost_equal, BaseTest

import numpy as np

class TestFamaMacBeth(BaseTest):
    def testFamaMacBethRolling(self):
        # self.checkFamaMacBethExtended('rolling', self.panel_x, self.panel_y,
        #                               nw_lags_beta=2)

        # df = DataFrame(np.random.randn(50, 10))
        x = dict((k, DataFrame(np.random.randn(50, 10))) for k in 'abcdefg')
        x = Panel.from_dict(x)
        y = (DataFrame(np.random.randn(50, 10)) +
             DataFrame(0.01 * np.random.randn(50, 10)))
        self.checkFamaMacBethExtended('rolling', x, y, nw_lags_beta=2)
        self.checkFamaMacBethExtended('expanding', x, y, nw_lags_beta=2)

    def checkFamaMacBethExtended(self, window_type, x, y, **kwds):
        window = 25

        result = fama_macbeth(y=y, x=x, window_type=window_type, window=window,
                              **kwds)
        self._check_stuff_works(result)

        index =  result._index
        time = len(index)

        for i in xrange(time - window + 1):
            if window_type == 'rolling':
                start = index[i]
            else:
                start = index[0]

            end = index[i + window - 1]

            x2 = {}
            for k, v in x.iteritems():
                x2[k] = v.truncate(start, end)
            y2 = y.truncate(start, end)

            reference = fama_macbeth(y=y2, x=x2, **kwds)
            assert_almost_equal(reference._stats, result._stats[:, i])

        static = fama_macbeth(y=y2, x=x2, **kwds)
        self._check_stuff_works(static)

    def _check_stuff_works(self, result):
        # does it work?
        attrs = ['mean_beta', 'std_beta', 't_stat']
        for attr in attrs:
            getattr(result, attr)

        # does it work?
        result.summary

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)
