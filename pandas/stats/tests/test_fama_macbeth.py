from pandas.stats.api import fama_macbeth
from common import assert_almost_equal, BaseTest

class TestFamaMacBeth(BaseTest):
    def testFamaMacBethRolling(self):
        self.checkFamaMacBethExtended('rolling', self.panel_x, self.panel_y)

    def checkFamaMacBethExtended(self, window_type, x, y, **kwds):
        window = 25

        result = fama_macbeth(y=y, x=x, window_type=window_type, window=window,
                              **kwds)

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
