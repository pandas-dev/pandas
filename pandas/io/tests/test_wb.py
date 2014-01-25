import nose

import pandas
from pandas.compat import u
from pandas.util.testing import network
from pandas.util.testing import assert_frame_equal
from numpy.testing.decorators import slow
from pandas.io.wb import search, download, get_countries
import pandas.util.testing as tm


class TestWB(tm.TestCase):

    @slow
    @network
    def test_wdi_search(self):
        raise nose.SkipTest

        expected = {u('id'): {2634: u('GDPPCKD'),
                              4649: u('NY.GDP.PCAP.KD'),
                              4651: u('NY.GDP.PCAP.KN'),
                              4653: u('NY.GDP.PCAP.PP.KD')},
                    u('name'): {2634: u('GDP per Capita, constant US$, '
                                        'millions'),
                                4649: u('GDP per capita (constant 2000 US$)'),
                                4651: u('GDP per capita (constant LCU)'),
                                4653: u('GDP per capita, PPP (constant 2005 '
                                        'international $)')}}
        result = search('gdp.*capita.*constant').ix[:, :2]
        expected = pandas.DataFrame(expected)
        expected.index = result.index
        assert_frame_equal(result, expected)

    @slow
    @network
    def test_wdi_download(self):
        raise nose.SkipTest

        expected = {'GDPPCKN': {(u('United States'), u('2003')): u('40800.0735367688'), (u('Canada'), u('2004')): u('37857.1261134552'), (u('United States'), u('2005')): u('42714.8594790102'), (u('Canada'), u('2003')): u('37081.4575704003'), (u('United States'), u('2004')): u('41826.1728310667'), (u('Mexico'), u('2003')): u('72720.0691255285'), (u('Mexico'), u('2004')): u('74751.6003347038'), (u('Mexico'), u('2005')): u('76200.2154469437'), (u('Canada'), u('2005')): u('38617.4563629611')}, 'GDPPCKD': {(u('United States'), u('2003')): u('40800.0735367688'), (u('Canada'), u('2004')): u('34397.055116118'), (u('United States'), u('2005')): u('42714.8594790102'), (u('Canada'), u('2003')): u('33692.2812368928'), (u('United States'), u('2004')): u('41826.1728310667'), (u('Mexico'), u('2003')): u('7608.43848670658'), (u('Mexico'), u('2004')): u('7820.99026814334'), (u('Mexico'), u('2005')): u('7972.55364129367'), (u('Canada'), u('2005')): u('35087.8925933298')}}
        expected = pandas.DataFrame(expected)
        result = download(country=['CA', 'MX', 'US', 'junk'], indicator=['GDPPCKD',
                                                                         'GDPPCKN', 'junk'], start=2003, end=2005)
        expected.index = result.index
        assert_frame_equal(result, pandas.DataFrame(expected))

    @slow
    @network
    def test_wdi_get_countries(self):
        result = get_countries()
        self.assertTrue('Zimbabwe' in list(result['name']))


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
