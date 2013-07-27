import nose

import pandas
from pandas.util.testing import network
from pandas.util.testing import assert_frame_equal
from numpy.testing.decorators import slow
from pandas.io.wb import search, download
import six


@slow
@network
def test_wdi_search():
    raise nose.SkipTest
    expected = {six.u('id'): {2634: six.u('GDPPCKD'),
                        4649: six.u('NY.GDP.PCAP.KD'),
                        4651: six.u('NY.GDP.PCAP.KN'),
                        4653: six.u('NY.GDP.PCAP.PP.KD')},
                six.u('name'): {2634: six.u('GDP per Capita, constant US$, '
                                             'millions'),
                          4649: six.u('GDP per capita (constant 2000 US$)'),
                          4651: six.u('GDP per capita (constant LCU)'),
                          4653: six.u('GDP per capita, PPP (constant 2005 '
                                      'international $)')}}
    result = search('gdp.*capita.*constant').ix[:, :2]
    expected = pandas.DataFrame(expected)
    expected.index = result.index
    assert_frame_equal(result, expected)


@slow
@network
def test_wdi_download():
    raise nose.SkipTest
    expected = {'GDPPCKN': {(six.u('United States'), six.u('2003')): six.u('40800.0735367688'), (six.u('Canada'), six.u('2004')): six.u('37857.1261134552'), (six.u('United States'), six.u('2005')): six.u('42714.8594790102'), (six.u('Canada'), six.u('2003')): six.u('37081.4575704003'), (six.u('United States'), six.u('2004')): six.u('41826.1728310667'), (six.u('Mexico'), six.u('2003')): six.u('72720.0691255285'), (six.u('Mexico'), six.u('2004')): six.u('74751.6003347038'), (six.u('Mexico'), six.u('2005')): six.u('76200.2154469437'), (six.u('Canada'), six.u('2005')): six.u('38617.4563629611')}, 'GDPPCKD': {(six.u('United States'), six.u('2003')): six.u('40800.0735367688'), (six.u('Canada'), six.u('2004')): six.u('34397.055116118'), (six.u('United States'), six.u('2005')): six.u('42714.8594790102'), (six.u('Canada'), six.u('2003')): six.u('33692.2812368928'), (six.u('United States'), six.u('2004')): six.u('41826.1728310667'), (six.u('Mexico'), six.u('2003')): six.u('7608.43848670658'), (six.u('Mexico'), six.u('2004')): six.u('7820.99026814334'), (six.u('Mexico'), six.u('2005')): six.u('7972.55364129367'), (six.u('Canada'), six.u('2005')): six.u('35087.8925933298')}}
    expected = pandas.DataFrame(expected)
    result = download(country=['CA', 'MX', 'US', 'junk'], indicator=['GDPPCKD',
                                                                     'GDPPCKN', 'junk'], start=2003, end=2005)
    expected.index = result.index
    assert_frame_equal(result, pandas.DataFrame(expected))


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
