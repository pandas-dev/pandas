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

        # Test that a name column exists, and that some results were returned
        # ...without being too strict about what the actual contents of the
        # results actually are.  The fact that there are some, is good enough.

        result = search('gdp.*capita.*constant')
        self.assertTrue(result.name.str.contains('GDP').any())

    @slow
    @network
    def test_wdi_download(self):

        # Test a bad indicator with double (US), triple (USA),
        # standard (CA, MX), non standard (KSV),
        # duplicated (US, US, USA), and unknown (BLA) country codes

        # ...but NOT a crash inducing country code (World bank strips pandas
        #    users of the luxury of laziness, because they create their
        #    own exceptions, and don't clean up legacy country codes.
        # ...but NOT a retired indicator (User should want it to error.)

        cntry_codes = ['CA', 'MX', 'USA', 'US', 'US', 'KSV', 'BLA']
        inds = ['NY.GDP.PCAP.CD','BAD.INDICATOR']

        expected = {'NY.GDP.PCAP.CD': {('Canada', '2003'): 28026.006013044702, ('Mexico', '2003'): 6601.0420648056606, ('Canada', '2004'): 31829.522562759001, ('Kosovo', '2003'): 1969.56271307405, ('Mexico', '2004'): 7042.0247834044303, ('United States', '2004'): 41928.886136479705, ('United States', '2003'): 39682.472247320402, ('Kosovo', '2004'): 2135.3328465238301}}
        expected = pandas.DataFrame(expected)
        #Round, to ignore revisions to data.
        expected = pandas.np.round(expected,decimals=-3)
        expected.sort(inplace=True)
        result = download(country=cntry_codes, indicator=inds,
                          start=2003, end=2004, errors='ignore')
        result.sort(inplace=True)
        #Round, to ignore revisions to data.
        result = pandas.np.round(result,decimals=-3)
        expected.index = result.index
        assert_frame_equal(result, pandas.DataFrame(expected))

    @slow
    @network
    def test_wdi_download_w_retired_indicator(self):

        cntry_codes = ['CA', 'MX', 'US']
        # Despite showing up in the search feature, and being listed online,
        # the api calls to GDPPCKD don't work in their own query builder, nor
        # pandas module.  GDPPCKD used to be a common symbol.
        # This test is written to ensure that error messages to pandas users
        # continue to make sense, rather than a user getting some missing
        # key error, cause their JSON message format changed.  If
        # World bank ever finishes the deprecation of this symbol,
        # this nose test should still pass.

        inds = ['GDPPCKD']

        try:
            result = download(country=cntry_codes, indicator=inds,
                              start=2003, end=2004, errors='ignore')
        # If for some reason result actually ever has data, it's cause WB
        # fixed the issue with this ticker.  Find another bad one.
        except ValueError as e:
            raise nose.SkipTest("No indicators returned data: {0}".format(e))

        # if it ever gets here, it means WB unretired the indicator.
        # even if they dropped it completely, it would still get caught above
        # or the WB API changed somehow in a really unexpected way.
        if len(result) > 0:
            raise nose.SkipTest("Invalid results")

    @slow
    @network
    def test_wdi_download_w_crash_inducing_countrycode(self):

        cntry_codes = ['CA', 'MX', 'US', 'XXX']
        inds = ['NY.GDP.PCAP.CD']

        try:
            result = download(country=cntry_codes, indicator=inds,
                              start=2003, end=2004, errors='ignore')
        except ValueError as e:
            raise nose.SkipTest("No indicators returned data: {0}".format(e))

        # if it ever gets here, it means the country code XXX got used by WB
        # or the WB API changed somehow in a really unexpected way.
        if len(result) > 0:
            raise nose.SkipTest("Invalid results")

    @slow
    @network
    def test_wdi_get_countries(self):
        result = get_countries()
        self.assertTrue('Zimbabwe' in list(result['name']))
        self.assertTrue(len(result) > 100)

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
