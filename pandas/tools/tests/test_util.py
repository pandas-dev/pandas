import os
import locale
import codecs
import nose

import numpy as np
from numpy.testing import assert_equal

import pandas.util.testing as tm
from pandas.tools.util import cartesian_product


CURRENT_LOCALE = locale.getlocale()
LOCALE_OVERRIDE = os.environ.get('LOCALE_OVERRIDE', None)


class TestCartesianProduct(tm.TestCase):

    def test_simple(self):
        x, y = list('ABC'), [1, 22]
        result = cartesian_product([x, y])
        expected = [np.array(['A', 'A', 'B', 'B', 'C', 'C']),
                    np.array([ 1, 22,  1, 22,  1, 22])]
        assert_equal(result, expected)


class TestLocaleUtils(tm.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestLocaleUtils, cls).setUpClass()
        cls.locales = tm.get_locales()

        if not cls.locales:
            raise nose.SkipTest("No locales found")

        if os.name == 'nt':  # we're on windows
            raise nose.SkipTest("Running on Windows")

    @classmethod
    def tearDownClass(cls):
        super(TestLocaleUtils, cls).tearDownClass()
        del cls.locales

    def test_get_locales(self):
        # all systems should have at least a single locale
        assert len(tm.get_locales()) > 0

    def test_get_locales_prefix(self):
        if len(self.locales) == 1:
            raise nose.SkipTest("Only a single locale found, no point in "
                                "trying to test filtering locale prefixes")
        first_locale = self.locales[0]
        assert len(tm.get_locales(prefix=first_locale[:2])) > 0

    def test_set_locale(self):
        if len(self.locales) == 1:
            raise nose.SkipTest("Only a single locale found, no point in "
                                "trying to test setting another locale")

        if LOCALE_OVERRIDE is not None:
            lang, enc = LOCALE_OVERRIDE.split('.')
        else:
            lang, enc = 'it_CH', 'UTF-8'

        enc = codecs.lookup(enc).name
        new_locale = lang, enc

        if not tm._can_set_locale(new_locale):
            with tm.assertRaises(locale.Error):
                with tm.set_locale(new_locale):
                    pass
        else:
            with tm.set_locale(new_locale) as normalized_locale:
                new_lang, new_enc = normalized_locale.split('.')
                new_enc = codecs.lookup(enc).name
                normalized_locale = new_lang, new_enc
                self.assertEqual(normalized_locale, new_locale)

        current_locale = locale.getlocale()
        self.assertEqual(current_locale, CURRENT_LOCALE)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
