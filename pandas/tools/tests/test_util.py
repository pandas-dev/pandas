import os
import locale
import codecs
import unittest

import nose

import numpy as np
from numpy.testing import assert_equal

from pandas.tools.util import cartesian_product
from pandas.util.testing import get_locales, set_locale


CURRENT_LOCALE = locale.getlocale()
LOCALE_OVERRIDE = os.environ.get('LOCALE_OVERRIDE', None)


class TestCartesianProduct(unittest.TestCase):

    def test_simple(self):
        x, y = list('ABC'), [1, 22]
        result = cartesian_product([x, y])
        expected = [np.array(['A', 'A', 'B', 'B', 'C', 'C']),
                    np.array([ 1, 22,  1, 22,  1, 22])]
        assert_equal(result, expected)


class TestLocaleUtils(unittest.TestCase):
    def test_get_locales(self):
        assert len(get_locales(prefix='en')) > 0

    def test_get_locales_prefix(self):
        with set_locale('en_US.UTF-8'):
            assert len(get_locales(prefix='en')) > 0

    def test_set_locale(self):
        old_locale = CURRENT_LOCALE

        if LOCALE_OVERRIDE is not None:
            lang, enc = LOCALE_OVERRIDE.split('.')
        else:
            lang, enc = 'it_CH', 'UTF-8'

        enc = codecs.lookup(enc).name
        new_locale = lang, enc

        with set_locale(new_locale) as normalized_locale:
            new_lang, new_enc = normalized_locale.split('.')
            new_enc = codecs.lookup(enc).name
            normalized_locale = new_lang, new_enc
            self.assertEqual(normalized_locale, new_locale)

        current_locale = locale.getlocale()
        self.assertEqual(current_locale, old_locale)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
