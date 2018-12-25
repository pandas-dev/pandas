# -*- coding: utf-8 -*-
import codecs
import locale
import os

import pytest

from pandas.compat import is_platform_windows

import pandas.core.common as com
import pandas.util.testing as tm

_all_locales = tm.get_locales() or []
_current_locale = locale.getlocale()

# Don't run any of these tests if we are on Windows or have no locales.
pytestmark = pytest.mark.skipif(is_platform_windows() or not _all_locales,
                                reason="Need non-Windows and locales")

_skip_if_only_one_locale = pytest.mark.skipif(
    len(_all_locales) <= 1, reason="Need multiple locales for meaningful test")


def test_can_set_locale_valid_set():
    # Can set the default locale.
    assert tm.can_set_locale("")


def test_can_set_locale_invalid_set():
    # Cannot set an invalid locale.
    assert not tm.can_set_locale("non-existent_locale")


def test_can_set_locale_invalid_get(monkeypatch):
    # see gh-22129
    #
    # In some cases, an invalid locale can be set,
    # but a subsequent getlocale() raises a ValueError.

    def mock_get_locale():
        raise ValueError()

    with monkeypatch.context() as m:
        m.setattr(locale, "getlocale", mock_get_locale)
        assert not tm.can_set_locale("")


def test_get_locales_at_least_one():
    # see gh-9744
    assert len(_all_locales) > 0


@_skip_if_only_one_locale
def test_get_locales_prefix():
    first_locale = _all_locales[0]
    assert len(tm.get_locales(prefix=first_locale[:2])) > 0


@_skip_if_only_one_locale
def test_set_locale():
    if com._all_none(_current_locale):
        # Not sure why, but on some Travis runs with pytest,
        # getlocale() returned (None, None).
        pytest.skip("Current locale is not set.")

    locale_override = os.environ.get("LOCALE_OVERRIDE", None)

    if locale_override is None:
        lang, enc = "it_CH", "UTF-8"
    elif locale_override == "C":
        lang, enc = "en_US", "ascii"
    else:
        lang, enc = locale_override.split(".")

    enc = codecs.lookup(enc).name
    new_locale = lang, enc

    if not tm.can_set_locale(new_locale):
        msg = "unsupported locale setting"

        with pytest.raises(locale.Error, match=msg):
            with tm.set_locale(new_locale):
                pass
    else:
        with tm.set_locale(new_locale) as normalized_locale:
            new_lang, new_enc = normalized_locale.split(".")
            new_enc = codecs.lookup(enc).name

            normalized_locale = new_lang, new_enc
            assert normalized_locale == new_locale

    # Once we exit the "with" statement, locale should be back to what it was.
    current_locale = locale.getlocale()
    assert current_locale == _current_locale
