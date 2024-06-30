"""
Test datetime formatting low-level routines
"""

from contextlib import nullcontext
from datetime import time
import locale

import pytest

from pandas._libs.tslibs.strftime import (
    UnsupportedStrFmtDirective,
    _create_escape_sequence,
    get_current_locale_specific_string,
)

import pandas._testing as tm

from pandas.tseries.api import convert_strftime_format


def get_local_am_pm():
    """Return the AM and PM strings returned by strftime in current locale."""
    am_local = time(1).strftime("%p")
    pm_local = time(13).strftime("%p")
    return am_local, pm_local


@pytest.mark.parametrize(
    "locale_str",
    [
        pytest.param(None, id=str(locale.getlocale())),
        "it_IT.utf8",
        "it_IT",  # Note: encoding will be 'ISO8859-1'
        "zh_CN.utf8",
        "zh_CN",  # Note: encoding will be 'gb2312'
    ],
)
def test_get_current_locale_specific_string(locale_str):
    """Test that `get_current_locale_specific_string` relies on runtime locale."""

    # Skip if locale cannot be set
    if locale_str is not None and not tm.can_set_locale(locale_str, locale.LC_ALL):
        pytest.skip(f"Skipping as locale '{locale_str}' cannot be set on host.")

    # Change locale temporarily for this test.
    with tm.set_locale(locale_str, locale.LC_ALL) if locale_str else nullcontext():
        # Get locale-specific reference
        am_local, pm_local = get_local_am_pm()

        # Test that the function returns the correct ones
        res = get_current_locale_specific_string()
        assert res.am == am_local
        assert res.pm == pm_local


class TestConvertStrftimeFormat:
    """Tests for `convert_strftime_format`."""

    @pytest.mark.parametrize(
        "strftime_fmt,res_fmt_old,res_fmt_new",
        (
            ("%p", "%(ampm)s", "{ampm:s}"),
            (
                "%m-%d-%Y",
                "%(month)02d-%(day)02d-%(year)04d",
                "{month:02d}-{day:02d}-{year:04d}",
            ),
            (
                "20%y-%m-%d__foo__%I:%M:%S%p",
                "20%(shortyear)02d-%(month)02d-%(day)02d__foo__"
                "%(hour12)02d:%(min)02d:%(sec)02d%(ampm)s",
                "20{shortyear:02d}-{month:02d}-{day:02d}__foo__"
                "{hour12:02d}:{min:02d}:{sec:02d}{ampm:s}",
            ),
        ),
    )
    def test_format_datetime(self, strftime_fmt, res_fmt_old, res_fmt_new):
        """Test that `convert_strftime_format` returns the correct template"""
        str_tmp, loc_s = convert_strftime_format(
            strftime_fmt, target="datetime", new_style_fmt=False
        )
        assert str_tmp == res_fmt_old

        str_tmp_new, loc_s2 = convert_strftime_format(
            strftime_fmt, target="datetime", new_style_fmt=True
        )
        assert loc_s2 == loc_s
        assert str_tmp_new == res_fmt_new

    @pytest.mark.parametrize(
        "strftime_fmt,res_fmt_old,res_fmt_new",
        (
            ("%p", "%(ampm)s", "{ampm:s}"),
            (
                "%m-%d-%Y",
                "%(month)02d-%(day)02d-%(year)04d",
                "{month:02d}-{day:02d}-{year:04d}",
            ),
            (
                "%y %I:%M:%S%p (ms=%l us=%u ns=%n)",
                "%(shortyear)02d %(hour12)02d:%(min)02d:%(sec)02d%(ampm)s "
                "(ms=%(ms)03d us=%(us)06d ns=%(ns)09d)",
                "{shortyear:02d} {hour12:02d}:{min:02d}:{sec:02d}{ampm:s} "
                "(ms={ms:03d} us={us:06d} ns={ns:09d})",
            ),
            (
                "20%y-%m-%d__f{o}o__%I:%M:%S%%%p",
                "20%(shortyear)02d-%(month)02d-%(day)02d__f{o}o__"
                "%(hour12)02d:%(min)02d:%(sec)02d%%%(ampm)s",
                "20{shortyear:02d}-{month:02d}-{day:02d}__f{{o}}o__"
                "{hour12:02d}:{min:02d}:{sec:02d}%%{ampm:s}",
            ),
        ),
    )
    def test_format_period(self, strftime_fmt, res_fmt_old, res_fmt_new):
        """Test that `convert_strftime_format` returns the correct template"""
        str_tmp, loc_s = convert_strftime_format(
            strftime_fmt, target="period", new_style_fmt=False
        )
        assert str_tmp == res_fmt_old

        str_tmp_new, loc_s2 = convert_strftime_format(
            strftime_fmt, target="period", new_style_fmt=True
        )
        assert loc_s2 == loc_s
        assert str_tmp_new == res_fmt_new

    @pytest.mark.parametrize(
        "locale_str",
        [
            pytest.param(None, id=str(locale.getlocale())),
            "it_IT.utf8",
            "it_IT",  # Note: encoding will be 'ISO8859-1'
            "zh_CN.utf8",
            "zh_CN",  # Note: encoding will be 'gb2312'
        ],
    )
    @pytest.mark.parametrize("target", ("datetime", "date", "time", "period"))
    def test_format_non_ascii(self, locale_str, target):
        """Test that `convert_strftime_format` is robust to locale and fmt encoding"""

        # Skip if locale cannot be set
        if locale_str is not None and not tm.can_set_locale(locale_str, locale.LC_ALL):
            pytest.skip(f"Skipping as locale '{locale_str}' cannot be set on host.")

        # Change locale temporarily for this test.
        with tm.set_locale(locale_str, locale.LC_ALL) if locale_str else nullcontext():
            strftime_fmt = "%y é"

            str_tmp, _ = convert_strftime_format(
                strftime_fmt, target="datetime", new_style_fmt=False
            )
            assert str_tmp == "%(shortyear)02d é"

            str_tmp_new, _ = convert_strftime_format(
                strftime_fmt, target="datetime", new_style_fmt=True
            )
            assert str_tmp_new == "{shortyear:02d} é"

    def test_invalid_datetime_directive(self):
        """Test that using invalid strftime directives for datetime raises an error"""
        with pytest.raises(UnsupportedStrFmtDirective, match="Unsupported directive"):
            convert_strftime_format("%F", target="datetime")

        # Make sure that the same directive is valid for periods
        assert convert_strftime_format("%F", target="period")[0] == "%(Fyear)d"

    def test_invalid_period_directive(self):
        """Test that using invalid strftime directives for period raises an error"""
        with pytest.raises(UnsupportedStrFmtDirective, match="Unsupported directive"):
            convert_strftime_format("%j", target="period")

    def test_unknown_directive(self):
        """Test that unknown/not available strftime directives lead to an error."""
        with pytest.raises(ValueError, match="Unsupported directive"):
            convert_strftime_format("%O", target="datetime")

        with pytest.raises(ValueError, match="Unsupported directive"):
            convert_strftime_format("%O", target="datetime", new_style_fmt=True)


def test_create_escape_sequence():
    txt = "-*"
    esc = _create_escape_sequence(txt, init_esc="*", prefix="-")
    assert esc not in txt
    assert esc == "--*"
