"""
Helpers for configuring locale settings.

Name `localization` is chosen to avoid overlap with builtin `locale` module.
"""
from contextlib import contextmanager
import locale


@contextmanager
def set_locale(new_locale, lc_var=locale.LC_ALL):
    """Context manager for temporarily setting a locale.

    Parameters
    ----------
    new_locale : str or tuple
        A string of the form <language_country>.<encoding>. For example to set
        the current locale to US English with a UTF8 encoding, you would pass
        "en_US.UTF-8".
    lc_var : int, default `locale.LC_ALL`
        The category of the locale being set.

    Notes
    -----
    This is useful when you want to run a particular block of code under a
    particular locale, without globally setting the locale. This probably isn't
    thread-safe.
    """
    current_locale = locale.getlocale()

    try:
        locale.setlocale(lc_var, new_locale)
        normalized_locale = locale.getlocale()
        if all(x is not None for x in normalized_locale):
            yield '.'.join(normalized_locale)
        else:
            yield new_locale
    finally:
        locale.setlocale(lc_var, current_locale)
