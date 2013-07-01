import numpy as np

from itertools import izip
from pandas.core.common import isnull
from pandas.core.series import Series
import re
import pandas.lib as lib


def _get_array_list(arr, others):
    if isinstance(others[0], (list, np.ndarray)):
        arrays = [arr] + list(others)
    else:
        arrays = [arr, others]

    return [np.asarray(x, dtype=object) for x in arrays]


def str_cat(arr, others=None, sep=None, na_rep=None):
    """
    Concatenate arrays of strings with given separator

    Parameters
    ----------
    arr : list or array-like
    others : list or array, or list of arrays
    sep : string or None, default None
    na_rep : string or None, default None
        If None, an NA in any array will propagate

    Returns
    -------
    concat : array
    """
    if sep is None:
        sep = ''

    if others is not None:
        arrays = _get_array_list(arr, others)

        n = _length_check(arrays)
        masks = np.array([isnull(x) for x in arrays])
        cats = None

        if na_rep is None:
            na_mask = np.logical_or.reduce(masks, axis=0)

            result = np.empty(n, dtype=object)
            np.putmask(result, na_mask, np.nan)

            notmask = -na_mask

            tuples = izip(*[x[notmask] for x in arrays])
            cats = [sep.join(tup) for tup in tuples]

            result[notmask] = cats
        else:
            for i, x in enumerate(arrays):
                x = np.where(masks[i], na_rep, x)
                if cats is None:
                    cats = x
                else:
                    cats = cats + sep + x

            result = cats

        return result
    else:
        arr = np.asarray(arr, dtype=object)
        mask = isnull(arr)
        if na_rep is None and mask.any():
            return np.nan
        return sep.join(np.where(mask, na_rep, arr))


def _length_check(others):
    n = None
    for x in others:
        if n is None:
            n = len(x)
        elif len(x) != n:
            raise ValueError('All arrays must be same length')

    return n


def _na_map(f, arr, na_result=np.nan):
    # should really _check_ for NA
    return _map(f, arr, na_mask=True, na_value=na_result)


def _map(f, arr, na_mask=False, na_value=np.nan):
    if not isinstance(arr, np.ndarray):
        arr = np.asarray(arr, dtype=object)
    if na_mask:
        mask = isnull(arr)
        try:
            result = lib.map_infer_mask(arr, f, mask.view(np.uint8))
        except (TypeError, AttributeError):
            def g(x):
                try:
                    return f(x)
                except (TypeError, AttributeError):
                    return na_value
            return _map(g, arr)
        if na_value is not np.nan:
            np.putmask(result, mask, na_value)
            if result.dtype == object:
                result = lib.maybe_convert_objects(result)
        return result
    else:
        return lib.map_infer(arr, f)

def str_title(arr):
    """
    Convert strings to titlecased version

    Returns
    -------
    titled : array
    """
    return _na_map(lambda x: x.title(), arr)


def str_count(arr, pat, flags=0):
    """
    Count occurrences of pattern in each string

    Parameters
    ----------
    arr : list or array-like
    pat : string, valid regular expression
    flags : int, default 0 (no flags)
        re module flags, e.g. re.IGNORECASE

    Returns
    -------
    counts : arrays
    """
    regex = re.compile(pat, flags=flags)
    f = lambda x: len(regex.findall(x))
    return _na_map(f, arr)


def str_contains(arr, pat, case=True, flags=0, na=np.nan):
    """
    Check whether given pattern is contained in each string in the array

    Parameters
    ----------
    pat : string
        Character sequence or regular expression
    case : boolean, default True
        If True, case sensitive
    flags : int, default 0 (no flags)
        re module flags, e.g. re.IGNORECASE
    na : default NaN, fill value for missing values.

    Returns
    -------

    """
    if not case:
        flags |= re.IGNORECASE

    regex = re.compile(pat, flags=flags)

    f = lambda x: bool(regex.search(x))
    return _na_map(f, arr, na)


def str_startswith(arr, pat, na=np.nan):
    """
    Return boolean array indicating whether each string starts with passed
    pattern

    Parameters
    ----------
    pat : string
        Character sequence
    na : bool, default NaN

    Returns
    -------
    startswith : array (boolean)
    """
    f = lambda x: x.startswith(pat)
    return _na_map(f, arr, na)


def str_endswith(arr, pat, na=np.nan):
    """
    Return boolean array indicating whether each string ends with passed
    pattern

    Parameters
    ----------
    pat : string
        Character sequence
    na : bool, default NaN

    Returns
    -------
    endswith : array (boolean)
    """
    f = lambda x: x.endswith(pat)
    return _na_map(f, arr, na)


def str_lower(arr):
    """
    Convert strings in array to lowercase

    Returns
    -------
    lowercase : array
    """
    return _na_map(lambda x: x.lower(), arr)


def str_upper(arr):
    """
    Convert strings in array to uppercase

    Returns
    -------
    uppercase : array
    """
    return _na_map(lambda x: x.upper(), arr)


def str_replace(arr, pat, repl, n=-1, case=True, flags=0):
    """
    Replace

    Parameters
    ----------
    pat : string
        Character sequence or regular expression
    repl : string
        Replacement sequence
    n : int, default -1 (all)
        Number of replacements to make from start
    case : boolean, default True
        If True, case sensitive
    flags : int, default 0 (no flags)
        re module flags, e.g. re.IGNORECASE

    Returns
    -------
    replaced : array
    """
    use_re = not case or len(pat) > 1 or flags

    if use_re:
        if not case:
            flags |= re.IGNORECASE
        regex = re.compile(pat, flags=flags)
        n = n if n >= 0 else 0

        def f(x):
            return regex.sub(repl, x, count=n)
    else:
        f = lambda x: x.replace(pat, repl, n)

    return _na_map(f, arr)


def str_repeat(arr, repeats):
    """
    Duplicate each string in the array by indicated number of times

    Parameters
    ----------
    repeats : int or array
        Same value for all (int) or different value per (array)

    Returns
    -------
    repeated : array
    """
    if np.isscalar(repeats):
        def rep(x):
            try:
                return str.__mul__(x, repeats)
            except TypeError:
                return unicode.__mul__(x, repeats)
        return _na_map(rep, arr)
    else:
        def rep(x, r):
            try:
                return str.__mul__(x, r)
            except TypeError:
                return unicode.__mul__(x, r)
        repeats = np.asarray(repeats, dtype=object)
        result = lib.vec_binop(arr, repeats, rep)
        return result


def str_match(arr, pat, flags=0):
    """
    Find groups in each string (from beginning) using passed regular expression

    Parameters
    ----------
    pat : string
        Pattern or regular expression
    flags : int, default 0 (no flags)
        re module flags, e.g. re.IGNORECASE

    Returns
    -------
    matches : array
    """
    regex = re.compile(pat, flags=flags)

    def f(x):
        m = regex.match(x)
        if m:
            return m.groups()
        else:
            return []

    return _na_map(f, arr)


def str_join(arr, sep):
    """
    Join lists contained as elements in array, a la str.join

    Parameters
    ----------
    sep : string
        Delimiter

    Returns
    -------
    joined : array
    """
    return _na_map(sep.join, arr)


def str_len(arr):
    """
    Compute length of each string in array.

    Returns
    -------
    lengths : array
    """
    return _na_map(len, arr)


def str_findall(arr, pat, flags=0):
    """
    Find all occurrences of pattern or regular expression

    Parameters
    ----------
    pat : string
        Pattern or regular expression
    flags : int, default 0 (no flags)
        re module flags, e.g. re.IGNORECASE

    Returns
    -------
    matches : array
    """
    regex = re.compile(pat, flags=flags)
    return _na_map(regex.findall, arr)


def str_pad(arr, width, side='left'):
    """
    Pad strings with whitespace

    Parameters
    ----------
    arr : list or array-like
    width : int
        Minimum width of resulting string; additional characters will be filled
        with spaces
    side : {'left', 'right', 'both'}, default 'left'

    Returns
    -------
    padded : array
    """
    if side == 'left':
        f = lambda x: x.rjust(width)
    elif side == 'right':
        f = lambda x: x.ljust(width)
    elif side == 'both':
        f = lambda x: x.center(width)
    else:  # pragma: no cover
        raise ValueError('Invalid side')

    return _na_map(f, arr)


def str_center(arr, width):
    """
    "Center" strings, filling left and right side with additional whitespace

    Parameters
    ----------
    width : int
        Minimum width of resulting string; additional characters will be filled
        with spaces

    Returns
    -------
    centered : array
    """
    return str_pad(arr, width, side='both')


def str_split(arr, pat=None, n=None):
    """
    Split each string (a la re.split) in array by given pattern, propagating NA
    values

    Parameters
    ----------
    pat : string, default None
        String or regular expression to split on. If None, splits on whitespace
    n : int, default None (all)

    Notes
    -----
    Both 0 and -1 will be interpreted as return all splits

    Returns
    -------
    split : array
    """
    if pat is None:
        if n is None or n == 0:
            n = -1
        f = lambda x: x.split()
    else:
        if len(pat) == 1:
            if n is None or n == 0:
                n = -1
            f = lambda x: x.split(pat, n)
        else:
            if n is None or n == -1:
                n = 0
            regex = re.compile(pat)
            f = lambda x: regex.split(x, maxsplit=n)

    return _na_map(f, arr)


def str_slice(arr, start=None, stop=None, step=1):
    """
    Slice substrings from each element in array

    Parameters
    ----------
    start : int or None
    stop : int or None

    Returns
    -------
    sliced : array
    """
    obj = slice(start, stop, step)
    f = lambda x: x[obj]
    return _na_map(f, arr)


def str_slice_replace(arr, start=None, stop=None, repl=None):
    """

    Parameters
    ----------

    Returns
    -------
    replaced : array
    """
    raise NotImplementedError


def str_strip(arr, to_strip=None):
    """
    Strip whitespace (including newlines) from each string in the array

    Parameters
    ----------
    to_strip : str or unicode

    Returns
    -------
    stripped : array
    """
    return _na_map(lambda x: x.strip(to_strip), arr)


def str_lstrip(arr, to_strip=None):
    """
    Strip whitespace (including newlines) from left side of each string in the
    array

    Parameters
    ----------
    to_strip : str or unicode

    Returns
    -------
    stripped : array
    """
    return _na_map(lambda x: x.lstrip(to_strip), arr)


def str_rstrip(arr, to_strip=None):
    """
    Strip whitespace (including newlines) from right side of each string in the
    array

    Parameters
    ----------
    to_strip : str or unicode

    Returns
    -------
    stripped : array
    """
    return _na_map(lambda x: x.rstrip(to_strip), arr)


def str_wrap(arr, width=80):
    """
    Wrap long strings to be formatted in paragraphs

    Parameters
    ----------
    width : int
        Maximum line-width

    Returns
    -------
    wrapped : array
    """
    raise NotImplementedError


def str_get(arr, i):
    """
    Extract element from lists, tuples, or strings in each element in the array

    Parameters
    ----------
    i : int
        Integer index (location)

    Returns
    -------
    items : array
    """
    f = lambda x: x[i] if len(x) > i else np.nan
    return _na_map(f, arr)


def str_decode(arr, encoding, errors="strict"):
    """
    Decode character string to unicode using indicated encoding

    Parameters
    ----------
    encoding : string
    errors : string

    Returns
    -------
    decoded : array
    """
    f = lambda x: x.decode(encoding, errors)
    return _na_map(f, arr)


def str_encode(arr, encoding, errors="strict"):
    """
    Encode character string to some other encoding using indicated encoding

    Parameters
    ----------
    encoding : string
    errors : string

    Returns
    -------
    encoded : array
    """
    f = lambda x: x.encode(encoding, errors)
    return _na_map(f, arr)


def _noarg_wrapper(f):
    def wrapper(self):
        result = f(self.series)
        return self._wrap_result(result)

    wrapper.__name__ = f.__name__
    if f.__doc__:
        wrapper.__doc__ = f.__doc__

    return wrapper


def _pat_wrapper(f, flags=False, na=False):
    def wrapper1(self, pat):
        result = f(self.series, pat)
        return self._wrap_result(result)

    def wrapper2(self, pat, flags=0):
        result = f(self.series, pat, flags=flags)
        return self._wrap_result(result)

    def wrapper3(self, pat, na=np.nan):
        result = f(self.series, pat, na=na)
        return self._wrap_result(result)

    wrapper = wrapper3 if na else wrapper2 if flags else wrapper1

    wrapper.__name__ = f.__name__
    if f.__doc__:
        wrapper.__doc__ = f.__doc__

    return wrapper


def copy(source):
    "Copy a docstring from another source function (if present)"
    def do_copy(target):
        if source.__doc__:
            target.__doc__ = source.__doc__
        return target
    return do_copy


class StringMethods(object):
    """
    Vectorized string functions for Series. NAs stay NA unless handled
    otherwise by a particular method. Patterned after Python's string methods,
    with some inspiration from R's stringr package.

    Examples
    --------
    >>> s.str.split('_')
    >>> s.str.replace('_', '')
    """
    def __init__(self, series):
        self.series = series

    def __getitem__(self, key):
        if isinstance(key, slice):
            return self.slice(start=key.start, stop=key.stop,
                              step=key.step)
        else:
            return self.get(key)

    def __iter__(self):
        i = 0
        g = self.get(i)
        while g.notnull().any():
            yield g
            i += 1
            g = self.get(i)

    def _wrap_result(self, result):
        return Series(result, index=self.series.index,
                      name=self.series.name)

    @copy(str_cat)
    def cat(self, others=None, sep=None, na_rep=None):
        result = str_cat(self.series, others=others, sep=sep, na_rep=na_rep)
        return self._wrap_result(result)

    @copy(str_split)
    def split(self, pat=None, n=-1):
        result = str_split(self.series, pat, n=n)
        return self._wrap_result(result)

    @copy(str_get)
    def get(self, i):
        result = str_get(self.series, i)
        return self._wrap_result(result)

    @copy(str_join)
    def join(self, sep):
        result = str_join(self.series, sep)
        return self._wrap_result(result)

    @copy(str_contains)
    def contains(self, pat, case=True, flags=0, na=np.nan):
        result = str_contains(self.series, pat, case=case, flags=flags,
                              na=na)
        return self._wrap_result(result)

    @copy(str_replace)
    def replace(self, pat, repl, n=-1, case=True, flags=0):
        result = str_replace(self.series, pat, repl, n=n, case=case,
                             flags=flags)
        return self._wrap_result(result)

    @copy(str_repeat)
    def repeat(self, repeats):
        result = str_repeat(self.series, repeats)
        return self._wrap_result(result)

    @copy(str_pad)
    def pad(self, width, side='left'):
        result = str_pad(self.series, width, side=side)
        return self._wrap_result(result)

    @copy(str_center)
    def center(self, width):
        result = str_center(self.series, width)
        return self._wrap_result(result)

    @copy(str_slice)
    def slice(self, start=None, stop=None, step=1):
        result = str_slice(self.series, start, stop)
        return self._wrap_result(result)

    @copy(str_slice)
    def slice_replace(self, i=None, j=None):
        raise NotImplementedError

    @copy(str_decode)
    def decode(self, encoding, errors="strict"):
        result = str_decode(self.series, encoding, errors)
        return self._wrap_result(result)

    @copy(str_encode)
    def encode(self, encoding, errors="strict"):
        result = str_encode(self.series, encoding, errors)
        return self._wrap_result(result)

    @copy(str_strip)
    def strip(self, to_strip=None):
        result = str_strip(self.series, to_strip)
        return self._wrap_result(result)

    @copy(str_lstrip)
    def lstrip(self, to_strip=None):
        result = str_lstrip(self.series, to_strip)
        return self._wrap_result(result)

    @copy(str_rstrip)
    def rstrip(self, to_strip=None):
        result = str_rstrip(self.series, to_strip)
        return self._wrap_result(result)

    count = _pat_wrapper(str_count, flags=True)
    startswith = _pat_wrapper(str_startswith, na=True)
    endswith = _pat_wrapper(str_endswith, na=True)
    findall = _pat_wrapper(str_findall, flags=True)
    match = _pat_wrapper(str_match, flags=True)

    len = _noarg_wrapper(str_len)
    lower = _noarg_wrapper(str_lower)
    upper = _noarg_wrapper(str_upper)
    title = _noarg_wrapper(str_title)
