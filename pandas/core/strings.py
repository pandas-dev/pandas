import numpy as np

from pandas.compat import zip
from pandas.core.common import isnull, _values_from_object, is_bool_dtype
import pandas.compat as compat
from pandas.util.decorators import Appender
import re
import pandas.lib as lib
import warnings
import textwrap


_shared_docs = dict()


def _get_array_list(arr, others):
    from pandas.core.series import Series

    if len(others) and isinstance(_values_from_object(others)[0],
                                  (list, np.ndarray, Series)):
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

            notmask = ~na_mask

            tuples = zip(*[x[notmask] for x in arrays])
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


def _na_map(f, arr, na_result=np.nan, dtype=object):
    # should really _check_ for NA
    return _map(f, arr, na_mask=True, na_value=na_result, dtype=dtype)


def _map(f, arr, na_mask=False, na_value=np.nan, dtype=object):
    from pandas.core.series import Series

    if not len(arr):
        return np.ndarray(0, dtype=dtype)

    if isinstance(arr, Series):
        arr = arr.values
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
            return _map(g, arr, dtype=dtype)
        if na_value is not np.nan:
            np.putmask(result, mask, na_value)
            if result.dtype == object:
                result = lib.maybe_convert_objects(result)
        return result
    else:
        return lib.map_infer(arr, f)


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
    return _na_map(f, arr, dtype=int)


def str_contains(arr, pat, case=True, flags=0, na=np.nan, regex=True):
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
    regex : bool, default True
        If True use re.search, otherwise use Python in operator

    Returns
    -------
    Series of boolean values

    See Also
    --------
    match : analagous, but stricter, relying on re.match instead of re.search

    """
    if regex:
        if not case:
            flags |= re.IGNORECASE

        regex = re.compile(pat, flags=flags)

        if regex.groups > 0:
            warnings.warn("This pattern has match groups. To actually get the"
                          " groups, use str.extract.", UserWarning)

        f = lambda x: bool(regex.search(x))
    else:
        if case:
            f = lambda x: pat in x
        else:
            upper_pat = pat.upper()
            f = lambda x: upper_pat in x
            uppered = _na_map(lambda x: x.upper(), arr)
            return _na_map(f, uppered, na, dtype=bool)
    return _na_map(f, arr, na, dtype=bool)


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
    return _na_map(f, arr, na, dtype=bool)


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
    return _na_map(f, arr, na, dtype=bool)


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
                return compat.binary_type.__mul__(x, repeats)
            except TypeError:
                return compat.text_type.__mul__(x, repeats)

        return _na_map(rep, arr)
    else:
        def rep(x, r):
            try:
                return compat.binary_type.__mul__(x, r)
            except TypeError:
                return compat.text_type.__mul__(x, r)

        repeats = np.asarray(repeats, dtype=object)
        result = lib.vec_binop(_values_from_object(arr), repeats, rep)
        return result


def str_match(arr, pat, case=True, flags=0, na=np.nan, as_indexer=False):
    """
    Deprecated: Find groups in each string using passed regular expression.
    If as_indexer=True, determine if each string matches a regular expression.

    Parameters
    ----------
    pat : string
        Character sequence or regular expression
    case : boolean, default True
        If True, case sensitive
    flags : int, default 0 (no flags)
        re module flags, e.g. re.IGNORECASE
    na : default NaN, fill value for missing values.
    as_indexer : False, by default, gives deprecated behavior better achieved
        using str_extract. True return boolean indexer.

    Returns
    -------
    Series of boolean values
        if as_indexer=True
    Series of tuples
        if as_indexer=False, default but deprecated

    See Also
    --------
    contains : analagous, but less strict, relying on re.search instead of
        re.match
    extract : now preferred to the deprecated usage of match (as_indexer=False)

    Notes
    -----
    To extract matched groups, which is the deprecated behavior of match, use
    str.extract.
    """

    if not case:
        flags |= re.IGNORECASE

    regex = re.compile(pat, flags=flags)

    if (not as_indexer) and regex.groups > 0:
        # Do this first, to make sure it happens even if the re.compile
        # raises below.
        warnings.warn("In future versions of pandas, match will change to"
                      " always return a bool indexer.", UserWarning)

    if as_indexer and regex.groups > 0:
        warnings.warn("This pattern has match groups. To actually get the"
                      " groups, use str.extract.", UserWarning)

    # If not as_indexer and regex.groups == 0, this returns empty lists
    # and is basically useless, so we will not warn.

    if (not as_indexer) and regex.groups > 0:
        dtype = object
        def f(x):
            m = regex.match(x)
            if m:
                return m.groups()
            else:
                return []
    else:
        # This is the new behavior of str_match.
        dtype = bool
        f = lambda x: bool(regex.match(x))

    return _na_map(f, arr, na, dtype=dtype)


def _get_single_group_name(rx):
    try:
        return list(rx.groupindex.keys()).pop()
    except IndexError:
        return None


def str_extract(arr, pat, flags=0):
    """
    Find groups in each string using passed regular expression

    Parameters
    ----------
    pat : string
        Pattern or regular expression
    flags : int, default 0 (no flags)
        re module flags, e.g. re.IGNORECASE

    Returns
    -------
    extracted groups : Series (one group) or DataFrame (multiple groups)
        Note that dtype of the result is always object, even when no match is
        found and the result is a Series or DataFrame containing only NaN
        values.

    Examples
    --------
    A pattern with one group will return a Series. Non-matches will be NaN.

    >>> Series(['a1', 'b2', 'c3']).str.extract('[ab](\d)')
    0      1
    1      2
    2    NaN
    dtype: object

    A pattern with more than one group will return a DataFrame.

    >>> Series(['a1', 'b2', 'c3']).str.extract('([ab])(\d)')
         0    1
    0    a    1
    1    b    2
    2  NaN  NaN

    A pattern may contain optional groups.

    >>> Series(['a1', 'b2', 'c3']).str.extract('([ab])?(\d)')
         0  1
    0    a  1
    1    b  2
    2  NaN  3

    Named groups will become column names in the result.

    >>> Series(['a1', 'b2', 'c3']).str.extract('(?P<letter>[ab])(?P<digit>\d)')
      letter digit
    0      a     1
    1      b     2
    2    NaN   NaN

    """
    from pandas.core.series import Series
    from pandas.core.frame import DataFrame

    regex = re.compile(pat, flags=flags)
    # just to be safe, check this
    if regex.groups == 0:
        raise ValueError("This pattern contains no groups to capture.")
    empty_row = [np.nan]*regex.groups
    def f(x):
        if not isinstance(x, compat.string_types):
            return empty_row
        m = regex.search(x)
        if m:
            return [np.nan if item is None else item for item in m.groups()]
        else:
            return empty_row
    if regex.groups == 1:
        result = Series([f(val)[0] for val in arr],
                        name=_get_single_group_name(regex),
                        index=arr.index, dtype=object)
    else:
        names = dict(zip(regex.groupindex.values(), regex.groupindex.keys()))
        columns = [names.get(1 + i, i) for i in range(regex.groups)]
        if arr.empty:
            result = DataFrame(columns=columns, dtype=object)
        else:
            result = DataFrame([f(val) for val in arr],
                               columns=columns,
                               index=arr.index,
                               dtype=object)
    return result


def str_get_dummies(arr, sep='|'):
    """
    Split each string by sep and return a frame of dummy/indicator variables.

    Examples
    --------
    >>> Series(['a|b', 'a', 'a|c']).str.get_dummies()
       a  b  c
    0  1  1  0
    1  1  0  0
    2  1  0  1

    >>> pd.Series(['a|b', np.nan, 'a|c']).str.get_dummies()
       a  b  c
    0  1  1  0
    1  0  0  0
    2  1  0  1

    See also ``pd.get_dummies``.

    """
    from pandas.core.frame import DataFrame

    # TODO remove this hack?
    arr = arr.fillna('')
    try:
        arr = sep + arr + sep
    except TypeError:
        arr = sep + arr.astype(str) + sep

    tags = set()
    for ts in arr.str.split(sep):
        tags.update(ts)
    tags = sorted(tags - set([""]))

    dummies = np.empty((len(arr), len(tags)), dtype=np.int64)

    for i, t in enumerate(tags):
        pat = sep + t + sep
        dummies[:, i] = lib.map_infer(arr.values, lambda x: pat in x)
    return DataFrame(dummies, arr.index, tags)


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


def str_find(arr, sub, start=0, end=None, side='left'):
    """
    Return indexes in each strings where the substring is
    fully contained between [start:end]. Return -1 on failure.

    Parameters
    ----------
    sub : str
        Substring being searched
    start : int
        Left edge index
    end : int
        Right edge index
    side : {'left', 'right'}, default 'left'
        Specifies a starting side, equivalent to ``find`` or ``rfind``

    Returns
    -------
    found : array
    """

    if not isinstance(sub, compat.string_types):
        msg = 'expected a string object, not {0}'
        raise TypeError(msg.format(type(sub).__name__))

    if side == 'left':
        method = 'find'
    elif side == 'right':
        method = 'rfind'
    else:  # pragma: no cover
        raise ValueError('Invalid side')

    if end is None:
        f = lambda x: getattr(x, method)(sub, start)
    else:
        f = lambda x: getattr(x, method)(sub, start, end)

    return _na_map(f, arr, dtype=int)


def str_pad(arr, width, side='left', fillchar=' '):
    """
    Pad strings with an additional character

    Parameters
    ----------
    arr : list or array-like
    width : int
        Minimum width of resulting string; additional characters will be filled
        with spaces
    side : {'left', 'right', 'both'}, default 'left'
    fillchar : str
        Additional character for filling, default is whitespace

    Returns
    -------
    padded : array
    """

    if not isinstance(fillchar, compat.string_types):
        msg = 'fillchar must be a character, not {0}'
        raise TypeError(msg.format(type(fillchar).__name__))

    if len(fillchar) != 1:
        raise TypeError('fillchar must be a character, not str')

    if side == 'left':
        f = lambda x: x.rjust(width, fillchar)
    elif side == 'right':
        f = lambda x: x.ljust(width, fillchar)
    elif side == 'both':
        f = lambda x: x.center(width, fillchar)
    else:  # pragma: no cover
        raise ValueError('Invalid side')

    return _na_map(f, arr)


def str_split(arr, pat=None, n=None, return_type='series'):
    """
    Split each string (a la re.split) in array by given pattern, propagating NA
    values

    Parameters
    ----------
    pat : string, default None
        String or regular expression to split on. If None, splits on whitespace
    n : int, default None (all)
    return_type : {'series', 'index', 'frame'}, default 'series'
        If frame, returns a DataFrame (elements are strings)
        If series or index, returns the same type as the original object
        (elements are lists of strings).

    Notes
    -----
    Both 0 and -1 will be interpreted as return all splits

    Returns
    -------
    split : array
    """
    from pandas.core.series import Series
    from pandas.core.frame import DataFrame
    from pandas.core.index import Index

    if return_type not in ('series', 'index', 'frame'):
        raise ValueError("return_type must be {'series', 'index', 'frame'}")
    if return_type == 'frame' and isinstance(arr, Index):
        raise ValueError("return_type='frame' is not supported for string "
                         "methods on Index")
    if pat is None:
        if n is None or n == 0:
            n = -1
        f = lambda x: x.split(pat, n)
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
    if return_type == 'frame':
        res = DataFrame((Series(x) for x in _na_map(f, arr)), index=arr.index)
    else:
        res = _na_map(f, arr)
    return res


def str_slice(arr, start=None, stop=None, step=None):
    """
    Slice substrings from each element in array

    Parameters
    ----------
    start : int or None
    stop : int or None
    step : int or None

    Returns
    -------
    sliced : array
    """
    obj = slice(start, stop, step)
    f = lambda x: x[obj]
    return _na_map(f, arr)


def str_slice_replace(arr, start=None, stop=None, repl=None):
    """
    Replace a slice of each string with another string.

    Parameters
    ----------
    start : int or None
    stop : int or None
    repl : str or None

    Returns
    -------
    replaced : array
    """
    if repl is None:
        repl = ''

    def f(x):
        if x[start:stop] == '':
            local_stop = start
        else:
            local_stop = stop
        y = ''
        if start is not None:
            y += x[:start]
        y += repl
        if stop is not None:
            y += x[local_stop:]
        return y
    return _na_map(f, arr)


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


def str_wrap(arr, width, **kwargs):
    r"""
    Wrap long strings to be formatted in paragraphs.

    This method has the same keyword parameters and defaults as
    :class:`textwrap.TextWrapper`.

    Parameters
    ----------
    width : int
        Maximum line-width
    expand_tabs : bool, optional
        If true, tab characters will be expanded to spaces (default: True)
    replace_whitespace : bool, optional
        If true, each whitespace character (as defined by string.whitespace) remaining
        after tab expansion will be replaced by a single space (default: True)
    drop_whitespace : bool, optional
        If true, whitespace that, after wrapping, happens to end up at the beginning
        or end of a line is dropped (default: True)
    break_long_words : bool, optional
        If true, then words longer than width will be broken in order to ensure that
        no lines are longer than width. If it is false, long words will not be broken,
        and some lines may be longer than width. (default: True)
    break_on_hyphens : bool, optional
        If true, wrapping will occur preferably on whitespace and right after hyphens
        in compound words, as it is customary in English. If false, only whitespaces
        will be considered as potentially good places for line breaks, but you need
        to set break_long_words to false if you want truly insecable words.
        (default: True)

    Returns
    -------
    wrapped : array

    Notes
    -----
    Internally, this method uses a :class:`textwrap.TextWrapper` instance with default
    settings. To achieve behavior matching R's stringr library str_wrap function, use
    the arguments:

    - expand_tabs = False
    - replace_whitespace = True
    - drop_whitespace = True
    - break_long_words = False
    - break_on_hyphens = False

    Examples
    --------

    >>> s = pd.Series(['line to be wrapped', 'another line to be wrapped'])
    >>> s.str.wrap(12)
    0             line to be\nwrapped
    1    another line\nto be\nwrapped
    """
    kwargs['width'] = width

    tw = textwrap.TextWrapper(**kwargs)

    return _na_map(lambda s: '\n'.join(tw.wrap(s)), arr)


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


def _noarg_wrapper(f, docstring=None, **kargs):
    def wrapper(self):
        result = _na_map(f, self.series, **kargs)
        return self._wrap_result(result)

    wrapper.__name__ = f.__name__
    if docstring is not None:
        wrapper.__doc__ = docstring
    else:
        raise ValueError('Provide docstring')

    return wrapper


def _pat_wrapper(f, flags=False, na=False, **kwargs):
    def wrapper1(self, pat):
        result = f(self.series, pat)
        return self._wrap_result(result)

    def wrapper2(self, pat, flags=0, **kwargs):
        result = f(self.series, pat, flags=flags, **kwargs)
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
    Vectorized string functions for Series and Index. NAs stay NA unless
    handled otherwise by a particular method. Patterned after Python's string
    methods, with some inspiration from R's stringr package.

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
        from pandas.core.series import Series
        from pandas.core.frame import DataFrame
        from pandas.core.index import Index

        if not hasattr(result, 'ndim'):
            return result
        elif result.ndim == 1:
            name = getattr(result, 'name', None)
            if isinstance(self.series, Index):
                # if result is a boolean np.array, return the np.array
                # instead of wrapping it into a boolean Index (GH 8875)
                if is_bool_dtype(result):
                    return result
                return Index(result, name=name or self.series.name)
            return Series(result, index=self.series.index,
                          name=name or self.series.name)
        else:
            assert result.ndim < 3
            return DataFrame(result, index=self.series.index)

    @copy(str_cat)
    def cat(self, others=None, sep=None, na_rep=None):
        result = str_cat(self.series, others=others, sep=sep, na_rep=na_rep)
        return self._wrap_result(result)

    @copy(str_split)
    def split(self, pat=None, n=-1, return_type='series'):
        result = str_split(self.series, pat, n=n, return_type=return_type)
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
    def contains(self, pat, case=True, flags=0, na=np.nan, regex=True):
        result = str_contains(self.series, pat, case=case, flags=flags,
                              na=na, regex=regex)
        return self._wrap_result(result)

    @copy(str_match)
    def match(self, pat, case=True, flags=0, na=np.nan, as_indexer=False):
        result = str_match(self.series, pat, case=case, flags=flags,
                              na=na, as_indexer=as_indexer)
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
    def pad(self, width, side='left', fillchar=' '):
        result = str_pad(self.series, width, side=side, fillchar=fillchar)
        return self._wrap_result(result)

    _shared_docs['str_pad'] = ("""
    Filling %s side of strings with an additional character

    Parameters
    ----------
    width : int
        Minimum width of resulting string; additional characters will be filled
        with ``fillchar``
    fillchar : str
        Additional character for filling, default is whitespace

    Returns
    -------
    filled : array
    """)

    @Appender(_shared_docs['str_pad'] % 'left and right')
    def center(self, width, fillchar=' '):
        return self.pad(width, side='both', fillchar=fillchar)

    @Appender(_shared_docs['str_pad'] % 'right')
    def ljust(self, width, fillchar=' '):
        return self.pad(width, side='right', fillchar=fillchar)

    @Appender(_shared_docs['str_pad'] % 'left')
    def rjust(self, width, fillchar=' '):
        return self.pad(width, side='left', fillchar=fillchar)

    def zfill(self, width):
        """"
        Filling left side with 0

        Parameters
        ----------
        width : int
            Minimum width of resulting string; additional characters will be filled
            with 0

        Returns
        -------
        filled : array
        """
        result = str_pad(self.series, width, side='left', fillchar='0')
        return self._wrap_result(result)

    @copy(str_slice)
    def slice(self, start=None, stop=None, step=None):
        result = str_slice(self.series, start, stop, step)
        return self._wrap_result(result)

    @copy(str_slice_replace)
    def slice_replace(self, start=None, stop=None, repl=None):
        result = str_slice_replace(self.series, start, stop, repl)
        return self._wrap_result(result)

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

    @copy(str_wrap)
    def wrap(self, width, **kwargs):
        result = str_wrap(self.series, width, **kwargs)
        return self._wrap_result(result)

    @copy(str_get_dummies)
    def get_dummies(self, sep='|'):
        result = str_get_dummies(self.series, sep)
        return self._wrap_result(result)

    count = _pat_wrapper(str_count, flags=True)
    startswith = _pat_wrapper(str_startswith, na=True)
    endswith = _pat_wrapper(str_endswith, na=True)
    findall = _pat_wrapper(str_findall, flags=True)
    extract = _pat_wrapper(str_extract, flags=True)

    _shared_docs['find'] = ("""
    Return %(side)s indexes in each strings where the substring is
    fully contained between [start:end]. Return -1 on failure.
    Equivalent to standard ``str.%(method)s``.

    Parameters
    ----------
    sub : str
        Substring being searched
    start : int
        Left edge index
    end : int
        Right edge index

    Returns
    -------
    found : array

    See Also
    --------
    %(also)s
    """)

    @Appender(_shared_docs['find'] % dict(side='lowest', method='find',
              also='rfind : Return highest indexes in each strings'))
    def find(self, sub, start=0, end=None):
        result = str_find(self.series, sub, start=start, end=end, side='left')
        return self._wrap_result(result)

    @Appender(_shared_docs['find'] % dict(side='highest', method='rfind',
              also='find : Return lowest indexes in each strings'))
    def rfind(self, sub, start=0, end=None):
        result = str_find(self.series, sub, start=start, end=end, side='right')
        return self._wrap_result(result)

    _shared_docs['len'] = ("""
    Compute length of each string in array.

    Returns
    -------
    lengths : array
    """)
    len = _noarg_wrapper(len, docstring=_shared_docs['len'], dtype=int)

    _shared_docs['casemethods'] = ("""
    Convert strings in array to %(type)s.
    Equivalent to ``str.%(method)s``.

    Returns
    -------
    converted : array
    """)
    _shared_docs['lower'] = dict(type='lowercase', method='lower')
    _shared_docs['upper'] = dict(type='uppercase', method='upper')
    _shared_docs['title'] = dict(type='titlecase', method='title')
    _shared_docs['capitalize'] = dict(type='be capitalized', method='capitalize')
    _shared_docs['swapcase'] = dict(type='be swapcased', method='swapcase')
    lower = _noarg_wrapper(lambda x: x.lower(),
                           docstring=_shared_docs['casemethods'] % _shared_docs['lower'])
    upper = _noarg_wrapper(lambda x: x.upper(),
                           docstring=_shared_docs['casemethods'] % _shared_docs['upper'])
    title = _noarg_wrapper(lambda x: x.title(),
                           docstring=_shared_docs['casemethods'] % _shared_docs['title'])
    capitalize = _noarg_wrapper(lambda x: x.capitalize(),
                                docstring=_shared_docs['casemethods'] % _shared_docs['capitalize'])
    swapcase = _noarg_wrapper(lambda x: x.swapcase(),
                              docstring=_shared_docs['casemethods'] % _shared_docs['swapcase'])

    _shared_docs['ismethods'] = ("""
    Check whether all characters in each string in the array are %(type)s.
    Equivalent to ``str.%(method)s``.

    Returns
    -------
    Series of boolean values
    """)
    _shared_docs['isalnum'] = dict(type='alphanumeric', method='isalnum')
    _shared_docs['isalpha'] = dict(type='alphabetic', method='isalpha')
    _shared_docs['isdigit'] = dict(type='digits', method='isdigit')
    _shared_docs['isspace'] = dict(type='whitespace', method='isspace')
    _shared_docs['islower'] = dict(type='lowercase', method='islower')
    _shared_docs['isupper'] = dict(type='uppercase', method='isupper')
    _shared_docs['istitle'] = dict(type='titlecase', method='istitle')
    _shared_docs['isnumeric'] = dict(type='numeric', method='isnumeric')
    _shared_docs['isdecimal'] = dict(type='decimal', method='isdecimal')
    isalnum = _noarg_wrapper(lambda x: x.isalnum(),
                             docstring=_shared_docs['ismethods'] % _shared_docs['isalnum'])
    isalpha = _noarg_wrapper(lambda x: x.isalpha(),
                             docstring=_shared_docs['ismethods'] % _shared_docs['isalpha'])
    isdigit = _noarg_wrapper(lambda x: x.isdigit(),
                             docstring=_shared_docs['ismethods'] % _shared_docs['isdigit'])
    isspace = _noarg_wrapper(lambda x: x.isspace(),
                             docstring=_shared_docs['ismethods'] % _shared_docs['isspace'])
    islower = _noarg_wrapper(lambda x: x.islower(),
                             docstring=_shared_docs['ismethods'] % _shared_docs['islower'])
    isupper = _noarg_wrapper(lambda x: x.isupper(),
                             docstring=_shared_docs['ismethods'] % _shared_docs['isupper'])
    istitle = _noarg_wrapper(lambda x: x.istitle(),
                             docstring=_shared_docs['ismethods'] % _shared_docs['istitle'])
    isnumeric = _noarg_wrapper(lambda x: compat.u_safe(x).isnumeric(),
                             docstring=_shared_docs['ismethods'] % _shared_docs['isnumeric'])
    isdecimal = _noarg_wrapper(lambda x: compat.u_safe(x).isdecimal(),
                             docstring=_shared_docs['ismethods'] % _shared_docs['isdecimal'])
