""":mod:`pandas.io.html` is a module containing functionality for dealing with
HTML IO.

"""

import os
import re
import numbers
import urllib2
import contextlib
import collections
import urlparse

try:
    from importlib import import_module
except ImportError:
    import_module = __import__

from pandas import DataFrame, MultiIndex
from pandas.io.parsers import _is_url


#############
# READ HTML #
#############
_RE_WHITESPACE = re.compile(r'([\r\n]+|\s{2,})')


def _remove_whitespace(s, regex=_RE_WHITESPACE):
    """Replace extra whitespace inside of a string with a single space.

    Parameters
    ----------
    s : str or unicode
        The string from which to remove extra whitespace.

    regex : regex
        The regular expression to use to remove extra whitespace.

    Returns
    -------
    subd : str or unicode
        `s` with all extra whitespace replaced with a single space.
    """
    return regex.sub(' ', s.strip())


def _get_skiprows_iter(skiprows):
    """Get an iterator given an integer, slice or container.

    Parameters
    ----------
    skiprows : int, slice, container
        The iterator to use to skip rows; can also be a slice.

    Raises
    ------
    TypeError
        * If `skiprows` is not a slice, integer, or Container

    Raises
    ------
    TypeError
        * If `skiprows` is not a slice, integer, or Container

    Returns
    -------
    it : iterable
        A proper iterator to use to skip rows of a DataFrame.
    """
    if isinstance(skiprows, slice):
        return range(skiprows.start or 0, skiprows.stop, skiprows.step or 1)
    elif isinstance(skiprows, numbers.Integral):
        return range(skiprows)
    elif isinstance(skiprows, collections.Container):
        return skiprows
    else:
        raise TypeError('{0} is not a valid type for skipping'
                        ' rows'.format(type(skiprows)))

    def _parse_columns(self, row):
        return row.xpath('.//td|.//th')


class _HtmlFrameParser(object):
    """Base class for parsers that parse HTML into DataFrames.

    Parameters
    ----------
    io : str or file-like
        This can be either a string of raw HTML, a valid URL using the HTTP,
        FTP, or FILE protocols or a file-like object.

    match : str or regex
        The text to match in the document.

    attrs : dict
        List of HTML <table> element attributes to match.

    Attributes
    ----------
    io : str or file-like
        raw HTML, URL, or file-like object

    match : regex
        The text to match in the raw HTML

    attrs : dict-like
        A dictionary of valid table attributes to use to search for table
        elements.

    Notes
    -----
    To subclass this class effectively you must override the following methods:
        * :func:`_build_doc`
        * :func:`_text_getter`
        * :func:`_parse_columns`
        * :func:`_parse_table`
        * :func:`_parse_rows`
    See each method's respective documentation for details on their
    functionality.
    """
    def __init__(self, io, match, attrs):
        self.io = io
        self.match = match
        self.attrs = attrs

    def parse_rows(self):
        """Return a list of list of each table's rows.

        Returns
        -------
        row_list : list of list of node-like
            A list of each table's rows, which are DOM nodes (usually <th> or
            <tr> elements).
        """
        tables = self._parse_tables(self._build_doc(), self.match, self.attrs)
        assert tables, 'No tables found'
        return (self._parse_rows(table) for table in tables)

    def parse_raw_data(self):
        """Return a list of the raw data from each table.

        Returns
        -------
        data : list of list of lists of str or unicode
            Each table's data is contained in a list of lists of str or
            unicode.
        """
        return [self._parse_raw_data(rows, self._text_getter,
                                     self._parse_columns)
                for rows in self.parse_rows()]

    def _parse_raw_data(self, rows, text_getter, column_finder):
        """Parse the raw data into a list of lists.

        Parameters
        ----------
        rows : iterable of node-like
            A list of row elements.

        text_getter : callable
            A callable that gets the text from an individual node. This must be
            defined by subclasses.

        column_finder : callable
            A callable that takes a row node as input and returns a list of the
            column node in that row. This must be defined by subclasses.

        Raises
        ------
        AssertionError
            * If `text_getter` is not callable
            * If `column_finder` is not callable

        Returns
        -------
        data : list of list of strings
        """
        # callable is back in Python 3.2
        assert callable(text_getter), '"text_getter" must be callable'
        assert callable(column_finder), '"column_finder" must be callable'

        data = []

        for row in rows:
            if _remove_whitespace(text_getter(row)):
                col = []

                for el in column_finder(row):
                    t = _remove_whitespace(text_getter(el))

                    if t:
                        col.append(t)
                data.append(col)

        return data

    def _text_getter(self, obj):
        """Return the text of an individual DOM node.

        Parameters
        ----------
        obj : node-like
            A DOM node.

        Returns
        -------
        text : str or unicode
            The text from an individual DOM node.
        """
        raise NotImplementedError

    def _parse_columns(self, obj):
        """Return the column elements from a row element.

        Parameters
        ----------
        obj : node-like

        Returns
        -------
        columns : list of node-like
            These are the elements of each row, i.e., the columns.
        """
        raise NotImplementedError

    def _parse_tables(self, doc, match, attrs):
        """Return all tables from the parsed DOM.

        Parameters
        ----------
        doc : tree-like
            The DOM from which to parse the table element.

        match : str or regular expression
            The text to search for in the DOM tree.

        attrs : dict
            A dictionary of table attributes that can be used to disambiguate
            mutliple tables on a page.

        Raises
        ------
        AssertionError
            * If `match` does not match any text in the document.

        Returns
        -------
        tables : list of node-like
            A list of <table> elements to be parsed into raw data.
        """
        raise NotImplementedError

    def _parse_rows(self, table):
        """Return the list of row elements from the parsed table element.

        Parameters
        ----------
        table : node-like
            A table element that contains row elements.

        Returns
        -------
        rows : list of node-like
            A list row elements of a table, usually <tr> or <th> elements.
        """
        raise NotImplementedError

    def _build_doc(self):
        """Return a tree-like object that can be used to iterate over the DOM.

        Returns
        -------
        obj : tree-like
        """
        raise NotImplementedError


class _BeautifulSoupFrameParser(_HtmlFrameParser):
    """HTML to DataFrame parser that uses BeautifulSoup under the hood.

    See Also
    --------
    pandas.io.html._HtmlFrameParser
    pandas.io.html._LxmlFrameParser

    Notes
    -----
    Documentation strings for this class are in the base class
    :class:`pandas.io.html._HtmlFrameParser`.
    """
    def __init__(self, *args, **kwargs):
        super(_BeautifulSoupFrameParser, self).__init__(*args, **kwargs)

    def _text_getter(self, obj):
        return obj.text

    def _parse_columns(self, row):
        return row.find_all(('td', 'th'))

    def _parse_rows(self, table):
        return table.find_all(('tr', 'thead', 'tfoot'))

    def _parse_tables(self, doc, match, attrs):
        tables = doc.find_all('table', attrs=attrs)
        assert tables, 'No tables found'

        tables = [table for table in tables
                  if table.find(text=match) is not None]
        assert tables, "No tables found matching '{0}'".format(match.pattern)
        return tables

    def _build_doc(self):
        if _is_url(self.io):
            try:
                with contextlib.closing(urllib2.urlopen(self.io)) as url:
                    raw_text = url.read()
            except urllib2.URLError:
                raise ValueError('Invalid URL: "{0}"'.format(self.io))
        elif hasattr(self.io, 'read'):
            raw_text = self.io.read()
        elif os.path.isfile(self.io):
            with open(self.io) as f:
                raw_text = f.read()
        elif isinstance(self.io, basestring):
            raw_text = self.io
        else:
            raise ValueError("Cannot read object of"
                             " type '{0}'".format(type(self.io)))
        assert raw_text, 'No text parsed from document'

        from bs4 import BeautifulSoup, SoupStrainer
        strainer = SoupStrainer('table')
        return BeautifulSoup(raw_text, parse_only=strainer)


def _build_node_xpath_expr(attrs):
    """Build an xpath expression to simulate bs4's ability to pass in kwargs to
    search for attributes when using the lxml parser.

    Parameters
    ----------
    attrs : dict
        A dict of HTML attributes. These are NOT checked for validity.

    Returns
    -------
    expr : unicode
        An XPath expression that checks for the given HTML attributes.
    """
    # give class attribute as class_ because class is a python keyword
    if 'class_' in attrs:
        attrs['class'] = attrs.pop('class_')

    s = (u"@{k}='{v}'".format(k=k, v=v) for k, v in attrs.iteritems())
    return u'[{0}]'.format(' and '.join(s))


_re_namespace = {'re': 'http://exslt.org/regular-expressions'}


class _LxmlFrameParser(_HtmlFrameParser):
    """HTML to DataFrame parser that uses lxml under the hood.

    Warning
    -------
    This parser can only handle HTTP, FTP, and FILE urls.

    See Also
    --------
    _HtmlFrameParser
    _BeautifulSoupFrameParser

    Notes
    -----
    Documentation strings for this class are in the base class
    :class:`_HtmlFrameParser`.
    """
    def __init__(self, *args, **kwargs):
        super(_LxmlFrameParser, self).__init__(*args, **kwargs)

    def _text_getter(self, obj):
        return obj.text_content()

    def _parse_columns(self, row):
        return row.xpath('.//td|.//th')

    def _parse_rows(self, table):
        return table.xpath('(.//tr|.//thead|.//tfoot)[normalize-space()]')

    def _parse_tables(self, doc, match, kwargs):
        pattern = match.pattern

        # check all descendants for the given pattern
        check_all_expr = u'//*'
        if pattern:
            check_all_expr += u"[re:test(text(), '{0}')]".format(pattern)

        # go up the tree until we find a table
        check_table_expr = '/ancestor::table'
        xpath_expr = check_all_expr + check_table_expr

        # if any table attributes were given build an xpath expression to
        # search for them
        if kwargs:
            xpath_expr += _build_node_xpath_expr(kwargs)
        tables = doc.xpath(xpath_expr, namespaces=_re_namespace)
        assert tables, "No tables found matching regex '{0}'".format(pattern)
        return tables

    def _build_doc(self):
        """
        Raises
        ------
        IOError
            * If a valid URL is detected, but for some reason cannot be parsed.
              This is probably due to a faulty or non-existent internet
              connection.
        ValueError
            * If a URL that lxml cannot parse is passed.

        See Also
        --------
        pandas.io.html._HtmlFrameParser._build_doc
        """
        from lxml.html import parse, fromstring

        try:
            # try to parse the input in the simplest way
            return parse(self.io)
        except (UnicodeDecodeError, IOError):
            # something went wrong, check for not-a-url because it's probably a
            # huge string blob
            if not _is_url(self.io):
                return fromstring(self.io)
            elif urlparse.urlparse(self.io).scheme not in ('http', 'ftp',
                                                           'file'):
                raise ValueError('"{0}" does not have a valid URL'
                                 ' protocol'.format(self.io))
            else:
                raise IOError('"{0}" is a valid URL, so you probably are not'
                              ' properly connected to the'
                              ' internet'.format(self.io))


def _data_to_frame(data, header, index_col, infer_types, skiprows):
    """Parse a BeautifulSoup table into a DataFrame.

    Parameters
    ----------
    data : list of lists of str or unicode
        The raw data to be placed into a DataFrame. This is a list of lists of
        strings or unicode. If it helps, it can be thought of as a matrix of
        strings instead.

    header : int or None
        An integer indicating the row to use for the column header or None
        indicating no header will be used.

    index_col : int or None
        An integer indicating the column to use for the index or None
        indicating no column will be used.

    infer_types : bool
        Whether to convert numbers and dates.

    skiprows : collections.Container or int or slice
        Iterable used to skip rows.

    Returns
    -------
    df : DataFrame
        A DataFrame containing the data from `data`

    Raises
    ------
    ValueError
        * If `skiprows` is not found in the rows of the parsed DataFrame.

    Raises
    ------
    ValueError
        * If `skiprows` is not found in the rows of the parsed DataFrame.

    See Also
    --------
    read_html

    Notes
    -----
    The `data` parameter is guaranteed not to be a list of empty lists.
    """
    df = DataFrame(data)

    if skiprows is not None:
        it = _get_skiprows_iter(skiprows)

        try:
            df = df.drop(it)
        except ValueError:
            raise ValueError('Labels {0} not found when trying to skip'
                             ' rows'.format(it))

    if header is not None:
        header_rows = df.iloc[header]

        if header_rows.ndim == 2:
            names = header_rows.index
            df.columns = MultiIndex.from_arrays(header_rows.values,
                                                names=names)
        else:
            df.columns = header_rows

        df = df.drop(df.index[header])

    # convert to numbers/dates where possible
    # must be sequential since dates trump numbers if both args are given
    if infer_types:
        df = df.convert_objects(convert_numeric=True)
        df = df.convert_objects(convert_dates='coerce')

    if index_col is not None:
        cols = df.columns[index_col]

        try:
            cols = cols.tolist()
        except AttributeError:
            pass

        # drop by default
        df.set_index(cols, inplace=True)

    return df


_possible_parsers = {'lxml': _LxmlFrameParser,
                     'bs4': _BeautifulSoupFrameParser}


def read_html(io, match='.+', flavor='bs4', header=None, index_col=None,
              skiprows=None, infer_types=True, attrs=None):
    r"""Read an HTML table into a DataFrame.

    Parameters
    ----------
    io : str or file-like
        A string or file like object that can be either a url, a file-like
        object, or a raw string containing HTML.  Note that lxml only accepts
        the http, ftp and file url protocols.

    match : str or regex, optional
        The set of tables containing text matching this regex or string will be
        returned. Unless the HTML is extremely simple you will probably need to
        pass a non-empty string here. Defaults to '.+' (match any non-empty
        string). The default value will return all tables contained on a page.
        This value is converted to a regular expression so that there is
        consistent behavior between Beautiful Soup and lxml.

    flavor : str, {'lxml', 'bs4'}
        The parsing engine to use under the hood. lxml is faster and bs4
        (Beautiful Soup 4) is better at parsing nested tags, which are not
        uncommon when parsing tables. Defaults to 'bs4'.

    header : int or array-like or None, optional
        The row (or rows for a MultiIndex) to use to make the columns headers.
        Note that this row will be removed from the data. Defaults to None.

    index_col : int or array-like or None, optional
        The column to use to make the index. Note that this column will be
        removed from the data. Defaults to None.

    skiprows : int or collections.Container or slice or None, optional
        If an integer is given then skip this many rows after parsing the
        column header. If a sequence of integers is given skip those specific
        rows (0-based). Defaults to None, i.e., no rows are skipped. Note that

        .. code-block:: python

           skiprows == 0

        yields the same result as

        .. code-block:: python

           skiprows is None

        If `skiprows` is a positive integer, say :math:`n`, then
        it is treated as "skip :math:`n` rows", *not* as "skip the
        :math:`n^\textrm{th}` row".

    infer_types : bool, optional
        Whether to convert numeric types and date-appearing strings to numbers
        and dates, respectively. Defaults to True.

    attrs : dict or None, optional
        This is a dictionary of attributes that you can pass to use to identify
        the table in the HTML. These are not checked for validity before being
        passed to lxml or Beautiful Soup. However, these attributes must be
        valid HTML table attributes to work correctly. Defaults to None. For
        example,

        .. code-block:: python

           attrs = {'id': 'table'}

        is a valid attribute dictionary because the 'id' HTML tag attribute is
        a valid HTML attribute for *any* HTML tag as per `this document
        <http://www.w3.org/TR/html-markup/global-attributes.html>`__.

        .. code-block:: python

           attrs = {'asdf': 'table'}

        is *not* a valid attribute dictionary because 'asdf' is not a valid
        HTML attribute even if it is a valid XML attribute.  Valid HTML 4.01
        table attributes can be found `here
        <http://www.w3.org/TR/REC-html40/struct/tables.html#h-11.2>`__. A
        working draft of the HTML 5 spec can be found `here
        <http://www.w3.org/TR/html-markup/table.html>`__. It contains the
        latest information on table attributes for the modern web.

    Returns
    -------
    dfs : list of DataFrames
        A list of DataFrames, each of which is the parsed data from each of the
        tables on the page.

    Notes
    -----
    There's as little cleaning of the data as possible due to the heterogeneity
    and general disorder of HTML on the web.

    Expect some cleanup after you call this function. For example,
    you might need to pass `infer_types=False` and perform manual conversion if
    the column names are converted to NaN when you pass the `header=0`
    argument. We try to assume as little as possible about the structure of the
    table and push the idiosyncrasies of the HTML contained in the table to
    you, the user.

    This function only searches for <table> elements and only for <tr> and <th>
    rows and <td> elements within those rows. This could be extended by
    subclassing one of the parser classes contained in :mod:`pandas.io.html`.

    Similar to :func:`read_csv` the `header` argument is applied **after**
    `skiprows` is applied.

    This function will *always* return a list of :class:`DataFrame` *or*
    it will fail, e.g., it will *not* return an empty list.

    Examples
    --------
    Parse a table from a list of failed banks from the FDIC:

    >>> from pandas import read_html, DataFrame
    >>> url = 'http://www.fdic.gov/bank/individual/failed/banklist.html'
    >>> dfs = read_html(url, match='Florida', attrs={'id': 'table'})
    >>> assert dfs  # will not be empty if the call to read_html doesn't fail
    >>> assert isinstance(dfs, list)  # read_html returns a list of DataFrames
    >>> assert all(map(lambda x: isinstance(x, DataFrame), dfs))

    Parse some spam infomation from the USDA:

    >>> url = ('http://ndb.nal.usda.gov/ndb/foods/show/1732?fg=&man=&'
    ...        'lfacet=&format=&count=&max=25&offset=&sort=&qlookup=spam')
    >>> dfs = read_html(url, match='Water', header=0)
    >>> assert dfs
    >>> assert isinstance(dfs, list)
    >>> assert all(map(lambda x: isinstance(x, DataFrame), dfs))

    You can pass nothing to the `match` argument:

    >>> url = 'http://www.fdic.gov/bank/individual/failed/banklist.html'
    >>> dfs = read_html(url)
    >>> print(len(dfs))  # this will most likely be greater than 1

    Try a different parser:

    >>> url = 'http://www.fdic.gov/bank/individual/failed/banklist.html'
    >>> dfs = read_html(url, 'Florida', flavor='lxml', attrs={'id': 'table'})
    >>> assert dfs
    >>> assert isinstance(dfs, list)
    >>> assert all(map(lambda x: isinstance(x, DataFrame), dfs))
    """
    # annoying type check here because we don't want to spend time parsing HTML
    # only to end up failing because of an invalid value of skiprows
    if isinstance(skiprows, numbers.Integral):
        assert skiprows >= 0, ('cannot skip rows starting from the end of the '
                               'data (you passed a negative value)')

    valid_backends = _possible_parsers.keys()
    assert flavor in valid_backends, ("'{0}' is not a valid backend, the valid"
                                      " backends are "
                                      "{1}".format(flavor, valid_backends))
    parser = _possible_parsers[flavor]

    # bonus: re.compile is idempotent under function iteration so you can pass
    # a compiled regex to it and it will return itself
    p = parser(io, re.compile(match), attrs)
    return [_data_to_frame(data, header, index_col, infer_types, skiprows)
            for data in p.parse_raw_data()]
