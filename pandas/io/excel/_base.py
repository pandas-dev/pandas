import abc
import datetime
from io import BytesIO
import os
from textwrap import fill

from pandas._config import config

from pandas._libs.parsers import STR_NA_VALUES
from pandas.errors import EmptyDataError
from pandas.util._decorators import Appender

from pandas.core.dtypes.common import is_bool, is_float, is_integer, is_list_like

from pandas.core.frame import DataFrame

from pandas.io.common import (
    get_filepath_or_buffer,
    is_url,
    stringify_path,
    urlopen,
    validate_header_arg,
)
from pandas.io.excel._util import (
    _fill_mi_header,
    _get_default_writer,
    _maybe_convert_usecols,
    _pop_header_name,
    get_writer,
)
from pandas.io.parsers import TextParser

_read_excel_doc = (
    """
Read an Excel file into a pandas DataFrame.

Supports `xls`, `xlsx`, `xlsm`, `xlsb`, and `odf` file extensions
read from a local filesystem or URL. Supports an option to read
a single sheet or a list of sheets.

Parameters
----------
io : str, bytes, ExcelFile, xlrd.Book, path object, or file-like object
    Any valid string path is acceptable. The string could be a URL. Valid
    URL schemes include http, ftp, s3, and file. For file URLs, a host is
    expected. A local file could be: ``file://localhost/path/to/table.xlsx``.

    If you want to pass in a path object, pandas accepts any ``os.PathLike``.

    By file-like object, we refer to objects with a ``read()`` method,
    such as a file handler (e.g. via builtin ``open`` function)
    or ``StringIO``.
sheet_name : str, int, list, or None, default 0
    Strings are used for sheet names. Integers are used in zero-indexed
    sheet positions. Lists of strings/integers are used to request
    multiple sheets. Specify None to get all sheets.

    Available cases:

    * Defaults to ``0``: 1st sheet as a `DataFrame`
    * ``1``: 2nd sheet as a `DataFrame`
    * ``"Sheet1"``: Load sheet with name "Sheet1"
    * ``[0, 1, "Sheet5"]``: Load first, second and sheet named "Sheet5"
      as a dict of `DataFrame`
    * None: All sheets.

header : int, list of int, default 0
    Row (0-indexed) to use for the column labels of the parsed
    DataFrame. If a list of integers is passed those row positions will
    be combined into a ``MultiIndex``. Use None if there is no header.
names : array-like, default None
    List of column names to use. If file contains no header row,
    then you should explicitly pass header=None.
index_col : int, list of int, default None
    Column (0-indexed) to use as the row labels of the DataFrame.
    Pass None if there is no such column.  If a list is passed,
    those columns will be combined into a ``MultiIndex``.  If a
    subset of data is selected with ``usecols``, index_col
    is based on the subset.
usecols : int, str, list-like, or callable default None
    * If None, then parse all columns.
    * If str, then indicates comma separated list of Excel column letters
      and column ranges (e.g. "A:E" or "A,C,E:F"). Ranges are inclusive of
      both sides.
    * If list of int, then indicates list of column numbers to be parsed.
    * If list of string, then indicates list of column names to be parsed.

      .. versionadded:: 0.24.0

    * If callable, then evaluate each column name against it and parse the
      column if the callable returns ``True``.

    Returns a subset of the columns according to behavior above.

      .. versionadded:: 0.24.0

squeeze : bool, default False
    If the parsed data only contains one column then return a Series.
dtype : Type name or dict of column -> type, default None
    Data type for data or columns. E.g. {'a': np.float64, 'b': np.int32}
    Use `object` to preserve data as stored in Excel and not interpret dtype.
    If converters are specified, they will be applied INSTEAD
    of dtype conversion.
engine : str, default None
    If io is not a buffer or path, this must be set to identify io.
    Acceptable values are None, "xlrd", "openpyxl" or "odf".
converters : dict, default None
    Dict of functions for converting values in certain columns. Keys can
    either be integers or column labels, values are functions that take one
    input argument, the Excel cell content, and return the transformed
    content.
true_values : list, default None
    Values to consider as True.
false_values : list, default None
    Values to consider as False.
skiprows : list-like
    Rows to skip at the beginning (0-indexed).
nrows : int, default None
    Number of rows to parse.

    .. versionadded:: 0.23.0

na_values : scalar, str, list-like, or dict, default None
    Additional strings to recognize as NA/NaN. If dict passed, specific
    per-column NA values. By default the following values are interpreted
    as NaN: '"""
    + fill("', '".join(sorted(STR_NA_VALUES)), 70, subsequent_indent="    ")
    + """'.
keep_default_na : bool, default True
    Whether or not to include the default NaN values when parsing the data.
    Depending on whether `na_values` is passed in, the behavior is as follows:

    * If `keep_default_na` is True, and `na_values` are specified, `na_values`
      is appended to the default NaN values used for parsing.
    * If `keep_default_na` is True, and `na_values` are not specified, only
      the default NaN values are used for parsing.
    * If `keep_default_na` is False, and `na_values` are specified, only
      the NaN values specified `na_values` are used for parsing.
    * If `keep_default_na` is False, and `na_values` are not specified, no
      strings will be parsed as NaN.

    Note that if `na_filter` is passed in as False, the `keep_default_na` and
    `na_values` parameters will be ignored.
na_filter : bool, default True
    Detect missing value markers (empty strings and the value of na_values). In
    data without any NAs, passing na_filter=False can improve the performance
    of reading a large file.
verbose : bool, default False
    Indicate number of NA values placed in non-numeric columns.
parse_dates : bool, list-like, or dict, default False
    The behavior is as follows:

    * bool. If True -> try parsing the index.
    * list of int or names. e.g. If [1, 2, 3] -> try parsing columns 1, 2, 3
      each as a separate date column.
    * list of lists. e.g.  If [[1, 3]] -> combine columns 1 and 3 and parse as
      a single date column.
    * dict, e.g. {'foo' : [1, 3]} -> parse columns 1, 3 as date and call
      result 'foo'

    If a column or index contains an unparseable date, the entire column or
    index will be returned unaltered as an object data type. If you don`t want to
    parse some cells as date just change their type in Excel to "Text".
    For non-standard datetime parsing, use ``pd.to_datetime`` after ``pd.read_excel``.

    Note: A fast-path exists for iso8601-formatted dates.
date_parser : function, optional
    Function to use for converting a sequence of string columns to an array of
    datetime instances. The default uses ``dateutil.parser.parser`` to do the
    conversion. Pandas will try to call `date_parser` in three different ways,
    advancing to the next if an exception occurs: 1) Pass one or more arrays
    (as defined by `parse_dates`) as arguments; 2) concatenate (row-wise) the
    string values from the columns defined by `parse_dates` into a single array
    and pass that; and 3) call `date_parser` once for each row using one or
    more strings (corresponding to the columns defined by `parse_dates`) as
    arguments.
thousands : str, default None
    Thousands separator for parsing string columns to numeric.  Note that
    this parameter is only necessary for columns stored as TEXT in Excel,
    any numeric columns will automatically be parsed, regardless of display
    format.
comment : str, default None
    Comments out remainder of line. Pass a character or characters to this
    argument to indicate comments in the input file. Any data between the
    comment string and the end of the current line is ignored.
skipfooter : int, default 0
    Rows at the end to skip (0-indexed).
convert_float : bool, default True
    Convert integral floats to int (i.e., 1.0 --> 1). If False, all numeric
    data will be read in as floats: Excel stores all numbers as floats
    internally.
mangle_dupe_cols : bool, default True
    Duplicate columns will be specified as 'X', 'X.1', ...'X.N', rather than
    'X'...'X'. Passing in False will cause data to be overwritten if there
    are duplicate names in the columns.
**kwds : optional
        Optional keyword arguments can be passed to ``TextFileReader``.

Returns
-------
DataFrame or dict of DataFrames
    DataFrame from the passed in Excel file. See notes in sheet_name
    argument for more information on when a dict of DataFrames is returned.

See Also
--------
DataFrame.to_excel : Write DataFrame to an Excel file.
DataFrame.to_csv : Write DataFrame to a comma-separated values (csv) file.
read_csv : Read a comma-separated values (csv) file into DataFrame.
read_fwf : Read a table of fixed-width formatted lines into DataFrame.

Examples
--------
The file can be read using the file name as string or an open file object:

>>> pd.read_excel('tmp.xlsx', index_col=0)  # doctest: +SKIP
       Name  Value
0   string1      1
1   string2      2
2  #Comment      3

>>> pd.read_excel(open('tmp.xlsx', 'rb'),
...               sheet_name='Sheet3')  # doctest: +SKIP
   Unnamed: 0      Name  Value
0           0   string1      1
1           1   string2      2
2           2  #Comment      3

Index and header can be specified via the `index_col` and `header` arguments

>>> pd.read_excel('tmp.xlsx', index_col=None, header=None)  # doctest: +SKIP
     0         1      2
0  NaN      Name  Value
1  0.0   string1      1
2  1.0   string2      2
3  2.0  #Comment      3

Column types are inferred but can be explicitly specified

>>> pd.read_excel('tmp.xlsx', index_col=0,
...               dtype={'Name': str, 'Value': float})  # doctest: +SKIP
       Name  Value
0   string1    1.0
1   string2    2.0
2  #Comment    3.0

True, False, and NA values, and thousands separators have defaults,
but can be explicitly specified, too. Supply the values you would like
as strings or lists of strings!

>>> pd.read_excel('tmp.xlsx', index_col=0,
...               na_values=['string1', 'string2'])  # doctest: +SKIP
       Name  Value
0       NaN      1
1       NaN      2
2  #Comment      3

Comment lines in the excel input file can be skipped using the `comment` kwarg

>>> pd.read_excel('tmp.xlsx', index_col=0, comment='#')  # doctest: +SKIP
      Name  Value
0  string1    1.0
1  string2    2.0
2     None    NaN
"""
)


@Appender(_read_excel_doc)
def read_excel(
    io,
    sheet_name=0,
    header=0,
    names=None,
    index_col=None,
    usecols=None,
    squeeze=False,
    dtype=None,
    engine=None,
    converters=None,
    true_values=None,
    false_values=None,
    skiprows=None,
    nrows=None,
    na_values=None,
    keep_default_na=True,
    verbose=False,
    parse_dates=False,
    date_parser=None,
    thousands=None,
    comment=None,
    skipfooter=0,
    convert_float=True,
    mangle_dupe_cols=True,
    **kwds,
):

    for arg in ("sheet", "sheetname", "parse_cols"):
        if arg in kwds:
            raise TypeError(f"read_excel() got an unexpected keyword argument `{arg}`")

    if not isinstance(io, ExcelFile):
        io = ExcelFile(io, engine=engine)
    elif engine and engine != io.engine:
        raise ValueError(
            "Engine should not be specified when passing "
            "an ExcelFile - ExcelFile already has the engine set"
        )

    return io.parse(
        sheet_name=sheet_name,
        header=header,
        names=names,
        index_col=index_col,
        usecols=usecols,
        squeeze=squeeze,
        dtype=dtype,
        converters=converters,
        true_values=true_values,
        false_values=false_values,
        skiprows=skiprows,
        nrows=nrows,
        na_values=na_values,
        keep_default_na=keep_default_na,
        verbose=verbose,
        parse_dates=parse_dates,
        date_parser=date_parser,
        thousands=thousands,
        comment=comment,
        skipfooter=skipfooter,
        convert_float=convert_float,
        mangle_dupe_cols=mangle_dupe_cols,
        **kwds,
    )


class _BaseExcelReader(metaclass=abc.ABCMeta):
    def __init__(self, filepath_or_buffer):
        # If filepath_or_buffer is a url, load the data into a BytesIO
        if is_url(filepath_or_buffer):
            filepath_or_buffer = BytesIO(urlopen(filepath_or_buffer).read())
        elif not isinstance(filepath_or_buffer, (ExcelFile, self._workbook_class)):
            filepath_or_buffer, _, _, _ = get_filepath_or_buffer(filepath_or_buffer)

        if isinstance(filepath_or_buffer, self._workbook_class):
            self.book = filepath_or_buffer
        elif hasattr(filepath_or_buffer, "read"):
            # N.B. xlrd.Book has a read attribute too
            filepath_or_buffer.seek(0)
            self.book = self.load_workbook(filepath_or_buffer)
        elif isinstance(filepath_or_buffer, str):
            self.book = self.load_workbook(filepath_or_buffer)
        elif isinstance(filepath_or_buffer, bytes):
            self.book = self.load_workbook(BytesIO(filepath_or_buffer))
        else:
            raise ValueError(
                "Must explicitly set engine if not passing in buffer or path for io."
            )

    @property
    @abc.abstractmethod
    def _workbook_class(self):
        pass

    @abc.abstractmethod
    def load_workbook(self, filepath_or_buffer):
        pass

    @property
    @abc.abstractmethod
    def sheet_names(self):
        pass

    @abc.abstractmethod
    def get_sheet_by_name(self, name):
        pass

    @abc.abstractmethod
    def get_sheet_by_index(self, index):
        pass

    @abc.abstractmethod
    def get_sheet_data(self, sheet, convert_float):
        pass

    def parse(
        self,
        sheet_name=0,
        header=0,
        names=None,
        index_col=None,
        usecols=None,
        squeeze=False,
        dtype=None,
        true_values=None,
        false_values=None,
        skiprows=None,
        nrows=None,
        na_values=None,
        verbose=False,
        parse_dates=False,
        date_parser=None,
        thousands=None,
        comment=None,
        skipfooter=0,
        convert_float=True,
        mangle_dupe_cols=True,
        **kwds,
    ):

        validate_header_arg(header)

        ret_dict = False

        # Keep sheetname to maintain backwards compatibility.
        if isinstance(sheet_name, list):
            sheets = sheet_name
            ret_dict = True
        elif sheet_name is None:
            sheets = self.sheet_names
            ret_dict = True
        else:
            sheets = [sheet_name]

        # handle same-type duplicates.
        sheets = list(dict.fromkeys(sheets).keys())

        output = {}

        for asheetname in sheets:
            if verbose:
                print(f"Reading sheet {asheetname}")

            if isinstance(asheetname, str):
                sheet = self.get_sheet_by_name(asheetname)
            else:  # assume an integer if not a string
                sheet = self.get_sheet_by_index(asheetname)

            data = self.get_sheet_data(sheet, convert_float)
            usecols = _maybe_convert_usecols(usecols)

            if not data:
                output[asheetname] = DataFrame()
                continue

            if is_list_like(header) and len(header) == 1:
                header = header[0]

            # forward fill and pull out names for MultiIndex column
            header_names = None
            if header is not None and is_list_like(header):
                header_names = []
                control_row = [True] * len(data[0])

                for row in header:
                    if is_integer(skiprows):
                        row += skiprows

                    data[row], control_row = _fill_mi_header(data[row], control_row)

                    if index_col is not None:
                        header_name, _ = _pop_header_name(data[row], index_col)
                        header_names.append(header_name)

            if is_list_like(index_col):
                # Forward fill values for MultiIndex index.
                if not is_list_like(header):
                    offset = 1 + header
                else:
                    offset = 1 + max(header)

                # Check if we have an empty dataset
                # before trying to collect data.
                if offset < len(data):
                    for col in index_col:
                        last = data[offset][col]

                        for row in range(offset + 1, len(data)):
                            if data[row][col] == "" or data[row][col] is None:
                                data[row][col] = last
                            else:
                                last = data[row][col]

            has_index_names = is_list_like(header) and len(header) > 1

            # GH 12292 : error when read one empty column from excel file
            try:
                parser = TextParser(
                    data,
                    names=names,
                    header=header,
                    index_col=index_col,
                    has_index_names=has_index_names,
                    squeeze=squeeze,
                    dtype=dtype,
                    true_values=true_values,
                    false_values=false_values,
                    skiprows=skiprows,
                    nrows=nrows,
                    na_values=na_values,
                    parse_dates=parse_dates,
                    date_parser=date_parser,
                    thousands=thousands,
                    comment=comment,
                    skipfooter=skipfooter,
                    usecols=usecols,
                    mangle_dupe_cols=mangle_dupe_cols,
                    **kwds,
                )

                output[asheetname] = parser.read(nrows=nrows)

                if not squeeze or isinstance(output[asheetname], DataFrame):
                    if header_names:
                        output[asheetname].columns = output[
                            asheetname
                        ].columns.set_names(header_names)

            except EmptyDataError:
                # No Data, return an empty DataFrame
                output[asheetname] = DataFrame()

        if ret_dict:
            return output
        else:
            return output[asheetname]


class ExcelWriter(metaclass=abc.ABCMeta):
    """
    Class for writing DataFrame objects into excel sheets.

    Default is to use xlwt for xls, openpyxl for xlsx.
    See DataFrame.to_excel for typical usage.

    Parameters
    ----------
    path : str
        Path to xls or xlsx file.
    engine : str (optional)
        Engine to use for writing. If None, defaults to
        ``io.excel.<extension>.writer``.  NOTE: can only be passed as a keyword
        argument.
    date_format : str, default None
        Format string for dates written into Excel files (e.g. 'YYYY-MM-DD').
    datetime_format : str, default None
        Format string for datetime objects written into Excel files.
        (e.g. 'YYYY-MM-DD HH:MM:SS').
    mode : {'w', 'a'}, default 'w'
        File mode to use (write or append).

        .. versionadded:: 0.24.0

    Attributes
    ----------
    None

    Methods
    -------
    None

    Notes
    -----
    None of the methods and properties are considered public.

    For compatibility with CSV writers, ExcelWriter serializes lists
    and dicts to strings before writing.

    Examples
    --------
    Default usage:

    >>> with ExcelWriter('path_to_file.xlsx') as writer:
    ...     df.to_excel(writer)

    To write to separate sheets in a single file:

    >>> with ExcelWriter('path_to_file.xlsx') as writer:
    ...     df1.to_excel(writer, sheet_name='Sheet1')
    ...     df2.to_excel(writer, sheet_name='Sheet2')

    You can set the date format or datetime format:

    >>> with ExcelWriter('path_to_file.xlsx',
                          date_format='YYYY-MM-DD',
                          datetime_format='YYYY-MM-DD HH:MM:SS') as writer:
    ...     df.to_excel(writer)

    You can also append to an existing Excel file:

    >>> with ExcelWriter('path_to_file.xlsx', mode='a') as writer:
    ...     df.to_excel(writer, sheet_name='Sheet3')
    """

    # Defining an ExcelWriter implementation (see abstract methods for more...)

    # - Mandatory
    #   - ``write_cells(self, cells, sheet_name=None, startrow=0, startcol=0)``
    #     --> called to write additional DataFrames to disk
    #   - ``supported_extensions`` (tuple of supported extensions), used to
    #      check that engine supports the given extension.
    #   - ``engine`` - string that gives the engine name. Necessary to
    #     instantiate class directly and bypass ``ExcelWriterMeta`` engine
    #     lookup.
    #   - ``save(self)`` --> called to save file to disk
    # - Mostly mandatory (i.e. should at least exist)
    #   - book, cur_sheet, path

    # - Optional:
    #   - ``__init__(self, path, engine=None, **kwargs)`` --> always called
    #     with path as first argument.

    # You also need to register the class with ``register_writer()``.
    # Technically, ExcelWriter implementations don't need to subclass
    # ExcelWriter.
    def __new__(cls, path, engine=None, **kwargs):
        # only switch class if generic(ExcelWriter)

        if cls is ExcelWriter:
            if engine is None or (isinstance(engine, str) and engine == "auto"):
                if isinstance(path, str):
                    ext = os.path.splitext(path)[-1][1:]
                else:
                    ext = "xlsx"

                try:
                    engine = config.get_option(f"io.excel.{ext}.writer")
                    if engine == "auto":
                        engine = _get_default_writer(ext)
                except KeyError as err:
                    raise ValueError(f"No engine for filetype: '{ext}'") from err
            cls = get_writer(engine)

        return object.__new__(cls)

    # declare external properties you can count on
    book = None
    curr_sheet = None
    path = None

    @property
    @abc.abstractmethod
    def supported_extensions(self):
        """Extensions that writer engine supports."""
        pass

    @property
    @abc.abstractmethod
    def engine(self):
        """Name of engine."""
        pass

    @abc.abstractmethod
    def write_cells(
        self, cells, sheet_name=None, startrow=0, startcol=0, freeze_panes=None
    ):
        """
        Write given formatted cells into Excel an excel sheet

        Parameters
        ----------
        cells : generator
            cell of formatted data to save to Excel sheet
        sheet_name : str, default None
            Name of Excel sheet, if None, then use self.cur_sheet
        startrow : upper left cell row to dump data frame
        startcol : upper left cell column to dump data frame
        freeze_panes: int tuple of length 2
            contains the bottom-most row and right-most column to freeze
        """
        pass

    @abc.abstractmethod
    def save(self):
        """
        Save workbook to disk.
        """
        pass

    def __init__(
        self,
        path,
        engine=None,
        date_format=None,
        datetime_format=None,
        mode="w",
        **engine_kwargs,
    ):
        # validate that this engine can handle the extension
        if isinstance(path, str):
            ext = os.path.splitext(path)[-1]
        else:
            ext = "xls" if engine == "xlwt" else "xlsx"

        self.check_extension(ext)

        self.path = path
        self.sheets = {}
        self.cur_sheet = None

        if date_format is None:
            self.date_format = "YYYY-MM-DD"
        else:
            self.date_format = date_format
        if datetime_format is None:
            self.datetime_format = "YYYY-MM-DD HH:MM:SS"
        else:
            self.datetime_format = datetime_format

        self.mode = mode

    def __fspath__(self):
        return stringify_path(self.path)

    def _get_sheet_name(self, sheet_name):
        if sheet_name is None:
            sheet_name = self.cur_sheet
        if sheet_name is None:  # pragma: no cover
            raise ValueError("Must pass explicit sheet_name or set cur_sheet property")
        return sheet_name

    def _value_with_fmt(self, val):
        """
        Convert numpy types to Python types for the Excel writers.

        Parameters
        ----------
        val : object
            Value to be written into cells

        Returns
        -------
        Tuple with the first element being the converted value and the second
            being an optional format
        """
        fmt = None

        if is_integer(val):
            val = int(val)
        elif is_float(val):
            val = float(val)
        elif is_bool(val):
            val = bool(val)
        elif isinstance(val, datetime.datetime):
            fmt = self.datetime_format
        elif isinstance(val, datetime.date):
            fmt = self.date_format
        elif isinstance(val, datetime.timedelta):
            val = val.total_seconds() / float(86400)
            fmt = "0"
        else:
            val = str(val)

        return val, fmt

    @classmethod
    def check_extension(cls, ext):
        """
        checks that path's extension against the Writer's supported
        extensions.  If it isn't supported, raises UnsupportedFiletypeError.
        """
        if ext.startswith("."):
            ext = ext[1:]
        if not any(ext in extension for extension in cls.supported_extensions):
            raise ValueError(f"Invalid extension for engine '{cls.engine}': '{ext}'")
        else:
            return True

    # Allow use as a contextmanager
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        """synonym for save, to make it more file-like"""
        return self.save()


class ExcelFile:
    """
    Class for parsing tabular excel sheets into DataFrame objects.
    Uses xlrd. See read_excel for more documentation

    Parameters
    ----------
    io : str, path object (pathlib.Path or py._path.local.LocalPath),
        a file-like object, xlrd workbook or openpypl workbook.
        If a string or path object, expected to be a path to xls, xlsx or odf file.
    engine : str, default None
        If io is not a buffer or path, this must be set to identify io.
        Acceptable values are None, ``xlrd``, ``openpyxl``,  ``odf``, or ``pyxlsb``.
        Note that ``odf`` reads tables out of OpenDocument formatted files.
    """

    from pandas.io.excel._odfreader import _ODFReader
    from pandas.io.excel._openpyxl import _OpenpyxlReader
    from pandas.io.excel._xlrd import _XlrdReader
    from pandas.io.excel._pyxlsb import _PyxlsbReader

    _engines = {
        "xlrd": _XlrdReader,
        "openpyxl": _OpenpyxlReader,
        "odf": _ODFReader,
        "pyxlsb": _PyxlsbReader,
    }

    def __init__(self, io, engine=None):
        if engine is None:
            engine = "xlrd"
        if engine not in self._engines:
            raise ValueError(f"Unknown engine: {engine}")

        self.engine = engine
        # could be a str, ExcelFile, Book, etc.
        self.io = io
        # Always a string
        self._io = stringify_path(io)

        self._reader = self._engines[engine](self._io)

    def __fspath__(self):
        return self._io

    def parse(
        self,
        sheet_name=0,
        header=0,
        names=None,
        index_col=None,
        usecols=None,
        squeeze=False,
        converters=None,
        true_values=None,
        false_values=None,
        skiprows=None,
        nrows=None,
        na_values=None,
        parse_dates=False,
        date_parser=None,
        thousands=None,
        comment=None,
        skipfooter=0,
        convert_float=True,
        mangle_dupe_cols=True,
        **kwds,
    ):
        """
        Parse specified sheet(s) into a DataFrame.

        Equivalent to read_excel(ExcelFile, ...)  See the read_excel
        docstring for more info on accepted parameters.

        Returns
        -------
        DataFrame or dict of DataFrames
            DataFrame from the passed in Excel file.
        """
        if "chunksize" in kwds:
            raise NotImplementedError(
                "chunksize keyword of read_excel is not implemented"
            )

        return self._reader.parse(
            sheet_name=sheet_name,
            header=header,
            names=names,
            index_col=index_col,
            usecols=usecols,
            squeeze=squeeze,
            converters=converters,
            true_values=true_values,
            false_values=false_values,
            skiprows=skiprows,
            nrows=nrows,
            na_values=na_values,
            parse_dates=parse_dates,
            date_parser=date_parser,
            thousands=thousands,
            comment=comment,
            skipfooter=skipfooter,
            convert_float=convert_float,
            mangle_dupe_cols=mangle_dupe_cols,
            **kwds,
        )

    @property
    def book(self):
        return self._reader.book

    @property
    def sheet_names(self):
        return self._reader.sheet_names

    def close(self):
        """close io if necessary"""
        if self.engine == "openpyxl":
            # https://stackoverflow.com/questions/31416842/
            #  openpyxl-does-not-close-excel-workbook-in-read-only-mode
            wb = self.book
            wb._archive.close()

        if hasattr(self.io, "close"):
            self.io.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __del__(self):
        # Ensure we don't leak file descriptors, but put in try/except in case
        # attributes are already deleted
        try:
            self.close()
        except AttributeError:
            pass
