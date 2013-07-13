"""
Module contains tools for processing Stata files into DataFrames

The StataReader below was originally written by Joe Presbrey as part of PyDTA.
It has been extended and improved by Skipper Seabold from the Statsmodels project
who also developed the StataWriter and was finally added to pandas in an once again
improved version.

You can find more information on http://presbrey.mit.edu/PyDTA and
http://statsmodels.sourceforge.net/devel/
"""

from StringIO import StringIO
import numpy as np

import sys
import struct
from pandas.core.base import StringMixin
from pandas.core.frame import DataFrame
from pandas.core.series import Series
from pandas.core.categorical import Categorical
import datetime
from pandas.util import py3compat
from pandas import isnull
from pandas.io.parsers import _parser_params, Appender
from pandas.io.common import get_filepath_or_buffer


_read_stata_doc = """
Read Stata file into DataFrame

%s
""" % (_parser_params)


@Appender(_read_stata_doc)
def read_stata(filepath_or_buffer, convert_dates=True, convert_categoricals=True, encoding=None, index=None):
    reader = StataReader(filepath_or_buffer, encoding)

    return reader.data(convert_dates, convert_categoricals, index)

_date_formats = ["%tc", "%tC", "%td", "%tw", "%tm", "%tq", "%th", "%ty"]

def _stata_elapsed_date_to_datetime(date, fmt):
    """
    Convert from SIF to datetime. http://www.stata.com/help.cgi?datetime

    Parameters
    ----------
    date : int
        The Stata Internal Format date to convert to datetime according to fmt
    fmt : str
        The format to convert to. Can be, tc, td, tw, tm, tq, th, ty

    Examples
    --------
    >>> _stata_elapsed_date_to_datetime(52, "%tw")                                datetime.datetime(1961, 1, 1, 0, 0)

    Notes
    -----
    datetime/c - tc
        milliseconds since 01jan1960 00:00:00.000, assuming 86,400 s/day
    datetime/C - tC - NOT IMPLEMENTED
        milliseconds since 01jan1960 00:00:00.000, adjusted for leap seconds
    date - td
        days since 01jan1960 (01jan1960 = 0)
    weekly date - tw
        weeks since 1960w1
        This assumes 52 weeks in a year, then adds 7 * remainder of the weeks.
        The datetime value is the start of the week in terms of days in the
        year, not ISO calendar weeks.
    monthly date - tm
        months since 1960m1
    quarterly date - tq
        quarters since 1960q1
    half-yearly date - th
        half-years since 1960h1 yearly
    date - ty
        years since 0000

    If you don't have pandas with datetime support, then you can't do
    milliseconds accurately.
    """
    #NOTE: we could run into overflow / loss of precision situations here
    # casting to int, but I'm not sure what to do. datetime won't deal with
    # numpy types and numpy datetime isn't mature enough / we can't rely on
    # pandas version > 0.7.1
    #TODO: IIRC relative delta doesn't play well with np.datetime?
    if np.isnan(date):
        return np.datetime64('nat')

    date = int(date)
    stata_epoch = datetime.datetime(1960, 1, 1)
    if fmt in ["%tc", "tc"]:
        from dateutil.relativedelta import relativedelta
        return stata_epoch + relativedelta(microseconds=date * 1000)
    elif fmt in ["%tC", "tC"]:
        from warnings import warn
        warn("Encountered %tC format. Leaving in Stata Internal Format.")
        return date
    elif fmt in ["%td", "td"]:
        return stata_epoch + datetime.timedelta(int(date))
    elif fmt in ["%tw", "tw"]:  # does not count leap days - 7 days is a week
        year = datetime.datetime(stata_epoch.year + date // 52, 1, 1)
        day_delta = (date % 52) * 7
        return year + datetime.timedelta(int(day_delta))
    elif fmt in ["%tm", "tm"]:
        year = stata_epoch.year + date // 12
        month_delta = (date % 12) + 1
        return datetime.datetime(year, month_delta, 1)
    elif fmt in ["%tq", "tq"]:
        year = stata_epoch.year + date // 4
        month_delta = (date % 4) * 3 + 1
        return datetime.datetime(year, month_delta, 1)
    elif fmt in ["%th", "th"]:
        year = stata_epoch.year + date // 2
        month_delta = (date % 2) * 6 + 1
        return datetime.datetime(year, month_delta, 1)
    elif fmt in ["%ty", "ty"]:
        if date > 0:
            return datetime.datetime(date, 1, 1)
        else:  # don't do negative years bc can't mix dtypes in column
            raise ValueError("Year 0 and before not implemented")
    else:
        raise ValueError("Date fmt %s not understood" % fmt)


def _datetime_to_stata_elapsed(date, fmt):
    """
    Convert from datetime to SIF. http://www.stata.com/help.cgi?datetime

    Parameters
    ----------
    date : datetime.datetime
        The date to convert to the Stata Internal Format given by fmt
    fmt : str
        The format to convert to. Can be, tc, td, tw, tm, tq, th, ty
    """
    if not isinstance(date, datetime.datetime):
        raise ValueError("date should be datetime.datetime format")
    stata_epoch = datetime.datetime(1960, 1, 1)
    if fmt in ["%tc", "tc"]:
        delta = date - stata_epoch
        return (delta.days * 86400000 + delta.seconds*1000 +
                delta.microseconds/1000)
    elif fmt in ["%tC", "tC"]:
        from warnings import warn
        warn("Stata Internal Format tC not supported.")
        return date
    elif fmt in ["%td", "td"]:
        return (date - stata_epoch).days
    elif fmt in ["%tw", "tw"]:
        return (52*(date.year-stata_epoch.year) +
                (date - datetime.datetime(date.year, 1, 1)).days / 7)
    elif fmt in ["%tm", "tm"]:
        return (12 * (date.year - stata_epoch.year) + date.month - 1)
    elif fmt in ["%tq", "tq"]:
        return 4*(date.year-stata_epoch.year) + int((date.month - 1)/3)
    elif fmt in ["%th", "th"]:
        return 2 * (date.year - stata_epoch.year) + int(date.month > 6)
    elif fmt in ["%ty", "ty"]:
        return date.year
    else:
        raise ValueError("fmt %s not understood" % fmt)


class StataMissingValue(StringMixin):
    """
    An observation's missing value.

    Parameters
    -----------
    offset
    value

    Attributes
    ----------
    string
    value

    Notes
    -----
    More information: <http://www.stata.com/help.cgi?missing>
    """

    def __init__(self, offset, value):
        self._value = value
        if type(value) is int or type(value) is long:
            self._str = value - offset is 1 and \
                '.' or ('.' + chr(value - offset + 96))
        else:
            self._str = '.'
    string = property(lambda self: self._str, doc="The Stata representation of the missing value: '.', '.a'..'.z'")
    value = property(lambda self: self._value, doc='The binary representation of the missing value.')

    def __unicode__(self):
        return self.string

    def __repr__(self):
        # not perfect :-/
        return "%s(%s)" % (self.__class__, self)


class StataParser(object):
    def __init__(self, encoding):
        if(encoding is None):
            self._encoding = 'cp1252'
        else:
            self._encoding = encoding

        #type          code.
        #--------------------
        #str1        1 = 0x01
        #str2        2 = 0x02
        #...
        #str244    244 = 0xf4
        #byte      251 = 0xfb  (sic)
        #int       252 = 0xfc
        #long      253 = 0xfd
        #float     254 = 0xfe
        #double    255 = 0xff
        #--------------------
        #NOTE: the byte type seems to be reserved for categorical variables
        # with a label, but the underlying variable is -127 to 100
        # we're going to drop the label and cast to int
        self.DTYPE_MAP = \
            dict(
                zip(range(1, 245), ['a' + str(i) for i in range(1, 245)]) +
                [
                    (251, np.int16),
                    (252, np.int32),
                    (253, np.int64),
                    (254, np.float32),
                    (255, np.float64)
                ]
            )
        self.TYPE_MAP = range(251) + list('bhlfd')
        #NOTE: technically, some of these are wrong. there are more numbers
        # that can be represented. it's the 27 ABOVE and BELOW the max listed
        # numeric data type in [U] 12.2.2 of the 11.2 manual
        self.MISSING_VALUES = \
            {
                'b': (-127, 100),
                'h': (-32767, 32740),
                'l': (-2147483647, 2147483620),
                'f': (-1.701e+38, +1.701e+38),
                'd': (-1.798e+308, +8.988e+307)
            }

        self.OLD_TYPE_MAPPING = \
            {
                'i': 252,
                'f': 254,
                'b': 251
            }

    def _decode_bytes(self, str, errors=None):
        if py3compat.PY3:
            return str.decode(self._encoding, errors)
        else:
            return str


class StataReader(StataParser):
    """
    Class for working with a Stata dataset. There are two possibilities for usage:

     * The from_dta() method on the DataFrame class.
       This will return a DataFrame with the Stata dataset. Note that when using the
       from_dta() method, you will not have access to meta-information like variable
       labels or the data label.

     * Work with this object directly. Upon instantiation, the header of the Stata data
       file is read, giving you access to attributes like variable_labels(), data_label(),
       nobs(), ... A DataFrame with the data is returned by the read() method; this will
       also fill up the value_labels. Note that calling the value_labels() method will
       result in an error if the read() method has not been called yet. This is because
       the value labels are stored at the end of a Stata dataset, after the data.

    Parameters
    ----------
    path_or_buf : string or file-like object
        Path to .dta file or object implementing a binary read() functions
    encoding : string, None or encoding
        Encoding used to parse the files. Note that Stata doesn't
        support unicode. None defaults to cp1252.
    """
    def __init__(self, path_or_buf, encoding=None):
        super(StataReader, self).__init__(encoding)
        self.col_sizes = ()
        self._has_string_data = False
        self._missing_values = False
        self._data_read = False
        self._value_labels_read = False
        if isinstance(path_or_buf, str):
            path_or_buf, encoding = get_filepath_or_buffer(path_or_buf, encoding='cp1252')
            if encoding is not None:
                self._encoding = encoding

        if type(path_or_buf) is str:
            self.path_or_buf = open(path_or_buf, 'rb')
        else:
            self.path_or_buf = path_or_buf

        self._read_header()

    def _read_header(self):
        # header
        self.format_version = struct.unpack('b', self.path_or_buf.read(1))[0]
        if self.format_version not in [104, 105, 108, 113, 114, 115]:
            raise ValueError("Version of given Stata file is not 104, 105, 108, 113 (Stata 8/9), 114 (Stata 10/11) or 115 (Stata 12)")
        self.byteorder = self.path_or_buf.read(1) == 0x1 and '>' or '<'
        self.filetype = struct.unpack('b', self.path_or_buf.read(1))[0]
        self.path_or_buf.read(1)  # unused

        self.nvar = struct.unpack(self.byteorder + 'H', self.path_or_buf.read(2))[0]
        self.nobs = struct.unpack(self.byteorder + 'I', self.path_or_buf.read(4))[0]
        if self.format_version > 105:
            self.data_label = self.path_or_buf.read(81)
        else:
            self.data_label = self.path_or_buf.read(32)
        if self.format_version > 104:
            self.time_stamp = self.path_or_buf.read(18)

        # descriptors
        if self.format_version > 108:
            typlist = [ord(self.path_or_buf.read(1)) for i in range(self.nvar)]
        else:
            typlist = [self.OLD_TYPE_MAPPING[self._decode_bytes(self.path_or_buf.read(1))] for i in range(self.nvar)]

        try:
            self.typlist = [self.TYPE_MAP[typ] for typ in typlist]
        except:
            raise ValueError("cannot convert stata types [{0}]".format(','.join(typlist)))
        try:
            self.dtyplist = [self.DTYPE_MAP[typ] for typ in typlist]
        except:
            raise ValueError("cannot convert stata dtypes [{0}]".format(','.join(typlist)))

        if self.format_version > 108:
            self.varlist = [self._null_terminate(self.path_or_buf.read(33)) for i in range(self.nvar)]
        else:
            self.varlist = [self._null_terminate(self.path_or_buf.read(9)) for i in range(self.nvar)]
        self.srtlist = struct.unpack(self.byteorder + ('h' * (self.nvar + 1)), self.path_or_buf.read(2 * (self.nvar + 1)))[:-1]
        if self.format_version > 113:
            self.fmtlist = [self._null_terminate(self.path_or_buf.read(49)) for i in range(self.nvar)]
        elif self.format_version > 104:
            self.fmtlist = [self._null_terminate(self.path_or_buf.read(12)) for i in range(self.nvar)]
        else:
            self.fmtlist = [self._null_terminate(self.path_or_buf.read(7)) for i in range(self.nvar)]
        if self.format_version > 108:
            self.lbllist = [self._null_terminate(self.path_or_buf.read(33)) for i in range(self.nvar)]
        else:
            self.lbllist = [self._null_terminate(self.path_or_buf.read(9)) for i in range(self.nvar)]
        if self.format_version > 105:
            self.vlblist = [self._null_terminate(self.path_or_buf.read(81)) for i in range(self.nvar)]
        else:
            self.vlblist = [self._null_terminate(self.path_or_buf.read(32)) for i in range(self.nvar)]

        # ignore expansion fields (Format 105 and later)
        # When reading, read five bytes; the last four bytes now tell you the
        # size of the next read, which you discard.  You then continue like
        # this until you read 5 bytes of zeros.

        if self.format_version > 104:
            while True:
                data_type = struct.unpack(self.byteorder + 'b', self.path_or_buf.read(1))[0]
                if self.format_version > 108:
                    data_len = struct.unpack(self.byteorder + 'i', self.path_or_buf.read(4))[0]
                else:
                    data_len = struct.unpack(self.byteorder + 'h', self.path_or_buf.read(2))[0]
                if data_type == 0:
                    break
                self.path_or_buf.read(data_len)

        # necessary data to continue parsing
        self.data_location = self.path_or_buf.tell()
        self.has_string_data = len([x for x in self.typlist if type(x) is int]) > 0
        self._col_size()

    def _calcsize(self, fmt):
        return type(fmt) is int and fmt or struct.calcsize(self.byteorder + fmt)

    def _col_size(self, k=None):
        """Calculate size of a data record."""
        if len(self.col_sizes) == 0:
            self.col_sizes = map(lambda x: self._calcsize(x), self.typlist)
        if k is None:
            return self.col_sizes
        else:
            return self.col_sizes[k]

    def _unpack(self, fmt, byt):
        d = struct.unpack(self.byteorder + fmt, byt)[0]
        if fmt[-1] in self.MISSING_VALUES:
            nmin, nmax = self.MISSING_VALUES[fmt[-1]]
            if d < nmin or d > nmax:
                if self._missing_values:
                    return StataMissingValue(nmax, d)
                else:
                    return None
        return d

    def _null_terminate(self, s):
        if py3compat.PY3:  # have bytes not strings, so must decode
            null_byte = b"\0"
            try:
                s = s[:s.index(null_byte)]
            except:
                pass
            return s.decode(self._encoding)
        else:
            null_byte = "\0"
            try:
                return s.lstrip(null_byte)[:s.index(null_byte)]
            except:
                return s

    def _next(self):
        typlist = self.typlist
        if self.has_string_data:
            data = [None] * self.nvar
            for i in range(len(data)):
                if type(typlist[i]) is int:
                    data[i] = self._null_terminate(self.path_or_buf.read(typlist[i]))
                else:
                    data[i] = self._unpack(typlist[i], self.path_or_buf.read(self._col_size(i)))
            return data
        else:
            return map(lambda i: self._unpack(typlist[i],
                                              self.path_or_buf.read(self._col_size(i))),
                       range(self.nvar))

    def _dataset(self):
        """
        Returns a Python generator object for iterating over the dataset.


        Parameters
        ----------

        Returns
        -------
        Generator object for iterating over the dataset.  Yields each row of
        observations as a list by default.

        Notes
        -----
        If missing_values is True during instantiation of StataReader then
        observations with _StataMissingValue(s) are not filtered and should
        be handled by your applcation.
        """

        try:
            self._file.seek(self._data_location)
        except Exception:
            pass

        for i in range(self.nobs):
            yield self._next()

    def _read_value_labels(self):
        if not self._data_read:
            raise Exception("Data has not been read. Because of the layout of Stata files, this is necessary before reading value labels.")
        if self._value_labels_read:
            raise Exception("Value labels have already been read.")

        self.value_label_dict = dict()

        if self.format_version <= 108:
            return  # Value labels are not supported in version 108 and earlier.

        while True:
            slength = self.path_or_buf.read(4)
            if not slength:
                break  # end of variable lable table
            labname = self._null_terminate(self.path_or_buf.read(33))
            self.path_or_buf.read(3)  # padding

            n = struct.unpack(self.byteorder + 'I', self.path_or_buf.read(4))[0]
            txtlen = struct.unpack(self.byteorder + 'I', self.path_or_buf.read(4))[0]
            off = []
            for i in range(n):
                off.append(struct.unpack(self.byteorder + 'I', self.path_or_buf.read(4))[0])
            val = []
            for i in range(n):
                val.append(struct.unpack(self.byteorder + 'I', self.path_or_buf.read(4))[0])
            txt = self.path_or_buf.read(txtlen)
            self.value_label_dict[labname] = dict()
            for i in range(n):
                self.value_label_dict[labname][val[i]] = self._null_terminate(txt[off[i]:])
        self._value_labels_read = True

    def data(self, convert_dates=True, convert_categoricals=True, index=None):
        """
        Reads observations from Stata file, converting them into a dataframe

        Parameters
        ----------
        convert_dates : boolean, defaults to True
            Convert date variables to DataFrame time values
        convert_categoricals : boolean, defaults to True
            Read value labels and convert columns to Categorical/Factor variables
        index : identifier of index column
            identifier of column that should be used as index of the DataFrame

        Returns
        -------
        y : DataFrame instance
        """
        if self._data_read:
            raise Exception("Data has already been read.")
        self._data_read = True

        stata_dta = self._dataset()

        data = []
        for rownum, line in enumerate(stata_dta):
            # doesn't handle missing value objects, just casts
            # None will only work without missing value object.
            for i, val in enumerate(line):
                #NOTE: This will only be scalar types because missing strings
                # are empty not None in Stata
                if val is None:
                    line[i] = np.nan
            data.append(tuple(line))

        if convert_categoricals:
            self._read_value_labels()

        data = DataFrame(data, columns=self.varlist, index=index)

        cols_ = np.where(self.dtyplist)[0]
        for i in cols_:
            if self.dtyplist[i] is not None:
                col = data.columns[i]
                if data[col].dtype is not np.dtype(object):
                    data[col] = Series(data[col], data[col].index, self.dtyplist[i])

        if convert_dates:
            cols = np.where(map(lambda x: x in _date_formats, self.fmtlist))[0]
            for i in cols:
                col = data.columns[i]
                data[col] = data[col].apply(_stata_elapsed_date_to_datetime, args=(self.fmtlist[i],))

        if convert_categoricals:
            cols = np.where(map(lambda x: x in self.value_label_dict.iterkeys(), self.lbllist))[0]
            for i in cols:
                col = data.columns[i]
                labeled_data = np.copy(data[col])
                labeled_data = labeled_data.astype(object)
                for k, v in self.value_label_dict[self.lbllist[i]].iteritems():
                    labeled_data[data[col] == k] = v
                data[col] = Categorical.from_array(labeled_data)

        return data

    def data_label(self):
        """Returns data label of Stata file"""
        return self.data_label

    def variable_labels(self):
        """Returns variable labels as a dict, associating each variable name with corresponding label"""
        return dict(zip(self.varlist, self.vlblist))

    def value_labels(self):
        """Returns a dict, associating each variable name a dict, associating each value its corresponding label"""
        if not self._value_labels_read:
            self._read_value_labels()

        return self.value_label_dict


def _open_file_binary_write(fname, encoding):
    if hasattr(fname, 'write'):
        #if 'b' not in fname.mode:
        return fname
    return open(fname, "wb")


def _set_endianness(endianness):
    if endianness.lower() in ["<", "little"]:
        return "<"
    elif endianness.lower() in [">", "big"]:
        return ">"
    else:  # pragma : no cover
        raise ValueError("Endianness %s not understood" % endianness)


def _pad_bytes(name, length):
    """
    Takes a char string and pads it wih null bytes until it's length chars
    """
    return name + "\x00" * (length - len(name))


def _default_names(nvar):
    """
    Returns default Stata names v1, v2, ... vnvar
    """
    return ["v%d" % i for i in range(1, nvar+1)]


def _convert_datetime_to_stata_type(fmt):
    """
    Converts from one of the stata date formats to a type in TYPE_MAP
    """
    if fmt in ["tc", "%tc", "td", "%td", "tw", "%tw", "tm", "%tm", "tq",
               "%tq", "th", "%th", "ty", "%ty"]:
        return np.float64  # Stata expects doubles for SIFs
    else:
        raise ValueError("fmt %s not understood" % fmt)


def _maybe_convert_to_int_keys(convert_dates, varlist):
    new_dict = {}
    for key in convert_dates:
        if not convert_dates[key].startswith("%"):  # make sure proper fmts
            convert_dates[key] = "%" + convert_dates[key]
        if key in varlist:
            new_dict.update({varlist.index(key): convert_dates[key]})
        else:
            if not isinstance(key, int):
                raise ValueError("convery_dates key is not in varlist and is not an int")
            new_dict.update({key: convert_dates[key]})
    return new_dict


def _dtype_to_stata_type(dtype):
    """
    Converts dtype types to stata types. Returns the byte of the given ordinal.
    See TYPE_MAP and comments for an explanation. This is also explained in
    the dta spec.
    1 - 244 are strings of this length
    251 - chr(251) - for int8 and int16, byte
    252 - chr(252) - for int32, int
    253 - chr(253) - for int64, long
    254 - chr(254) - for float32, float
    255 - chr(255) - double, double

    If there are dates to convert, then dtype will already have the correct
    type inserted.
    """
    #TODO: expand to handle datetime to integer conversion
    if dtype.type == np.string_:
        return chr(dtype.itemsize)
    elif dtype.type == np.object_:  # try to coerce it to the biggest string
                                    # not memory efficient, what else could we do?
        return chr(244)
    elif dtype == np.float64:
        return chr(255)
    elif dtype == np.float32:
        return chr(254)
    elif dtype == np.int64:
        return chr(253)
    elif dtype == np.int32:
        return chr(252)
    elif dtype == np.int8 or dtype == np.int16:
        return chr(251)
    else:  # pragma : no cover
        raise ValueError("Data type %s not currently understood. "
                         "Please report an error to the developers." % dtype)


def _dtype_to_default_stata_fmt(dtype):
    """
    Maps numpy dtype to stata's default format for this type. Not terribly
    important since users can change this in Stata. Semantics are

    string  -> "%DDs" where DD is the length of the string
    float64 -> "%10.0g"
    float32 -> "%9.0g"
    int64   -> "%9.0g"
    int32   -> "%12.0g"
    int16   -> "%8.0g"
    int8    -> "%8.0g"
    """
    #TODO: expand this to handle a default datetime format?
    if dtype.type == np.string_:
        return "%" + str(dtype.itemsize) + "s"
    elif dtype.type == np.object_:
        return "%244s"
    elif dtype == np.float64:
        return "%10.0g"
    elif dtype == np.float32:
        return "%9.0g"
    elif dtype == np.int64:
        return "%9.0g"
    elif dtype == np.int32:
        return "%12.0g"
    elif dtype == np.int8 or dtype == np.int16:
        return "%8.0g"
    else:  # pragma : no cover
        raise ValueError("Data type %s not currently understood. "
                         "Please report an error to the developers." % dtype)


class StataWriter(StataParser):
    """
    A class for writing Stata binary dta files from array-like objects

    Parameters
    ----------
    fname : file path or buffer
        Where to save the dta file.
    data : array-like
        Array-like input to save. Pandas objects are also accepted.
    convert_dates : dict
        Dictionary mapping column of datetime types to the stata internal
        format that you want to use for the dates. Options are
        'tc', 'td', 'tm', 'tw', 'th', 'tq', 'ty'. Column can be either a
        number or a name.
    encoding : str
        Default is latin-1. Note that Stata does not support unicode.
    byteorder : str
        Can be ">", "<", "little", or "big". The default is None which uses
        `sys.byteorder`

    Returns
    -------
    writer : StataWriter instance
        The StataWriter instance has a write_file method, which will
        write the file to the given `fname`.

    Examples
    --------
    >>> writer = StataWriter('./data_file.dta', data)
    >>> writer.write_file()

    Or with dates

    >>> writer = StataWriter('./date_data_file.dta', date, {2 : 'tw'})
    >>> writer.write_file()
    """
    def __init__(self, fname, data, convert_dates=None, write_index=True, encoding="latin-1",
                 byteorder=None):
        super(StataWriter, self).__init__(encoding)
        self._convert_dates = convert_dates
        self._write_index = write_index
        # attach nobs, nvars, data, varlist, typlist
        self._prepare_pandas(data)

        if byteorder is None:
            byteorder = sys.byteorder
        self._byteorder = _set_endianness(byteorder)
        self._file = _open_file_binary_write(fname, self._encoding)
        self.type_converters = {253: np.long, 252: int}

    def _write(self, to_write):
        """
        Helper to call encode before writing to file for Python 3 compat.
        """
        if py3compat.PY3:
            self._file.write(to_write.encode(self._encoding))
        else:
            self._file.write(to_write)

    def _prepare_pandas(self, data):
        #NOTE: we might need a different API / class for pandas objects so
        # we can set different semantics - handle this with a PR to pandas.io
        class DataFrameRowIter(object):
            def __init__(self, data):
                self.data = data

            def __iter__(self):
                for i, row in data.iterrows():
                    yield row

        if self._write_index:
            data = data.reset_index()
        self.datarows = DataFrameRowIter(data)
        self.nobs, self.nvar = data.shape
        self.data = data
        self.varlist = data.columns.tolist()
        dtypes = data.dtypes
        if self._convert_dates is not None:
            self._convert_dates = _maybe_convert_to_int_keys(self._convert_dates, self.varlist)
            for key in self._convert_dates:
                new_type = _convert_datetime_to_stata_type(self._convert_dates[key])
                dtypes[key] = np.dtype(new_type)
        self.typlist = [_dtype_to_stata_type(dt) for dt in dtypes]
        self.fmtlist = [_dtype_to_default_stata_fmt(dt) for dt in dtypes]
        # set the given format for the datetime cols
        if self._convert_dates is not None:
            for key in self._convert_dates:
                self.fmtlist[key] = self._convert_dates[key]

    def write_file(self):
        self._write_header()
        self._write_descriptors()
        self._write_variable_labels()
        # write 5 zeros for expansion fields
        self._write(_pad_bytes("", 5))
        if self._convert_dates is None:
            self._write_data_nodates()
        else:
            self._write_data_dates()
        #self._write_value_labels()
        self._file.close()

    def _write_header(self, data_label=None, time_stamp=None):
        byteorder = self._byteorder
        # ds_format - just use 114
        self._file.write(struct.pack("b", 114))
        # byteorder
        self._write(byteorder == ">" and "\x01" or "\x02")
        # filetype
        self._write("\x01")
        # unused
        self._write("\x00")
        # number of vars, 2 bytes
        self._file.write(struct.pack(byteorder+"h", self.nvar)[:2])
        # number of obs, 4 bytes
        self._file.write(struct.pack(byteorder+"i", self.nobs)[:4])
        # data label 81 bytes, char, null terminated
        if data_label is None:
            self._file.write(self._null_terminate(_pad_bytes("", 80)))
        else:
            self._file.write(self._null_terminate(_pad_bytes(data_label[:80], 80)))
        # time stamp, 18 bytes, char, null terminated
        # format dd Mon yyyy hh:mm
        if time_stamp is None:
            time_stamp = datetime.datetime.now()
        elif not isinstance(time_stamp, datetime):
            raise ValueError("time_stamp should be datetime type")
        self._file.write(self._null_terminate(time_stamp.strftime("%d %b %Y %H:%M")))

    def _write_descriptors(self, typlist=None, varlist=None, srtlist=None,
                           fmtlist=None, lbllist=None):
        nvar = self.nvar
        # typlist, length nvar, format byte array
        for typ in self.typlist:
            self._write(typ)

        # varlist, length 33*nvar, char array, null terminated
        for name in self.varlist:
            name = self._null_terminate(name, True)
            name = _pad_bytes(name[:32], 33)
            self._write(name)

        # srtlist, 2*(nvar+1), int array, encoded by byteorder
        srtlist = _pad_bytes("", (2*(nvar+1)))
        self._write(srtlist)

        # fmtlist, 49*nvar, char array
        for fmt in self.fmtlist:
            self._write(_pad_bytes(fmt, 49))

        # lbllist, 33*nvar, char array
        #NOTE: this is where you could get fancy with pandas categorical type
        for i in range(nvar):
            self._write(_pad_bytes("", 33))

    def _write_variable_labels(self, labels=None):
        nvar = self.nvar
        if labels is None:
            for i in range(nvar):
                self._write(_pad_bytes("", 81))

    def _write_data_nodates(self):
        data = self.datarows
        byteorder = self._byteorder
        TYPE_MAP = self.TYPE_MAP
        typlist = self.typlist
        for row in data:
            #row = row.squeeze().tolist() # needed for structured arrays
            for i, var in enumerate(row):
                typ = ord(typlist[i])
                if typ <= 244:  # we've got a string
                    if len(var) < typ:
                        var = _pad_bytes(var, typ)
                    self._write(var)
                else:
                    try:
                        self._file.write(struct.pack(byteorder + TYPE_MAP[typ], var))
                    except struct.error:
                        # have to be strict about type pack won't do any
                        # kind of casting
                        self._file.write(struct.pack(byteorder+TYPE_MAP[typ],
                                         self.type_converters[typ](var)))

    def _write_data_dates(self):
        convert_dates = self._convert_dates
        data = self.datarows
        byteorder = self._byteorder
        TYPE_MAP = self.TYPE_MAP
        MISSING_VALUES = self.MISSING_VALUES
        typlist = self.typlist
        for row in data:
            #row = row.squeeze().tolist() # needed for structured arrays
            for i, var in enumerate(row):
                typ = ord(typlist[i])
                #NOTE: If anyone finds this terribly slow, there is
                # a vectorized way to convert dates, see genfromdta for going
                # from int to datetime and reverse it. will copy data though
                if i in convert_dates:
                    var = _datetime_to_stata_elapsed(var, self.fmtlist[i])
                if typ <= 244:  # we've got a string
                    if len(var) < typ:
                        var = _pad_bytes(var, typ)
                    self._write(var)
                else:
                    if isnull(var):  # this only matters for floats
                        var = MISSING_VALUES[typ]
                    self._file.write(struct.pack(byteorder+TYPE_MAP[typ], var))

    def _null_terminate(self, s, as_string=False):
        null_byte = '\x00'
        if py3compat.PY3 and not as_string:
            s += null_byte
            return s.encode(self._encoding)
        else:
            s += null_byte
            return s
