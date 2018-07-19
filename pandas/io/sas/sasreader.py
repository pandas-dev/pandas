"""
Read SAS sas7bdat or xport files.
"""
from pandas import compat
from pandas.io.common import _stringify_path


def read_sas(filepath_or_buffer, format=None, index=None, encoding=None,
             chunksize=None, iterator=False, convert_dates=True):
    """
    Read SAS files stored as either XPORT or SAS7BDAT format files.

    Parameters
    ----------
    filepath_or_buffer : string or file-like object
        Path to the SAS file.
    format : string {'xport', 'sas7bdat'} or None
        If None, file format is inferred.  If 'xport' or 'sas7bdat',
        uses the corresponding format.
    index : identifier of index column, defaults to None
        Identifier of column that should be used as index of the DataFrame.
    encoding : string, default is None
        Encoding for text data.  If None, text data are stored as raw bytes.
    chunksize : int
        Read file `chunksize` lines at a time, returns iterator.
    iterator : bool, defaults to False
        If True, returns an iterator for reading the file incrementally.
    convert_dates: bool, default to True
        If True convert SAS date and datetime columns to Pandas datetime
        NB.  For datetimes larger than pd.Timestamp.max
        '2262-04-11 23:47:16.854775807' an exception
        pandas._libs.tslibs.np_datetime.OutOfBoundsDatetime is thrown
        If False SAS date and datetime columns are read as their native
        float64 and can be converted after the import to
        datetime.datetime or datetime.date values (or high values capped
        to pd.Timestamp.max and converted with pandas.to_datetime)

    Returns
    -------
    DataFrame if iterator=False and chunksize=None, else SAS7BDATReader
    or XportReader
    """
    if format is None:
        buffer_error_msg = ("If this is a buffer object rather "
                            "than a string name, you must specify "
                            "a format string")
        filepath_or_buffer = _stringify_path(filepath_or_buffer)
        if not isinstance(filepath_or_buffer, compat.string_types):
            raise ValueError(buffer_error_msg)
        try:
            fname = filepath_or_buffer.lower()
            if fname.endswith(".xpt"):
                format = "xport"
            elif fname.endswith(".sas7bdat"):
                format = "sas7bdat"
            else:
                raise ValueError("unable to infer format of SAS file")
        except:
            pass

    if format.lower() == 'xport':
        from pandas.io.sas.sas_xport import XportReader
        reader = XportReader(filepath_or_buffer, index=index,
                             encoding=encoding,
                             chunksize=chunksize)
    elif format.lower() == 'sas7bdat':
        from pandas.io.sas.sas7bdat import SAS7BDATReader
        reader = SAS7BDATReader(filepath_or_buffer, index=index,
                                encoding=encoding,
                                chunksize=chunksize,
                                convert_dates=convert_dates)
    else:
        raise ValueError('unknown SAS format')

    if iterator or chunksize:
        return reader

    data = reader.read()
    reader.close()
    return data
