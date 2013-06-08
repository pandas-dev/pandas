
# pylint: disable-msg=E1101,W0613,W0603
from pandas import Series, DataFrame
from pandas.io.common import get_filepath_or_buffer
from StringIO import StringIO

import pandas.json as _json
loads = _json.loads
dumps = _json.dumps

### interface to/from ###

def to_json(path_or_buf, obj, orient=None, double_precision=10,
            force_ascii=True):
        """
        Convert the object to a JSON string

        Note NaN's and None will be converted to null and datetime objects
        will be converted to UNIX timestamps.

        Parameters
        ----------
        path_or_buf : the pathname or buffer to write the output
            if this is None, return a StringIO of the converted string
        orient : {'split', 'records', 'index', 'columns', 'values'},
            default is 'index' for Series, 'columns' for DataFrame

            The format of the JSON string
            split : dict like
                {index -> [index], columns -> [columns], data -> [values]}
            records : list like [{column -> value}, ... , {column -> value}]
            index : dict like {index -> {column -> value}}
            columns : dict like {column -> {index -> value}}
            values : just the values array
        double_precision : The number of decimal places to use when encoding
            floating point values, default 10.
        force_ascii : force encoded string to be ASCII, default True.

        Returns
        -------
        result : a JSON compatible string written to the path_or_buf;
                 if the path_or_buf is none, return a StringIO of the result

        """
        
        if orient is None:
            if isinstance(obj, Series):
                orient = 'index'
            elif isinstance(obj, DataFrame):
                orient = 'columns'

        s = dumps(obj, orient=orient, double_precision=double_precision,
                          ensure_ascii=force_ascii)
        if isinstance(path_or_buf, basestring):
            with open(path_or_buf,'w') as fh:
                fh.write(s)
        elif path_or_buf is None:
            return StringIO(s)
        else:
            path_or_buf.write(s)

def read_json(filepath_or_buffer=None, json=None, typ='frame', orient=None, dtype=None, numpy=True):
    """
    Convert JSON string to pandas object

    Parameters
    ----------
    filepath_or_buffer : a VALID JSON StringIO or file handle / StringIO. The string could be
        a URL. Valid URL schemes include http, ftp, s3, and file. For file URLs, a host
        is expected. For instance, a local file could be
        file ://localhost/path/to/table.json
    json : a VALID JSON string, optional, used if filepath_or_buffer is not provided
    typ : type of object to recover (series or frame), default 'frame'
    orient : {'split', 'records', 'index'}, default 'index'
        The format of the JSON string
        split : dict like
            {index -> [index], name -> name, data -> [values]}
        records : list like [value, ... , value]
        index : dict like {index -> value}
    dtype : dtype of the resulting Series
    numpy: direct decoding to numpy arrays. default True but falls back
        to standard decoding if a problem occurs.

    Returns
    -------
    result : Series or DataFrame
    """

    if json is None:
        filepath_or_buffer,_ = get_filepath_or_buffer(filepath_or_buffer)
        if isinstance(filepath_or_buffer, basestring):
                with open(filepath_or_buffer,'r') as fh:
                        json = fh.read()
        elif hasattr(filepath_or_buffer, 'read'):
                json = filepath_or_buffer.read()
        else:
                json = filepath_or_buffer

    obj = None
    if typ == 'frame':
        if orient is None:
            orient = 'columns'
        obj = load_frame(json, orient, dtype, numpy)

    if typ == 'series' or obj is None:
        if orient == 'columns':
            orient = 'index'
        obj = load_series(json, orient, dtype, numpy)

    return obj

def load_series(json, orient, dtype, numpy):
    s = None

    if dtype is not None and orient == "split":
        numpy = False

    if numpy:
        try:
            if orient == "split":
                decoded = loads(json, dtype=dtype, numpy=True)
                decoded = dict((str(k), v) for k, v in decoded.iteritems())
                s = Series(**decoded)
            elif orient == "columns" or orient == "index":
                s = Series(*loads(json, dtype=dtype, numpy=True,
                                  labelled=True))
            else:
                s = Series(loads(json, dtype=dtype, numpy=True))
        except ValueError:
            numpy = False

    if not numpy:
        if orient == "split":
            decoded = dict((str(k), v)
                           for k, v in loads(json).iteritems())
            s = Series(dtype=dtype, **decoded)
        else:
            s = Series(loads(json), dtype=dtype)

    return s

        
def load_frame(json, orient, dtype, numpy):
    """ try to recover a frame, return None if we didn't get anything """

    if dtype is not None and orient == "split":
        numpy = False

    if numpy:
        try:
            if orient == "columns":
                args = loads(json, dtype=dtype, numpy=True, labelled=True)
                if args:
                    args = (args[0].T, args[2], args[1])
                df = DataFrame(*args)
            elif orient == "split":
                decoded = loads(json, dtype=dtype, numpy=True)
                decoded = dict((str(k), v) for k, v in decoded.iteritems())
                df = DataFrame(**decoded)
            elif orient == "values":
                df = DataFrame(loads(json, dtype=dtype, numpy=True))
            else:
                df = DataFrame(*loads(json, dtype=dtype, numpy=True,
                                     labelled=True))
        except ValueError:
            numpy = False

    if not numpy:
        if orient == "columns":
            df = DataFrame(loads(json), dtype=dtype)
        elif orient == "split":
            decoded = dict((str(k), v)
                           for k, v in loads(json).iteritems())
            df = DataFrame(dtype=dtype, **decoded)
        elif orient == "index":
            df = DataFrame(loads(json), dtype=dtype).T
        else:
            df = DataFrame(loads(json), dtype=dtype)

    return df
