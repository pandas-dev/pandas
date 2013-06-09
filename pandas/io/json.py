
# pylint: disable-msg=E1101,W0613,W0603
from StringIO import StringIO
import os

from pandas import Series, DataFrame, to_datetime
from pandas.io.common import get_filepath_or_buffer
import pandas.json as _json
loads = _json.loads
dumps = _json.dumps

import numpy as np

### interface to/from ###

def to_json(path_or_buf, obj, orient=None, double_precision=10, force_ascii=True):
        
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
            return s
        else:
            path_or_buf.write(s)

def read_json(path_or_buf=None, orient=None, typ='frame', dtype=None, numpy=True,
              parse_dates=False, keep_default_dates=True):
    """
    Convert JSON string to pandas object

    Parameters
    ----------
    filepath_or_buffer : a VALID JSON string or file handle / StringIO. The string could be
        a URL. Valid URL schemes include http, ftp, s3, and file. For file URLs, a host
        is expected. For instance, a local file could be
        file ://localhost/path/to/table.json
    orient : {'split', 'records', 'index'}, default 'index'
        The format of the JSON string
        split : dict like
            {index -> [index], name -> name, data -> [values]}
        records : list like [value, ... , value]
        index : dict like {index -> value}
    typ : type of object to recover (series or frame), default 'frame'
    dtype : dtype of the resulting object
    numpy: direct decoding to numpy arrays. default True but falls back
        to standard decoding if a problem occurs.
    parse_dates : a list of columns to parse for dates; If True, then try to parse datelike columns
        default is False
    keep_default_dates : boolean, default True. If parsing dates,
        then parse the default datelike columns

    Returns
    -------
    result : Series or DataFrame
    """

    filepath_or_buffer,_ = get_filepath_or_buffer(path_or_buf)
    if isinstance(filepath_or_buffer, basestring):
        if os.path.exists(filepath_or_buffer):
            with open(filepath_or_buffer,'r') as fh:
                json = fh.read()
        else:
            json = filepath_or_buffer
    elif hasattr(filepath_or_buffer, 'read'):
        json = filepath_or_buffer.read()
    else:
        json = filepath_or_buffer

    obj = None
    if typ == 'frame':
        obj = FrameParser(json, orient, dtype, numpy, parse_dates, keep_default_dates).parse()

    if typ == 'series' or obj is None:
        obj = SeriesParser(json, orient, dtype, numpy).parse()

    return obj

class Parser(object):
    _min_date = 31536000000000000L
    
    def __init__(self, json, orient, dtype, numpy, parse_dates=False, keep_default_dates=False):
        self.json = json

        if orient is None:
            orient = self._default_orient
            
        self.orient = orient
        self.dtype = dtype

        if dtype is not None and orient == "split":
            numpy = False

        self.numpy = numpy
        self.parse_dates = parse_dates
        self.keep_default_dates = keep_default_dates
        self.obj = None

    def parse(self):
        self._parse()
        if self.obj is not None:
            self.convert_axes()
            if self.parse_dates:
                self.try_parse_dates()
        return self.obj

    def try_parse_dates(self):
        raise NotImplementedError

class SeriesParser(Parser):
    _default_orient = 'index'

    def _parse(self):

        json = self.json
        dtype = self.dtype
        orient = self.orient
        numpy = self.numpy
    
        if numpy:
            try:
                if orient == "split":
                    decoded = loads(json, dtype=dtype, numpy=True)
                    decoded = dict((str(k), v) for k, v in decoded.iteritems())
                    self.obj = Series(**decoded)
                elif orient == "columns" or orient == "index":
                    self.obj = Series(*loads(json, dtype=dtype, numpy=True,
                                             labelled=True))
                else:
                    self.obj = Series(loads(json, dtype=dtype, numpy=True))
            except ValueError:
                numpy = False

        if not numpy:
            if orient == "split":
                decoded = dict((str(k), v)
                               for k, v in loads(json).iteritems())
                self.obj = Series(dtype=dtype, **decoded)
            else:
                self.obj = Series(loads(json), dtype=dtype)

    def convert_axes(self):
        """ try to axes if they are datelike """
        if self.obj is None: return

        try:
           self.obj.index = to_datetime(self.obj.index.astype('int64'))
        except:
           pass

class FrameParser(Parser):
    _default_orient = 'columns'

    def _parse(self):

        json = self.json
        dtype = self.dtype
        orient = self.orient
        numpy = self.numpy

        if numpy:
            try:
                if orient == "columns":
                    args = loads(json, dtype=dtype, numpy=True, labelled=True)
                    if args:
                        args = (args[0].T, args[2], args[1])
                    self.obj = DataFrame(*args)
                elif orient == "split":
                    decoded = loads(json, dtype=dtype, numpy=True)
                    decoded = dict((str(k), v) for k, v in decoded.iteritems())
                    self.obj = DataFrame(**decoded)
                elif orient == "values":
                    self.obj = DataFrame(loads(json, dtype=dtype, numpy=True))
                else:
                    self.obj = DataFrame(*loads(json, dtype=dtype, numpy=True,
                                         labelled=True))
            except ValueError:
                numpy = False

        if not numpy:
            if orient == "columns":
                self.obj = DataFrame(loads(json), dtype=dtype)
            elif orient == "split":
                decoded = dict((str(k), v)
                               for k, v in loads(json).iteritems())
                self.obj = DataFrame(dtype=dtype, **decoded)
            elif orient == "index":
                self.obj = DataFrame(loads(json), dtype=dtype).T
            else:
                self.obj = DataFrame(loads(json), dtype=dtype)

    def convert_axes(self):
        """ try to axes if they are datelike """
        if self.obj is None: return

        if self.orient == 'columns':
            axis = 'index'
        elif self.orient == 'index':
            axis = 'columns'
        else:
            return

        try:
            a = getattr(self.obj,axis).astype('int64')
            if (a>self._min_date).all():
                setattr(self.obj,axis,to_datetime(a))
        except:
            pass

    def try_parse_dates(self):
        """
        try to parse out dates
        these are only in in64 columns
        """

        if self.obj is None: return

        # our columns to parse
        parse_dates = self.parse_dates
        if parse_dates is True:
            parse_dates = []
        parse_dates = set(parse_dates)

        def is_ok(col, c):
            """ return if this col is ok to try for a date parse """
            if not isinstance(col, basestring): return False

            if issubclass(c.dtype.type,np.number) and (c<self._min_date).all():
                return False
                    
            if (col.endswith('_at') or
                col.endswith('_time') or
                col.lower() == 'modified' or
                col.lower() == 'date' or
                col.lower() == 'datetime'):
                    return True
            return False


        for col, c in self.obj.iteritems():
            if (self.keep_default_dates and is_ok(col, c)) or col in parse_dates:
                try:
                    self.obj[col] = to_datetime(c)
                except:
                    try:
                        self.obj[col] = to_datetime(c.astype('int64'))
                    except:
                        pass
