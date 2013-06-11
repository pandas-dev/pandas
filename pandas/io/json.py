
# pylint: disable-msg=E1101,W0613,W0603
from StringIO import StringIO
import os

from pandas import Series, DataFrame, to_datetime
from pandas.io.common import get_filepath_or_buffer
import pandas.json as _json
loads = _json.loads
dumps = _json.dumps

import numpy as np
from pandas.tslib import iNaT

### interface to/from ###

def to_json(path_or_buf, obj, orient=None, date_format='epoch', double_precision=10, force_ascii=True):
        
    if isinstance(obj, Series):
        s = SeriesWriter(obj, orient=orient, date_format=date_format, double_precision=double_precision, 
                         ensure_ascii=force_ascii).write()
    elif isinstance(obj, DataFrame):
        s = FrameWriter(obj, orient=orient, date_format=date_format, double_precision=double_precision,
                        ensure_ascii=force_ascii).write()
    else:
        raise NotImplementedError

    if isinstance(path_or_buf, basestring):
        with open(path_or_buf,'w') as fh:
            fh.write(s)
    elif path_or_buf is None:
        return s
    else:
        path_or_buf.write(s)

class Writer(object):

    def __init__(self, obj, orient, date_format, double_precision, ensure_ascii):
        self.obj = obj

        if orient is None:
            orient = self._default_orient
            
        self.orient = orient
        self.date_format = date_format
        self.double_precision = double_precision
        self.ensure_ascii = ensure_ascii

        self.is_copy = False
        self._format_axes()
        self._format_dates()

    def _format_dates(self):
        raise NotImplementedError

    def _format_axes(self):
        raise NotImplementedError

    def _needs_to_date(self, data):
        return self.date_format == 'iso' and data.dtype == 'datetime64[ns]'

    def _format_to_date(self, data):
        if self._needs_to_date(data):
            return data.apply(lambda x: x.isoformat())
        return data
    
    def copy_if_needed(self):
        """ copy myself if necessary """
        if not self.is_copy:
            self.obj = self.obj.copy()
            self.is_copy = True

    def write(self):
        return dumps(self.obj, orient=self.orient, double_precision=self.double_precision, ensure_ascii=self.ensure_ascii)

class SeriesWriter(Writer):
    _default_orient = 'index'

    def _format_axes(self):
        if self._needs_to_date(self.obj.index):
            self.copy_if_needed()
            self.obj.index = self._format_to_date(self.obj.index.to_series())

    def _format_dates(self):
        if self._needs_to_date(self.obj):
            self.copy_if_needed()
            self.obj = self._format_to_date(self.obj)

class FrameWriter(Writer):
    _default_orient = 'columns'

    def _format_axes(self):
        """ try to axes if they are datelike """
        if self.orient == 'columns':
            axis = 'index'
        elif self.orient == 'index':
            axis = 'columns'
        else:
            return

        a = getattr(self.obj,axis)
        if self._needs_to_date(a):
            self.copy_if_needed()
            setattr(self.obj,axis,self._format_to_date(a.to_series()))

    def _format_dates(self):
        if self.date_format == 'iso':
            dtypes = self.obj.dtypes
            dtypes = dtypes[dtypes == 'datetime64[ns]']
            if len(dtypes):
                self.copy_if_needed()
                for c in dtypes.index:
                    self.obj[c] = self._format_to_date(self.obj[c])

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
        obj = SeriesParser(json, orient, dtype, numpy, parse_dates, keep_default_dates).parse()

    return obj

class Parser(object):
    
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
            self._convert_axes()
            if self.parse_dates:
                self._try_parse_dates()
        return self.obj


    def _try_parse_to_date(self, data):
        """ try to parse a ndarray like into a date column
            try to coerce object in epoch/iso formats and
            integer/float in epcoh formats """

        new_data = data
        if new_data.dtype == 'object':
            try:
                new_data = data.astype('int64')
            except:
                pass


        # ignore numbers that are out of range
        if issubclass(new_data.dtype.type,np.number):
            if not ((new_data == iNaT) | (new_data > 31536000000000000L)).all():
                return data
                
        try:
            new_data = to_datetime(new_data)
        except:
            try:
                new_data = to_datetime(new_data.astype('int64'))
            except:

                # return old, noting more we can do
                new_data = data

        return new_data

    def _try_parse_dates(self):
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

    def _convert_axes(self):
        """ try to axes if they are datelike """
        try:
           self.obj.index = self._try_parse_to_date(self.obj.index)
        except:
           pass

    def _try_parse_dates(self):
        if self.obj is None: return

        if self.parse_dates:
            self.obj = self._try_parse_to_date(self.obj)

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

    def _convert_axes(self):
        """ try to axes if they are datelike """
        if self.orient == 'columns':
            axis = 'index'
        elif self.orient == 'index':
            axis = 'columns'
        else:
            return

        try:
            a = getattr(self.obj,axis)
            setattr(self.obj,axis,self._try_parse_to_date(a))
        except:
            pass

    def _try_parse_dates(self):
        if self.obj is None: return

        # our columns to parse
        parse_dates = self.parse_dates
        if parse_dates is True:
            parse_dates = []
        parse_dates = set(parse_dates)

        def is_ok(col):
            """ return if this col is ok to try for a date parse """
            if not isinstance(col, basestring): return False

            if (col.endswith('_at') or
                col.endswith('_time') or
                col.lower() == 'modified' or
                col.lower() == 'date' or
                col.lower() == 'datetime'):
                    return True
            return False


        for col, c in self.obj.iteritems():
            if (self.keep_default_dates and is_ok(col)) or col in parse_dates:
                self.obj[col] = self._try_parse_to_date(c)
