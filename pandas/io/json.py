
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
import pandas.lib as lib

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

    def _format_bools(self):
        if self._needs_to_bool(self.obj):
            self.copy_if_needed()
            self.obj = self._format_to_bool(self.obj)

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

def read_json(path_or_buf=None, orient=None, typ='frame', dtype=True,
              convert_axes=True, convert_dates=True, keep_default_dates=True,
              numpy=False, precise_float=False):
    """
    Convert JSON string to pandas object

    Parameters
    ----------
    filepath_or_buffer : a VALID JSON string or file handle / StringIO. The string could be
        a URL. Valid URL schemes include http, ftp, s3, and file. For file URLs, a host
        is expected. For instance, a local file could be
        file ://localhost/path/to/table.json
    orient :
        Series :
          default is 'index'
          allowed values are: {'split','records','index'}

        DataFrame :
          default is 'columns'
          allowed values are: {'split','records','index','columns','values'}

        The format of the JSON string
          split : dict like {index -> [index], columns -> [columns], data -> [values]}
          records : list like [{column -> value}, ... , {column -> value}]
          index : dict like {index -> {column -> value}}
          columns : dict like {column -> {index -> value}}
          values : just the values array

    typ : type of object to recover (series or frame), default 'frame'
    dtype : if True, infer dtypes, if a dict of column to dtype, then use those,
        if False, then don't infer dtypes at all, default is True,
        apply only to the data
    convert_axes : boolean, try to convert the axes to the proper dtypes, default is True
    convert_dates : a list of columns to parse for dates; If True, then try to parse datelike columns
        default is True
    keep_default_dates : boolean, default True. If parsing dates,
        then parse the default datelike columns
    numpy : direct decoding to numpy arrays. default is False.Note that the JSON ordering MUST be the same
        for each term if numpy=True.
    precise_float : boolean, default False. Set to enable usage of higher precision (strtod) function
        when decoding string to double values. Default (False) is to use fast but less precise builtin functionality

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
        obj = FrameParser(json, orient, dtype, convert_axes, convert_dates, keep_default_dates, numpy).parse()

    if typ == 'series' or obj is None:
        if not isinstance(dtype,bool):
            dtype = dict(data = dtype)
        obj = SeriesParser(json, orient, dtype, convert_axes, convert_dates, keep_default_dates, numpy).parse()

    return obj

class Parser(object):

    def __init__(self, json, orient, dtype=True, convert_axes=True,
                 convert_dates=True, keep_default_dates=False, numpy=False,
                 precise_float=False):
        self.json = json

        if orient is None:
            orient = self._default_orient

        self.orient = orient
        self.dtype = dtype

        if orient == "split":
            numpy = False

        self.numpy = numpy
        self.precise_float = precise_float
        self.convert_axes  = convert_axes
        self.convert_dates = convert_dates
        self.keep_default_dates = keep_default_dates
        self.obj = None

    def parse(self):

        # try numpy
        numpy = self.numpy
        if numpy:
            self._parse_numpy()

        else:
            self._parse_no_numpy()

        if self.obj is None: return None
        if self.convert_axes:
            self._convert_axes()
        self._try_convert_types()
        return self.obj

    def _convert_axes(self):
        """ try to convert axes """
        for axis in self.obj._AXIS_NUMBERS.keys():
            new_axis, result = self._try_convert_data(axis, self.obj._get_axis(axis), use_dtypes=False, convert_dates=True)
            if result:
                setattr(self.obj,axis,new_axis)

    def _try_convert_types(self):
        raise NotImplementedError

    def _try_convert_data(self, name, data, use_dtypes=True, convert_dates=True):
        """ try to parse a ndarray like into a column by inferring dtype """

        # don't try to coerce, unless a force conversion
        if use_dtypes:
            if self.dtype is False:
                return data, False
            elif self.dtype is True:
                pass

            else:

                # dtype to force
                dtype = self.dtype.get(name) if isinstance(self.dtype,dict) else self.dtype
                if dtype is not None:
                    try:
                        dtype = np.dtype(dtype)
                        return data.astype(dtype), True
                    except:
                        return data, False

        if convert_dates:
            new_data, result = self._try_convert_to_date(data)
            if result:
                return new_data, True

        result = False

        if data.dtype == 'object':

            # try float
            try:
                data = data.astype('float64')
                result = True
            except:
                pass

        if data.dtype == 'float':

            # coerce floats to 64
            try:
                data = data.astype('float64')
                result = True
            except:
                pass

        # do't coerce 0-len data
        if len(data) and (data.dtype == 'float' or data.dtype == 'object'):

            # coerce ints if we can
            try:
                new_data = data.astype('int64')
                if (new_data == data).all():
                    data = new_data
                    result = True
            except:
                pass

        # coerce ints to 64
        if data.dtype == 'int':

            # coerce floats to 64
            try:
                data = data.astype('int64')
                result = True
            except:
                pass

        return data, result

    def _try_convert_to_date(self, data):
        """ try to parse a ndarray like into a date column
            try to coerce object in epoch/iso formats and
            integer/float in epcoh formats, return a boolean if parsing
            was successful """

        # no conversion on empty
        if not len(data): return data, False

        new_data = data
        if new_data.dtype == 'object':
            try:
                new_data = data.astype('int64')
            except:
                pass


        # ignore numbers that are out of range
        if issubclass(new_data.dtype.type,np.number):
            if not ((new_data == iNaT) | (new_data > 31536000000000000L)).all():
                return data, False

        try:
            new_data = to_datetime(new_data)
        except:
            try:
                new_data = to_datetime(new_data.astype('int64'))
            except:

                # return old, noting more we can do
                return data, False

        return new_data, True

    def _try_convert_dates(self):
        raise NotImplementedError

class SeriesParser(Parser):
    _default_orient = 'index'

    def _parse_no_numpy(self):

        json = self.json
        orient = self.orient
        if orient == "split":
            decoded = dict((str(k), v)
                           for k, v in loads(
                               json,
                               precise_float=self.precise_float).iteritems())
            self.obj = Series(dtype=None, **decoded)
        else:
            self.obj = Series(
                loads(json, precise_float=self.precise_float), dtype=None)

    def _parse_numpy(self):

        json = self.json
        orient = self.orient
        if orient == "split":
            decoded = loads(json, dtype=None, numpy=True,
                            precise_float=self.precise_float)
            decoded = dict((str(k), v) for k, v in decoded.iteritems())
            self.obj = Series(**decoded)
        elif orient == "columns" or orient == "index":
            self.obj = Series(*loads(json, dtype=None, numpy=True,
                                     labelled=True,
                                     precise_float=self.precise_float))
        else:
            self.obj = Series(loads(json, dtype=None, numpy=True,
                                    precise_float=self.precise_float))

    def _try_convert_types(self):
        if self.obj is None: return
        obj, result = self._try_convert_data('data', self.obj, convert_dates=self.convert_dates)
        if result:
            self.obj = obj

class FrameParser(Parser):
    _default_orient = 'columns'

    def _parse_numpy(self):

        json = self.json
        orient = self.orient

        if orient == "columns":
            args = loads(json, dtype=None, numpy=True, labelled=True,
                         precise_float=self.precise_float)
            if args:
                args = (args[0].T, args[2], args[1])
            self.obj = DataFrame(*args)
        elif orient == "split":
            decoded = loads(json, dtype=None, numpy=True,
                            precise_float=self.precise_float)
            decoded = dict((str(k), v) for k, v in decoded.iteritems())
            self.obj = DataFrame(**decoded)
        elif orient == "values":
            self.obj = DataFrame(loads(json, dtype=None, numpy=True,
                                       precise_float=self.precise_float))
        else:
            self.obj = DataFrame(*loads(json, dtype=None, numpy=True, labelled=True,
                                        precise_float=self.precise_float))

    def _parse_no_numpy(self):

        json = self.json
        orient = self.orient

        if orient == "columns":
            self.obj = DataFrame(
                loads(json, precise_float=self.precise_float), dtype=None)
        elif orient == "split":
            decoded = dict((str(k), v)
                           for k, v in loads(
                               json,
                               precise_float=self.precise_float).iteritems())
            self.obj = DataFrame(dtype=None, **decoded)
        elif orient == "index":
            self.obj = DataFrame(
                loads(json, precise_float=self.precise_float), dtype=None).T
        else:
            self.obj = DataFrame(
                loads(json, precise_float=self.precise_float), dtype=None)

    def _try_convert_types(self):
        if self.obj is None: return
        if self.convert_dates:
            self._try_convert_dates()
        for col in self.obj.columns:
            new_data, result = self._try_convert_data(col, self.obj[col], convert_dates=False)
            if result:
                self.obj[col] = new_data

    def _try_convert_dates(self):
        if self.obj is None: return

        # our columns to parse
        convert_dates = self.convert_dates
        if convert_dates is True:
            convert_dates = []
        convert_dates = set(convert_dates)

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


        for col in self.obj.columns:
            if (self.keep_default_dates and is_ok(col)) or col in convert_dates:
                new_data, result = self._try_convert_to_date(self.obj[col])
                if result:
                    self.obj[col] = new_data
