# pylint: disable-msg=E1101,W0613,W0603

import os
import copy
from collections import defaultdict
import numpy as np

import pandas.json as _json
from pandas.tslib import iNaT
from pandas.compat import long, u
from pandas import compat, isnull
from pandas import Series, DataFrame, to_datetime
from pandas.io.common import get_filepath_or_buffer
from pandas.core.common import AbstractMethodError
from pandas.formats.printing import pprint_thing

loads = _json.loads
dumps = _json.dumps

# interface to/from


def to_json(path_or_buf, obj, orient=None, date_format='epoch',
            double_precision=10, force_ascii=True, date_unit='ms',
            default_handler=None):

    if isinstance(obj, Series):
        s = SeriesWriter(
            obj, orient=orient, date_format=date_format,
            double_precision=double_precision, ensure_ascii=force_ascii,
            date_unit=date_unit, default_handler=default_handler).write()
    elif isinstance(obj, DataFrame):
        s = FrameWriter(
            obj, orient=orient, date_format=date_format,
            double_precision=double_precision, ensure_ascii=force_ascii,
            date_unit=date_unit, default_handler=default_handler).write()
    else:
        raise NotImplementedError("'obj' should be a Series or a DataFrame")

    if isinstance(path_or_buf, compat.string_types):
        with open(path_or_buf, 'w') as fh:
            fh.write(s)
    elif path_or_buf is None:
        return s
    else:
        path_or_buf.write(s)


class Writer(object):

    def __init__(self, obj, orient, date_format, double_precision,
                 ensure_ascii, date_unit, default_handler=None):
        self.obj = obj

        if orient is None:
            orient = self._default_orient

        self.orient = orient
        self.date_format = date_format
        self.double_precision = double_precision
        self.ensure_ascii = ensure_ascii
        self.date_unit = date_unit
        self.default_handler = default_handler

        self.is_copy = None
        self._format_axes()

    def _format_axes(self):
        raise AbstractMethodError(self)

    def write(self):
        return dumps(
            self.obj,
            orient=self.orient,
            double_precision=self.double_precision,
            ensure_ascii=self.ensure_ascii,
            date_unit=self.date_unit,
            iso_dates=self.date_format == 'iso',
            default_handler=self.default_handler)


class SeriesWriter(Writer):
    _default_orient = 'index'

    def _format_axes(self):
        if not self.obj.index.is_unique and self.orient == 'index':
            raise ValueError("Series index must be unique for orient="
                             "'%s'" % self.orient)


class FrameWriter(Writer):
    _default_orient = 'columns'

    def _format_axes(self):
        """ try to axes if they are datelike """
        if not self.obj.index.is_unique and self.orient in (
                'index', 'columns'):
            raise ValueError("DataFrame index must be unique for orient="
                             "'%s'." % self.orient)
        if not self.obj.columns.is_unique and self.orient in (
                'index', 'columns', 'records'):
            raise ValueError("DataFrame columns must be unique for orient="
                             "'%s'." % self.orient)


def read_json(path_or_buf=None, orient=None, typ='frame', dtype=True,
              convert_axes=True, convert_dates=True, keep_default_dates=True,
              numpy=False, precise_float=False, date_unit=None):
    """
    Convert a JSON string to pandas object

    Parameters
    ----------
    path_or_buf : a valid JSON string or file-like, default: None
        The string could be a URL. Valid URL schemes include http, ftp, s3, and
        file. For file URLs, a host is expected. For instance, a local file
        could be ``file://localhost/path/to/table.json``

    orient

        * `Series`

          - default is ``'index'``
          - allowed values are: ``{'split','records','index'}``
          - The Series index must be unique for orient ``'index'``.

        * `DataFrame`

          - default is ``'columns'``
          - allowed values are: {'split','records','index','columns','values'}
          - The DataFrame index must be unique for orients 'index' and
            'columns'.
          - The DataFrame columns must be unique for orients 'index',
            'columns', and 'records'.

        * The format of the JSON string

          - split : dict like
            ``{index -> [index], columns -> [columns], data -> [values]}``
          - records : list like
            ``[{column -> value}, ... , {column -> value}]``
          - index : dict like ``{index -> {column -> value}}``
          - columns : dict like ``{column -> {index -> value}}``
          - values : just the values array

    typ : type of object to recover (series or frame), default 'frame'
    dtype : boolean or dict, default True
        If True, infer dtypes, if a dict of column to dtype, then use those,
        if False, then don't infer dtypes at all, applies only to the data.
    convert_axes : boolean, default True
        Try to convert the axes to the proper dtypes.
    convert_dates : boolean, default True
        List of columns to parse for dates; If True, then try to parse
        datelike columns default is True; a column label is datelike if

        * it ends with ``'_at'``,

        * it ends with ``'_time'``,

        * it begins with ``'timestamp'``,

        * it is ``'modified'``, or

        * it is ``'date'``

    keep_default_dates : boolean, default True
        If parsing dates, then parse the default datelike columns
    numpy : boolean, default False
        Direct decoding to numpy arrays. Supports numeric data only, but
        non-numeric column and index labels are supported. Note also that the
        JSON ordering MUST be the same for each term if numpy=True.
    precise_float : boolean, default False
        Set to enable usage of higher precision (strtod) function when
        decoding string to double values. Default (False) is to use fast but
        less precise builtin functionality
    date_unit : string, default None
        The timestamp unit to detect if converting dates. The default behaviour
        is to try and detect the correct precision, but if this is not desired
        then pass one of 's', 'ms', 'us' or 'ns' to force parsing only seconds,
        milliseconds, microseconds or nanoseconds respectively.

    Returns
    -------
    result : Series or DataFrame
    """

    filepath_or_buffer, _, _ = get_filepath_or_buffer(path_or_buf)
    if isinstance(filepath_or_buffer, compat.string_types):
        try:
            exists = os.path.exists(filepath_or_buffer)

        # if the filepath is too long will raise here
        # 5874
        except (TypeError, ValueError):
            exists = False

        if exists:
            with open(filepath_or_buffer, 'r') as fh:
                json = fh.read()
        else:
            json = filepath_or_buffer
    elif hasattr(filepath_or_buffer, 'read'):
        json = filepath_or_buffer.read()
    else:
        json = filepath_or_buffer

    obj = None
    if typ == 'frame':
        obj = FrameParser(json, orient, dtype, convert_axes, convert_dates,
                          keep_default_dates, numpy, precise_float,
                          date_unit).parse()

    if typ == 'series' or obj is None:
        if not isinstance(dtype, bool):
            dtype = dict(data=dtype)
        obj = SeriesParser(json, orient, dtype, convert_axes, convert_dates,
                           keep_default_dates, numpy, precise_float,
                           date_unit).parse()

    return obj


class Parser(object):

    _STAMP_UNITS = ('s', 'ms', 'us', 'ns')
    _MIN_STAMPS = {
        's': long(31536000),
        'ms': long(31536000000),
        'us': long(31536000000000),
        'ns': long(31536000000000000)}

    def __init__(self, json, orient, dtype=True, convert_axes=True,
                 convert_dates=True, keep_default_dates=False, numpy=False,
                 precise_float=False, date_unit=None):
        self.json = json

        if orient is None:
            orient = self._default_orient

        self.orient = orient
        self.dtype = dtype

        if orient == "split":
            numpy = False

        if date_unit is not None:
            date_unit = date_unit.lower()
            if date_unit not in self._STAMP_UNITS:
                raise ValueError('date_unit must be one of %s' %
                                 (self._STAMP_UNITS,))
            self.min_stamp = self._MIN_STAMPS[date_unit]
        else:
            self.min_stamp = self._MIN_STAMPS['s']

        self.numpy = numpy
        self.precise_float = precise_float
        self.convert_axes = convert_axes
        self.convert_dates = convert_dates
        self.date_unit = date_unit
        self.keep_default_dates = keep_default_dates
        self.obj = None

    def check_keys_split(self, decoded):
        "checks that dict has only the appropriate keys for orient='split'"
        bad_keys = set(decoded.keys()).difference(set(self._split_keys))
        if bad_keys:
            bad_keys = ", ".join(bad_keys)
            raise ValueError(u("JSON data had unexpected key(s): %s") %
                             pprint_thing(bad_keys))

    def parse(self):

        # try numpy
        numpy = self.numpy
        if numpy:
            self._parse_numpy()

        else:
            self._parse_no_numpy()

        if self.obj is None:
            return None
        if self.convert_axes:
            self._convert_axes()
        self._try_convert_types()
        return self.obj

    def _convert_axes(self):
        """ try to convert axes """
        for axis in self.obj._AXIS_NUMBERS.keys():
            new_axis, result = self._try_convert_data(
                axis, self.obj._get_axis(axis), use_dtypes=False,
                convert_dates=True)
            if result:
                setattr(self.obj, axis, new_axis)

    def _try_convert_types(self):
        raise AbstractMethodError(self)

    def _try_convert_data(self, name, data, use_dtypes=True,
                          convert_dates=True):
        """ try to parse a ndarray like into a column by inferring dtype """

        # don't try to coerce, unless a force conversion
        if use_dtypes:
            if self.dtype is False:
                return data, False
            elif self.dtype is True:
                pass

            else:

                # dtype to force
                dtype = (self.dtype.get(name)
                         if isinstance(self.dtype, dict) else self.dtype)
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

        if data.dtype.kind == 'f':

            if data.dtype != 'float64':

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
        if not len(data):
            return data, False

        new_data = data
        if new_data.dtype == 'object':
            try:
                new_data = data.astype('int64')
            except:
                pass

        # ignore numbers that are out of range
        if issubclass(new_data.dtype.type, np.number):
            in_range = (isnull(new_data.values) | (new_data > self.min_stamp) |
                        (new_data.values == iNaT))
            if not in_range.all():
                return data, False

        date_units = (self.date_unit,) if self.date_unit else self._STAMP_UNITS
        for date_unit in date_units:
            try:
                new_data = to_datetime(new_data, errors='raise',
                                       unit=date_unit)
            except ValueError:
                continue
            except:
                break
            return new_data, True
        return data, False

    def _try_convert_dates(self):
        raise AbstractMethodError(self)


class SeriesParser(Parser):
    _default_orient = 'index'
    _split_keys = ('name', 'index', 'data')

    def _parse_no_numpy(self):

        json = self.json
        orient = self.orient
        if orient == "split":
            decoded = dict((str(k), v)
                           for k, v in compat.iteritems(loads(
                               json,
                               precise_float=self.precise_float)))
            self.check_keys_split(decoded)
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
            decoded = dict((str(k), v) for k, v in compat.iteritems(decoded))
            self.check_keys_split(decoded)
            self.obj = Series(**decoded)
        elif orient == "columns" or orient == "index":
            self.obj = Series(*loads(json, dtype=None, numpy=True,
                                     labelled=True,
                                     precise_float=self.precise_float))
        else:
            self.obj = Series(loads(json, dtype=None, numpy=True,
                                    precise_float=self.precise_float))

    def _try_convert_types(self):
        if self.obj is None:
            return
        obj, result = self._try_convert_data(
            'data', self.obj, convert_dates=self.convert_dates)
        if result:
            self.obj = obj


class FrameParser(Parser):
    _default_orient = 'columns'
    _split_keys = ('columns', 'index', 'data')

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
            decoded = dict((str(k), v) for k, v in compat.iteritems(decoded))
            self.check_keys_split(decoded)
            self.obj = DataFrame(**decoded)
        elif orient == "values":
            self.obj = DataFrame(loads(json, dtype=None, numpy=True,
                                       precise_float=self.precise_float))
        else:
            self.obj = DataFrame(*loads(json, dtype=None, numpy=True,
                                        labelled=True,
                                        precise_float=self.precise_float))

    def _parse_no_numpy(self):

        json = self.json
        orient = self.orient

        if orient == "columns":
            self.obj = DataFrame(
                loads(json, precise_float=self.precise_float), dtype=None)
        elif orient == "split":
            decoded = dict((str(k), v)
                           for k, v in compat.iteritems(loads(
                               json,
                               precise_float=self.precise_float)))
            self.check_keys_split(decoded)
            self.obj = DataFrame(dtype=None, **decoded)
        elif orient == "index":
            self.obj = DataFrame(
                loads(json, precise_float=self.precise_float), dtype=None).T
        else:
            self.obj = DataFrame(
                loads(json, precise_float=self.precise_float), dtype=None)

    def _process_converter(self, f, filt=None):
        """ take a conversion function and possibly recreate the frame """

        if filt is None:
            filt = lambda col, c: True

        needs_new_obj = False
        new_obj = dict()
        for i, (col, c) in enumerate(self.obj.iteritems()):
            if filt(col, c):
                new_data, result = f(col, c)
                if result:
                    c = new_data
                    needs_new_obj = True
            new_obj[i] = c

        if needs_new_obj:

            # possibly handle dup columns
            new_obj = DataFrame(new_obj, index=self.obj.index)
            new_obj.columns = self.obj.columns
            self.obj = new_obj

    def _try_convert_types(self):
        if self.obj is None:
            return
        if self.convert_dates:
            self._try_convert_dates()

        self._process_converter(
            lambda col, c: self._try_convert_data(col, c, convert_dates=False))

    def _try_convert_dates(self):
        if self.obj is None:
            return

        # our columns to parse
        convert_dates = self.convert_dates
        if convert_dates is True:
            convert_dates = []
        convert_dates = set(convert_dates)

        def is_ok(col):
            """ return if this col is ok to try for a date parse """
            if not isinstance(col, compat.string_types):
                return False

            col_lower = col.lower()
            if (col_lower.endswith('_at') or
                    col_lower.endswith('_time') or
                    col_lower == 'modified' or
                    col_lower == 'date' or
                    col_lower == 'datetime' or
                    col_lower.startswith('timestamp')):
                return True
            return False

        self._process_converter(
            lambda col, c: self._try_convert_to_date(c),
            lambda col, c: ((self.keep_default_dates and is_ok(col)) or
                            col in convert_dates))

# ---------------------------------------------------------------------
# JSON normalization routines


def nested_to_record(ds, prefix="", level=0):
    """a simplified json_normalize

    converts a nested dict into a flat dict ("record"), unlike json_normalize,
    it does not attempt to extract a subset of the data.

    Parameters
    ----------
    ds : dict or list of dicts
    prefix: the prefix, optional, default: ""
    level: the number of levels in the jason string, optional, default: 0

    Returns
    -------
    d - dict or list of dicts, matching `ds`

    Examples
    --------

    IN[52]: nested_to_record(dict(flat1=1,dict1=dict(c=1,d=2),
                                  nested=dict(e=dict(c=1,d=2),d=2)))
    Out[52]:
    {'dict1.c': 1,
     'dict1.d': 2,
     'flat1': 1,
     'nested.d': 2,
     'nested.e.c': 1,
     'nested.e.d': 2}
    """
    singleton = False
    if isinstance(ds, dict):
        ds = [ds]
        singleton = True

    new_ds = []
    for d in ds:

        new_d = copy.deepcopy(d)
        for k, v in d.items():
            # each key gets renamed with prefix
            if level == 0:
                newkey = str(k)
            else:
                newkey = prefix + '.' + str(k)

            # only dicts gets recurse-flattend
            # only at level>1 do we rename the rest of the keys
            if not isinstance(v, dict):
                if level != 0:  # so we skip copying for top level, common case
                    v = new_d.pop(k)
                    new_d[newkey] = v
                continue
            else:
                v = new_d.pop(k)
                new_d.update(nested_to_record(v, newkey, level + 1))
        new_ds.append(new_d)

    if singleton:
        return new_ds[0]
    return new_ds


def json_normalize(data, record_path=None, meta=None,
                   meta_prefix=None,
                   record_prefix=None):
    """
    "Normalize" semi-structured JSON data into a flat table

    Parameters
    ----------
    data : dict or list of dicts
        Unserialized JSON objects
    record_path : string or list of strings, default None
        Path in each object to list of records. If not passed, data will be
        assumed to be an array of records
    meta : list of paths (string or list of strings), default None
        Fields to use as metadata for each record in resulting table
    record_prefix : string, default None
        If True, prefix records with dotted (?) path, e.g. foo.bar.field if
        path to records is ['foo', 'bar']
    meta_prefix : string, default None

    Returns
    -------
    frame : DataFrame

    Examples
    --------

    >>> data = [{'state': 'Florida',
    ...          'shortname': 'FL',
    ...          'info': {
    ...               'governor': 'Rick Scott'
    ...          },
    ...          'counties': [{'name': 'Dade', 'population': 12345},
    ...                      {'name': 'Broward', 'population': 40000},
    ...                      {'name': 'Palm Beach', 'population': 60000}]},
    ...         {'state': 'Ohio',
    ...          'shortname': 'OH',
    ...          'info': {
    ...               'governor': 'John Kasich'
    ...          },
    ...          'counties': [{'name': 'Summit', 'population': 1234},
    ...                       {'name': 'Cuyahoga', 'population': 1337}]}]
    >>> from pandas.io.json import json_normalize
    >>> result = json_normalize(data, 'counties', ['state', 'shortname',
    ...                                           ['info', 'governor']])
    >>> result
             name  population info.governor    state shortname
    0        Dade       12345    Rick Scott  Florida        FL
    1     Broward       40000    Rick Scott  Florida        FL
    2  Palm Beach       60000    Rick Scott  Florida        FL
    3      Summit        1234   John Kasich     Ohio        OH
    4    Cuyahoga        1337   John Kasich     Ohio        OH

    """
    def _pull_field(js, spec):
        result = js
        if isinstance(spec, list):
            for field in spec:
                result = result[field]
        else:
            result = result[spec]

        return result

    # A bit of a hackjob
    if isinstance(data, dict):
        data = [data]

    if record_path is None:
        if any([isinstance(x, dict) for x in compat.itervalues(data[0])]):
            # naive normalization, this is idempotent for flat records
            # and potentially will inflate the data considerably for
            # deeply nested structures:
            #  {VeryLong: { b: 1,c:2}} -> {VeryLong.b:1 ,VeryLong.c:@}
            #
            # TODO: handle record value which are lists, at least error
            #       reasonably
            data = nested_to_record(data)
        return DataFrame(data)
    elif not isinstance(record_path, list):
        record_path = [record_path]

    if meta is None:
        meta = []
    elif not isinstance(meta, list):
        meta = [meta]

    for i, x in enumerate(meta):
        if not isinstance(x, list):
            meta[i] = [x]

    # Disastrously inefficient for now
    records = []
    lengths = []

    meta_vals = defaultdict(list)
    meta_keys = ['.'.join(val) for val in meta]

    def _recursive_extract(data, path, seen_meta, level=0):
        if len(path) > 1:
            for obj in data:
                for val, key in zip(meta, meta_keys):
                    if level + 1 == len(val):
                        seen_meta[key] = _pull_field(obj, val[-1])

                _recursive_extract(obj[path[0]], path[1:],
                                   seen_meta, level=level + 1)
        else:
            for obj in data:
                recs = _pull_field(obj, path[0])

                # For repeating the metadata later
                lengths.append(len(recs))

                for val, key in zip(meta, meta_keys):
                    if level + 1 > len(val):
                        meta_val = seen_meta[key]
                    else:
                        meta_val = _pull_field(obj, val[level:])
                    meta_vals[key].append(meta_val)

                records.extend(recs)

    _recursive_extract(data, record_path, {}, level=0)

    result = DataFrame(records)

    if record_prefix is not None:
        result.rename(columns=lambda x: record_prefix + x, inplace=True)

    # Data types, a problem
    for k, v in compat.iteritems(meta_vals):
        if meta_prefix is not None:
            k = meta_prefix + k

        if k in result:
            raise ValueError('Conflicting metadata name %s, '
                             'need distinguishing prefix ' % k)

        result[k] = np.array(v).repeat(lengths)

    return result
