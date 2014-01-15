"""
Msgpack serializer support for reading and writing pandas data structures
to disk
"""

# portions of msgpack_numpy package, by Lev Givon were incorporated
# into this module (and tests_packers.py)

"""
License
=======

Copyright (c) 2013, Lev Givon.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above
  copyright notice, this list of conditions and the following
  disclaimer in the documentation and/or other materials provided
  with the distribution.
* Neither the name of Lev Givon nor the names of any
  contributors may be used to endorse or promote products derived
  from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import os
from datetime import datetime, date, timedelta
from dateutil.parser import parse

import numpy as np
from pandas import compat
from pandas.compat import u, PY3
from pandas import (
    Timestamp, Period, Series, DataFrame, Panel, Panel4D,
    Index, MultiIndex, Int64Index, PeriodIndex, DatetimeIndex, Float64Index,
    NaT
)
from pandas.sparse.api import SparseSeries, SparseDataFrame, SparsePanel
from pandas.sparse.array import BlockIndex, IntIndex
from pandas.core.generic import NDFrame
from pandas.core.common import needs_i8_conversion
from pandas.io.common import get_filepath_or_buffer
from pandas.core.internals import BlockManager, make_block
import pandas.core.internals as internals

from pandas.msgpack import Unpacker as _Unpacker, Packer as _Packer
import zlib

try:
    import blosc
    _BLOSC = True
except:
    _BLOSC = False

# until we can pass this into our conversion functions,
# this is pretty hacky
compressor = None


def to_msgpack(path_or_buf, *args, **kwargs):
    """
    msgpack (serialize) object to input file path

    THIS IS AN EXPERIMENTAL LIBRARY and the storage format
    may not be stable until a future release.

    Parameters
    ----------
    path_or_buf : string File path, buffer-like, or None
                  if None, return generated string
    args : an object or objects to serialize
    append : boolean whether to append to an existing msgpack
             (default is False)
    compress : type of compressor (zlib or blosc), default to None (no
               compression)
    """
    global compressor
    compressor = kwargs.pop('compress', None)
    append = kwargs.pop('append', None)
    if append:
        mode = 'a+b'
    else:
        mode = 'wb'

    def writer(fh):
        for a in args:
            fh.write(pack(a, **kwargs))

    if isinstance(path_or_buf, compat.string_types):
        with open(path_or_buf, mode) as fh:
            writer(fh)
    elif path_or_buf is None:
        buf = compat.BytesIO()
        writer(buf)
        return buf.getvalue()
    else:
        writer(path_or_buf)


def read_msgpack(path_or_buf, iterator=False, **kwargs):
    """
    Load msgpack pandas object from the specified
    file path

    THIS IS AN EXPERIMENTAL LIBRARY and the storage format
    may not be stable until a future release.

    Parameters
    ----------
    path_or_buf : string File path, BytesIO like or string
    iterator : boolean, if True, return an iterator to the unpacker
               (default is False)

    Returns
    -------
    obj : type of object stored in file

    """
    path_or_buf, _ = get_filepath_or_buffer(path_or_buf)
    if iterator:
        return Iterator(path_or_buf)

    def read(fh):
        l = list(unpack(fh))
        if len(l) == 1:
            return l[0]
        return l

    # see if we have an actual file
    if isinstance(path_or_buf, compat.string_types):

        try:
            exists = os.path.exists(path_or_buf)
        except (TypeError,ValueError):
            exists = False

        if exists:
            with open(path_or_buf, 'rb') as fh:
                return read(fh)

    # treat as a string-like
    if not hasattr(path_or_buf, 'read'):

        try:
            fh = compat.BytesIO(path_or_buf)
            return read(fh)
        finally:
            fh.close()

    # a buffer like
    return read(path_or_buf)

dtype_dict = {21: np.dtype('M8[ns]'),
              u('datetime64[ns]'): np.dtype('M8[ns]'),
              u('datetime64[us]'): np.dtype('M8[us]'),
              22: np.dtype('m8[ns]'),
              u('timedelta64[ns]'): np.dtype('m8[ns]'),
              u('timedelta64[us]'): np.dtype('m8[us]')}


def dtype_for(t):
    if t in dtype_dict:
        return dtype_dict[t]
    return np.typeDict[t]

c2f_dict = {'complex':    np.float64,
            'complex128': np.float64,
            'complex64':  np.float32}

# numpy 1.6.1 compat
if hasattr(np, 'float128'):
    c2f_dict['complex256'] = np.float128


def c2f(r, i, ctype_name):
    """
    Convert strings to complex number instance with specified numpy type.
    """

    ftype = c2f_dict[ctype_name]
    return np.typeDict[ctype_name](ftype(r) + 1j * ftype(i))


def convert(values):
    """ convert the numpy values to a list """

    dtype = values.dtype
    if needs_i8_conversion(dtype):
        values = values.view('i8')
    v = values.ravel()

    # convert object
    if dtype == np.object_:
        return v.tolist()

    if compressor == 'zlib':

        # return string arrays like they are
        if dtype == np.object_:
            return v.tolist()

        # convert to a bytes array
        v = v.tostring()
        return zlib.compress(v)

    elif compressor == 'blosc' and _BLOSC:

        # return string arrays like they are
        if dtype == np.object_:
            return v.tolist()

        # convert to a bytes array
        v = v.tostring()
        return blosc.compress(v, typesize=dtype.itemsize)

    # ndarray (on original dtype)
    return v.tostring()


def unconvert(values, dtype, compress=None):

    if dtype == np.object_:
        return np.array(values, dtype=object)

    if compress == 'zlib':

        values = zlib.decompress(values)
        return np.frombuffer(values, dtype=dtype)

    elif compress == 'blosc':

        if not _BLOSC:
            raise Exception("cannot uncompress w/o blosc")

        # decompress
        values = blosc.decompress(values)

        return np.frombuffer(values, dtype=dtype)

    # from a string
    return np.fromstring(values.encode('latin1'), dtype=dtype)


def encode(obj):
    """
    Data encoder
    """

    tobj = type(obj)
    if isinstance(obj, Index):
        if isinstance(obj, PeriodIndex):
            return {'typ': 'period_index',
                    'klass': obj.__class__.__name__,
                    'name': getattr(obj, 'name', None),
                    'freq': getattr(obj, 'freqstr', None),
                    'dtype': obj.dtype.num,
                    'data': convert(obj.asi8)}
        elif isinstance(obj, DatetimeIndex):
            tz = getattr(obj, 'tz', None)

            # store tz info and data as UTC
            if tz is not None:
                tz = tz.zone
                obj = obj.tz_convert('UTC')
            return {'typ': 'datetime_index',
                    'klass': obj.__class__.__name__,
                    'name': getattr(obj, 'name', None),
                    'dtype': obj.dtype.num,
                    'data': convert(obj.asi8),
                    'freq': getattr(obj, 'freqstr', None),
                    'tz': tz}
        elif isinstance(obj, MultiIndex):
            return {'typ': 'multi_index',
                    'klass': obj.__class__.__name__,
                    'names': getattr(obj, 'names', None),
                    'dtype': obj.dtype.num,
                    'data': convert(obj.values)}
        else:
            return {'typ': 'index',
                    'klass': obj.__class__.__name__,
                    'name': getattr(obj, 'name', None),
                    'dtype': obj.dtype.num,
                    'data': convert(obj.values)}
    elif isinstance(obj, Series):
        if isinstance(obj, SparseSeries):
            raise NotImplementedError(
                'msgpack sparse series is not implemented'
            )
            #d = {'typ': 'sparse_series',
            #     'klass': obj.__class__.__name__,
            #     'dtype': obj.dtype.num,
            #     'index': obj.index,
            #     'sp_index': obj.sp_index,
            #     'sp_values': convert(obj.sp_values),
            #     'compress': compressor}
            #for f in ['name', 'fill_value', 'kind']:
            #    d[f] = getattr(obj, f, None)
            #return d
        else:
            return {'typ': 'series',
                    'klass': obj.__class__.__name__,
                    'name': getattr(obj, 'name', None),
                    'index': obj.index,
                    'dtype': obj.dtype.num,
                    'data': convert(obj.values),
                    'compress': compressor}
    elif issubclass(tobj, NDFrame):
        if isinstance(obj, SparseDataFrame):
            raise NotImplementedError(
                'msgpack sparse frame is not implemented'
            )
            #d = {'typ': 'sparse_dataframe',
            #     'klass': obj.__class__.__name__,
            #     'columns': obj.columns}
            #for f in ['default_fill_value', 'default_kind']:
            #    d[f] = getattr(obj, f, None)
            #d['data'] = dict([(name, ss)
            #                 for name, ss in compat.iteritems(obj)])
            #return d
        elif isinstance(obj, SparsePanel):
            raise NotImplementedError(
                'msgpack sparse frame is not implemented'
            )
            #d = {'typ': 'sparse_panel',
            #     'klass': obj.__class__.__name__,
            #     'items': obj.items}
            #for f in ['default_fill_value', 'default_kind']:
            #    d[f] = getattr(obj, f, None)
            #d['data'] = dict([(name, df)
            #                 for name, df in compat.iteritems(obj)])
            #return d
        else:

            data = obj._data
            if not data.is_consolidated():
                data = data.consolidate()

           # the block manager
            return {'typ': 'block_manager',
                    'klass': obj.__class__.__name__,
                    'axes': data.axes,
                    'blocks': [{'items': b.items,
                                'values': convert(b.values),
                                'shape': b.values.shape,
                                'dtype': b.dtype.num,
                                'klass': b.__class__.__name__,
                                'compress': compressor
                                } for b in data.blocks]}

    elif isinstance(obj, (datetime, date, np.datetime64, timedelta,
                          np.timedelta64)):
        if isinstance(obj, Timestamp):
            tz = obj.tzinfo
            if tz is not None:
                tz = tz.zone
            offset = obj.offset
            if offset is not None:
                offset = offset.freqstr
            return {'typ': 'timestamp',
                    'value': obj.value,
                    'offset': offset,
                    'tz': tz}
        elif isinstance(obj, np.timedelta64):
            return {'typ': 'timedelta64',
                    'data': obj.view('i8')}
        elif isinstance(obj, timedelta):
            return {'typ': 'timedelta',
                    'data': (obj.days, obj.seconds, obj.microseconds)}
        elif isinstance(obj, np.datetime64):
            return {'typ': 'datetime64',
                    'data': str(obj)}
        elif isinstance(obj, datetime):
            return {'typ': 'datetime',
                    'data': obj.isoformat()}
        elif isinstance(obj, date):
            return {'typ': 'date',
                    'data': obj.isoformat()}
        raise Exception("cannot encode this datetimelike object: %s" % obj)
    elif isinstance(obj, Period):
        return {'typ': 'period',
                'ordinal': obj.ordinal,
                'freq': obj.freq}
    elif isinstance(obj, BlockIndex):
        return {'typ': 'block_index',
                'klass': obj.__class__.__name__,
                'blocs': obj.blocs,
                'blengths': obj.blengths,
                'length': obj.length}
    elif isinstance(obj, IntIndex):
        return {'typ': 'int_index',
                'klass': obj.__class__.__name__,
                'indices': obj.indices,
                'length': obj.length}
    elif isinstance(obj, np.ndarray):
        return {'typ': 'ndarray',
                'shape': obj.shape,
                'ndim': obj.ndim,
                'dtype': obj.dtype.num,
                'data': convert(obj),
                'compress': compressor}
    elif isinstance(obj, np.number):
        if np.iscomplexobj(obj):
            return {'typ': 'np_scalar',
                    'sub_typ': 'np_complex',
                    'dtype': obj.dtype.name,
                    'real': obj.real.__repr__(),
                    'imag': obj.imag.__repr__()}
        else:
            return {'typ': 'np_scalar',
                    'dtype': obj.dtype.name,
                    'data': obj.__repr__()}
    elif isinstance(obj, complex):
        return {'typ': 'np_complex',
                'real': obj.real.__repr__(),
                'imag': obj.imag.__repr__()}

    return obj


def decode(obj):
    """
    Decoder for deserializing numpy data types.
    """

    typ = obj.get('typ')
    if typ is None:
        return obj
    elif typ == 'timestamp':
        return Timestamp(obj['value'], tz=obj['tz'], offset=obj['offset'])
    elif typ == 'period':
        return Period(ordinal=obj['ordinal'], freq=obj['freq'])
    elif typ == 'index':
        dtype = dtype_for(obj['dtype'])
        data = unconvert(obj['data'], np.typeDict[obj['dtype']],
                         obj.get('compress'))
        return globals()[obj['klass']](data, dtype=dtype, name=obj['name'])
    elif typ == 'multi_index':
        data = unconvert(obj['data'], np.typeDict[obj['dtype']],
                         obj.get('compress'))
        data = [tuple(x) for x in data]
        return globals()[obj['klass']].from_tuples(data, names=obj['names'])
    elif typ == 'period_index':
        data = unconvert(obj['data'], np.int64, obj.get('compress'))
        d = dict(name=obj['name'], freq=obj['freq'])
        return globals()[obj['klass']](data, **d)
    elif typ == 'datetime_index':
        data = unconvert(obj['data'], np.int64, obj.get('compress'))
        d = dict(name=obj['name'], freq=obj['freq'], verify_integrity=False)
        result = globals()[obj['klass']](data, **d)
        tz = obj['tz']

        # reverse tz conversion
        if tz is not None:
            result = result.tz_localize('UTC').tz_convert(tz)
        return result

    elif typ == 'series':
        dtype = dtype_for(obj['dtype'])
        index = obj['index']
        return globals()[obj['klass']](unconvert(obj['data'], dtype,
                                                 obj['compress']),
                                       index=index, name=obj['name'])
    elif typ == 'block_manager':
        axes = obj['axes']

        def create_block(b):
            dtype = dtype_for(b['dtype'])
            return make_block(unconvert(b['values'], dtype, b['compress'])
                              .reshape(b['shape']), b['items'], axes[0],
                              klass=getattr(internals, b['klass']))

        blocks = [create_block(b) for b in obj['blocks']]
        return globals()[obj['klass']](BlockManager(blocks, axes))
    elif typ == 'datetime':
        return parse(obj['data'])
    elif typ == 'datetime64':
        return np.datetime64(parse(obj['data']))
    elif typ == 'date':
        return parse(obj['data']).date()
    elif typ == 'timedelta':
        return timedelta(*obj['data'])
    elif typ == 'timedelta64':
        return np.timedelta64(int(obj['data']))
    #elif typ == 'sparse_series':
    #    dtype = dtype_for(obj['dtype'])
    #    return globals()[obj['klass']](
    #        unconvert(obj['sp_values'], dtype, obj['compress']),
    #        sparse_index=obj['sp_index'], index=obj['index'],
    #        fill_value=obj['fill_value'], kind=obj['kind'], name=obj['name'])
    #elif typ == 'sparse_dataframe':
    #    return globals()[obj['klass']](
    #        obj['data'], columns=obj['columns'],
    #        default_fill_value=obj['default_fill_value'],
    #        default_kind=obj['default_kind']
    #    )
    #elif typ == 'sparse_panel':
    #    return globals()[obj['klass']](
    #        obj['data'], items=obj['items'],
    #        default_fill_value=obj['default_fill_value'],
    #        default_kind=obj['default_kind'])
    elif typ == 'block_index':
        return globals()[obj['klass']](obj['length'], obj['blocs'],
                                       obj['blengths'])
    elif typ == 'int_index':
        return globals()[obj['klass']](obj['length'], obj['indices'])
    elif typ == 'ndarray':
        return unconvert(obj['data'], np.typeDict[obj['dtype']],
                         obj.get('compress')).reshape(obj['shape'])
    elif typ == 'np_scalar':
        if obj.get('sub_typ') == 'np_complex':
            return c2f(obj['real'], obj['imag'], obj['dtype'])
        else:
            dtype = dtype_for(obj['dtype'])
            try:
                return dtype(obj['data'])
            except:
                return dtype.type(obj['data'])
    elif typ == 'np_complex':
        return complex(obj['real'] + '+' + obj['imag'] + 'j')
    elif isinstance(obj, (dict, list, set)):
        return obj
    else:
        return obj


def pack(o, default=encode,
         encoding='latin1', unicode_errors='strict', use_single_float=False):
    """
    Pack an object and return the packed bytes.
    """

    return Packer(default=default, encoding=encoding,
                  unicode_errors=unicode_errors,
                  use_single_float=use_single_float).pack(o)


def unpack(packed, object_hook=decode,
           list_hook=None, use_list=False, encoding='latin1',
           unicode_errors='strict', object_pairs_hook=None):
    """
    Unpack a packed object, return an iterator
    Note: packed lists will be returned as tuples
    """

    return Unpacker(packed, object_hook=object_hook,
                    list_hook=list_hook,
                    use_list=use_list, encoding=encoding,
                    unicode_errors=unicode_errors,
                    object_pairs_hook=object_pairs_hook)


class Packer(_Packer):

    def __init__(self, default=encode,
                 encoding='latin1',
                 unicode_errors='strict',
                 use_single_float=False):
        super(Packer, self).__init__(default=default,
                                     encoding=encoding,
                                     unicode_errors=unicode_errors,
                                     use_single_float=use_single_float)


class Unpacker(_Unpacker):

    def __init__(self, file_like=None, read_size=0, use_list=False,
                 object_hook=decode,
                 object_pairs_hook=None, list_hook=None, encoding='latin1',
                 unicode_errors='strict', max_buffer_size=0):
        super(Unpacker, self).__init__(file_like=file_like,
                                       read_size=read_size,
                                       use_list=use_list,
                                       object_hook=object_hook,
                                       object_pairs_hook=object_pairs_hook,
                                       list_hook=list_hook,
                                       encoding=encoding,
                                       unicode_errors=unicode_errors,
                                       max_buffer_size=max_buffer_size)


class Iterator(object):

    """ manage the unpacking iteration,
        close the file on completion """

    def __init__(self, path, **kwargs):
        self.path = path
        self.kwargs = kwargs

    def __iter__(self):

        needs_closing = True
        try:

            # see if we have an actual file
            if isinstance(self.path, compat.string_types):

                try:
                    path_exists = os.path.exists(self.path)
                except TypeError:
                    path_exists = False

                if path_exists:
                    fh = open(self.path, 'rb')
                else:
                    fh = compat.BytesIO(self.path)

            else:

                if not hasattr(self.path, 'read'):
                    fh = compat.BytesIO(self.path)

                else:

                    # a file-like
                    needs_closing = False
                    fh = self.path

            unpacker = unpack(fh)
            for o in unpacker:
                yield o
        finally:
            if needs_closing:
                fh.close()
