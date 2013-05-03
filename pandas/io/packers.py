"""
Msgpack serializer support for reading and writing pandas data structures
to disk
"""

# porfions of msgpack_numpy package, by Lev Givon were incorporated
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

from datetime import datetime
import time
import re
import copy
import itertools
import warnings

import numpy as np
from pandas import (
    Timestamp, Period, Series, TimeSeries, DataFrame, Panel, Panel4D,
    Index, MultiIndex, Int64Index, PeriodIndex, DatetimeIndex, NaT
)
from pandas.sparse.api import SparseSeries, SparseDataFrame, SparsePanel
from pandas.sparse.array import BlockIndex, IntIndex
from pandas.tseries.api import PeriodIndex, DatetimeIndex
from pandas.core.index import Int64Index, _ensure_index
import pandas.core.common as com
from pandas.core.common import needs_i8_conversion
from pandas.core.internals import BlockManager, make_block
import pandas.core.internals as internals

try:
    import msgpack
    from msgpack import _packer, _unpacker
    _USE_MSGPACK = True
except:
    _USE_MSGPACK = False

def to_msgpack(path, *args, **kwargs):
    """
    msgpack (serialize) object to input file path

    Parameters
    ----------
    path : string
        File path
    args : an object or objects to serialize

    append : boolean whether to append to an existing msgpack
             (default is False)
    """
    if not _USE_MSGPACK:
        raise Exception("please install msgpack to create msgpack stores!")

    append = kwargs.get('append')
    if append:
        f = open(path, 'a+b')
    else:
        f = open(path, 'wb')
    try:
        if len(args) == 1:
            f.write(pack(args[0]))
        else:
            for a in args:
                f.write(pack(a))
    finally:
        f.close()


def read_msgpack(path, iterator=False, **kwargs):
    """
    Load msgpack pandas object from the specified
    file path

    Parameters
    ----------
    path : string
        File path
    iterator : boolean, if True, return an iterator to the unpacker
               (default is False)

    Returns
    -------
    obj : type of object stored in file

    """
    if not _USE_MSGPACK:
        raise Exception("please install msgpack to read msgpack stores!")
    if iterator:
        return Iterator(path)

    with open(path,'rb') as fh:
        l = list(unpack(fh))
        if len(l) == 1:
            return l[0]
        return l

dtype_dict = { 'datetime64[ns]'  : np.dtype('M8[ns]'),
               'timedelta64[ns]' : np.dtype('m8[ns]') }

def dtype_for(t):
    if t in dtype_dict:
        return dtype_dict[t]
    return np.typeDict[t]

c2f_dict = {'complex':    np.float64,
            'complex128': np.float64,
            'complex256': np.float128,
            'complex64':  np.float32}

def c2f(r, i, ctype_name):
    """
    Convert strings to complex number instance with specified numpy type.
    """
    
    ftype = c2f_dict[ctype_name]
    return np.typeDict[ctype_name](ftype(r)+1j*ftype(i))

def convert(values):
    """ convert the numpy values to a list """

    dtype = values.dtype
    if needs_i8_conversion(dtype):
        values = values.view('i8')
    return values.ravel().tolist()


def encode(obj):
    """
    Data encoder
    """
        
    if isinstance(obj, Index):
        if isinstance(obj, PeriodIndex):
            return {'typ' : 'period_index',
                    'klass' : obj.__class__.__name__,
                    'name' : getattr(obj,'name',None),
                    'dtype': obj.dtype.name,
                    'data': obj.tolist() }
        elif isinstance(obj, DatetimeIndex):
            return {'typ' : 'datetime_index',
                    'klass' : obj.__class__.__name__,
                    'name' : getattr(obj,'name',None),
                    'dtype': obj.dtype.name,
                    'data': obj.values.view('i8').tolist(),
                    'freq' : obj.freqstr,
                    'tz'   : obj.tz}
        elif isinstance(obj, MultiIndex):
            return {'typ' : 'multi_index',
                    'klass' : obj.__class__.__name__,
                    'names' : getattr(obj,'names',None),
                    'dtype': obj.dtype.name,
                    'data': obj.tolist() }
        else:
            return {'typ' : 'index',
                    'klass' : obj.__class__.__name__,
                    'name' : getattr(obj,'name',None),
                    'dtype': obj.dtype.name,
                    'data': obj.tolist() }
    elif isinstance(obj, Series):
        if isinstance(obj, SparseSeries):
            import pdb; pdb.set_trace()
        else:
            return {'typ' : 'series',
                    'klass' : obj.__class__.__name__,
                    'name' : getattr(obj,'name',None),
                    'index' : obj.index,
                    'dtype': obj.dtype.name,
                    'data': convert(obj.values) }
    elif isinstance(obj, DataFrame):
        if isinstance(obj, SparseDataFrame):
            import pdb; pdb.set_trace()
        else:

            data = obj._data
            if not data.is_consolidated():
                data = data.consolidate()

           # the block manager
            return {'typ' : 'dataframe',
                    'klass'  : obj.__class__.__name__,
                    'axes'   : data.axes,
                    'blocks' : [ { 'items'  : b.items, 
                                   'values' : convert(b.values), 
                                   'shape'  : b.values.shape,
                                   'dtype'  : b.dtype.name,
                                   'klass' : b.__class__.__name__ 
                                   } for b in data.blocks ] }

    elif isinstance(obj, datetime):
        if isinstance(obj, Timestamp):
            tz = obj.tzinfo
            if tz is not None:
                tz = tz.zone
            offset = obj.offset
            if offset is not None:
                offset = offset.freqstr
            return {'typ' : 'timestamp',
                    'value': obj.value,
                    'offset' : offset,
                    'tz' : tz}
        return { 'typ' : 'datetime',
                 'data' : obj.isoformat() }
    elif isinstance(obj, Period):
        return {'typ' : 'period',
                'ordinal' : obj.ordinal,
                'freq' : obj.freq }
    elif isinstance(obj, np.ndarray):
        return {'typ' : 'ndarray',
                'shape': obj.shape,
                'ndim': obj.ndim,
                'dtype': obj.dtype.name,
                'data': convert(obj)}
    elif isinstance(obj, np.number):
        if np.iscomplexobj(obj):
            return {'typ' : 'np_scalar',
                    'sub_typ' : 'np_complex',
                    'dtype': obj.dtype.name,
                    'real': obj.real.__repr__(),
                    'imag': obj.imag.__repr__()}
        else:
            return {'typ' : 'np_scalar',
                    'dtype': obj.dtype.name,
                    'data': obj.__repr__()}
    elif isinstance(obj, complex):
        return {'typ' : 'np_complex',
                'real': obj.real.__repr__(),
                'imag': obj.imag.__repr__()}
    else:
        import pdb; pdb.set_trace()
        return obj

def decode(obj):
    """
    Decoder for deserializing numpy data types.
    """
    
    typ = obj.get('typ')
    if typ is None:
        return obj
    elif typ == 'timestamp':
        return Timestamp(obj['value'],tz=obj['tz'],offset=obj['offset'])
    elif typ == 'period':
        return Period(ordinal=obj['ordinal'],freq=obj['freq'])
    elif typ == 'index':
        dtype = dtype_for(obj['dtype'])
        data = obj['data']
        return globals()[obj['klass']](data,dtype=dtype,name=obj['name'])
    elif typ == 'multi_index':
        return globals()[obj['klass']].from_tuples(obj['data'],names=obj['names'])
    elif typ == 'period_index':
        return globals()[obj['klass']](obj['data'],name=obj['name'])
    elif typ == 'datetime_index':
        return globals()[obj['klass']](obj['data'],freq=obj['freq'],tz=obj['tz'],name=obj['name'])
    elif typ == 'series':
        dtype = dtype_for(obj['dtype'])
        index = obj['index']
        return globals()[obj['klass']](obj['data'],index=index,dtype=dtype,name=obj['name'])
    elif typ == 'dataframe':
        axes = obj['axes']

        def create_block(b):
            dtype = dtype_for(b['dtype'])
            return make_block(np.array(b['values'],dtype=dtype).reshape(b['shape']),b['items'],axes[0],klass=getattr(internals,b['klass'])) 

        blocks = [ create_block(b) for b in obj['blocks'] ]
        return globals()[obj['klass']](BlockManager(blocks, axes))
    elif typ == 'datetime':
        import pdb; pdb.set_trace()
        return datetime.fromtimestamp(obj['data'])
    elif typ == 'ndarray':
        return np.array(obj['data'],
                        dtype=np.typeDict[obj['dtype']],
                        ndmin=obj['ndim']).reshape(obj['shape'])
    elif typ == 'np_scalar':
        if obj.get('sub_typ') == 'np_complex':
            return c2f(obj['real'], obj['imag'], obj['dtype'])
        else:
            return np.typeDict[obj['dtype']](obj['data'])
    elif typ == 'np_complex':
        return complex(obj['real']+'+'+obj['imag']+'j')
    elif isinstance(obj, (dict,list,set)):
        return obj
    else:
        import pdb; pdb.set_trace()
        return obj

def pack(o, default=encode, 
         encoding='utf-8', unicode_errors='strict', use_single_float=False):
    """
    Pack an object and return the packed bytes.
    """

    return Packer(default=default, encoding=encoding,
           unicode_errors=unicode_errors, 
           use_single_float=use_single_float).pack(o)

def unpack(packed, object_hook=decode, 
           list_hook=None, use_list=False, encoding='utf-8',
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

if _USE_MSGPACK:

    class Packer(_packer.Packer):
        def __init__(self, default=encode, 
                     encoding='utf-8',
                     unicode_errors='strict',
                     use_single_float=False):
            super(Packer, self).__init__(default=default, 
                                         encoding=encoding,
                                         unicode_errors=unicode_errors,
                                         use_single_float=use_single_float)

    class Unpacker(_unpacker.Unpacker):
        def __init__(self, file_like=None, read_size=0, use_list=False,
                     object_hook=decode,
                     object_pairs_hook=None, list_hook=None, encoding='utf-8',
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

        try:
            fh   = open(self.path,'rb')
            unpacker = unpack(fh)
            for o in unpacker:
                yield o
        finally:
            fh.close()
