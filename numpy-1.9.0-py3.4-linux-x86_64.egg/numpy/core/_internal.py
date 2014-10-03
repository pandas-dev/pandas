"""
A place for code to be called from core C-code.

Some things are more easily handled Python.

"""
from __future__ import division, absolute_import, print_function

import re
import sys
import warnings

from numpy.compat import asbytes, bytes

if (sys.byteorder == 'little'):
    _nbo = asbytes('<')
else:
    _nbo = asbytes('>')

def _makenames_list(adict, align):
    from .multiarray import dtype
    allfields = []
    fnames = list(adict.keys())
    for fname in fnames:
        obj = adict[fname]
        n = len(obj)
        if not isinstance(obj, tuple) or n not in [2, 3]:
            raise ValueError("entry not a 2- or 3- tuple")
        if (n > 2) and (obj[2] == fname):
            continue
        num = int(obj[1])
        if (num < 0):
            raise ValueError("invalid offset.")
        format = dtype(obj[0], align=align)
        if (format.itemsize == 0):
            raise ValueError("all itemsizes must be fixed.")
        if (n > 2):
            title = obj[2]
        else:
            title = None
        allfields.append((fname, format, num, title))
    # sort by offsets
    allfields.sort(key=lambda x: x[2])
    names = [x[0] for x in allfields]
    formats = [x[1] for x in allfields]
    offsets = [x[2] for x in allfields]
    titles = [x[3] for x in allfields]

    return names, formats, offsets, titles

# Called in PyArray_DescrConverter function when
#  a dictionary without "names" and "formats"
#  fields is used as a data-type descriptor.
def _usefields(adict, align):
    from .multiarray import dtype
    try:
        names = adict[-1]
    except KeyError:
        names = None
    if names is None:
        names, formats, offsets, titles = _makenames_list(adict, align)
    else:
        formats = []
        offsets = []
        titles = []
        for name in names:
            res = adict[name]
            formats.append(res[0])
            offsets.append(res[1])
            if (len(res) > 2):
                titles.append(res[2])
            else:
                titles.append(None)

    return dtype({"names" : names,
                  "formats" : formats,
                  "offsets" : offsets,
                  "titles" : titles}, align)


# construct an array_protocol descriptor list
#  from the fields attribute of a descriptor
# This calls itself recursively but should eventually hit
#  a descriptor that has no fields and then return
#  a simple typestring

def _array_descr(descriptor):
    fields = descriptor.fields
    if fields is None:
        subdtype = descriptor.subdtype
        if subdtype is None:
            if descriptor.metadata is None:
                return descriptor.str
            else:
                new = descriptor.metadata.copy()
                if new:
                    return (descriptor.str, new)
                else:
                    return descriptor.str
        else:
            return (_array_descr(subdtype[0]), subdtype[1])


    names = descriptor.names
    ordered_fields = [fields[x] + (x,) for x in names]
    result = []
    offset = 0
    for field in ordered_fields:
        if field[1] > offset:
            num = field[1] - offset
            result.append(('', '|V%d' % num))
            offset += num
        if len(field) > 3:
            name = (field[2], field[3])
        else:
            name = field[2]
        if field[0].subdtype:
            tup = (name, _array_descr(field[0].subdtype[0]),
                   field[0].subdtype[1])
        else:
            tup = (name, _array_descr(field[0]))
        offset += field[0].itemsize
        result.append(tup)

    return result

# Build a new array from the information in a pickle.
# Note that the name numpy.core._internal._reconstruct is embedded in
# pickles of ndarrays made with NumPy before release 1.0
# so don't remove the name here, or you'll
# break backward compatibilty.
def _reconstruct(subtype, shape, dtype):
    from .multiarray import ndarray
    return ndarray.__new__(subtype, shape, dtype)


# format_re was originally from numarray by J. Todd Miller

format_re = re.compile(asbytes(
                           r'(?P<order1>[<>|=]?)'
                           r'(?P<repeats> *[(]?[ ,0-9L]*[)]? *)'
                           r'(?P<order2>[<>|=]?)'
                           r'(?P<dtype>[A-Za-z0-9.]*(?:\[[a-zA-Z0-9,.]+\])?)'))
sep_re = re.compile(asbytes(r'\s*,\s*'))
space_re = re.compile(asbytes(r'\s+$'))

# astr is a string (perhaps comma separated)

_convorder = {asbytes('='): _nbo}

def _commastring(astr):
    startindex = 0
    result = []
    while startindex < len(astr):
        mo = format_re.match(astr, pos=startindex)
        try:
            (order1, repeats, order2, dtype) = mo.groups()
        except (TypeError, AttributeError):
            raise ValueError('format number %d of "%s" is not recognized' %
                                            (len(result)+1, astr))
        startindex = mo.end()
        # Separator or ending padding
        if startindex < len(astr):
            if space_re.match(astr, pos=startindex):
                startindex = len(astr)
            else:
                mo = sep_re.match(astr, pos=startindex)
                if not mo:
                    raise ValueError(
                            'format number %d of "%s" is not recognized' %
                                            (len(result)+1, astr))
                startindex = mo.end()

        if order2 == asbytes(''):
            order = order1
        elif order1 == asbytes(''):
            order = order2
        else:
            order1 = _convorder.get(order1, order1)
            order2 = _convorder.get(order2, order2)
            if (order1 != order2):
                raise ValueError('inconsistent byte-order specification %s and %s' % (order1, order2))
            order = order1

        if order in [asbytes('|'), asbytes('='), _nbo]:
            order = asbytes('')
        dtype = order + dtype
        if (repeats == asbytes('')):
            newitem = dtype
        else:
            newitem = (dtype, eval(repeats))
        result.append(newitem)

    return result

def _getintp_ctype():
    from .multiarray import dtype
    val = _getintp_ctype.cache
    if val is not None:
        return val
    char = dtype('p').char
    import ctypes
    if (char == 'i'):
        val = ctypes.c_int
    elif char == 'l':
        val = ctypes.c_long
    elif char == 'q':
        val = ctypes.c_longlong
    else:
        val = ctypes.c_long
    _getintp_ctype.cache = val
    return val
_getintp_ctype.cache = None

# Used for .ctypes attribute of ndarray

class _missing_ctypes(object):
    def cast(self, num, obj):
        return num

    def c_void_p(self, num):
        return num

class _ctypes(object):
    def __init__(self, array, ptr=None):
        try:
            import ctypes
            self._ctypes = ctypes
        except ImportError:
            self._ctypes = _missing_ctypes()
        self._arr = array
        self._data = ptr
        if self._arr.ndim == 0:
            self._zerod = True
        else:
            self._zerod = False

    def data_as(self, obj):
        return self._ctypes.cast(self._data, obj)

    def shape_as(self, obj):
        if self._zerod:
            return None
        return (obj*self._arr.ndim)(*self._arr.shape)

    def strides_as(self, obj):
        if self._zerod:
            return None
        return (obj*self._arr.ndim)(*self._arr.strides)

    def get_data(self):
        return self._data

    def get_shape(self):
        if self._zerod:
            return None
        return (_getintp_ctype()*self._arr.ndim)(*self._arr.shape)

    def get_strides(self):
        if self._zerod:
            return None
        return (_getintp_ctype()*self._arr.ndim)(*self._arr.strides)

    def get_as_parameter(self):
        return self._ctypes.c_void_p(self._data)

    data = property(get_data, None, doc="c-types data")
    shape = property(get_shape, None, doc="c-types shape")
    strides = property(get_strides, None, doc="c-types strides")
    _as_parameter_ = property(get_as_parameter, None, doc="_as parameter_")


# Given a datatype and an order object
#  return a new names tuple
#  with the order indicated
def _newnames(datatype, order):
    oldnames = datatype.names
    nameslist = list(oldnames)
    if isinstance(order, str):
        order = [order]
    if isinstance(order, (list, tuple)):
        for name in order:
            try:
                nameslist.remove(name)
            except ValueError:
                raise ValueError("unknown field name: %s" % (name,))
        return tuple(list(order) + nameslist)
    raise ValueError("unsupported order value: %s" % (order,))

# Given an array with fields and a sequence of field names
# construct a new array with just those fields copied over
def _index_fields(ary, fields):
    from .multiarray import empty, dtype, array
    dt = ary.dtype

    names = [name for name in fields if name in dt.names]
    formats = [dt.fields[name][0] for name in fields if name in dt.names]
    offsets = [dt.fields[name][1] for name in fields if name in dt.names]

    view_dtype = {'names':names, 'formats':formats, 'offsets':offsets, 'itemsize':dt.itemsize}
    view = ary.view(dtype=view_dtype)

    # Return a copy for now until behavior is fully deprecated
    # in favor of returning view
    copy_dtype = {'names':view_dtype['names'], 'formats':view_dtype['formats']}
    return array(view, dtype=copy_dtype, copy=True)

# Given a string containing a PEP 3118 format specifier,
# construct a Numpy dtype

_pep3118_native_map = {
    '?': '?',
    'b': 'b',
    'B': 'B',
    'h': 'h',
    'H': 'H',
    'i': 'i',
    'I': 'I',
    'l': 'l',
    'L': 'L',
    'q': 'q',
    'Q': 'Q',
    'e': 'e',
    'f': 'f',
    'd': 'd',
    'g': 'g',
    'Zf': 'F',
    'Zd': 'D',
    'Zg': 'G',
    's': 'S',
    'w': 'U',
    'O': 'O',
    'x': 'V', # padding
}
_pep3118_native_typechars = ''.join(_pep3118_native_map.keys())

_pep3118_standard_map = {
    '?': '?',
    'b': 'b',
    'B': 'B',
    'h': 'i2',
    'H': 'u2',
    'i': 'i4',
    'I': 'u4',
    'l': 'i4',
    'L': 'u4',
    'q': 'i8',
    'Q': 'u8',
    'e': 'f2',
    'f': 'f',
    'd': 'd',
    'Zf': 'F',
    'Zd': 'D',
    's': 'S',
    'w': 'U',
    'O': 'O',
    'x': 'V', # padding
}
_pep3118_standard_typechars = ''.join(_pep3118_standard_map.keys())

def _dtype_from_pep3118(spec, byteorder='@', is_subdtype=False):
    from numpy.core.multiarray import dtype

    fields = {}
    offset = 0
    explicit_name = False
    this_explicit_name = False
    common_alignment = 1
    is_padding = False
    last_offset = 0

    dummy_name_index = [0]
    def next_dummy_name():
        dummy_name_index[0] += 1
    def get_dummy_name():
        while True:
            name = 'f%d' % dummy_name_index[0]
            if name not in fields:
                return name
            next_dummy_name()

    # Parse spec
    while spec:
        value = None

        # End of structure, bail out to upper level
        if spec[0] == '}':
            spec = spec[1:]
            break

        # Sub-arrays (1)
        shape = None
        if spec[0] == '(':
            j = spec.index(')')
            shape = tuple(map(int, spec[1:j].split(',')))
            spec = spec[j+1:]

        # Byte order
        if spec[0] in ('@', '=', '<', '>', '^', '!'):
            byteorder = spec[0]
            if byteorder == '!':
                byteorder = '>'
            spec = spec[1:]

        # Byte order characters also control native vs. standard type sizes
        if byteorder in ('@', '^'):
            type_map = _pep3118_native_map
            type_map_chars = _pep3118_native_typechars
        else:
            type_map = _pep3118_standard_map
            type_map_chars = _pep3118_standard_typechars

        # Item sizes
        itemsize = 1
        if spec[0].isdigit():
            j = 1
            for j in range(1, len(spec)):
                if not spec[j].isdigit():
                    break
            itemsize = int(spec[:j])
            spec = spec[j:]

        # Data types
        is_padding = False

        if spec[:2] == 'T{':
            value, spec, align, next_byteorder = _dtype_from_pep3118(
                spec[2:], byteorder=byteorder, is_subdtype=True)
        elif spec[0] in type_map_chars:
            next_byteorder = byteorder
            if spec[0] == 'Z':
                j = 2
            else:
                j = 1
            typechar = spec[:j]
            spec = spec[j:]
            is_padding = (typechar == 'x')
            dtypechar = type_map[typechar]
            if dtypechar in 'USV':
                dtypechar += '%d' % itemsize
                itemsize = 1
            numpy_byteorder = {'@': '=', '^': '='}.get(byteorder, byteorder)
            value = dtype(numpy_byteorder + dtypechar)
            align = value.alignment
        else:
            raise ValueError("Unknown PEP 3118 data type specifier %r" % spec)

        #
        # Native alignment may require padding
        #
        # Here we assume that the presence of a '@' character implicitly implies
        # that the start of the array is *already* aligned.
        #
        extra_offset = 0
        if byteorder == '@':
            start_padding = (-offset) % align
            intra_padding = (-value.itemsize) % align

            offset += start_padding

            if intra_padding != 0:
                if itemsize > 1 or (shape is not None and _prod(shape) > 1):
                    # Inject internal padding to the end of the sub-item
                    value = _add_trailing_padding(value, intra_padding)
                else:
                    # We can postpone the injection of internal padding,
                    # as the item appears at most once
                    extra_offset += intra_padding

            # Update common alignment
            common_alignment = (align*common_alignment
                                / _gcd(align, common_alignment))

        # Convert itemsize to sub-array
        if itemsize != 1:
            value = dtype((value, (itemsize,)))

        # Sub-arrays (2)
        if shape is not None:
            value = dtype((value, shape))

        # Field name
        this_explicit_name = False
        if spec and spec.startswith(':'):
            i = spec[1:].index(':') + 1
            name = spec[1:i]
            spec = spec[i+1:]
            explicit_name = True
            this_explicit_name = True
        else:
            name = get_dummy_name()

        if not is_padding or this_explicit_name:
            if name in fields:
                raise RuntimeError("Duplicate field name '%s' in PEP3118 format"
                                   % name)
            fields[name] = (value, offset)
            last_offset = offset
            if not this_explicit_name:
                next_dummy_name()

        byteorder = next_byteorder

        offset += value.itemsize
        offset += extra_offset

    # Check if this was a simple 1-item type
    if len(fields) == 1 and not explicit_name and fields['f0'][1] == 0 \
           and not is_subdtype:
        ret = fields['f0'][0]
    else:
        ret = dtype(fields)

    # Trailing padding must be explicitly added
    padding = offset - ret.itemsize
    if byteorder == '@':
        padding += (-offset) % common_alignment
    if is_padding and not this_explicit_name:
        ret = _add_trailing_padding(ret, padding)

    # Finished
    if is_subdtype:
        return ret, spec, common_alignment, byteorder
    else:
        return ret

def _add_trailing_padding(value, padding):
    """Inject the specified number of padding bytes at the end of a dtype"""
    from numpy.core.multiarray import dtype

    if value.fields is None:
        vfields = {'f0': (value, 0)}
    else:
        vfields = dict(value.fields)

    if value.names and value.names[-1] == '' and \
           value[''].char == 'V':
        # A trailing padding field is already present
        vfields[''] = ('V%d' % (vfields[''][0].itemsize + padding),
                       vfields[''][1])
        value = dtype(vfields)
    else:
        # Get a free name for the padding field
        j = 0
        while True:
            name = 'pad%d' % j
            if name not in vfields:
                vfields[name] = ('V%d' % padding, value.itemsize)
                break
            j += 1

        value = dtype(vfields)
        if '' not in vfields:
            # Strip out the name of the padding field
            names = list(value.names)
            names[-1] = ''
            value.names = tuple(names)
    return value

def _prod(a):
    p = 1
    for x in a:
        p *= x
    return p

def _gcd(a, b):
    """Calculate the greatest common divisor of a and b"""
    while b:
        a, b = b, a%b
    return a
