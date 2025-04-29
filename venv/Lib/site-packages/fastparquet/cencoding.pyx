# https://cython.readthedocs.io/en/latest/src/userguide/
#   source_files_and_compilation.html#compiler-directives
# cython: profile=False
# cython: linetrace=False
# cython: binding=False
# cython: language_level=3
# cython: initializedcheck=False
# cython: boundscheck=False
# cython: wraparound=False
# cython: overflowcheck=False
# cython: cdivision=True
# cython: always_allow_keywords=False

import cython
import numpy as np
cdef extern from "string.h":
    void *memcpy(void *dest, const void *src, size_t n)
from cpython cimport (
    PyBytes_FromStringAndSize, PyBytes_GET_SIZE, PyUnicode_DecodeUTF8,
)
from libc.stdint cimport int8_t, uint8_t, uint32_t, int32_t, uint64_t, int64_t


cpdef void read_rle(NumpyIO file_obj, int32_t header, int32_t bit_width, NumpyIO o, int32_t itemsize=4):
    """Read a run-length encoded run from the given fo with the given header and bit_width.

    The count is determined from the header and the width is used to grab the
    value that's repeated. Yields the value repeated count times.
    """
    cdef:
        uint32_t count, width, i, vals_left
        int32_t data = 0
        char * inptr = file_obj.get_pointer()
        char * outptr = o.get_pointer()
    count = header >> 1
    width = (bit_width + 7) // 8
    for i in range(width):
        data |= (inptr[0] & 0xff) << (i * 8)
        inptr += 1
    vals_left = (o.nbytes - o.loc) // itemsize
    if count > vals_left:
        count = vals_left
    if itemsize == 4:
        for i in range(count):
            (<int32_t*>outptr)[0] = data
            outptr += 4
    else:
        for i in range(count):
            outptr[0] = data & 0xff
            outptr += 1
    o.loc += outptr - o.get_pointer()
    file_obj.loc += inptr - file_obj.get_pointer()


cpdef int32_t width_from_max_int(int64_t value):
    """Convert the value specified to a bit_width."""
    cdef int32_t i
    for i in range(0, 64):
        if value == 0:
            return i
        value >>= 1


cdef int32_t _mask_for_bits(int32_t i):
    """Generate a mask to grab `i` bits from an int value."""
    return (1 << i) - 1


cpdef void read_bitpacked1(NumpyIO file_obj, int32_t count, NumpyIO o):
    # implementation of np.unpackbits with output array. Output is int8 array
    cdef:
        char * inptr = file_obj.get_pointer()
        char * outptr = o.get_pointer()
        char * endptr
        unsigned char data
        int32_t counter, i, startcount=count
    if count > o.nbytes - o.loc:
        count = o.nbytes - o.loc
    for counter in range(count // 8):
        # whole bytes
        data = inptr[0]
        inptr += 1
        for i in range(8):
            outptr[0] = data & 1
            outptr += 1
            data >>= 1
    if count % 8:
        # remaining values in the last byte
        data = <int32_t>inptr[0]
        inptr += 1
        for i in range(count % 8):
            outptr[0] = data & 1
            outptr += 1
            data >>= 1
    file_obj.loc += (startcount + 7) // 8
    o.loc += count


cpdef void write_bitpacked1(NumpyIO file_obj, int32_t count, NumpyIO o):
    # implementation of np.packbits with output array. Input is int8 array
    cdef char * inptr
    cdef char * outptr
    cdef char data = 0
    cdef int32_t counter, i
    cdef int64_t indata
    outptr = o.get_pointer()
    inptr = file_obj.get_pointer()
    for counter in range(count // 8):
        # fetch a long in one op, instead of byte by byte
        indata = (<int64_t*>inptr)[0]
        inptr += 8
        for i in range(8):
            data = data << 1 | (indata & 1)
            indata >>= 8
        outptr[0] = data
        outptr += 1
    if count % 8:
        # leftover partial byte
        data = 0
        for i in range(count % 8):
            data = data << 1 | (inptr[0] != 0)
            inptr += 1
        outptr[0] = data
        outptr += 1
    file_obj.loc += count * 4
    o.loc += (count + 7) // 8


cpdef void read_bitpacked(NumpyIO file_obj, int32_t header, int32_t width, NumpyIO o, int32_t itemsize=4):
    """
    Read values packed into width-bits each (which can be >8)
    """
    cdef:
        uint32_t count, mask, data, vals_left
        unsigned char left = 8, right = 0
        char * inptr = file_obj.get_pointer()
        char * outptr = o.get_pointer()
        char * endptr

    count = (header >> 1) * 8
    # TODO: special case for width=1, 2, 4, 8
    if width == 1 and itemsize == 1:
        read_bitpacked1(file_obj, count, o)
        return
    endptr = (o.nbytes - o.loc) + outptr - itemsize
    mask = _mask_for_bits(width)
    data = 0xff & <int32_t>inptr[0]
    inptr += 1
    while count:
        if right > 8:
            data >>= 8
            left -= 8
            right -= 8
        elif left - right < width:
            data |= (inptr[0] & 0xff) << left
            inptr += 1
            left += 8
        else:
            if outptr <= endptr:
                if itemsize == 4:
                    (<int32_t*>outptr)[0] = <int32_t>(data >> right & mask)
                    outptr += 4
                else:
                    outptr[0] = data >> right & mask
                    outptr += 1
            count -= 1
            right += width
    o.loc = o.loc + outptr - o.get_pointer()
    file_obj.loc += inptr - file_obj.get_pointer()


cpdef uint64_t read_unsigned_var_int(NumpyIO file_obj):
    """Read a value using the unsigned, variable int encoding.
    file-obj is a NumpyIO of bytes; avoids struct to allow numba-jit
    """
    cdef uint64_t result = 0
    cdef int32_t shift = 0
    cdef char byte
    cdef char * inptr = file_obj.get_pointer()

    while True:
        byte = inptr[0]
        inptr += 1
        result |= (<int64_t>(byte & 0x7F) << shift)
        if (byte & 0x80) == 0:
            break
        shift += 7
    file_obj.loc += inptr - file_obj.get_pointer()
    return result


cpdef void read_rle_bit_packed_hybrid(NumpyIO io_obj, int32_t width, uint32_t length, NumpyIO o,
                                      int32_t itemsize=4):
    """Read values from `io_obj` using the rel/bit-packed hybrid encoding.

    If length is not specified, then a 32-bit int is read first to grab the
    length of the encoded data.

    file-obj is a NumpyIO of bytes; o if an output NumpyIO of int32 or int8/bool

    The caller can tell the number of elements in the output by looking
    at .tell().
    """
    cdef int32_t start, header
    if length is False:
        length = <uint32_t>io_obj.read_int()
    start = io_obj.loc
    while io_obj.loc - start < length and o.loc < o.nbytes:
        header = <int32_t>read_unsigned_var_int(io_obj)
        if header & 1 == 0:
            read_rle(io_obj, header, width, o, itemsize)
        else:
            read_bitpacked(io_obj, header, width, o, itemsize)


cdef void delta_read_bitpacked(NumpyIO file_obj, uint8_t bitwidth,
                               NumpyIO o, uint64_t count, uint8_t longval=0):
    cdef:
        uint64_t data = 0
        int8_t left = 0
        int8_t right = 0
        uint64_t mask = 0XFFFFFFFFFFFFFFFF >> (64 - bitwidth)
    while count > 0:
        if (left - right) < bitwidth:
            data = data | (<uint64_t>file_obj.read_byte() << left)
            left += 8
        elif right > 8:
            data >>= 8
            left -= 8
            right -= 8
        else:
            if longval:
                o.write_long((data >> right) & mask)
            else:
                o.write_int((data >> right) & mask)
            right += bitwidth
            count -= 1


cpdef void delta_binary_unpack(NumpyIO file_obj, NumpyIO o, uint8_t longval=0):
    cdef:
        uint64_t block_size = read_unsigned_var_int(file_obj)
        uint64_t miniblock_per_block = read_unsigned_var_int(file_obj)
        int64_t count = read_unsigned_var_int(file_obj)
        int64_t value = zigzag_long(read_unsigned_var_int(file_obj))
        int64_t block, min_delta, i, j, values_per_miniblock, temp
        const uint8_t[:] bitwidths
        uint8_t bitwidth
    values_per_miniblock = block_size // miniblock_per_block
    while True:
        min_delta = zigzag_long(read_unsigned_var_int(file_obj))
        bitwidths = file_obj.read(miniblock_per_block)
        for i in range(miniblock_per_block):
            bitwidth = bitwidths[i]
            if bitwidth:
                temp = o.loc
                if count > 1:
                    # no more diffs if on last value
                    delta_read_bitpacked(file_obj, bitwidth, o, values_per_miniblock, longval)
                o.loc = temp
                for j in range(values_per_miniblock):
                    if longval:
                        temp = o.read_long()
                        o.loc -= 8
                        o.write_long(value)
                    else:
                        temp = o.read_int()
                        o.loc -= 4
                        o.write_int(value)
                    value += min_delta + temp
                    count -= 1
                    if count <= 0:
                        return
            else:
                for j in range(values_per_miniblock):
                    if longval:
                        o.write_long(value)
                    else:
                        o.write_int(value)
                    value += min_delta
                    count -= 1
                    if count <= 0:
                        return


cpdef void encode_unsigned_varint(uint64_t x, NumpyIO o):  # pragma: no cover
    while x > 127:
        o.write_byte((x & 0x7F) | 0x80)
        x >>= 7
    o.write_byte(x)


cpdef encode_bitpacked(int32_t[:] values, int32_t width, NumpyIO o):
    """
    Write values packed into width-bits each (which can be >8)
    """

    cdef int32_t bit_packed_count = (values.shape[0] + 7) // 8
    encode_unsigned_varint(bit_packed_count << 1 | 1, o)  # write run header
    cdef int32_t bit=0, bits=0, v, counter
    for counter in range(values.shape[0]):
        v = values[counter]
        bits |= v << bit
        bit += width
        while bit >= 8:
            o.write_byte(bits & 0xff)
            bit -= 8
            bits >>= 8
    if bit:
        o.write_byte(bits)


cpdef void encode_rle_bp(int32_t[:] data, int32_t width, NumpyIO o, int32_t withlength = 0):
    cdef uint32_t start, end
    if withlength:
        start = o.tell()
        o.seek(4, 1)
    encode_bitpacked(data, width, o)
    if withlength:
        end = o.tell()
        o.seek(start)
        o.write_int(end - start - 4)
        o.seek(end)


@cython.freelist(100)
@cython.final
cdef class NumpyIO(object):
    """
    Read or write from a numpy array like a file object

    The main purpose is to keep track of the current location in the memory
    """
    cdef const uint8_t[:] data
    cdef uint32_t loc, nbytes
    cdef char* ptr
    cdef char writable

    def __cinit__(self, const uint8_t[::1] data):
        self.data = data
        self.loc = 0
        self.ptr = <char*>&data[0]
        self.nbytes = data.shape[0]

    cdef char* get_pointer(self):
        return self.ptr + self.loc

    @property
    def len(self):
        return self.nbytes

    cpdef const uint8_t[:] read(self, int32_t x=-1):
        cdef const uint8_t[:] out
        if x < 1:
            x = self.nbytes - self.loc
        out = self.data[self.loc:self.loc + x]
        self.loc += x
        return out

    cpdef uint8_t read_byte(self):
        cdef char out
        out = self.ptr[self.loc]
        self.loc += 1
        return out

    cpdef int32_t read_int(self):
        cdef int32_t i
        if self.nbytes - self.loc < 4:
            return 0
        i = (<int32_t*> self.get_pointer())[0]
        self.loc += 4
        return i

    cpdef void write(self, const char[::1] d):
        memcpy(<void*>self.ptr[self.loc], <void*>&d[0], d.shape[0])
        self.loc += d.shape[0]

    cpdef void write_byte(self, uint8_t b):
        if self.loc >= self.nbytes:
            # ignore attempt to write past end of buffer
            return
        self.ptr[self.loc] = b
        self.loc += 1

    cpdef void write_int(self, int32_t i):
        if self.nbytes - self.loc < 4:
            return
        (<int32_t*> self.get_pointer())[0] = i
        self.loc += 4

    cdef void write_long(self, int64_t i):
        if self.nbytes - self.loc < 8:
            return
        (<int64_t*> self.get_pointer())[0] = i
        self.loc += 8

    cdef int64_t read_long(self):
        cdef int64_t i
        if self.nbytes - self.loc < 8:
            return 0
        i = (<int64_t*> self.get_pointer())[0]
        self.loc += 8
        return i

    cdef void write_many(self, char b, int32_t count):
        cdef int32_t i
        for i in range(count):
            self.write_byte(b)

    cpdef int32_t tell(self):
        return self.loc

    cpdef uint32_t seek(self, int32_t loc, int32_t whence=0):
        if whence == 0:
            self.loc = loc
        elif whence == 1:
            self.loc += loc
        elif whence == 2:
            self.loc = self.nbytes + loc
        if self.loc > self.nbytes:
            self.loc = self.nbytes
        return self.loc

    @cython.wraparound(False)
    cpdef const uint8_t[:] so_far(self):
        """ In write mode, the data we have gathered until now
        """
        return self.data[:self.loc]


def _assemble_objects(object[:] assign, const uint8_t[:] defi, const uint8_t[:] rep,
                      val, dic, d,
                      char null, null_val, int32_t max_defi, int32_t prev_i):
    """Dremel-assembly of arrays of values into lists

    Parameters
    ----------
    assign: array dtype O
        To insert lists into
    defi: int array
        Definition levels, max 3
    rep: int array
        Repetition levels, max 1
    dic: array of labels or None
        Applied if d is True
    d: bool
        Whether to dereference dict values
    null: bool
        Can an entry be None?
    null_val: bool
        can list elements be None
    max_defi: int
        value of definition level that corresponds to non-null
    prev_i: int
        1 + index where the last row in the previous page was inserted (0 if first page)
    """
    cdef int32_t counter, i, re, de
    cdef int32_t vali = 0
    cdef char started = False, have_null = False
    if d:
        # dereference dict values
        val = dic[val]
    i = prev_i
    part = []
    for counter in range(rep.shape[0]):
        de = defi[counter] if defi is not None else max_defi
        re = rep[counter]
        if not re:
            # new row - save what we have
            if started:
                assign[i] = None if have_null else part
                part = []
                i += 1
            else:
                # first time: no row to save yet, unless it's a row continued from previous page
                if vali > 0:
                    assign[i - 1].extend(part) # add the items to previous row
                    part = []
                    # don't increment i since we only filled i-1
                started = True
        if de == max_defi:
            # append real value to current item
            part.append(val[vali])
            vali += 1
        elif de > null:
            # append null to current item
            part.append(None)
        # next object is None as opposed to an object
        have_null = de == 0 and null
    if started: # normal case - add the leftovers to the next row
        assign[i] = None if have_null else part
    else: # can only happen if the only elements in this page are the continuation of the last row from previous page
        assign[i - 1].extend(part)
    return i


cdef int64_t nat = -9223372036854775808


cpdef void time_shift(const int64_t[::1] data, int32_t factor=1000):
    cdef int32_t i
    cdef int64_t * ptr
    cdef int64_t value
    ptr = <int64_t*>&data[0]
    for i in range(data.shape[0]):
        if ptr[0] != nat:
            ptr[0] *= factor
        ptr += 1


cdef int32_t zigzag_int(uint64_t n):
    return (n >> 1) ^ -(n & 1)


cdef int64_t zigzag_long(uint64_t n):
    return (n >> 1) ^ -(n & 1)


cdef uint64_t long_zigzag(int64_t n):
    return (n << 1) ^ (n >> 63)


cpdef dict read_thrift(NumpyIO data):
    cdef char byte, id = 0, bit
    cdef int32_t size
    cdef dict out = {}
    cdef bint hasi64 = 0
    cdef bint hasi32 = 0
    cdef list i32 = None
    while True:
        byte = data.read_byte()
        if byte == 0:
            break
        id += (byte & 0b11110000) >> 4
        bit = byte & 0b00001111
        if bit == 5:
            out[id] = zigzag_long(read_unsigned_var_int(data))
            hasi32 = True
            if i32 is None:
                i32 = list()
            i32.append(id)
        elif bit == 6:
            out[id] = zigzag_long(read_unsigned_var_int(data))
            hasi64 = True
        elif bit == 7:
            out[id] = <double>data.get_pointer()[0]
            data.seek(8, 1)
        elif bit == 8:
            size = read_unsigned_var_int(data)
            out[id] = PyBytes_FromStringAndSize(data.get_pointer(), size)
            data.seek(size, 1)
        elif bit == 9:
            out[id] = read_list(data)
        elif bit == 12:
            out[id] = read_thrift(data)
        elif bit == 1:
            out[id] = True
        elif bit == 2:
            out[id] = False
        elif bit == 4:
            # I16
            out[id] = zigzag_long(read_unsigned_var_int(data))
        elif bit == 3:
            # I8
            out[id] = data.read_byte()
        else:
            print("Corrupted thrift data at ", data.tell(), ": ", id, bit)
    if hasi32:
        if hasi64:
            out["i32list"] = i32
        else:
            out["i32"] = 1
    return out


cdef list read_list(NumpyIO data):
    cdef unsigned char byte, typ
    cdef int32_t size, bsize, _
    byte = data.read_byte()
    if byte >= 0xf0:  # 0b11110000
        size = read_unsigned_var_int(data)
    else:
        size = ((byte & 0xf0) >> 4)
    out = []
    typ = byte & 0x0f # 0b00001111
    if typ == 5 or typ == 6:
        for _ in range(size):
            out.append(zigzag_long(read_unsigned_var_int(data)))
    elif typ == 8:
        for _ in range(size):
            # all parquet list types contain str, not bytes
            bsize = read_unsigned_var_int(data)
            out.append(PyUnicode_DecodeUTF8(data.get_pointer(), bsize, "ignore"))
            data.seek(bsize, 1)
    else:
        for _ in range(size):
            out.append(read_thrift(data))

    return out


cpdef void write_thrift(dict data, NumpyIO output):
    cdef int i, l, prev = 0
    cdef int delt = 0
    cdef double d
    cdef bytes b
    cdef char * c
    cdef int i32 = "i32" in data
    cdef list i32s
    if "i32list" in data:
        i32 = 2
        i32s = data['i32list']
    for i in range(1, 14):  # 14 is the max number of fields
        if i not in data:
            continue
        val = data.get(i)
        if val is None:
            # not defined - skip (None is default on load)
            continue
        delt = i - prev
        prev = i
        if isinstance(val, bool):
            if val is True:
                output.write_byte((delt << 4) | 1)
            else:
                output.write_byte((delt << 4) | 2)
        elif isinstance(val, int):
            if i32 == 1 or (i32 == 2 and i in i32s):
                output.write_byte((delt << 4) | 5)
            else:
                output.write_byte((delt << 4) | 6)
            encode_unsigned_varint(long_zigzag(<int64_t>val), output)
        elif isinstance(val, float):
            output.write_byte((delt << 4) | 7)
            d = val
            (<double*>output.get_pointer())[0] = d
            output.loc += 8
        elif isinstance(val, bytes):
            output.write_byte((delt << 4) | 8)
            l = PyBytes_GET_SIZE(<bytes>val)
            encode_unsigned_varint(l, output)
            c = val
            memcpy(<void*>output.get_pointer(), <void*>c, l)
            output.loc += l
        elif isinstance(val, str):
            output.write_byte((delt << 4) | 8)
            b = (<str>val).encode()
            l = PyBytes_GET_SIZE(b)
            encode_unsigned_varint(l, output)
            c = b
            memcpy(<void*>output.get_pointer(), <void*>c, l)
            output.loc += l
        elif isinstance(val, list):
            output.write_byte((delt << 4) | 9)
            write_list(<list>val, output)
        elif isinstance(val, ThriftObject):
            output.write_byte((delt << 4) | 12)
            write_thrift((<ThriftObject>val).data, output)
        else:
            output.write_byte((delt << 4) | 12)
            write_thrift(<dict>val, output)
    output.write_byte(0)


cdef void write_list(list data, NumpyIO output):
    cdef int l = len(data)
    cdef int i
    cdef ThriftObject dd
    cdef bytes b
    cdef str s
    cdef char * c
    if l:
        first = data[0]
        if isinstance(first, int):
            if l > 14:  # all lists are i64
                output.write_byte(5 | 0b11110000)
                encode_unsigned_varint(l, output)
            else:
                output.write_byte(5 | (l << 4))
            for i in data:
                encode_unsigned_varint(long_zigzag(i), output)
        elif isinstance(first, bytes):
            if l > 14:
                output.write_byte(8 | 0b11110000)
                encode_unsigned_varint(l, output)
            else:
                output.write_byte(8 | (l << 4))
            for b in data:
                i = PyBytes_GET_SIZE(b)
                encode_unsigned_varint(i, output)
                c = b
                memcpy(<void*>output.get_pointer(), <void*>c, i)
                output.loc += i
        elif isinstance(first, str):
            if l > 14:
                output.write_byte(8 | 0b11110000)
                encode_unsigned_varint(l, output)
            else:
                output.write_byte(8 | (l << 4))
            for s in data:
                b = s.encode("utf8", "ignore")
                i = PyBytes_GET_SIZE(b)
                encode_unsigned_varint(i, output)
                c = b
                memcpy(<void*>output.get_pointer(), <void*>c, i)
                output.loc += i
        else: # STRUCT
            if l > 14:
                output.write_byte(12 | 0b11110000)
                encode_unsigned_varint(l, output)
            else:
                output.write_byte(12 | (l << 4))
            for d in data:
                if isinstance(d, ThriftObject):
                    write_thrift((<ThriftObject>d).data, output)
                else:
                    write_thrift(d, output)
    else:
        # Not sure if zero-length list is allowed
        encode_unsigned_varint(0, output)


def from_buffer(buffer, name=None):
    cdef NumpyIO buf
    if isinstance(buffer, NumpyIO):
        buf = buffer
    else:
        buf = NumpyIO(buffer)
    cdef dict o = read_thrift(buf)
    if name is not None:
        return ThriftObject(name, o)
    return o


@cython.freelist(1000)
@cython.final
cdef class ThriftObject:

    cdef str name
    cdef dict spec
    cdef dict children
    cdef dict data

    def __init__(self, str name, dict indict):
        self.name = name
        self.spec = specs[name]
        self.children = children.get(name, {})
        self.data = indict

    def __getattr__(self, str item):
        cdef str ch
        if item in self.spec:
            out = self.get(self.spec[item], None)
            ch = self.children.get(item)
            if ch is not None and out is not None:
                if isinstance(out, list):
                    return [ThriftObject(ch, o) if isinstance(o, dict) else o for o in out]
                return ThriftObject(ch, out) if isinstance(out, dict) else out
            return out
        else:
            try:
                return self.data[item]
            except KeyError:
                raise AttributeError

    def __setitem__(self, key, value):
        self.data[key] = value

    def __getitem__(self, item):
        return self.data.get(item)

    def __delitem__(self, key):
        self.data.pop(key)

    def get(self, key, default=None):
        return self.data.get(key, default)

    def __setattr__(self, str item, value):
        cdef int i = self.spec[item]
        cdef int j
        if isinstance(value, ThriftObject):
            self.data[i] = value.data
        elif isinstance(value, list):
            self.data[i] = [(<ThriftObject>v).data for v in value]
        else:
            self.data[i] = value

    def __delattr__(self, item):
        cdef int i = self.spec[item]
        del self.data[i]

    cpdef const uint8_t[:] to_bytes(self):
        """raw serialise of internal state"""
        cdef int size = 0
        if self.name == "RowGroup":
            size = 1000 * len(self[1])  # num-columns
        elif self.name == "FileMetaData":
            # num-cols * num-rgs + size of key-values
            size = 1000 * len(self[4]) * len(self[2]) + len(str(self[5]))
        if size < 500000:
            size = 500000
        cdef uint8_t[::1] ser_buf = np.empty(size, dtype='uint8')
        cdef NumpyIO o = NumpyIO(ser_buf)
        write_thrift(self.data, o)
        return o.so_far()

    def __reduce_ex__(self, _):
        # TODO: to_bytes returns a memoryview, so could sideband for pickle 5
        return from_buffer, (bytes(self.to_bytes()), self.name)

    @property
    def thrift_name(self):
        return self.name

    @property
    def contents(self):
        return self.data

    from_buffer = from_buffer

    def copy(self):
        """shallow copy"""
        return type(self)(self.name, self.data.copy())

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memodict={}):
        import copy
        d = copy.deepcopy(self.data)
        return ThriftObject(self.name, d)

    cpdef _asdict(self):
        """Create dict version with field names instead of integers"""
        cdef str k
        cdef out = {}
        for k in self.spec:
            if k in self.children:
                lower = getattr(self, k)
                if lower is None:
                    out[k] = None
                elif isinstance(lower, list):
                    out[k] = [l._asdict() for l in lower]
                else:
                    out[k] = lower._asdict()
            else:
                lower = getattr(self, k)
                if isinstance(lower, bytes):
                    lower = str(lower)
                elif isinstance(lower, list) and lower and isinstance(lower[0], bytes):
                    lower = [str(l) for l in lower]
                out[k] = lower
        return out

    def __dir__(self):
        """Lists attributed"""
        return list(self.spec)

    def __repr__(self):
        alt = self._asdict()
        try:
            import yaml
            return yaml.dump(alt)
        except ImportError:
            return str(alt)

    def __eq__(self, other):
        if isinstance(other, ThriftObject):
            return dict_eq(self.contents, other.contents)
        elif isinstance(other, dict):
            return dict_eq(self.contents, other)
        return False

    @staticmethod
    def from_fields(thrift_name,bint i32=0, list i32list=None, **kwargs):
        cdef spec = specs[thrift_name]
        cdef int i
        cdef str k
        cdef dict out = {}
        for k, i in spec.items():  # ensure field index increases monotonically
            if k in kwargs:
                # missing fields are implicitly None
                v = kwargs[k]
                if isinstance(v, ThriftObject):
                    out[i] = (<ThriftObject>v).data
                elif isinstance(v, list) and v and isinstance(v[0], ThriftObject):
                    out[i] = [(<ThriftObject>it).data for it in v]
                else:
                    out[i] = v
        if i32:
            # integer fields are all 32-bit
            out['i32'] = 1
        if i32list:
            # given integer fields are 32-bit
            out['i32list'] = i32list
        return ThriftObject(thrift_name, out)


def dict_eq(d1, d2):
    """ dicts are equal if none-None keys match """
    if isinstance(d1, ThriftObject):
        d1 = d1.contents
    if isinstance(d2, ThriftObject):
        d2 = d2.contents
    for k in set(d1).union(d2):
        if not isinstance(k, int):
            # dynamic fields are immaterial
            continue
        if d1.get(k, None) is None:
            if d2.get(k, None) is None:
                continue
            return False
        if d2.get(k, None) is None:
            return False
        elif isinstance(d1[k], dict):
            if not dict_eq(d1[k], d2[k]):
                return False
        elif isinstance(d1[k], list):
            if len(d1[k]) != len(d2[k]):
                return False
            # Recursive call as described in
            # https://github.com/dask/fastparquet/pull/723#issuecomment-995147362
            if any(not dict_eq(a,b) if isinstance(a, dict) else (a != b)
                   for a, b in zip(d1[k], d2[k])):
                return False
        elif isinstance(d1[k], str):
            s = d2[k]
            if d1[k] != (s.decode() if isinstance(s, bytes) else s):
                return False
        else:
            if d1.get(k, None) != d2.get(k, None):
                return False
    return True


cdef dict specs = {
    'Statistics': {'max': 1,
                   'min': 2,
                   'null_count': 3,
                   'distinct_count': 4,
                   'max_value': 5,
                   'min_value': 6},
    'StringType': {},
    'UUIDType': {},
    'MapType': {},
    'ListType': {},
    'EnumType': {},
    'DateType': {},
    'NullType': {},
    'DecimalType': {'scale': 1, 'precision': 2},
    'MilliSeconds': {},
    'MicroSeconds': {},
    'NanoSeconds': {},
    'TimeUnit': {'MILLIS': 1, 'MICROS': 2, 'NANOS': 3},
    'TimestampType': {'isAdjustedToUTC': 1, 'unit': 2},
    'TimeType': {'isAdjustedToUTC': 1, 'unit': 2},
    'IntType': {'bitWidth': 1, 'isSigned': 2},
    'JsonType': {},
    'BsonType': {},
    'LogicalType': {'STRING': 1,
                    'MAP': 2,
                    'LIST': 3,
                    'ENUM': 4,
                    'DECIMAL': 5,
                    'DATE': 6,
                    'TIME': 7,
                    'TIMESTAMP': 8,
                    'INTEGER': 10,
                    'UNKNOWN': 11,
                    'JSON': 12,
                    'BSON': 13,
                    'UUID': 14},
    'SchemaElement': {'type': 1,
                      'type_length': 2,
                      'repetition_type': 3,
                      'name': 4,
                      'num_children': 5,
                      'converted_type': 6,
                      'scale': 7,
                      'precision': 8,
                      'field_id': 9,
                      'logicalType': 10},
    'DataPageHeader': {'num_values': 1,
                       'encoding': 2,
                       'definition_level_encoding': 3,
                       'repetition_level_encoding': 4,
                       'statistics': 5},
    'IndexPageHeader': {},
    'DictionaryPageHeader': {'num_values': 1, 'encoding': 2, 'is_sorted': 3},
    'DataPageHeaderV2': {'num_values': 1,
                         'num_nulls': 2,
                         'num_rows': 3,
                         'encoding': 4,
                         'definition_levels_byte_length': 5,
                         'repetition_levels_byte_length': 6,
                         'is_compressed': 7,
                         'statistics': 8},
    'SplitBlockAlgorithm': {},
    'BloomFilterAlgorithm': {'BLOCK': 1},
    'XxHash': {},
    'BloomFilterHash': {'XXHASH': 1},
    'Uncompressed': {},
    'PageHeader': {'type': 1,
                   'uncompressed_page_size': 2,
                   'compressed_page_size': 3,
                   'crc': 4,
                   'data_page_header': 5,
                   'index_page_header': 6,
                   'dictionary_page_header': 7,
                   'data_page_header_v2': 8},
    'KeyValue': {'key': 1, 'value': 2},
    'SortingColumn': {'column_idx': 1, 'descending': 2, 'nulls_first': 3},
    'PageEncodingStats': {'page_type': 1, 'encoding': 2, 'count': 3},
    'ColumnMetaData': {'type': 1,
                       'encodings': 2,
                       'path_in_schema': 3,
                       'codec': 4,
                       'num_values': 5,
                       'total_uncompressed_size': 6,
                       'total_compressed_size': 7,
                       'key_value_metadata': 8,
                       'data_page_offset': 9,
                       'index_page_offset': 10,
                       'dictionary_page_offset': 11,
                       'statistics': 12,
                       'encoding_stats': 13,
                       'bloom_filter_offset': 14},
    'ColumnChunk': {'file_path': 1,
                    'file_offset': 2,
                    'meta_data': 3,
                    'offset_index_offset': 4,
                    'offset_index_length': 5,
                    'column_index_offset': 6,
                    'column_index_length': 7,
                    'crypto_metadata': 8,
                    'encrypted_column_metadata': 9},
    'RowGroup': {'columns': 1,
                 'total_byte_size': 2,
                 'num_rows': 3,
                 'sorting_columns': 4,
                 'file_offset': 5,
                 'total_compressed_size': 6,
                 'ordinal': 7},
    'TypeDefinedOrder': {},
    'ColumnOrder': {'TYPE_ORDER': 1},
    'PageLocation': {'offset': 1,
                     'compressed_page_size': 2,
                     'first_row_index': 3},
    'OffsetIndex': {'page_locations': 1},
    'ColumnIndex': {'null_pages': 1,
                    'min_values': 2,
                    'max_values': 3,
                    'boundary_order': 4,
                    'null_counts': 5},
    'FileMetaData': {'version': 1,
                     'schema': 2,
                     'num_rows': 3,
                     'row_groups': 4,
                     'key_value_metadata': 5,
                     'created_by': 6,
                     'column_orders': 7,
                     'encryption_algorithm': 8,
                     'footer_signing_key_metadata': 9},
}

cdef dict children = {
    'TimeUnit': {'MILLIS': 'MilliSeconds',
                 'MICROS': 'MicroSeconds',
                 'NANOS': 'NanoSeconds'},
    'TimestampType': {'unit': 'TimeUnit'},
    'TimeType': {'unit': 'TimeUnit'},
    'LogicalType': {'STRING': 'StringType',
                    'MAP': 'MapType',
                    'LIST': 'ListType',
                    'ENUM': 'EnumType',
                    'DECIMAL': 'DecimalType',
                    'DATE': 'DateType',
                    'TIME': 'TimeType',
                    'TIMESTAMP': 'TimestampType',
                    'INTEGER': 'IntType',
                    'UNKNOWN': 'NullType',
                    'JSON': 'JsonType',
                    'BSON': 'BsonType',
                    'UUID': 'UUIDType'},
    'SchemaElement': {'logicalType': 'LogicalType'},
    'DataPageHeader': {'statistics': 'Statistics'},
    'DataPageHeaderV2': {'statistics': 'Statistics'},
    'PageHeader': {'data_page_header': 'DataPageHeader',
                   'index_page_header': 'IndexPageHeader',
                   'dictionary_page_header': 'DictionaryPageHeader',
                   'data_page_header_v2': 'DataPageHeaderV2'},
    'ColumnMetaData': {'key_value_metadata': 'KeyValue',
                       'statistics': 'Statistics',
                       'encoding_stats': 'PageEncodingStats'},
    'ColumnCryptoMetaData': {'ENCRYPTION_WITH_FOOTER_KEY': 'EncryptionWithFooterKey',
                             'ENCRYPTION_WITH_COLUMN_KEY': 'EncryptionWithColumnKey'},
    'ColumnChunk': {'meta_data': 'ColumnMetaData',
                    'crypto_metadata': 'ColumnCryptoMetaData'},
    'RowGroup': {'columns': 'ColumnChunk', 'sorting_columns': 'SortingColumn'},
    'ColumnOrder': {'TYPE_ORDER': 'TypeDefinedOrder'},
    'OffsetIndex': {'page_locations': 'PageLocation'},
    'FileMetaData': {'schema': 'SchemaElement',
                     'row_groups': 'RowGroup',
                     'key_value_metadata': 'KeyValue',
                     'column_orders': 'ColumnOrder',
                     'encryption_algorithm': 'EncryptionAlgorithm'},
}

# specs = {}
# for o in [o for o in fastparquet.parquet_thrift.__dict__.values() if isinstance(o, type)]:
#     if hasattr(o, "thrift_spec"):
#         specs[o.__name__] = {k[2]: k[0] for k in o.thrift_spec if k}
#
#
#
# children = {}
# for o in [o for o in fastparquet.parquet_thrift.__dict__.values() if isinstance(o, type)]:
#     if hasattr(o, "thrift_spec"):
#         bit = {}
#         for k in o.thrift_spec:
#             if k and k[1] == fastparquet.parquet_thrift.TType.STRUCT and hasattr(k[3][0], "thrift_spec"):
#                 bit[k[2]] =  k[3][0].__name__
#             elif k and k[1] == fastparquet.parquet_thrift.TType.LIST and k[3][0] == \
#                 fastparquet.parquet_thrift.TType.STRUCT:
#                 bit[k[2]] =  k[3][1][0].__name__
#         if bit:
#             children[o.__name__] = bit
#
