"""gdb printing extension for Numba types.
"""
import re

try:
    import gdb.printing
    import gdb
except ImportError:
    raise ImportError("GDB python support is not available.")


class NumbaArrayPrinter:

    def __init__(self, val):
        self.val = val

    def to_string(self):
        try:
            import numpy as np
            HAVE_NUMPY = True
        except ImportError:
            HAVE_NUMPY = False

        try:
            NULL = 0x0

            # Raw data references, these need unpacking/interpreting.

            # Member "data" is...
            # DW_TAG_member of DIDerivedType, tag of DW_TAG_pointer_type
            # encoding e.g. DW_ATE_float
            data = self.val["data"]

            # Member "itemsize" is...
            # DW_TAG_member of DIBasicType encoding DW_ATE_signed
            itemsize = self.val["itemsize"]

            # Members "shape" and "strides" are...
            # DW_TAG_member of DIDerivedType, the type is a DICompositeType
            # (it's a Numba UniTuple) with tag: DW_TAG_array_type, i.e. it's an
            # array repr, it has a basetype of e.g. DW_ATE_unsigned and also
            # "elements" which are referenced with a DISubrange(count: <const>)
            # to say how many elements are in the array.
            rshp = self.val["shape"]
            rstrides = self.val["strides"]

            # bool on whether the data is aligned.
            is_aligned = False

            # type information decode, simple type:
            ty_str = str(self.val.type)
            if HAVE_NUMPY and ('aligned' in ty_str or 'Record' in ty_str):
                ty_str = ty_str.replace('unaligned ','').strip()
                matcher = re.compile(r"array\((Record.*), (.*), (.*)\)\ \(.*")
                # NOTE: need to deal with "Alignment" else dtype size is wrong
                arr_info = [x.strip() for x in matcher.match(ty_str).groups()]
                dtype_str, ndim_str, order_str = arr_info
                rstr = r'Record\((.*\[.*\]);([0-9]+);(True|False)'
                rstr_match = re.match(rstr, dtype_str)
                # balign is unused, it's the alignment
                fields, balign, is_aligned_str = rstr_match.groups()
                is_aligned = is_aligned_str == 'True'
                field_dts = fields.split(',')
                struct_entries = []
                for f in field_dts:
                    splitted = f.split('[')
                    name = splitted[0]
                    dt_part = splitted[1:]
                    if len(dt_part) > 1:
                        raise TypeError('Unsupported sub-type: %s' % f)
                    else:
                        dt_part = dt_part[0]
                        if "nestedarray" in dt_part:
                            raise TypeError('Unsupported sub-type: %s' % f)
                        dt_as_str = dt_part.split(';')[0].split('=')[1]
                        dtype = np.dtype(dt_as_str)
                    struct_entries.append((name, dtype))
                    # The dtype is actually a record of some sort
                dtype_str = struct_entries
            else:  # simple type
                matcher = re.compile(r"array\((.*),(.*),(.*)\)\ \(.*")
                arr_info = [x.strip() for x in matcher.match(ty_str).groups()]
                dtype_str, ndim_str, order_str = arr_info
                # fix up unichr dtype
                if 'unichr x ' in dtype_str:
                    dtype_str = dtype_str[1:-1].replace('unichr x ', '<U')

            def dwarr2inttuple(dwarr):
                # Converts a gdb handle to a dwarf array to a tuple of ints
                fields = dwarr.type.fields()
                lo, hi = fields[0].type.range()
                return tuple([int(dwarr[x]) for x in range(lo, hi + 1)])

            # shape/strides extraction
            shape = dwarr2inttuple(rshp)
            strides = dwarr2inttuple(rstrides)

            # if data is not NULL
            if data != NULL:
                if HAVE_NUMPY:
                    # The data extent in bytes is:
                    # sum(shape * strides)
                    # get the data, then wire to as_strided
                    shp_arr = np.array([max(0, x - 1) for x in shape])
                    strd_arr = np.array(strides)
                    extent = np.sum(shp_arr * strd_arr)
                    extent += int(itemsize)
                    dtype_clazz = np.dtype(dtype_str, align=is_aligned)
                    dtype = dtype_clazz
                    this_proc = gdb.selected_inferior()
                    mem = this_proc.read_memory(int(data), extent)
                    arr_data = np.frombuffer(mem, dtype=dtype)
                    new_arr = np.lib.stride_tricks.as_strided(arr_data,
                                                              shape=shape,
                                                              strides=strides,)
                    return '\n' + str(new_arr)
                # Catch all for no NumPy
                return "array([...], dtype=%s, shape=%s)" % (dtype_str, shape)
            else:
                # Not yet initialized or NULLed out data
                buf = list(["NULL/Uninitialized"])
                return "array([" + ', '.join(buf) + "]" + ")"
        except Exception as e:
            return 'array[Exception: Failed to parse. %s]' % e


class NumbaComplexPrinter:

    def __init__(self, val):
        self.val = val

    def to_string(self):
        return "%s+%sj" % (self.val['real'], self.val['imag'])


class NumbaTuplePrinter:

    def __init__(self, val):
        self.val = val

    def to_string(self):
        buf = []
        fields = self.val.type.fields()
        for f in fields:
            buf.append(str(self.val[f.name]))
        return "(%s)" % ', '.join(buf)


class NumbaUniTuplePrinter:

    def __init__(self, val):
        self.val = val

    def to_string(self):
        # unituples are arrays
        fields = self.val.type.fields()
        lo, hi = fields[0].type.range()
        buf = []
        for i in range(lo, hi + 1):
            buf.append(str(self.val[i]))
        return "(%s)" % ', '.join(buf)


class NumbaUnicodeTypePrinter:

    def __init__(self, val):
        self.val = val

    def to_string(self):
        NULL = 0x0
        data = self.val["data"]
        nitems = self.val["length"]
        kind = self.val["kind"]
        if data != NULL:
            # This needs sorting out, encoding is wrong
            this_proc = gdb.selected_inferior()
            mem = this_proc.read_memory(int(data), nitems * kind)
            if isinstance(mem, memoryview):
                buf = bytes(mem).decode()
            else:
                buf = mem.decode('utf-8')
        else:
            buf = str(data)
        return "'%s'" % buf


def _create_printers():
    printer = gdb.printing.RegexpCollectionPrettyPrinter("Numba")
    printer.add_printer('Numba unaligned array printer', '^unaligned array\\(',
                        NumbaArrayPrinter)
    printer.add_printer('Numba array printer', '^array\\(', NumbaArrayPrinter)
    printer.add_printer('Numba complex printer', '^complex[0-9]+\\ ',
                        NumbaComplexPrinter)
    printer.add_printer('Numba Tuple printer', '^Tuple\\(',
                        NumbaTuplePrinter)
    printer.add_printer('Numba UniTuple printer', '^UniTuple\\(',
                        NumbaUniTuplePrinter)
    printer.add_printer('Numba unicode_type printer', '^unicode_type\\s+\\(',
                        NumbaUnicodeTypePrinter)
    return printer


# register the Numba pretty printers for the current object
gdb.printing.register_pretty_printer(gdb.current_objfile(), _create_printers())
