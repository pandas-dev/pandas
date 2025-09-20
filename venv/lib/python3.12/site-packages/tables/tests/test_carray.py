import sys
from pathlib import Path

import numpy as np

import tables as tb
from tables.tests import common


def foreign_byteorder():
    return {"little": "big", "big": "little"}[sys.byteorder]


class BasicTestCase(common.TempFileMixin, common.PyTablesTestCase):
    # Default values
    obj = None
    flavor = "numpy"
    type = "int32"
    shape = (2, 2)
    start = 0
    stop = 10
    step = 1
    length = 1
    chunkshape = (5, 5)
    byteorder = None
    compress = 0
    complib = "zlib"  # Default compression library
    shuffle = 0
    bitshuffle = 0
    fletcher32 = 0
    reopen = 1  # Tells whether the file has to be reopened on each test or not

    def setUp(self):
        super().setUp()
        # Create an instance of an HDF5 Table
        self.rootgroup = self.h5file.root
        self.populateFile()
        if self.reopen:
            # Close the file
            self.h5file.close()

    def populateFile(self):
        group = self.rootgroup
        obj = self.obj
        if obj is None:
            if self.type == "string":
                atom = tb.StringAtom(itemsize=self.length)
            else:
                atom = tb.Atom.from_type(self.type)
        else:
            atom = None
        title = self.__class__.__name__
        filters = tb.Filters(
            complevel=self.compress,
            complib=self.complib,
            shuffle=self.shuffle,
            bitshuffle=self.bitshuffle,
            fletcher32=self.fletcher32,
        )
        carray = self.h5file.create_carray(
            group,
            "carray1",
            atom=atom,
            shape=self.shape,
            title=title,
            filters=filters,
            chunkshape=self.chunkshape,
            byteorder=self.byteorder,
            obj=obj,
        )
        carray.flavor = self.flavor

        # Fill it with data
        self.rowshape = list(carray.shape)
        self.objsize = self.length * np.prod(carray.shape)

        if self.flavor == "numpy":
            if self.type == "string":
                object = np.ndarray(
                    buffer=b"a" * self.objsize,
                    shape=self.shape,
                    dtype="S%s" % carray.atom.itemsize,
                )
            else:
                object = np.arange(self.objsize, dtype=carray.atom.dtype)
                object.shape = carray.shape
        if common.verbose:
            print("Object to append -->", repr(object))

        carray[...] = object

    def _get_shape(self):
        if self.shape is not None:
            shape = self.shape
        else:
            shape = np.asarray(self.obj).shape

        return shape

    def test00_attributes(self):
        if self.reopen:
            self.h5file = tb.open_file(self.h5fname, "r")
        obj = self.h5file.get_node("/carray1")

        shape = self._get_shape()

        self.assertEqual(obj.flavor, self.flavor)
        self.assertEqual(obj.shape, shape)
        self.assertEqual(obj.ndim, len(shape))
        self.assertEqual(obj.chunkshape, self.chunkshape)
        self.assertEqual(obj.nrows, shape[0])
        self.assertEqual(obj.atom.type, self.type)

    def test01_readCArray(self):
        """Checking read() of chunked layout arrays."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01_readCArray..." % self.__class__.__name__)

        # Create an instance of an HDF5 Table
        if self.reopen:
            self.h5file = tb.open_file(self.h5fname, "r")
        carray = self.h5file.get_node("/carray1")

        # Choose a small value for buffer size
        carray.nrowsinbuf = 3
        if common.verbose:
            print("CArray descr:", repr(carray))
            print("shape of read array ==>", carray.shape)
            print("reopening?:", self.reopen)

        shape = self._get_shape()

        # Build the array to do comparisons
        if self.flavor == "numpy":
            if self.type == "string":
                object_ = np.ndarray(
                    buffer=b"a" * self.objsize,
                    shape=self.shape,
                    dtype=f"S{carray.atom.itemsize}",
                )
            else:
                object_ = np.arange(self.objsize, dtype=carray.atom.dtype)
                object_.shape = shape

        stop = self.stop
        # stop == None means read only the element designed by start
        # (in read() contexts)
        if self.stop is None:
            if self.start == -1:  # corner case
                stop = carray.nrows
            else:
                stop = self.start + 1
        # Protection against number of elements less than existing
        # if rowshape[self.extdim] < self.stop or self.stop == 0:
        if carray.nrows < stop:
            # self.stop == 0 means last row only in read()
            # and not in [::] slicing notation
            stop = int(carray.nrows)
        # do a copy() in order to ensure that len(object._data)
        # actually do a measure of its length
        obj = object_[self.start : stop : self.step].copy()

        # Read all the array
        try:
            data = carray.read(self.start, stop, self.step)
        except IndexError:
            if self.flavor == "numpy":
                data = np.empty(shape=self.shape, dtype=self.type)
            else:
                data = np.empty(shape=self.shape, dtype=self.type)

        if common.verbose:
            if hasattr(obj, "shape"):
                print("shape should look as:", obj.shape)
            print("Object read ==>", repr(data))
            print("Should look like ==>", repr(obj))

        if hasattr(data, "shape"):
            self.assertEqual(len(data.shape), len(shape))
        else:
            # Scalar case
            self.assertEqual(len(self.shape), 1)
        self.assertEqual(carray.chunkshape, self.chunkshape)
        self.assertTrue(common.allequal(data, obj, self.flavor))

    def test01_readCArray_out_argument(self):
        """Checking read() of chunked layout arrays."""

        # Create an instance of an HDF5 Table
        if self.reopen:
            self.h5file = tb.open_file(self.h5fname, "r")
        carray = self.h5file.get_node("/carray1")

        shape = self._get_shape()

        # Choose a small value for buffer size
        carray.nrowsinbuf = 3
        # Build the array to do comparisons
        if self.flavor == "numpy":
            if self.type == "string":
                object_ = np.ndarray(
                    buffer=b"a" * self.objsize,
                    shape=self.shape,
                    dtype=f"S{carray.atom.itemsize}",
                )
            else:
                object_ = np.arange(self.objsize, dtype=carray.atom.dtype)
                object_.shape = shape

        stop = self.stop
        # stop == None means read only the element designed by start
        # (in read() contexts)
        if self.stop is None:
            if self.start == -1:  # corner case
                stop = carray.nrows
            else:
                stop = self.start + 1
        # Protection against number of elements less than existing
        # if rowshape[self.extdim] < self.stop or self.stop == 0:
        if carray.nrows < stop:
            # self.stop == 0 means last row only in read()
            # and not in [::] slicing notation
            stop = int(carray.nrows)
        # do a copy() in order to ensure that len(object._data)
        # actually do a measure of its length
        obj = object_[self.start : stop : self.step].copy()

        # Read all the array
        try:
            data = np.empty(shape, dtype=carray.atom.dtype)
            data = data[self.start : stop : self.step].copy()
            carray.read(self.start, stop, self.step, out=data)
        except IndexError:
            if self.flavor == "numpy":
                data = np.empty(shape=shape, dtype=self.type)
            else:
                data = np.empty(shape=shape, dtype=self.type)

        if hasattr(data, "shape"):
            self.assertEqual(len(data.shape), len(shape))
        else:
            # Scalar case
            self.assertEqual(len(shape), 1)
        self.assertEqual(carray.chunkshape, self.chunkshape)
        self.assertTrue(common.allequal(data, obj, self.flavor))

    def test02_getitemCArray(self):
        """Checking chunked layout array __getitem__ special method."""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test02_getitemCArray..." % self.__class__.__name__
            )

        if not hasattr(self, "slices"):
            # If there is not a slices attribute, create it
            self.slices = (slice(self.start, self.stop, self.step),)

        # Create an instance of an HDF5 Table
        if self.reopen:
            self.h5file = tb.open_file(self.h5fname, "r")
        carray = self.h5file.get_node("/carray1")

        if common.verbose:
            print("CArray descr:", repr(carray))
            print("shape of read array ==>", carray.shape)
            print("reopening?:", self.reopen)

        shape = self._get_shape()

        # Build the array to do comparisons
        if self.type == "string":
            object_ = np.ndarray(
                buffer=b"a" * self.objsize,
                shape=self.shape,
                dtype=f"S{carray.atom.itemsize}",
            )
        else:
            object_ = np.arange(self.objsize, dtype=carray.atom.dtype)
            object_.shape = shape

        # do a copy() in order to ensure that len(object._data)
        # actually do a measure of its length
        obj = object_.__getitem__(self.slices).copy()

        # Read data from the array
        try:
            data = carray.__getitem__(self.slices)
        except IndexError:
            print("IndexError!")
            if self.flavor == "numpy":
                data = np.empty(shape=self.shape, dtype=self.type)
            else:
                data = np.empty(shape=self.shape, dtype=self.type)

        if common.verbose:
            print("Object read:\n", repr(data))  # , data.info()
            print("Should look like:\n", repr(obj))  # , object.info()
            if hasattr(obj, "shape"):
                print("Original object shape:", self.shape)
                print("Shape read:", data.shape)
                print("shape should look as:", obj.shape)

        if not hasattr(data, "shape"):
            # Scalar case
            self.assertEqual(len(self.shape), 1)
        self.assertEqual(carray.chunkshape, self.chunkshape)
        self.assertTrue(common.allequal(data, obj, self.flavor))

    def test03_setitemCArray(self):
        """Checking chunked layout array __setitem__ special method."""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test03_setitemCArray..." % self.__class__.__name__
            )

        if not hasattr(self, "slices"):
            # If there is not a slices attribute, create it
            self.slices = (slice(self.start, self.stop, self.step),)

        # Create an instance of an HDF5 Table
        if self.reopen:
            self.h5file = tb.open_file(self.h5fname, "a")
        carray = self.h5file.get_node("/carray1")

        if common.verbose:
            print("CArray descr:", repr(carray))
            print("shape of read array ==>", carray.shape)
            print("reopening?:", self.reopen)

        shape = self._get_shape()

        # Build the array to do comparisons
        if self.type == "string":
            object_ = np.ndarray(
                buffer=b"a" * self.objsize,
                shape=self.shape,
                dtype=f"S{carray.atom.itemsize}",
            )
        else:
            object_ = np.arange(self.objsize, dtype=carray.atom.dtype)
            object_.shape = shape

        # do a copy() in order to ensure that len(object._data)
        # actually do a measure of its length
        obj = object_.__getitem__(self.slices).copy()

        if self.type == "string":
            if hasattr(self, "wslice"):
                obj[self.wslize] = "xXx"
                carray[self.wslice] = "xXx"
            elif sum(obj[self.slices].shape) != 0:
                obj[:] = "xXx"
                if obj.size > 0:
                    carray[self.slices] = obj
        else:
            if hasattr(self, "wslice"):
                obj[self.wslice] = obj[self.wslice] * 2 + 3
                carray[self.wslice] = carray[self.wslice] * 2 + 3
            elif sum(obj[self.slices].shape) != 0:
                obj = obj * 2 + 3
                if np.prod(obj.shape) > 0:
                    carray[self.slices] = carray[self.slices] * 2 + 3
            # Cast again object to its original type
            obj = np.array(obj, dtype=carray.atom.dtype)
        # Read datafrom the array
        try:
            data = carray.__getitem__(self.slices)
        except IndexError:
            print("IndexError!")
            if self.flavor == "numpy":
                data = np.empty(shape=self.shape, dtype=self.type)
            else:
                data = np.empty(shape=self.shape, dtype=self.type)

        if common.verbose:
            print("Object read:\n", repr(data))  # , data.info()
            print("Should look like:\n", repr(obj))  # , object.info()
            if hasattr(obj, "shape"):
                print("Original object shape:", self.shape)
                print("Shape read:", data.shape)
                print("shape should look as:", obj.shape)

        if not hasattr(data, "shape"):
            # Scalar case
            self.assertEqual(len(self.shape), 1)
        self.assertEqual(carray.chunkshape, self.chunkshape)
        self.assertTrue(common.allequal(data, obj, self.flavor))


class BasicWriteTestCase(BasicTestCase):
    type = "int32"
    shape = (2,)
    chunkshape = (5,)
    step = 1
    wslice = 1  # single element case


class BasicWrite2TestCase(BasicTestCase):
    type = "int32"
    shape = (2,)
    chunkshape = (5,)
    step = 1
    wslice = slice(shape[0] - 2, shape[0], 2)  # range of elements
    reopen = 0  # This case does not reopen files


class BasicWrite3TestCase(BasicTestCase):
    obj = [1, 2]
    type = np.asarray(obj).dtype.name
    shape = None
    chunkshape = (5,)
    step = 1
    reopen = 0  # This case does not reopen files


class BasicWrite4TestCase(BasicTestCase):
    obj = np.array([1, 2])
    type = obj.dtype.name
    shape = None
    chunkshape = (5,)
    step = 1
    reopen = 0  # This case does not reopen files


class BasicWrite5TestCase(BasicTestCase):
    obj = [[1, 2], [3, 4]]
    type = np.asarray(obj).dtype.name
    shape = None
    chunkshape = (5, 1)
    step = 1
    reopen = 0  # This case does not reopen files


class BasicWrite6TestCase(BasicTestCase):
    obj = [1, 2]
    type = np.asarray(obj).dtype.name
    shape = None
    chunkshape = (5,)
    step = 1
    reopen = 1  # This case does reopen files


class BasicWrite7TestCase(BasicTestCase):
    obj = np.array([1, 2])
    type = obj.dtype.name
    shape = None
    chunkshape = (5,)
    step = 1
    reopen = 1  # This case does reopen files


class BasicWrite8TestCase(BasicTestCase):
    obj = [[1, 2], [3, 4]]
    type = np.asarray(obj).dtype.name
    shape = None
    chunkshape = (5, 1)
    step = 1
    reopen = 1  # This case does reopen files


class EmptyCArrayTestCase(BasicTestCase):
    type = "int32"
    shape = (2, 2)
    chunkshape = (5, 5)
    start = 0
    stop = 10
    step = 1


class EmptyCArray2TestCase(BasicTestCase):
    type = "int32"
    shape = (2, 2)
    chunkshape = (5, 5)
    start = 0
    stop = 10
    step = 1
    reopen = 0  # This case does not reopen files


@common.unittest.skipIf(
    not common.lzo_avail, "LZO compression library not available"
)
class SlicesCArrayTestCase(BasicTestCase):
    compress = 1
    complib = "lzo"
    type = "int32"
    shape = (2, 2)
    chunkshape = (5, 5)
    slices = (slice(1, 2, 1), slice(1, 3, 1))


class EllipsisCArrayTestCase(BasicTestCase):
    type = "int32"
    shape = (2, 2)
    chunkshape = (5, 5)
    # slices = (slice(1,2,1), Ellipsis)
    slices = (Ellipsis, slice(1, 2, 1))


@common.unittest.skipIf(
    not common.lzo_avail, "LZO compression library not available"
)
class Slices2CArrayTestCase(BasicTestCase):
    compress = 1
    complib = "lzo"
    type = "int32"
    shape = (2, 2, 4)
    chunkshape = (5, 5, 5)
    slices = (slice(1, 2, 1), slice(None, None, None), slice(1, 4, 2))


class Ellipsis2CArrayTestCase(BasicTestCase):
    type = "int32"
    shape = (2, 2, 4)
    chunkshape = (5, 5, 5)
    slices = (slice(1, 2, 1), Ellipsis, slice(1, 4, 2))


@common.unittest.skipIf(
    not common.lzo_avail, "LZO compression library not available"
)
class Slices3CArrayTestCase(BasicTestCase):
    compress = 1  # To show the chunks id DEBUG is on
    complib = "lzo"
    type = "int32"
    shape = (2, 3, 4, 2)
    chunkshape = (5, 5, 5, 5)
    slices = (
        slice(1, 2, 1),
        slice(0, None, None),
        slice(1, 4, 2),
    )  # Don't work
    # slices = (slice(None, None, None), slice(0, None, None),
    #           slice(1,4,1))  # W
    # slices = (slice(None, None, None), slice(None, None, None),
    #           slice(1,4,2))  # N
    # slices = (slice(1,2,1), slice(None, None, None), slice(1,4,2))  # N
    # Disable the failing test temporarily with a working test case
    slices = (slice(1, 2, 1), slice(1, 4, None), slice(1, 4, 2))  # Y
    # slices = (slice(1,2,1), slice(0, 4, None), slice(1,4,1))  # Y
    slices = (slice(1, 2, 1), slice(0, 4, None), slice(1, 4, 2))  # N
    # slices = (slice(1,2,1), slice(0, 4, None), slice(1,4,2),
    #           slice(0,100,1))  # N


class Slices4CArrayTestCase(BasicTestCase):
    type = "int32"
    shape = (2, 3, 4, 2, 5, 6)
    chunkshape = (5, 5, 5, 5, 5, 5)
    slices = (
        slice(1, 2, 1),
        slice(0, None, None),
        slice(1, 4, 2),
        slice(0, 4, 2),
        slice(3, 5, 2),
        slice(2, 7, 1),
    )


class Ellipsis3CArrayTestCase(BasicTestCase):
    type = "int32"
    shape = (2, 3, 4, 2)
    chunkshape = (5, 5, 5, 5)
    slices = (Ellipsis, slice(0, 4, None), slice(1, 4, 2))
    slices = (slice(1, 2, 1), slice(0, 4, None), slice(1, 4, 2), Ellipsis)


class Ellipsis4CArrayTestCase(BasicTestCase):
    type = "int32"
    shape = (2, 3, 4, 5)
    chunkshape = (5, 5, 5, 5)
    slices = (Ellipsis, slice(0, 4, None), slice(1, 4, 2))
    slices = (slice(1, 2, 1), Ellipsis, slice(1, 4, 2))


class Ellipsis5CArrayTestCase(BasicTestCase):
    type = "int32"
    shape = (2, 3, 4, 5)
    chunkshape = (5, 5, 5, 5)
    slices = (slice(1, 2, 1), slice(0, 4, None), Ellipsis)


class Ellipsis6CArrayTestCase(BasicTestCase):
    type = "int32"
    shape = (2, 3, 4, 5)
    chunkshape = (5, 5, 5, 5)
    # The next slices gives problems with setting values (test03)
    # This is a problem on the test design, not the Array.__setitem__
    # code, though. See # see test_earray.py Ellipsis6EArrayTestCase
    slices = (slice(1, 2, 1), slice(0, 4, None), 2, Ellipsis)


class Ellipsis7CArrayTestCase(BasicTestCase):
    type = "int32"
    shape = (2, 3, 4, 5)
    chunkshape = (5, 5, 5, 5)
    slices = (slice(1, 2, 1), slice(0, 4, None), slice(2, 3), Ellipsis)


class MD3WriteTestCase(BasicTestCase):
    type = "int32"
    shape = (2, 2, 3)
    chunkshape = (4, 4, 4)
    step = 2


class MD5WriteTestCase(BasicTestCase):
    type = "int32"
    shape = (2, 2, 3, 4, 5)  # ok
    # shape = (1, 1, 2, 1)  # Minimum shape that shows problems with HDF5 1.6.1
    # shape = (2, 3, 2, 4, 5)  # Floating point exception (HDF5 1.6.1)
    # shape = (2, 3, 3, 2, 5, 6) # Segmentation fault (HDF5 1.6.1)
    chunkshape = (1, 1, 1, 1, 1)
    start = 1
    stop = 10
    step = 10


class MD6WriteTestCase(BasicTestCase):
    type = "int32"
    shape = (2, 3, 3, 2, 5, 6)
    chunkshape = (1, 1, 1, 1, 5, 6)
    start = 1
    stop = 10
    step = 3


class MD6WriteTestCase__(BasicTestCase):
    type = "int32"
    shape = (2, 2)
    chunkshape = (1, 1)
    start = 1
    stop = 3
    step = 1


class MD7WriteTestCase(BasicTestCase):
    type = "int32"
    shape = (2, 3, 3, 4, 5, 2, 3)
    chunkshape = (10, 10, 10, 10, 10, 10, 10)
    start = 1
    stop = 10
    step = 2


class MD10WriteTestCase(BasicTestCase):
    type = "int32"
    shape = (1, 2, 3, 4, 5, 5, 4, 3, 2, 2)
    chunkshape = (5, 5, 5, 5, 5, 5, 5, 5, 5, 5)
    start = -1
    stop = -1
    step = 10


class ZlibComprTestCase(BasicTestCase):
    compress = 1
    complib = "zlib"
    start = 3
    # stop = 0   # means last row
    stop = None  # means last row from 0.8 on
    step = 10


class ZlibShuffleTestCase(BasicTestCase):
    shuffle = 1
    compress = 1
    complib = "zlib"
    # case start < stop , i.e. no rows read
    start = 3
    stop = 1
    step = 10


@common.unittest.skipIf(
    not common.blosc_avail, "BLOSC compression library not available"
)
class BloscComprTestCase(BasicTestCase):
    compress = 1  # sss
    complib = "blosc"
    chunkshape = (10, 10)
    start = 3
    stop = 10
    step = 3


@common.unittest.skipIf(
    not common.blosc_avail, "BLOSC compression library not available"
)
class BloscShuffleTestCase(BasicTestCase):
    shape = (20, 30)
    compress = 1
    shuffle = 1
    complib = "blosc"
    chunkshape = (100, 100)
    start = 3
    stop = 10
    step = 7


@common.unittest.skipIf(
    not common.blosc_avail, "BLOSC compression library not available"
)
class BloscBitShuffleTestCase(BasicTestCase):
    shape = (20, 30)
    compress = 1
    bitshuffle = 1
    complib = "blosc"
    chunkshape = (200, 100)
    start = 2
    stop = 11
    step = 7


@common.unittest.skipIf(
    not common.blosc_avail, "BLOSC compression library not available"
)
class BloscFletcherTestCase(BasicTestCase):
    # see gh-21
    shape = (200, 300)
    compress = 1
    shuffle = 1
    fletcher32 = 1
    complib = "blosc"
    chunkshape = (100, 100)
    start = 3
    stop = 10
    step = 7


@common.unittest.skipIf(
    not common.blosc_avail, "BLOSC compression library not available"
)
class BloscBloscLZTestCase(BasicTestCase):
    shape = (20, 30)
    compress = 1
    shuffle = 1
    complib = "blosc:blosclz"
    chunkshape = (200, 100)
    start = 2
    stop = 11
    step = 7


@common.unittest.skipIf(
    not common.blosc_avail, "BLOSC compression library not available"
)
@common.unittest.skipIf(
    "lz4" not in tb.blosc_compressor_list(), "lz4 required"
)
class BloscLZ4TestCase(BasicTestCase):
    shape = (20, 30)
    compress = 1
    shuffle = 1
    complib = "blosc:lz4"
    chunkshape = (100, 100)
    start = 3
    stop = 10
    step = 7


@common.unittest.skipIf(
    not common.blosc_avail, "BLOSC compression library not available"
)
@common.unittest.skipIf(
    "lz4" not in tb.blosc_compressor_list(), "lz4 required"
)
class BloscLZ4HCTestCase(BasicTestCase):
    shape = (20, 30)
    compress = 1
    shuffle = 1
    complib = "blosc:lz4hc"
    chunkshape = (100, 100)
    start = 3
    stop = 10
    step = 7


@common.unittest.skipIf(
    not common.blosc_avail, "BLOSC compression library not available"
)
@common.unittest.skipIf(
    "snappy" not in tb.blosc_compressor_list(), "snappy required"
)
class BloscSnappyTestCase(BasicTestCase):
    shape = (20, 30)
    compress = 1
    shuffle = 1
    complib = "blosc:snappy"
    chunkshape = (100, 100)
    start = 3
    stop = 10
    step = 7


@common.unittest.skipIf(
    not common.blosc_avail, "BLOSC compression library not available"
)
@common.unittest.skipIf(
    "zlib" not in tb.blosc_compressor_list(), "zlib required"
)
class BloscZlibTestCase(BasicTestCase):
    shape = (20, 30)
    compress = 1
    shuffle = 1
    complib = "blosc:zlib"
    chunkshape = (100, 100)
    start = 3
    stop = 10
    step = 7


@common.unittest.skipIf(
    not common.blosc_avail, "BLOSC compression library not available"
)
@common.unittest.skipIf(
    "zstd" not in tb.blosc_compressor_list(), "zstd required"
)
class BloscZstdTestCase(BasicTestCase):
    shape = (20, 30)
    compress = 1
    shuffle = 1
    complib = "blosc:zstd"
    chunkshape = (100, 100)
    start = 3
    stop = 10
    step = 7


@common.unittest.skipIf(
    not common.blosc2_avail, "BLOSC2 compression library not available"
)
class Blosc2ComprTestCase(BasicTestCase):
    compress = 1  # sss
    complib = "blosc2"
    chunkshape = (10, 10)
    start = 3
    stop = 10
    step = 3
    byteorder = foreign_byteorder()


@common.unittest.skipIf(
    not common.blosc2_avail, "BLOSC2 compression library not available"
)
class Blosc2FletcherTestCase(Blosc2ComprTestCase):
    fletcher32 = 1
    start = 0


@common.unittest.skipIf(
    not common.blosc2_avail, "BLOSC2 compression library not available"
)
class Blosc2CrossChunkTestCase(BasicTestCase):
    shape = (10, 10)
    compress = 1  # sss
    complib = "blosc2"
    chunkshape = (4, 4)
    start = 3
    stop = 6
    step = 3
    byteorder = foreign_byteorder()


@common.unittest.skipIf(
    not common.blosc2_avail, "BLOSC2 compression library not available"
)
class Blosc2CrossChunkOptTestCase(Blosc2CrossChunkTestCase):
    step = 1  # optimized
    byteorder = sys.byteorder


@common.unittest.skipIf(
    not common.blosc2_avail, "BLOSC2 compression library not available"
)
class Blosc2PastLastChunkTestCase(BasicTestCase):
    shape = (10, 10)
    compress = 1  # sss
    complib = "blosc2"
    chunkshape = (4, 4)
    start = 8
    stop = 100
    step = 3
    byteorder = foreign_byteorder()


@common.unittest.skipIf(
    not common.blosc2_avail, "BLOSC2 compression library not available"
)
class Blosc2PastLastChunkOptTestCase(Blosc2PastLastChunkTestCase):
    step = 1  # optimized
    byteorder = sys.byteorder


# Minimal test which can be figured out manually::
#
#     z  Data: 1   Chunk0:   Chunk1: 1   Slice:
#    /        /|\                    |\
#   |\       0 5 3       0           5 3        5
#   x y      |X X|       |\           \|       / \
#            4 2 7       4 2           7      4   7
#             \|/         \|                   \ /
#              6           6                    6
#
#                  Chunk0 & Slice: 4   Chunk1 & Slice: 5
#                                   \                   \
#                                    6                   7
@common.unittest.skipIf(
    not common.blosc2_avail, "BLOSC2 compression library not available"
)
class Blosc2Ndim3MinChunkOptTestCase(BasicTestCase):
    shape = (2, 2, 2)
    compress = 1
    complib = "blosc2"
    chunkshape = (2, 2, 1)
    byteorder = sys.byteorder
    type = "int8"
    slices = (slice(1, 2), slice(0, 2), slice(0, 2))


@common.unittest.skipIf(
    not common.blosc2_avail, "BLOSC2 compression library not available"
)
class Blosc2Ndim3ChunkOptTestCase(BasicTestCase):
    shape = (10, 10, 10)
    compress = 1
    complib = "blosc2"
    chunkshape = (7, 7, 7)
    byteorder = sys.byteorder
    type = "int32"
    slices = (slice(1, 2), Ellipsis, slice(1, 4))


@common.unittest.skipIf(
    not common.blosc2_avail, "BLOSC2 compression library not available"
)
class Blosc2Ndim4ChunkOptTestCase(BasicTestCase):
    shape = (13, 13, 13, 3)
    compress = 1
    complib = "blosc2"
    chunkshape = (5, 5, 5, 3)
    byteorder = sys.byteorder
    type = "int32"
    slices = (slice(0, 8), slice(7, 13), slice(3, 12), slice(1, 3))


# The file used in the test below is created with this script,
# producing a chunked array that lacks chunk rank/shape in filter args.
# It is a reduced version of ``examples/direct-chunk-shape.py``,
# check there for more info and the assemblage of the data array.
# An h5py release is used which contains a version of hdf5-blosc2
# that does not include chunk rank/shape in filter arguments.
#
# ::
#
#   import blosc2
#   import h5py
#   import hdf5plugin
#   import numpy
#
#   assert(hdf5plugin.version_info < (4, 2, 1))
#
#   fparams = hdf5plugin.Blosc2(cname='zstd', clevel=1,
#                               filters=hdf5plugin.Blosc2.SHUFFLE)
#   cparams = {
#       "codec": blosc2.Codec.ZSTD,
#       "clevel": 1,
#       "filters": [blosc2.Filter.SHUFFLE],
#   }
#
#   achunk = numpy.arange(4 * 4, dtype='int8').reshape((4, 4))
#   adata = numpy.zeros((6, 6), dtype=achunk.dtype)
#   adata[0:4, 0:4] = achunk[:, :]
#   adata[0:4, 4:6] = achunk[:, 0:2]
#   adata[4:6, 0:4] = achunk[0:2, :]
#   adata[4:6, 4:6] = achunk[0:2, 0:2]
#
#   h5f = h5py.File("b2nd-no-chunkshape.h5", "w")
#   dataset = h5f.create_dataset(
#       "data", adata.shape, dtype=adata.dtype, chunks=achunk.shape,
#       **fparams)
#   b2chunk = blosc2.asarray(achunk,
#                            chunks=achunk.shape, blocks=achunk.shape,
#                            cparams=cparams)
#   b2frame = b2chunk._schunk.to_cframe()
#   dataset.id.write_direct_chunk((0, 0), b2frame)
#   dataset.id.write_direct_chunk((0, 4), b2frame)
#   dataset.id.write_direct_chunk((4, 0), b2frame)
#   dataset.id.write_direct_chunk((4, 4), b2frame)
#   h5f.close()
@common.unittest.skipIf(
    not common.blosc2_avail, "BLOSC2 compression library not available"
)
class Blosc2NDNoChunkshape(common.TestFileMixin, common.PyTablesTestCase):
    h5fname = common.test_filename("b2nd-no-chunkshape.h5")

    adata = np.array(
        [
            [0, 1, 2, 3, 0, 1],
            [4, 5, 6, 7, 4, 5],
            [8, 9, 10, 11, 8, 9],
            [12, 13, 14, 15, 12, 13],
            [0, 1, 2, 3, 0, 1],
            [4, 5, 6, 7, 4, 5],
        ],
        dtype="int8",
    )

    def test_data_opt(self):
        array = self.h5file.get_node("/data")
        self.assertTrue(common.areArraysEqual(array[:], self.adata[:]))

    def test_data_filter(self):
        array = self.h5file.get_node("/data")
        self.assertTrue(common.areArraysEqual(array[::2], self.adata[::2]))


@common.unittest.skipIf(
    not common.lzo_avail, "LZO compression library not available"
)
class LZOComprTestCase(BasicTestCase):
    compress = 1  # sss
    complib = "lzo"
    chunkshape = (10, 10)
    start = 3
    stop = 10
    step = 3


@common.unittest.skipIf(
    not common.lzo_avail, "LZO compression library not available"
)
class LZOShuffleTestCase(BasicTestCase):
    shape = (20, 30)
    compress = 1
    shuffle = 1
    complib = "lzo"
    chunkshape = (100, 100)
    start = 3
    stop = 10
    step = 7


@common.unittest.skipIf(
    not common.bzip2_avail, "BZIP2 compression library not available"
)
class Bzip2ComprTestCase(BasicTestCase):
    shape = (20, 30)
    compress = 1
    complib = "bzip2"
    chunkshape = (100, 100)
    start = 3
    stop = 10
    step = 8


@common.unittest.skipIf(
    not common.bzip2_avail, "BZIP2 compression library not available"
)
class Bzip2ShuffleTestCase(BasicTestCase):
    shape = (20, 30)
    compress = 1
    shuffle = 1
    complib = "bzip2"
    chunkshape = (100, 100)
    start = 3
    stop = 10
    step = 6


class Fletcher32TestCase(BasicTestCase):
    shape = (60, 50)
    compress = 0
    fletcher32 = 1
    chunkshape = (50, 50)
    start = 4
    stop = 20
    step = 7


class AllFiltersTestCase(BasicTestCase):
    compress = 1
    shuffle = 1
    fletcher32 = 1
    complib = "zlib"
    chunkshape = (20, 20)  # sss
    start = 2
    stop = 99
    step = 6


class FloatTypeTestCase(BasicTestCase):
    type = "float64"
    shape = (2, 2)
    chunkshape = (5, 5)
    start = 3
    stop = 10
    step = 20


class ComplexTypeTestCase(BasicTestCase):
    type = "complex128"
    shape = (2, 2)
    chunkshape = (5, 5)
    start = 3
    stop = 10
    step = 20


class StringTestCase(BasicTestCase):
    type = "string"
    length = 20
    shape = (2, 2)
    # shape = (2,2,20)
    chunkshape = (5, 5)
    start = 3
    stop = 10
    step = 20
    slices = (slice(0, 1), slice(1, 2))


class String2TestCase(BasicTestCase):
    type = "string"
    length = 20
    shape = (2, 20)
    chunkshape = (5, 5)
    start = 1
    stop = 10
    step = 2


class StringComprTestCase(BasicTestCase):
    type = "string"
    length = 20
    shape = (20, 2, 10)
    # shape = (20,0,10,20)
    compr = 1
    # shuffle = 1  # this shouldn't do nothing on chars
    chunkshape = (50, 50, 2)
    start = -1
    stop = 100
    step = 20


class Int8TestCase(BasicTestCase):
    type = "int8"
    shape = (2, 2)
    compress = 1
    shuffle = 1
    chunkshape = (50, 50)
    start = -1
    stop = 100
    step = 20


class Int16TestCase(BasicTestCase):
    type = "int16"
    shape = (2, 2)
    compress = 1
    shuffle = 1
    chunkshape = (50, 50)
    start = 1
    stop = 100
    step = 1


class Int32TestCase(BasicTestCase):
    type = "int32"
    shape = (2, 2)
    compress = 1
    shuffle = 1
    chunkshape = (50, 50)
    start = -1
    stop = 100
    step = 20


@common.unittest.skipUnless(
    hasattr(tb, "Float16Atom"), "Float16Atom not available"
)
class Float16TestCase(BasicTestCase):
    type = "float16"
    shape = (200,)
    compress = 1
    shuffle = 1
    chunkshape = (20,)
    start = -1
    stop = 100
    step = 20


class Float32TestCase(BasicTestCase):
    type = "float32"
    shape = (200,)
    compress = 1
    shuffle = 1
    chunkshape = (20,)
    start = -1
    stop = 100
    step = 20


class Float64TestCase(BasicTestCase):
    type = "float64"
    shape = (200,)
    compress = 1
    shuffle = 1
    chunkshape = (20,)
    start = -1
    stop = 100
    step = 20


@common.unittest.skipUnless(
    hasattr(tb, "Float96Atom"), "Float96Atom not available"
)
class Float96TestCase(BasicTestCase):
    type = "float96"
    shape = (200,)
    compress = 1
    shuffle = 1
    chunkshape = (20,)
    start = -1
    stop = 100
    step = 20


@common.unittest.skipUnless(
    hasattr(tb, "Float128Atom"), "Float128Atom not available"
)
class Float128TestCase(BasicTestCase):
    type = "float128"
    shape = (200,)
    compress = 1
    shuffle = 1
    chunkshape = (20,)
    start = -1
    stop = 100
    step = 20


class Complex64TestCase(BasicTestCase):
    type = "complex64"
    shape = (4,)
    compress = 1
    shuffle = 1
    chunkshape = (2,)
    start = -1
    stop = 100
    step = 20


class Complex128TestCase(BasicTestCase):
    type = "complex128"
    shape = (20,)
    compress = 1
    shuffle = 1
    chunkshape = (2,)
    start = -1
    stop = 100
    step = 20


@common.unittest.skipUnless(
    hasattr(tb, "Complex192Atom"), "Complex192Atom not available"
)
class Complex192TestCase(BasicTestCase):
    type = "complex192"
    shape = (20,)
    compress = 1
    shuffle = 1
    chunkshape = (2,)
    start = -1
    stop = 100
    step = 20


@common.unittest.skipUnless(
    hasattr(tb, "Complex256Atom"), "Complex256Atom not available"
)
class Complex256TestCase(BasicTestCase):
    type = "complex256"
    shape = (20,)
    compress = 1
    shuffle = 1
    chunkshape = (2,)
    start = -1
    stop = 100
    step = 20


class ComprTestCase(BasicTestCase):
    type = "float64"
    compress = 1
    shuffle = 1
    shape = (200,)
    compr = 1
    chunkshape = (21,)
    start = 51
    stop = 100
    step = 7


# this is a subset of the tests in test_array.py, mostly to verify that errors
# are handled in the same way
class ReadOutArgumentTests(common.TempFileMixin, common.PyTablesTestCase):

    def setUp(self):
        super().setUp()
        self.size = 1000
        self.filters = tb.Filters(complevel=1, complib="blosc")

    def create_array(self):
        array = np.arange(self.size, dtype="i8")
        disk_array = self.h5file.create_carray(
            "/",
            "array",
            atom=tb.Int64Atom(),
            shape=(self.size,),
            filters=self.filters,
        )
        disk_array[:] = array
        return array, disk_array

    def test_read_entire_array(self):
        array, disk_array = self.create_array()
        out_buffer = np.empty((self.size,), "i8")
        disk_array.read(out=out_buffer)
        np.testing.assert_equal(out_buffer, array)

    def test_read_non_contiguous_buffer(self):
        array, disk_array = self.create_array()
        out_buffer = np.empty((self.size,), "i8")
        out_buffer_slice = out_buffer[0 : self.size : 2]

        with self.assertRaisesRegex(
            ValueError, "output array not C contiguous"
        ):
            disk_array.read(0, self.size, 2, out_buffer_slice)

    def test_buffer_too_small(self):
        array, disk_array = self.create_array()
        out_buffer = np.empty((self.size // 2,), "i8")
        self.assertRaises(
            ValueError, disk_array.read, 0, self.size, 1, out_buffer
        )
        try:
            disk_array.read(0, self.size, 1, out_buffer)
        except ValueError as exc:
            self.assertIn("output array size invalid, got", str(exc))


class SizeOnDiskInMemoryPropertyTestCase(
    common.TempFileMixin, common.PyTablesTestCase
):

    def setUp(self):
        super().setUp()
        self.array_size = (10_000, 10)
        # set chunkshape so it divides evenly into array_size, to avoid
        # partially filled chunks
        self.chunkshape = (1000, 10)
        # approximate size (in bytes) of non-data portion of hdf5 file
        self.hdf_overhead = 6000

    def create_array(self, complevel):
        filters = tb.Filters(complevel=complevel, complib="blosc")
        self.array = self.h5file.create_carray(
            "/",
            "somearray",
            atom=tb.Int16Atom(),
            shape=self.array_size,
            filters=filters,
            chunkshape=self.chunkshape,
        )

    def test_no_data(self):
        complevel = 0
        self.create_array(complevel)
        self.assertEqual(self.array.size_on_disk, 0)
        self.assertEqual(self.array.size_in_memory, 10_000 * 10 * 2)

    def test_data_no_compression(self):
        complevel = 0
        self.create_array(complevel)
        self.array[:] = 1
        self.assertEqual(self.array.size_on_disk, 10_000 * 10 * 2)
        self.assertEqual(self.array.size_in_memory, 10_000 * 10 * 2)

    def test_highly_compressible_data(self):
        complevel = 1
        self.create_array(complevel)
        self.array[:] = 1
        self.h5file.flush()
        file_size = Path(self.h5fname).stat().st_size
        self.assertTrue(
            abs(self.array.size_on_disk - file_size) <= self.hdf_overhead
        )
        self.assertTrue(self.array.size_on_disk < self.array.size_in_memory)
        self.assertEqual(self.array.size_in_memory, 10_000 * 10 * 2)

    # XXX
    def test_random_data(self):
        complevel = 1
        self.create_array(complevel)
        self.array[:] = np.random.randint(0, 1e6, self.array_size)
        self.h5file.flush()
        file_size = Path(self.h5fname).stat().st_size
        self.assertTrue(
            abs(self.array.size_on_disk - file_size) <= self.hdf_overhead
        )

        # XXX: check. The test fails if blosc is not available
        if tb.which_lib_version("blosc") is not None:
            self.assertAlmostEqual(self.array.size_on_disk, 10_000 * 10 * 2)
        else:
            self.assertTrue(
                abs(self.array.size_on_disk - 10_000 * 10 * 2) < 200
            )


class OffsetStrideTestCase(common.TempFileMixin, common.PyTablesTestCase):
    compress = 0
    complib = "zlib"  # Default compression library

    def setUp(self):
        super().setUp()
        # Create an instance of an HDF5 Table
        self.rootgroup = self.h5file.root

    def test01a_String(self):
        """Checking carray with offset NumPy strings appends."""

        root = self.rootgroup
        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01a_String..." % self.__class__.__name__)

        shape = (3, 2, 2)
        # Create a string atom
        carray = self.h5file.create_carray(
            root,
            "strings",
            atom=tb.StringAtom(itemsize=3),
            shape=shape,
            title="Array of strings",
            chunkshape=(1, 2, 2),
        )
        a = np.array([[["a", "b"], ["123", "45"], ["45", "123"]]], dtype="S3")
        carray[0] = a[0, 1:]
        a = np.array([[["s", "a"], ["ab", "f"], ["s", "abc"], ["abc", "f"]]])
        carray[1] = a[0, 2:]

        # Read all the data:
        data = carray.read()
        if common.verbose:
            print("Object read:", data)
            print("Nrows in", carray._v_pathname, ":", carray.nrows)
            print("Second row in carray ==>", data[1].tolist())

        self.assertEqual(carray.nrows, 3)
        self.assertEqual(data[0].tolist(), [[b"123", b"45"], [b"45", b"123"]])
        self.assertEqual(data[1].tolist(), [[b"s", b"abc"], [b"abc", b"f"]])
        self.assertEqual(len(data[0]), 2)
        self.assertEqual(len(data[1]), 2)

    def test01b_String(self):
        """Checking carray with strided NumPy strings appends."""

        root = self.rootgroup
        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01b_String..." % self.__class__.__name__)

        shape = (3, 2, 2)

        # Create a string atom
        carray = self.h5file.create_carray(
            root,
            "strings",
            atom=tb.StringAtom(itemsize=3),
            shape=shape,
            title="Array of strings",
            chunkshape=(1, 2, 2),
        )
        a = np.array([[["a", "b"], ["123", "45"], ["45", "123"]]], dtype="S3")
        carray[0] = a[0, ::2]
        a = np.array([[["s", "a"], ["ab", "f"], ["s", "abc"], ["abc", "f"]]])
        carray[1] = a[0, ::2]

        # Read all the rows:
        data = carray.read()
        if common.verbose:
            print("Object read:", data)
            print("Nrows in", carray._v_pathname, ":", carray.nrows)
            print("Second row in carray ==>", data[1].tolist())

        self.assertEqual(carray.nrows, 3)
        self.assertEqual(data[0].tolist(), [[b"a", b"b"], [b"45", b"123"]])
        self.assertEqual(data[1].tolist(), [[b"s", b"a"], [b"s", b"abc"]])
        self.assertEqual(len(data[0]), 2)
        self.assertEqual(len(data[1]), 2)

    def test02a_int(self):
        """Checking carray with offset NumPy ints appends."""

        root = self.rootgroup
        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02a_int..." % self.__class__.__name__)

        shape = (3, 3)

        # Create a string atom
        carray = self.h5file.create_carray(
            root,
            "CAtom",
            atom=tb.Int32Atom(),
            shape=shape,
            title="array of ints",
            chunkshape=(1, 3),
        )
        a = np.array(
            [(0, 0, 0), (1, 0, 3), (1, 1, 1), (0, 0, 0)], dtype="int32"
        )
        carray[0:2] = a[2:]  # Introduce an offset
        a = np.array([(1, 1, 1), (-1, 0, 0)], dtype="int32")
        carray[2:3] = a[1:]  # Introduce an offset

        # Read all the rows:
        data = carray.read()
        if common.verbose:
            print("Object read:", data)
            print("Nrows in", carray._v_pathname, ":", carray.nrows)
            print("Third row in carray ==>", data[2])

        self.assertEqual(carray.nrows, 3)
        self.assertTrue(
            common.allequal(data[0], np.array([1, 1, 1], dtype="int32"))
        )
        self.assertTrue(
            common.allequal(data[1], np.array([0, 0, 0], dtype="int32"))
        )
        self.assertTrue(
            common.allequal(data[2], np.array([-1, 0, 0], dtype="int32"))
        )

    def test02b_int(self):
        """Checking carray with strided NumPy ints appends."""

        root = self.rootgroup
        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02b_int..." % self.__class__.__name__)

        shape = (3, 3)

        # Create a string atom
        carray = self.h5file.create_carray(
            root,
            "CAtom",
            atom=tb.Int32Atom(),
            shape=shape,
            title="array of ints",
            chunkshape=(1, 3),
        )
        a = np.array(
            [(0, 0, 0), (1, 0, 3), (1, 1, 1), (3, 3, 3)], dtype="int32"
        )
        carray[0:2] = a[::3]  # Create an offset
        a = np.array([(1, 1, 1), (-1, 0, 0)], dtype="int32")
        carray[2:3] = a[::2]  # Create an offset

        # Read all the rows:
        data = carray.read()
        if common.verbose:
            print("Object read:", data)
            print("Nrows in", carray._v_pathname, ":", carray.nrows)
            print("Third row in carray ==>", data[2])

        self.assertEqual(carray.nrows, 3)
        self.assertTrue(
            common.allequal(data[0], np.array([0, 0, 0], dtype="int32"))
        )
        self.assertTrue(
            common.allequal(data[1], np.array([3, 3, 3], dtype="int32"))
        )
        self.assertTrue(
            common.allequal(data[2], np.array([1, 1, 1], dtype="int32"))
        )


class CopyTestCase(common.TempFileMixin, common.PyTablesTestCase):

    def test01a_copy(self):
        """Checking CArray.copy() method."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01a_copy..." % self.__class__.__name__)

        # Create an CArray
        shape = (2, 2)
        atom = tb.Int16Atom()
        array1 = self.h5file.create_carray(
            self.h5file.root,
            "array1",
            atom=atom,
            shape=shape,
            title="title array1",
            chunkshape=(2, 2),
        )
        array1[...] = np.array([[456, 2], [3, 457]], dtype="int16")

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="a")
            array1 = self.h5file.root.array1

        # Copy it to another location
        array2 = array1.copy("/", "array2")

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen()
            array1 = self.h5file.root.array1
            array2 = self.h5file.root.array2

        if common.verbose:
            print("array1-->", array1.read())
            print("array2-->", array2.read())
            # print("dirs-->", dir(array1), dir(array2))
            print("attrs array1-->", repr(array1.attrs))
            print("attrs array2-->", repr(array2.attrs))

        # Check that all the elements are equal
        self.assertTrue(common.allequal(array1.read(), array2.read()))

        # Assert other properties in array
        self.assertEqual(array1.nrows, array2.nrows)
        self.assertEqual(array1.shape, array2.shape)
        self.assertEqual(array1.extdim, array2.extdim)
        self.assertEqual(array1.flavor, array2.flavor)
        self.assertEqual(array1.atom.dtype, array2.atom.dtype)
        self.assertEqual(array1.atom.type, array2.atom.type)
        self.assertEqual(array1.title, array2.title)
        self.assertEqual(str(array1.atom), str(array2.atom))
        # The next line is commented out because a copy should not
        # keep the same chunkshape anymore.
        # F. Alted 2006-11-27
        # self.assertEqual(array1.chunkshape, array2.chunkshape)

    def test01b_copy(self):
        """Checking CArray.copy() method."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01b_copy..." % self.__class__.__name__)

        # Create an CArray
        shape = (2, 2)
        atom = tb.Int16Atom()
        array1 = self.h5file.create_carray(
            self.h5file.root,
            "array1",
            atom=atom,
            shape=shape,
            title="title array1",
            chunkshape=(5, 5),
        )
        array1[...] = np.array([[456, 2], [3, 457]], dtype="int16")

        if self.close:
            if common.verbose:
                print("(closing h5fname version)")
            self._reopen(mode="a")
            array1 = self.h5file.root.array1

        # Copy it to another location
        array2 = array1.copy("/", "array2")

        if self.close:
            if common.verbose:
                print("(closing h5fname version)")
            self._reopen()
            array1 = self.h5file.root.array1
            array2 = self.h5file.root.array2

        if common.verbose:
            print("array1-->", array1.read())
            print("array2-->", array2.read())
            # print("dirs-->", dir(array1), dir(array2))
            print("attrs array1-->", repr(array1.attrs))
            print("attrs array2-->", repr(array2.attrs))

        # Check that all the elements are equal
        self.assertTrue(common.allequal(array1.read(), array2.read()))

        # Assert other properties in array
        self.assertEqual(array1.nrows, array2.nrows)
        self.assertEqual(array1.shape, array2.shape)
        self.assertEqual(array1.extdim, array2.extdim)
        self.assertEqual(array1.flavor, array2.flavor)
        self.assertEqual(array1.atom.dtype, array2.atom.dtype)
        self.assertEqual(array1.atom.type, array2.atom.type)
        self.assertEqual(array1.title, array2.title)
        self.assertEqual(str(array1.atom), str(array2.atom))
        # By default, the chunkshape should be the same
        self.assertEqual(array1.chunkshape, array2.chunkshape)

    def test01c_copy(self):
        """Checking CArray.copy() method."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01c_copy..." % self.__class__.__name__)

        # Create an CArray
        shape = (5, 5)
        atom = tb.Int16Atom()
        array1 = self.h5file.create_carray(
            self.h5file.root,
            "array1",
            atom=atom,
            shape=shape,
            title="title array1",
            chunkshape=(2, 2),
        )
        array1[:2, :2] = np.array([[456, 2], [3, 457]], dtype="int16")

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="a")
            array1 = self.h5file.root.array1

        # Copy it to another location
        array2 = array1.copy("/", "array2")

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen()
            array1 = self.h5file.root.array1
            array2 = self.h5file.root.array2

        if common.verbose:
            print("array1-->", array1.read())
            print("array2-->", array2.read())
            # print("dirs-->", dir(array1), dir(array2))
            print("attrs array1-->", repr(array1.attrs))
            print("attrs array2-->", repr(array2.attrs))

        # Check that all the elements are equal
        self.assertTrue(common.allequal(array1.read(), array2.read()))

        # Assert other properties in array
        self.assertEqual(array1.nrows, array2.nrows)
        self.assertEqual(array1.shape, array2.shape)
        self.assertEqual(array1.extdim, array2.extdim)
        self.assertEqual(array1.flavor, array2.flavor)
        self.assertEqual(array1.atom.dtype, array2.atom.dtype)
        self.assertEqual(array1.atom.type, array2.atom.type)
        self.assertEqual(array1.title, array2.title)
        self.assertEqual(str(array1.atom), str(array2.atom))
        # The next line is commented out because a copy should not
        # keep the same chunkshape anymore.
        # F. Alted 2006-11-27
        # self.assertEqual(array1.chunkshape, array2.chunkshape)

    def test02_copy(self):
        """Checking CArray.copy() method (where specified)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02_copy..." % self.__class__.__name__)

        # Create an CArray
        shape = (5, 5)
        atom = tb.Int16Atom()
        array1 = self.h5file.create_carray(
            self.h5file.root,
            "array1",
            atom=atom,
            shape=shape,
            title="title array1",
            chunkshape=(2, 2),
        )
        array1[:2, :2] = np.array([[456, 2], [3, 457]], dtype="int16")

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="a")
            array1 = self.h5file.root.array1

        # Copy to another location
        group1 = self.h5file.create_group("/", "group1")
        array2 = array1.copy(group1, "array2")

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen()
            array1 = self.h5file.root.array1
            array2 = self.h5file.root.group1.array2

        if common.verbose:
            print("array1-->", array1.read())
            print("array2-->", array2.read())
            # print("dirs-->", dir(array1), dir(array2))
            print("attrs array1-->", repr(array1.attrs))
            print("attrs array2-->", repr(array2.attrs))

        # Check that all the elements are equal
        self.assertTrue(common.allequal(array1.read(), array2.read()))

        # Assert other properties in array
        self.assertEqual(array1.nrows, array2.nrows)
        self.assertEqual(array1.shape, array2.shape)
        self.assertEqual(array1.extdim, array2.extdim)
        self.assertEqual(array1.flavor, array2.flavor)
        self.assertEqual(array1.atom.dtype, array2.atom.dtype)
        self.assertEqual(array1.atom.type, array2.atom.type)
        self.assertEqual(array1.title, array2.title)
        self.assertEqual(str(array1.atom), str(array2.atom))
        # The next line is commented out because a copy should not
        # keep the same chunkshape anymore.
        # F. Alted 2006-11-27
        # self.assertEqual(array1.chunkshape, array2.chunkshape)

    def test03a_copy(self):
        """Checking CArray.copy() method (python flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03c_copy..." % self.__class__.__name__)

        shape = (2, 2)
        atom = tb.Int16Atom()
        array1 = self.h5file.create_carray(
            self.h5file.root,
            "array1",
            atom=atom,
            shape=shape,
            title="title array1",
            chunkshape=(2, 2),
        )
        array1.flavor = "python"
        array1[...] = [[456, 2], [3, 457]]

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="a")
            array1 = self.h5file.root.array1

        # Copy to another location
        array2 = array1.copy("/", "array2")

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen()
            array1 = self.h5file.root.array1
            array2 = self.h5file.root.array2

        if common.verbose:
            print("attrs array1-->", repr(array1.attrs))
            print("attrs array2-->", repr(array2.attrs))

        # Check that all elements are equal
        self.assertEqual(array1.read(), array2.read())
        # Assert other properties in array
        self.assertEqual(array1.nrows, array2.nrows)
        self.assertEqual(array1.shape, array2.shape)
        self.assertEqual(array1.extdim, array2.extdim)
        self.assertEqual(array1.flavor, array2.flavor)  # Very important here!
        self.assertEqual(array1.atom.dtype, array2.atom.dtype)
        self.assertEqual(array1.atom.type, array2.atom.type)
        self.assertEqual(array1.title, array2.title)
        self.assertEqual(str(array1.atom), str(array2.atom))
        # The next line is commented out because a copy should not
        # keep the same chunkshape anymore.
        # F. Alted 2006-11-27
        # self.assertEqual(array1.chunkshape, array2.chunkshape)

    def test03b_copy(self):
        """Checking CArray.copy() method (string python flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03d_copy..." % self.__class__.__name__)

        shape = (2, 2)
        atom = tb.StringAtom(itemsize=4)
        array1 = self.h5file.create_carray(
            self.h5file.root,
            "array1",
            atom=atom,
            shape=shape,
            title="title array1",
            chunkshape=(2, 2),
        )
        array1.flavor = "python"
        array1[...] = [["456", "2"], ["3", "457"]]

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="a")
            array1 = self.h5file.root.array1

        # Copy to another location
        array2 = array1.copy("/", "array2")

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen()
            array1 = self.h5file.root.array1
            array2 = self.h5file.root.array2

        if common.verbose:
            print("type value-->", type(array2[:][0][0]))
            print("value-->", array2[:])
            print("attrs array1-->", repr(array1.attrs))
            print("attrs array2-->", repr(array2.attrs))

        # Check that all elements are equal
        self.assertEqual(array1.read(), array2.read())

        # Assert other properties in array
        self.assertEqual(array1.nrows, array2.nrows)
        self.assertEqual(array1.shape, array2.shape)
        self.assertEqual(array1.extdim, array2.extdim)
        self.assertEqual(array1.flavor, array2.flavor)  # Very important here!
        self.assertEqual(array1.atom.dtype, array2.atom.dtype)
        self.assertEqual(array1.atom.type, array2.atom.type)
        self.assertEqual(array1.title, array2.title)
        self.assertEqual(str(array1.atom), str(array2.atom))
        # The next line is commented out because a copy should not
        # keep the same chunkshape anymore.
        # F. Alted 2006-11-27
        # self.assertEqual(array1.chunkshape, array2.chunkshape)

    def test03c_copy(self):
        """Checking CArray.copy() method (chararray flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03e_copy..." % self.__class__.__name__)

        shape = (2, 2)
        atom = tb.StringAtom(itemsize=4)
        array1 = self.h5file.create_carray(
            self.h5file.root,
            "array1",
            atom=atom,
            shape=shape,
            title="title array1",
            chunkshape=(2, 2),
        )
        array1[...] = np.array([["456", "2"], ["3", "457"]], dtype="S4")

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="a")
            array1 = self.h5file.root.array1

        # Copy to another location
        array2 = array1.copy("/", "array2")

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen()
            array1 = self.h5file.root.array1
            array2 = self.h5file.root.array2

        if common.verbose:
            print("attrs array1-->", repr(array1.attrs))
            print("attrs array2-->", repr(array2.attrs))

        # Check that all elements are equal
        self.assertTrue(common.allequal(array1.read(), array2.read()))
        # Assert other properties in array
        self.assertEqual(array1.nrows, array2.nrows)
        self.assertEqual(array1.shape, array2.shape)
        self.assertEqual(array1.extdim, array2.extdim)
        self.assertEqual(array1.flavor, array2.flavor)  # Very important here!
        self.assertEqual(array1.atom.dtype, array2.atom.dtype)
        self.assertEqual(array1.atom.type, array2.atom.type)
        self.assertEqual(array1.title, array2.title)
        self.assertEqual(str(array1.atom), str(array2.atom))
        # The next line is commented out because a copy should not
        # keep the same chunkshape anymore.
        # F. Alted 2006-11-27
        # self.assertEqual(array1.chunkshape, array2.chunkshape)

    def test04_copy(self):
        """Checking CArray.copy() method (checking title copying)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test04_copy..." % self.__class__.__name__)

        # Create an CArray
        shape = (2, 2)
        atom = tb.Int16Atom()
        array1 = self.h5file.create_carray(
            self.h5file.root,
            "array1",
            atom=atom,
            shape=shape,
            title="title array1",
            chunkshape=(2, 2),
        )
        array1[...] = np.array([[456, 2], [3, 457]], dtype="int16")

        # Append some user attrs
        array1.attrs.attr1 = "attr1"
        array1.attrs.attr2 = 2

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="a")
            array1 = self.h5file.root.array1

        # Copy it to another Array
        array2 = array1.copy("/", "array2", title="title array2")

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen()
            array1 = self.h5file.root.array1
            array2 = self.h5file.root.array2

        # Assert user attributes
        if common.verbose:
            print("title of destination array-->", array2.title)
        self.assertEqual(array2.title, "title array2")

    def test05_copy(self):
        """Checking CArray.copy() method (user attributes copied)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test05_copy..." % self.__class__.__name__)

        # Create an CArray
        shape = (2, 2)
        atom = tb.Int16Atom()
        array1 = self.h5file.create_carray(
            self.h5file.root,
            "array1",
            atom=atom,
            shape=shape,
            title="title array1",
            chunkshape=(2, 2),
        )
        array1[...] = np.array([[456, 2], [3, 457]], dtype="int16")

        # Append some user attrs
        array1.attrs.attr1 = "attr1"
        array1.attrs.attr2 = 2

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="a")
            array1 = self.h5file.root.array1

        # Copy it to another Array
        array2 = array1.copy("/", "array2", copyuserattrs=1)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen()
            array1 = self.h5file.root.array1
            array2 = self.h5file.root.array2

        if common.verbose:
            print("attrs array1-->", repr(array1.attrs))
            print("attrs array2-->", repr(array2.attrs))

        # Assert user attributes
        self.assertEqual(array2.attrs.attr1, "attr1")
        self.assertEqual(array2.attrs.attr2, 2)

    def test05b_copy(self):
        """Checking CArray.copy() method (user attributes not copied)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test05b_copy..." % self.__class__.__name__)

        # Create an Array
        shape = (2, 2)
        atom = tb.Int16Atom()
        array1 = self.h5file.create_carray(
            self.h5file.root,
            "array1",
            atom=atom,
            shape=shape,
            title="title array1",
            chunkshape=(2, 2),
        )
        array1[...] = np.array([[456, 2], [3, 457]], dtype="int16")

        # Append some user attrs
        array1.attrs.attr1 = "attr1"
        array1.attrs.attr2 = 2

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="a")
            array1 = self.h5file.root.array1

        # Copy it to another Array
        array2 = array1.copy("/", "array2", copyuserattrs=0)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen()
            array1 = self.h5file.root.array1
            array2 = self.h5file.root.array2

        if common.verbose:
            print("attrs array1-->", repr(array1.attrs))
            print("attrs array2-->", repr(array2.attrs))

        # Assert user attributes
        self.assertEqual(hasattr(array2.attrs, "attr1"), 0)
        self.assertEqual(hasattr(array2.attrs, "attr2"), 0)


class CloseCopyTestCase(CopyTestCase):
    close = 1


class OpenCopyTestCase(CopyTestCase):
    close = 0


class CopyIndexTestCase(common.TempFileMixin, common.PyTablesTestCase):
    nrowsinbuf = 2

    def test01_index(self):
        """Checking CArray.copy() method with indexes."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01_index..." % self.__class__.__name__)

        # Create an CArray
        shape = (100, 2)
        atom = tb.Int32Atom()
        array1 = self.h5file.create_carray(
            self.h5file.root,
            "array1",
            atom=atom,
            shape=shape,
            title="title array1",
            chunkshape=(2, 2),
        )
        r = np.arange(200, dtype="int32")
        r.shape = shape
        array1[...] = r

        # Select a different buffer size:
        array1.nrowsinbuf = self.nrowsinbuf

        # Copy to another array
        array2 = array1.copy(
            "/", "array2", start=self.start, stop=self.stop, step=self.step
        )
        if common.verbose:
            print("array1-->", array1.read())
            print("array2-->", array2.read())
            print("attrs array1-->", repr(array1.attrs))
            print("attrs array2-->", repr(array2.attrs))

        # Check that all the elements are equal
        r2 = r[self.start : self.stop : self.step]
        self.assertTrue(common.allequal(r2, array2.read()))

        # Assert the number of rows in array
        if common.verbose:
            print("nrows in array2-->", array2.nrows)
            print("and it should be-->", r2.shape[0])

        # The next line is commented out because a copy should not
        # keep the same chunkshape anymore.
        # F. Alted 2006-11-27
        # assert array1.chunkshape == array2.chunkshape
        self.assertEqual(r2.shape[0], array2.nrows)

    def _test02_indexclosef(self):
        """Checking CArray.copy() method with indexes (close file version)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02_indexclosef..." % self.__class__.__name__)

        # Create an CArray
        shape = (100, 2)
        atom = tb.Int32Atom()
        array1 = self.h5file.create_carray(
            self.h5file.root,
            "array1",
            atom=atom,
            shape=shape,
            title="title array1",
            chunkshape=(2, 2),
        )
        r = np.arange(200, dtype="int32")
        r.shape = shape
        array1[...] = r

        # Select a different buffer size:
        array1.nrowsinbuf = self.nrowsinbuf

        # Copy to another array
        array2 = array1.copy(
            "/", "array2", start=self.start, stop=self.stop, step=self.step
        )

        # Close and reopen the file
        self._reopen()
        array1 = self.h5file.root.array1
        array2 = self.h5file.root.array2

        if common.verbose:
            print("array1-->", array1.read())
            print("array2-->", array2.read())
            print("attrs array1-->", repr(array1.attrs))
            print("attrs array2-->", repr(array2.attrs))

        # Check that all the elements are equal
        r2 = r[self.start : self.stop : self.step]
        self.assertEqual(array1.chunkshape, array2.chunkshape)
        self.assertTrue(common.allequal(r2, array2.read()))

        # Assert the number of rows in array
        if common.verbose:
            print("nrows in array2-->", array2.nrows)
            print("and it should be-->", r2.shape[0])
        self.assertEqual(r2.shape[0], array2.nrows)


class CopyIndex1TestCase(CopyIndexTestCase):
    nrowsinbuf = 1
    start = 0
    stop = 7
    step = 1


class CopyIndex2TestCase(CopyIndexTestCase):
    nrowsinbuf = 2
    start = 0
    stop = -1
    step = 1


class CopyIndex3TestCase(CopyIndexTestCase):
    nrowsinbuf = 3
    start = 1
    stop = 7
    step = 1


class CopyIndex4TestCase(CopyIndexTestCase):
    nrowsinbuf = 4
    start = 0
    stop = 6
    step = 1


class CopyIndex5TestCase(CopyIndexTestCase):
    nrowsinbuf = 2
    start = 3
    stop = 7
    step = 1


class CopyIndex6TestCase(CopyIndexTestCase):
    nrowsinbuf = 2
    start = 3
    stop = 6
    step = 2


class CopyIndex7TestCase(CopyIndexTestCase):
    start = 0
    stop = 7
    step = 10


class CopyIndex8TestCase(CopyIndexTestCase):
    start = 6
    stop = -1  # Negative values means starting from the end
    step = 1


class CopyIndex9TestCase(CopyIndexTestCase):
    start = 3
    stop = 4
    step = 1


class CopyIndex10TestCase(CopyIndexTestCase):
    nrowsinbuf = 1
    start = 3
    stop = 4
    step = 2


class CopyIndex11TestCase(CopyIndexTestCase):
    start = -3
    stop = -1
    step = 2


class CopyIndex12TestCase(CopyIndexTestCase):
    start = -1  # Should point to the last element
    stop = None  # None should mean the last element (including it)
    step = 1


# The next test should be run only in **heavy** mode
class Rows64bitsTestCase(common.TempFileMixin, common.PyTablesTestCase):
    narows = 1000 * 1000  # each array will have 1 million entries
    # narows = 1000        # for testing only
    nanumber = 1000 * 3  # That should account for more than 2**31-1

    def setUp(self):
        super().setUp()

        # Create an CArray
        shape = (self.narows * self.nanumber,)
        array = self.h5file.create_carray(
            self.h5file.root,
            "array",
            atom=tb.Int8Atom(),
            shape=shape,
            filters=tb.Filters(complib="lzo", complevel=1),
        )

        # Fill the array
        na = np.arange(self.narows, dtype="int8")
        # for i in xrange(self.nanumber):
        #     s = slice(i * self.narows, (i + 1)*self.narows)
        #     array[s] = na
        s = slice(0, self.narows)
        array[s] = na
        s = slice(
            (self.nanumber - 1) * self.narows, self.nanumber * self.narows
        )
        array[s] = na

    def test01_basiccheck(self):
        """Some basic checks for carrays exceeding 2**31 rows"""

        array = self.h5file.root.array

        if self.close:
            if common.verbose:
                # Check how many entries there are in the array
                print("Before closing")
                print("Entries:", array.nrows, type(array.nrows))
                print("Entries:", array.nrows / (1000 * 1000), "Millions")
                print("Shape:", array.shape)

            # Re-open the file
            self._reopen()
            array = self.h5file.root.array
            if common.verbose:
                print("After re-open")

        # Check how many entries there are in the array
        if common.verbose:
            print("Entries:", array.nrows, type(array.nrows))
            print("Entries:", array.nrows / (1000 * 1000), "Millions")
            print("Shape:", array.shape)
            print("Last 10 elements-->", array[-10:])
            stop = self.narows % 256
            if stop > 127:
                stop -= 256
            start = stop - 10
            # print("start, stop-->", start, stop)
            print("Should look like:", np.arange(start, stop, dtype="int8"))

        nrows = self.narows * self.nanumber

        # check nrows
        self.assertEqual(array.nrows, nrows)

        # Check shape
        self.assertEqual(array.shape, (nrows,))

        # check the 10 first elements
        self.assertTrue(
            common.allequal(array[:10], np.arange(10, dtype="int8"))
        )

        # check the 10 last elements
        stop = self.narows % 256
        if stop > 127:
            stop -= 256
        start = stop - 10
        self.assertTrue(
            common.allequal(array[-10:], np.arange(start, stop, dtype="int8"))
        )


class Rows64bitsTestCase1(Rows64bitsTestCase):
    close = 0


class Rows64bitsTestCase2(Rows64bitsTestCase):
    close = 1


class BigArrayTestCase(common.TempFileMixin, common.PyTablesTestCase):
    shape = (3_000_000_000,)  # more than 2**31-1

    def setUp(self):
        super().setUp()
        # This should be fast since disk space isn't actually allocated,
        # so this case is OK for non-heavy test runs.
        self.h5file.create_carray(
            "/", "array", atom=tb.Int8Atom(), shape=self.shape
        )

    def test00_shape(self):
        """Check that the shape doesn't overflow."""
        # See ticket #147.
        self.assertEqual(self.h5file.root.array.shape, self.shape)
        try:
            self.assertEqual(len(self.h5file.root.array), self.shape[0])
        except OverflowError:
            # This can't be avoided in 32-bit platforms.
            self.assertTrue(
                self.shape[0] > np.iinfo(int).max,
                "Array length overflowed but ``int`` " "is wide enough.",
            )

    def test01_shape_reopen(self):
        """Check that the shape doesn't overflow after reopening."""
        self._reopen("r")
        self.test00_shape()


# Test for default values when creating arrays.
class DfltAtomTestCase(common.TempFileMixin, common.PyTablesTestCase):

    def test00_dflt(self):
        """Check that Atom.dflt is honored (string version)."""

        # Create a CArray with default values
        self.h5file.create_carray(
            "/",
            "bar",
            atom=tb.StringAtom(itemsize=5, dflt=b"abdef"),
            shape=(10, 10),
        )

        if self.reopen:
            self._reopen()

        # Check the values
        values = self.h5file.root.bar[:]
        if common.verbose:
            print("Read values:", values)
        self.assertTrue(
            common.allequal(
                values, np.array(["abdef"] * 100, "S5").reshape(10, 10)
            )
        )

    def test01_dflt(self):
        """Check that Atom.dflt is honored (int version)."""

        # Create a CArray with default values
        self.h5file.create_carray(
            "/", "bar", atom=tb.IntAtom(dflt=1), shape=(10, 10)
        )

        if self.reopen:
            self._reopen()

        # Check the values
        values = self.h5file.root.bar[:]
        if common.verbose:
            print("Read values:", values)
        self.assertTrue(common.allequal(values, np.ones((10, 10), "i4")))

    def test02_dflt(self):
        """Check that Atom.dflt is honored (float version)."""

        # Create a CArray with default values
        self.h5file.create_carray(
            "/", "bar", atom=tb.FloatAtom(dflt=1.134), shape=(10, 10)
        )

        if self.reopen:
            self._reopen()

        # Check the values
        values = self.h5file.root.bar[:]
        if common.verbose:
            print("Read values:", values)
        self.assertTrue(
            common.allequal(values, np.ones((10, 10), "f8") * 1.134)
        )


class DfltAtomNoReopen(DfltAtomTestCase):
    reopen = False


class DfltAtomReopen(DfltAtomTestCase):
    reopen = True


# Test for representation of defaults in atoms. Ticket #212.
class AtomDefaultReprTestCase(common.TempFileMixin, common.PyTablesTestCase):

    def test00a_zeros(self):
        """Testing default values.  Zeros (scalar)."""

        N = ()
        atom = tb.StringAtom(itemsize=3, shape=N, dflt=b"")
        ca = self.h5file.create_carray("/", "test", atom=atom, shape=(1,))
        if self.reopen:
            self._reopen("a")
            ca = self.h5file.root.test
        # Check the value
        if common.verbose:
            print("First row-->", repr(ca[0]))
            print("Defaults-->", repr(ca.atom.dflt))
        self.assertTrue(common.allequal(ca[0], np.zeros(N, "S3")))
        self.assertTrue(common.allequal(ca.atom.dflt, np.zeros(N, "S3")))

    def test00b_zeros(self):
        """Testing default values.  Zeros (array)."""

        N = 2
        atom = tb.StringAtom(itemsize=3, shape=N, dflt=b"")
        ca = self.h5file.create_carray("/", "test", atom=atom, shape=(1,))
        if self.reopen:
            self._reopen("a")
            ca = self.h5file.root.test
        # Check the value
        if common.verbose:
            print("First row-->", ca[0])
            print("Defaults-->", ca.atom.dflt)
        self.assertTrue(common.allequal(ca[0], np.zeros(N, "S3")))
        self.assertTrue(common.allequal(ca.atom.dflt, np.zeros(N, "S3")))

    def test01a_values(self):
        """Testing default values.  Ones."""

        N = 2
        atom = tb.Int32Atom(shape=N, dflt=1)
        ca = self.h5file.create_carray("/", "test", atom=atom, shape=(1,))
        if self.reopen:
            self._reopen("a")
            ca = self.h5file.root.test
        # Check the value
        if common.verbose:
            print("First row-->", ca[0])
            print("Defaults-->", ca.atom.dflt)
        self.assertTrue(common.allequal(ca[0], np.ones(N, "i4")))
        self.assertTrue(common.allequal(ca.atom.dflt, np.ones(N, "i4")))

    def test01b_values(self):
        """Testing default values.  Generic value."""

        N = 2
        generic = 112.32
        atom = tb.Float32Atom(shape=N, dflt=generic)
        ca = self.h5file.create_carray("/", "test", atom=atom, shape=(1,))
        if self.reopen:
            self._reopen("a")
            ca = self.h5file.root.test
        # Check the value
        if common.verbose:
            print("First row-->", ca[0])
            print("Defaults-->", ca.atom.dflt)
        self.assertTrue(common.allequal(ca[0], np.ones(N, "f4") * generic))
        self.assertTrue(
            common.allequal(ca.atom.dflt, np.ones(N, "f4") * generic)
        )

    def test02a_None(self):
        """Testing default values.  None (scalar)."""

        N = ()
        atom = tb.Int32Atom(shape=N, dflt=None)
        ca = self.h5file.create_carray("/", "test", atom=atom, shape=(1,))
        if self.reopen:
            self._reopen("a")
            ca = self.h5file.root.test
        # Check the value
        if common.verbose:
            print("First row-->", repr(ca[0]))
            print("Defaults-->", repr(ca.atom.dflt))
        self.assertTrue(common.allequal(ca.atom.dflt, np.zeros(N, "i4")))

    def test02b_None(self):
        """Testing default values.  None (array)."""

        N = 2
        atom = tb.Int32Atom(shape=N, dflt=None)
        ca = self.h5file.create_carray("/", "test", atom=atom, shape=(1,))
        if self.reopen:
            self._reopen("a")
            ca = self.h5file.root.test
        # Check the value
        if common.verbose:
            print("First row-->", ca[0])
            print("Defaults-->", ca.atom.dflt)
        self.assertTrue(common.allequal(ca.atom.dflt, np.zeros(N, "i4")))


class AtomDefaultReprNoReopen(AtomDefaultReprTestCase):
    reopen = False


class AtomDefaultReprReopen(AtomDefaultReprTestCase):
    reopen = True


class TruncateTestCase(common.TempFileMixin, common.PyTablesTestCase):
    def test(self):
        """Test for inability to truncate Array objects."""

        array1 = self.h5file.create_carray("/", "array1", tb.IntAtom(), [2, 2])
        self.assertRaises(TypeError, array1.truncate, 0)


# Test for dealing with multidimensional atoms
class MDAtomTestCase(common.TempFileMixin, common.PyTablesTestCase):

    def test01a_assign(self):
        """Assign a row to a (unidimensional) CArray with a MD atom."""

        # Create an CArray
        ca = self.h5file.create_carray(
            "/", "test", atom=tb.Int32Atom((2, 2)), shape=(1,)
        )
        if self.reopen:
            self._reopen("a")
            ca = self.h5file.root.test
        # Assign one row
        ca[0] = [[1, 3], [4, 5]]
        self.assertEqual(ca.nrows, 1)
        if common.verbose:
            print("First row-->", ca[0])
        self.assertTrue(
            common.allequal(ca[0], np.array([[1, 3], [4, 5]], "i4"))
        )

    def test01b_assign(self):
        """Assign several rows to a (unidimensional) CArray with a MD atom."""

        # Create an CArray
        ca = self.h5file.create_carray(
            "/", "test", atom=tb.Int32Atom((2, 2)), shape=(3,)
        )
        if self.reopen:
            self._reopen("a")
            ca = self.h5file.root.test
        # Assign three rows
        ca[:] = [[[1]], [[2]], [[3]]]  # Simple broadcast
        self.assertEqual(ca.nrows, 3)
        if common.verbose:
            print("Third row-->", ca[2])
        self.assertTrue(
            common.allequal(ca[2], np.array([[3, 3], [3, 3]], "i4"))
        )

    def test02a_assign(self):
        """Assign a row to a (multidimensional) CArray with a MD atom."""

        # Create an CArray
        ca = self.h5file.create_carray(
            "/", "test", atom=tb.Int32Atom((2,)), shape=(1, 3)
        )
        if self.reopen:
            self._reopen("a")
            ca = self.h5file.root.test
        # Assign one row
        ca[:] = [[[1, 3], [4, 5], [7, 9]]]
        self.assertEqual(ca.nrows, 1)
        if common.verbose:
            print("First row-->", ca[0])
        self.assertTrue(
            common.allequal(ca[0], np.array([[1, 3], [4, 5], [7, 9]], "i4"))
        )

    def test02b_assign(self):
        """Assign several rows to a (multidimensional) CArray with
        a MD atom."""

        # Create an CArray
        ca = self.h5file.create_carray(
            "/", "test", atom=tb.Int32Atom((2,)), shape=(3, 3)
        )
        if self.reopen:
            self._reopen("a")
            ca = self.h5file.root.test
        # Assign three rows
        ca[:] = [
            [[1, -3], [4, -5], [-7, 9]],
            [[-1, 3], [-4, 5], [7, -8]],
            [[-2, 3], [-5, 5], [7, -9]],
        ]
        self.assertEqual(ca.nrows, 3)
        if common.verbose:
            print("Third row-->", ca[2])
        self.assertTrue(
            common.allequal(ca[2], np.array([[-2, 3], [-5, 5], [7, -9]], "i4"))
        )

    def test03a_MDMDMD(self):
        """Complex assign of a MD array in a MD CArray with a MD atom."""

        # Create an CArray
        ca = self.h5file.create_carray(
            "/", "test", atom=tb.Int32Atom((2, 4)), shape=(3, 2, 3)
        )
        if self.reopen:
            self._reopen("a")
            ca = self.h5file.root.test

        # Assign values
        # The shape of the atom should be added at the end of the arrays
        a = np.arange(2 * 3 * 2 * 4, dtype="i4").reshape((2, 3, 2, 4))
        ca[:] = [a * 1, a * 2, a * 3]
        self.assertEqual(ca.nrows, 3)
        if common.verbose:
            print("Third row-->", ca[2])
        self.assertTrue(common.allequal(ca[2], a * 3))

    def test03b_MDMDMD(self):
        """Complex assign of a MD array in a MD CArray with a MD atom (II)."""

        # Create an CArray
        ca = self.h5file.create_carray(
            "/", "test", atom=tb.Int32Atom((2, 4)), shape=(2, 3, 3)
        )
        if self.reopen:
            self._reopen("a")
            ca = self.h5file.root.test

        # Assign values
        # The shape of the atom should be added at the end of the arrays
        a = np.arange(2 * 3 * 3 * 2 * 4, dtype="i4").reshape((2, 3, 3, 2, 4))
        ca[:] = a
        self.assertEqual(ca.nrows, 2)
        if common.verbose:
            print("Third row-->", ca[:, 2, ...])
        self.assertTrue(common.allequal(ca[:, 2, ...], a[:, 2, ...]))

    def test03c_MDMDMD(self):
        """Complex assign of a MD array in a MD CArray with a MD atom (III)."""

        # Create an CArray
        ca = self.h5file.create_carray(
            "/", "test", atom=tb.Int32Atom((2, 4)), shape=(3, 1, 2)
        )
        if self.reopen:
            self._reopen("a")
            ca = self.h5file.root.test

        # Assign values
        # The shape of the atom should be added at the end of the arrays
        a = np.arange(3 * 1 * 2 * 2 * 4, dtype="i4").reshape((3, 1, 2, 2, 4))
        ca[:] = a
        self.assertEqual(ca.nrows, 3)
        if common.verbose:
            print("Second row-->", ca[:, :, 1, ...])
        self.assertTrue(common.allequal(ca[:, :, 1, ...], a[:, :, 1, ...]))


class MDAtomNoReopen(MDAtomTestCase):
    reopen = False


class MDAtomReopen(MDAtomTestCase):
    reopen = True


# Test for building very large MD atoms without defaults.  Ticket #211.
class MDLargeAtomTestCase(common.TempFileMixin, common.PyTablesTestCase):

    def test01_create(self):
        """Create a CArray with a very large MD atom."""

        N = 2**16  # 4x larger than maximum object header size (64 KB)
        ca = self.h5file.create_carray(
            "/", "test", atom=tb.Int32Atom(shape=N), shape=(1,)
        )
        if self.reopen:
            self._reopen("a")
            ca = self.h5file.root.test

        # Check the value
        if common.verbose:
            print("First row-->", ca[0])
        self.assertTrue(common.allequal(ca[0], np.zeros(N, "i4")))


class MDLargeAtomNoReopen(MDLargeAtomTestCase):
    reopen = False


class MDLargeAtomReopen(MDLargeAtomTestCase):
    reopen = True


class AccessClosedTestCase(common.TempFileMixin, common.PyTablesTestCase):

    def setUp(self):
        super().setUp()
        self.array = self.h5file.create_carray(
            self.h5file.root, "array", atom=tb.Int32Atom(), shape=(10, 10)
        )
        self.array[...] = np.zeros((10, 10))

    def test_read(self):
        self.h5file.close()
        self.assertRaises(tb.ClosedNodeError, self.array.read)

    def test_getitem(self):
        self.h5file.close()
        self.assertRaises(tb.ClosedNodeError, self.array.__getitem__, 0)

    def test_setitem(self):
        self.h5file.close()
        self.assertRaises(tb.ClosedNodeError, self.array.__setitem__, 0, 0)


class TestCreateCArrayArgs(common.TempFileMixin, common.PyTablesTestCase):
    obj = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    where = "/"
    name = "carray"
    atom = tb.Atom.from_dtype(obj.dtype)
    shape = obj.shape
    title = "title"
    filters = None
    chunkshape = (1, 2)
    byteorder = None
    createparents = False

    def test_positional_args_01(self):
        self.h5file.create_carray(
            self.where,
            self.name,
            self.atom,
            self.shape,
            self.title,
            self.filters,
            self.chunkshape,
        )
        self.h5file.close()

        self.h5file = tb.open_file(self.h5fname)
        ptarr = self.h5file.get_node(self.where, self.name)
        nparr = ptarr.read()

        self.assertEqual(ptarr.title, self.title)
        self.assertEqual(ptarr.shape, self.shape)
        self.assertEqual(ptarr.atom, self.atom)
        self.assertEqual(ptarr.atom.dtype, self.atom.dtype)
        self.assertEqual(ptarr.chunkshape, self.chunkshape)
        self.assertTrue(common.allequal(np.zeros_like(self.obj), nparr))

    def test_positional_args_02(self):
        ptarr = self.h5file.create_carray(
            self.where,
            self.name,
            self.atom,
            self.shape,
            self.title,
            self.filters,
            self.chunkshape,
        )
        ptarr[...] = self.obj
        self.h5file.close()

        self.h5file = tb.open_file(self.h5fname)
        ptarr = self.h5file.get_node(self.where, self.name)
        nparr = ptarr.read()

        self.assertEqual(ptarr.title, self.title)
        self.assertEqual(ptarr.shape, self.shape)
        self.assertEqual(ptarr.atom, self.atom)
        self.assertEqual(ptarr.atom.dtype, self.atom.dtype)
        self.assertEqual(ptarr.chunkshape, self.chunkshape)
        self.assertTrue(common.allequal(self.obj, nparr))

    def test_positional_args_obj(self):
        self.h5file.create_carray(
            self.where,
            self.name,
            None,
            None,
            self.title,
            self.filters,
            self.chunkshape,
            self.byteorder,
            self.createparents,
            self.obj,
        )
        self.h5file.close()

        self.h5file = tb.open_file(self.h5fname)
        ptarr = self.h5file.get_node(self.where, self.name)
        nparr = ptarr.read()

        self.assertEqual(ptarr.title, self.title)
        self.assertEqual(ptarr.shape, self.shape)
        self.assertEqual(ptarr.atom, self.atom)
        self.assertEqual(ptarr.atom.dtype, self.atom.dtype)
        self.assertEqual(ptarr.chunkshape, self.chunkshape)
        self.assertTrue(common.allequal(self.obj, nparr))

    def test_kwargs_obj(self):
        self.h5file.create_carray(
            self.where,
            self.name,
            title=self.title,
            chunkshape=self.chunkshape,
            obj=self.obj,
        )
        self.h5file.close()

        self.h5file = tb.open_file(self.h5fname)
        ptarr = self.h5file.get_node(self.where, self.name)
        nparr = ptarr.read()

        self.assertEqual(ptarr.title, self.title)
        self.assertEqual(ptarr.shape, self.shape)
        self.assertEqual(ptarr.atom, self.atom)
        self.assertEqual(ptarr.atom.dtype, self.atom.dtype)
        self.assertEqual(ptarr.chunkshape, self.chunkshape)
        self.assertTrue(common.allequal(self.obj, nparr))

    def test_kwargs_atom_shape_01(self):
        ptarr = self.h5file.create_carray(
            self.where,
            self.name,
            title=self.title,
            chunkshape=self.chunkshape,
            atom=self.atom,
            shape=self.shape,
        )
        ptarr[...] = self.obj
        self.h5file.close()

        self.h5file = tb.open_file(self.h5fname)
        ptarr = self.h5file.get_node(self.where, self.name)
        nparr = ptarr.read()

        self.assertEqual(ptarr.title, self.title)
        self.assertEqual(ptarr.shape, self.shape)
        self.assertEqual(ptarr.atom, self.atom)
        self.assertEqual(ptarr.atom.dtype, self.atom.dtype)
        self.assertEqual(ptarr.chunkshape, self.chunkshape)
        self.assertTrue(common.allequal(self.obj, nparr))

    def test_kwargs_atom_shape_02(self):
        ptarr = self.h5file.create_carray(
            self.where,
            self.name,
            title=self.title,
            chunkshape=self.chunkshape,
            atom=self.atom,
            shape=self.shape,
        )
        # ptarr[...] = self.obj
        self.h5file.close()

        self.h5file = tb.open_file(self.h5fname)
        ptarr = self.h5file.get_node(self.where, self.name)
        nparr = ptarr.read()

        self.assertEqual(ptarr.title, self.title)
        self.assertEqual(ptarr.shape, self.shape)
        self.assertEqual(ptarr.atom, self.atom)
        self.assertEqual(ptarr.atom.dtype, self.atom.dtype)
        self.assertEqual(ptarr.chunkshape, self.chunkshape)
        self.assertTrue(common.allequal(np.zeros_like(self.obj), nparr))

    def test_kwargs_obj_atom(self):
        ptarr = self.h5file.create_carray(
            self.where,
            self.name,
            title=self.title,
            chunkshape=self.chunkshape,
            obj=self.obj,
            atom=self.atom,
        )
        self.h5file.close()

        self.h5file = tb.open_file(self.h5fname)
        ptarr = self.h5file.get_node(self.where, self.name)
        nparr = ptarr.read()

        self.assertEqual(ptarr.title, self.title)
        self.assertEqual(ptarr.shape, self.shape)
        self.assertEqual(ptarr.atom, self.atom)
        self.assertEqual(ptarr.atom.dtype, self.atom.dtype)
        self.assertEqual(ptarr.chunkshape, self.chunkshape)
        self.assertTrue(common.allequal(self.obj, nparr))

    def test_kwargs_obj_shape(self):
        ptarr = self.h5file.create_carray(
            self.where,
            self.name,
            title=self.title,
            chunkshape=self.chunkshape,
            obj=self.obj,
            shape=self.shape,
        )
        self.h5file.close()

        self.h5file = tb.open_file(self.h5fname)
        ptarr = self.h5file.get_node(self.where, self.name)
        nparr = ptarr.read()

        self.assertEqual(ptarr.title, self.title)
        self.assertEqual(ptarr.shape, self.shape)
        self.assertEqual(ptarr.atom, self.atom)
        self.assertEqual(ptarr.atom.dtype, self.atom.dtype)
        self.assertEqual(ptarr.chunkshape, self.chunkshape)
        self.assertTrue(common.allequal(self.obj, nparr))

    def test_kwargs_obj_atom_shape(self):
        ptarr = self.h5file.create_carray(
            self.where,
            self.name,
            title=self.title,
            chunkshape=self.chunkshape,
            obj=self.obj,
            atom=self.atom,
            shape=self.shape,
        )
        self.h5file.close()

        self.h5file = tb.open_file(self.h5fname)
        ptarr = self.h5file.get_node(self.where, self.name)
        nparr = ptarr.read()

        self.assertEqual(ptarr.title, self.title)
        self.assertEqual(ptarr.shape, self.shape)
        self.assertEqual(ptarr.atom, self.atom)
        self.assertEqual(ptarr.atom.dtype, self.atom.dtype)
        self.assertEqual(ptarr.chunkshape, self.chunkshape)
        self.assertTrue(common.allequal(self.obj, nparr))

    def test_kwargs_obj_atom_error(self):
        atom = tb.Atom.from_dtype(np.dtype("complex"))
        # shape = self.shape + self.shape
        self.assertRaises(
            TypeError,
            self.h5file.create_carray,
            self.where,
            self.name,
            title=self.title,
            obj=self.obj,
            atom=atom,
        )

    def test_kwargs_obj_shape_error(self):
        # atom = Atom.from_dtype(np.dtype('complex'))
        shape = self.shape + self.shape
        self.assertRaises(
            TypeError,
            self.h5file.create_carray,
            self.where,
            self.name,
            title=self.title,
            obj=self.obj,
            shape=shape,
        )

    def test_kwargs_obj_atom_shape_error_01(self):
        atom = tb.Atom.from_dtype(np.dtype("complex"))
        # shape = self.shape + self.shape
        self.assertRaises(
            TypeError,
            self.h5file.create_carray,
            self.where,
            self.name,
            title=self.title,
            obj=self.obj,
            atom=atom,
            shape=self.shape,
        )

    def test_kwargs_obj_atom_shape_error_02(self):
        # atom = Atom.from_dtype(np.dtype('complex'))
        shape = self.shape + self.shape
        self.assertRaises(
            TypeError,
            self.h5file.create_carray,
            self.where,
            self.name,
            title=self.title,
            obj=self.obj,
            atom=self.atom,
            shape=shape,
        )

    def test_kwargs_obj_atom_shape_error_03(self):
        atom = tb.Atom.from_dtype(np.dtype("complex"))
        shape = self.shape + self.shape
        self.assertRaises(
            TypeError,
            self.h5file.create_carray,
            self.where,
            self.name,
            title=self.title,
            obj=self.obj,
            atom=atom,
            shape=shape,
        )


def suite():
    theSuite = common.unittest.TestSuite()
    niter = 1
    # common.heavy = 1  # uncomment this only for testing purposes

    # theSuite.addTest(make_suite(BasicTestCase))
    for n in range(niter):
        theSuite.addTest(common.make_suite(BasicWriteTestCase))
        theSuite.addTest(common.make_suite(BasicWrite2TestCase))
        theSuite.addTest(common.make_suite(BasicWrite3TestCase))
        theSuite.addTest(common.make_suite(BasicWrite4TestCase))
        theSuite.addTest(common.make_suite(BasicWrite5TestCase))
        theSuite.addTest(common.make_suite(BasicWrite6TestCase))
        theSuite.addTest(common.make_suite(BasicWrite7TestCase))
        theSuite.addTest(common.make_suite(BasicWrite8TestCase))
        theSuite.addTest(common.make_suite(EmptyCArrayTestCase))
        theSuite.addTest(common.make_suite(EmptyCArray2TestCase))
        theSuite.addTest(common.make_suite(SlicesCArrayTestCase))
        theSuite.addTest(common.make_suite(Slices2CArrayTestCase))
        theSuite.addTest(common.make_suite(EllipsisCArrayTestCase))
        theSuite.addTest(common.make_suite(Ellipsis2CArrayTestCase))
        theSuite.addTest(common.make_suite(Ellipsis3CArrayTestCase))
        theSuite.addTest(common.make_suite(ZlibComprTestCase))
        theSuite.addTest(common.make_suite(ZlibShuffleTestCase))
        theSuite.addTest(common.make_suite(BloscComprTestCase))
        theSuite.addTest(common.make_suite(BloscShuffleTestCase))
        theSuite.addTest(common.make_suite(BloscBitShuffleTestCase))
        theSuite.addTest(common.make_suite(BloscFletcherTestCase))
        theSuite.addTest(common.make_suite(BloscBloscLZTestCase))
        theSuite.addTest(common.make_suite(BloscLZ4TestCase))
        theSuite.addTest(common.make_suite(BloscLZ4HCTestCase))
        theSuite.addTest(common.make_suite(BloscSnappyTestCase))
        theSuite.addTest(common.make_suite(BloscZlibTestCase))
        theSuite.addTest(common.make_suite(BloscZstdTestCase))
        theSuite.addTest(common.make_suite(Blosc2ComprTestCase))
        theSuite.addTest(common.make_suite(Blosc2FletcherTestCase))
        theSuite.addTest(common.make_suite(Blosc2CrossChunkTestCase))
        theSuite.addTest(common.make_suite(Blosc2CrossChunkOptTestCase))
        theSuite.addTest(common.make_suite(Blosc2PastLastChunkTestCase))
        theSuite.addTest(common.make_suite(Blosc2PastLastChunkOptTestCase))
        theSuite.addTest(common.make_suite(Blosc2Ndim3MinChunkOptTestCase))
        theSuite.addTest(common.make_suite(Blosc2Ndim3ChunkOptTestCase))
        theSuite.addTest(common.make_suite(Blosc2Ndim4ChunkOptTestCase))
        theSuite.addTest(common.make_suite(Blosc2NDNoChunkshape))
        theSuite.addTest(common.make_suite(LZOComprTestCase))
        theSuite.addTest(common.make_suite(LZOShuffleTestCase))
        theSuite.addTest(common.make_suite(Bzip2ComprTestCase))
        theSuite.addTest(common.make_suite(Bzip2ShuffleTestCase))
        theSuite.addTest(common.make_suite(FloatTypeTestCase))
        theSuite.addTest(common.make_suite(ComplexTypeTestCase))
        theSuite.addTest(common.make_suite(StringTestCase))
        theSuite.addTest(common.make_suite(String2TestCase))
        theSuite.addTest(common.make_suite(StringComprTestCase))
        theSuite.addTest(common.make_suite(Int8TestCase))
        theSuite.addTest(common.make_suite(Int16TestCase))
        theSuite.addTest(common.make_suite(Int32TestCase))
        theSuite.addTest(common.make_suite(Float16TestCase))
        theSuite.addTest(common.make_suite(Float32TestCase))
        theSuite.addTest(common.make_suite(Float64TestCase))
        theSuite.addTest(common.make_suite(Float96TestCase))
        theSuite.addTest(common.make_suite(Float128TestCase))
        theSuite.addTest(common.make_suite(Complex64TestCase))
        theSuite.addTest(common.make_suite(Complex128TestCase))
        theSuite.addTest(common.make_suite(Complex192TestCase))
        theSuite.addTest(common.make_suite(Complex256TestCase))
        theSuite.addTest(common.make_suite(ComprTestCase))
        theSuite.addTest(common.make_suite(OffsetStrideTestCase))
        theSuite.addTest(common.make_suite(Fletcher32TestCase))
        theSuite.addTest(common.make_suite(AllFiltersTestCase))
        theSuite.addTest(common.make_suite(ReadOutArgumentTests))
        theSuite.addTest(common.make_suite(SizeOnDiskInMemoryPropertyTestCase))
        theSuite.addTest(common.make_suite(CloseCopyTestCase))
        theSuite.addTest(common.make_suite(OpenCopyTestCase))
        theSuite.addTest(common.make_suite(CopyIndex1TestCase))
        theSuite.addTest(common.make_suite(CopyIndex2TestCase))
        theSuite.addTest(common.make_suite(CopyIndex3TestCase))
        theSuite.addTest(common.make_suite(CopyIndex4TestCase))
        theSuite.addTest(common.make_suite(CopyIndex5TestCase))
        theSuite.addTest(common.make_suite(BigArrayTestCase))
        theSuite.addTest(common.make_suite(DfltAtomNoReopen))
        theSuite.addTest(common.make_suite(DfltAtomReopen))
        theSuite.addTest(common.make_suite(AtomDefaultReprNoReopen))
        theSuite.addTest(common.make_suite(AtomDefaultReprReopen))
        theSuite.addTest(common.make_suite(TruncateTestCase))
        theSuite.addTest(common.make_suite(MDAtomNoReopen))
        theSuite.addTest(common.make_suite(MDAtomReopen))
        theSuite.addTest(common.make_suite(MDLargeAtomNoReopen))
        theSuite.addTest(common.make_suite(MDLargeAtomReopen))
        theSuite.addTest(common.make_suite(AccessClosedTestCase))
        theSuite.addTest(common.make_suite(TestCreateCArrayArgs))
    if common.heavy:
        theSuite.addTest(common.make_suite(Slices3CArrayTestCase))
        theSuite.addTest(common.make_suite(Slices4CArrayTestCase))
        theSuite.addTest(common.make_suite(Ellipsis4CArrayTestCase))
        theSuite.addTest(common.make_suite(Ellipsis5CArrayTestCase))
        theSuite.addTest(common.make_suite(Ellipsis6CArrayTestCase))
        theSuite.addTest(common.make_suite(Ellipsis7CArrayTestCase))
        theSuite.addTest(common.make_suite(MD3WriteTestCase))
        theSuite.addTest(common.make_suite(MD5WriteTestCase))
        theSuite.addTest(common.make_suite(MD6WriteTestCase))
        theSuite.addTest(common.make_suite(MD7WriteTestCase))
        theSuite.addTest(common.make_suite(MD10WriteTestCase))
        theSuite.addTest(common.make_suite(CopyIndex6TestCase))
        theSuite.addTest(common.make_suite(CopyIndex7TestCase))
        theSuite.addTest(common.make_suite(CopyIndex8TestCase))
        theSuite.addTest(common.make_suite(CopyIndex9TestCase))
        theSuite.addTest(common.make_suite(CopyIndex10TestCase))
        theSuite.addTest(common.make_suite(CopyIndex11TestCase))
        theSuite.addTest(common.make_suite(CopyIndex12TestCase))
        theSuite.addTest(common.make_suite(Rows64bitsTestCase1))
        theSuite.addTest(common.make_suite(Rows64bitsTestCase2))

    return theSuite


if __name__ == "__main__":
    common.parse_argv(sys.argv)
    common.print_versions()
    common.unittest.main(defaultTest="suite")
