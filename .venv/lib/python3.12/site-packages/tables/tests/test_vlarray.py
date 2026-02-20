import sys

import numpy as np

import tables as tb
from tables.tests import common


class C:
    c = (3, 4.5)


class BasicTestCase(common.TempFileMixin, common.PyTablesTestCase):
    compress = 0
    complib = "zlib"
    shuffle = 0
    bitshuffle = 0
    fletcher32 = 0
    flavor = "numpy"

    def setUp(self):
        super().setUp()

        # Create an instance of an HDF5 Table
        self.rootgroup = self.h5file.root
        self.populateFile()
        self.h5file.close()

    def populateFile(self):
        group = self.rootgroup
        filters = tb.Filters(
            complevel=self.compress,
            complib=self.complib,
            shuffle=self.shuffle,
            bitshuffle=self.bitshuffle,
            fletcher32=self.fletcher32,
        )
        vlarray = self.h5file.create_vlarray(
            group,
            "vlarray1",
            atom=tb.Int32Atom(),
            title="ragged array of ints",
            filters=filters,
            expectedrows=1000,
        )
        vlarray.flavor = self.flavor

        # Fill it with 5 rows
        vlarray.append([1, 2])
        if self.flavor == "numpy":
            vlarray.append(np.array([3, 4, 5], dtype="int32"))
            vlarray.append(np.array([], dtype="int32"))  # Empty entry
        elif self.flavor == "python":
            vlarray.append((3, 4, 5))
            vlarray.append(())  # Empty entry
        vlarray.append([6, 7, 8, 9])
        vlarray.append([10, 11, 12, 13, 14])

    def test00_attributes(self):
        self.h5file = tb.open_file(self.h5fname, "r")
        obj = self.h5file.get_node("/vlarray1")

        self.assertEqual(obj.flavor, self.flavor)
        self.assertEqual(obj.shape, (5,))
        self.assertEqual(obj.ndim, 1)
        self.assertEqual(obj.nrows, 5)
        self.assertEqual(obj.atom.type, "int32")

    def test01_read(self):
        """Checking vlarray read."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01_read..." % self.__class__.__name__)

        # Create an instance of an HDF5 Table
        self.h5file = tb.open_file(self.h5fname, "r")
        vlarray = self.h5file.get_node("/vlarray1")

        # Choose a small value for buffer size
        vlarray.nrowsinbuf = 3
        # Read some rows
        row = vlarray.read(0)[0]
        row2 = vlarray.read(2)[0]
        if common.verbose:
            print("Flavor:", vlarray.flavor)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("First row in vlarray ==>", row)

        nrows = 5
        self.assertEqual(nrows, vlarray.nrows)
        if self.flavor == "numpy":
            self.assertEqual(type(row), np.ndarray)
            self.assertTrue(
                common.allequal(
                    row, np.array([1, 2], dtype="int32"), self.flavor
                )
            )
            self.assertTrue(
                common.allequal(row2, np.array([], dtype="int32"), self.flavor)
            )
        elif self.flavor == "python":
            self.assertEqual(row, [1, 2])
            self.assertEqual(row2, [])
        self.assertEqual(len(row), 2)

        # Check filters:
        if self.compress != vlarray.filters.complevel and common.verbose:
            print("Error in compress. Class:", self.__class__.__name__)
            print("self, vlarray:", self.compress, vlarray.filters.complevel)
        self.assertEqual(vlarray.filters.complevel, self.compress)
        if self.compress > 0 and tb.which_lib_version(self.complib):
            self.assertEqual(vlarray.filters.complib, self.complib)
        if self.shuffle != vlarray.filters.shuffle and common.verbose:
            print("Error in shuffle. Class:", self.__class__.__name__)
            print("self, vlarray:", self.shuffle, vlarray.filters.shuffle)
        self.assertEqual(self.shuffle, vlarray.filters.shuffle)
        if self.bitshuffle != vlarray.filters.bitshuffle and common.verbose:
            print("Error in shuffle. Class:", self.__class__.__name__)
            print(
                "self, vlarray:", self.bitshuffle, vlarray.filters.bitshuffle
            )
        self.assertEqual(self.shuffle, vlarray.filters.shuffle)
        if self.fletcher32 != vlarray.filters.fletcher32 and common.verbose:
            print("Error in fletcher32. Class:", self.__class__.__name__)
            print(
                "self, vlarray:", self.fletcher32, vlarray.filters.fletcher32
            )
        self.assertEqual(self.fletcher32, vlarray.filters.fletcher32)

    def test02a_getitem(self):
        """Checking vlarray __getitem__ (slices)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02a_getitem..." % self.__class__.__name__)

        # Create an instance of an HDF5 Table
        self.h5file = tb.open_file(self.h5fname, "r")
        vlarray = self.h5file.get_node("/vlarray1")

        rows = [[1, 2], [3, 4, 5], [], [6, 7, 8, 9], [10, 11, 12, 13, 14]]

        slices = [
            slice(None, None, None),
            slice(1, 1, 1),
            slice(30, None, None),
            slice(0, None, None),
            slice(3, None, 1),
            slice(3, None, 2),
            slice(None, 1, None),
            slice(None, 2, 1),
            slice(None, 30, 2),
            slice(None, None, 1),
            slice(None, None, 2),
            slice(None, None, 3),
        ]
        for slc in slices:
            # Read the rows in slc
            rows2 = vlarray[slc]
            rows1 = rows[slc]
            rows1f = []
            if common.verbose:
                print("Flavor:", vlarray.flavor)
                print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
                print("Original rows ==>", rows1)
                print("Rows read in vlarray ==>", rows2)

            if self.flavor == "numpy":
                for val in rows1:
                    rows1f.append(np.array(val, dtype="int32"))
                for i in range(len(rows1f)):
                    self.assertTrue(
                        common.allequal(rows2[i], rows1f[i], self.flavor)
                    )
            elif self.flavor == "python":
                self.assertEqual(rows2, rows1)

    def test02b_getitem(self):
        """Checking vlarray __getitem__ (scalars)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02b_getitem..." % self.__class__.__name__)

        if self.flavor != "numpy":
            # This test is only valid for NumPy
            return

        # Create an instance of an HDF5 Table
        self.h5file = tb.open_file(self.h5fname, "r")
        vlarray = self.h5file.get_node("/vlarray1")

        # Get a numpy array of objects
        rows = np.array(vlarray[:], dtype=object)

        for slc in [0, np.array(1), 2, np.array([3]), [4]]:
            # Read the rows in slc
            rows2 = vlarray[slc]
            rows1 = rows[slc]
            if common.verbose:
                print("Flavor:", vlarray.flavor)
                print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
                print("Original rows ==>", rows1)
                print("Rows read in vlarray ==>", rows2)

            for i in range(len(rows1)):
                self.assertTrue(
                    common.allequal(rows2[i], rows1[i], self.flavor)
                )

    def test03_append(self):
        """Checking vlarray append."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03_append..." % self.__class__.__name__)

        # Create an instance of an HDF5 Table
        self.h5file = tb.open_file(self.h5fname, "a")
        vlarray = self.h5file.get_node("/vlarray1")

        # Append a new row
        vlarray.append([7, 8, 9, 10])

        # Choose a small value for buffer size
        vlarray.nrowsinbuf = 3

        # Read some rows:
        row1 = vlarray[0]
        row2 = vlarray[2]
        row3 = vlarray[-1]
        if common.verbose:
            print("Flavor:", vlarray.flavor)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("First row in vlarray ==>", row1)

        nrows = 6
        self.assertEqual(nrows, vlarray.nrows)
        if self.flavor == "numpy":
            self.assertEqual(type(row1), type(np.array([1, 2])))
            self.assertTrue(
                common.allequal(
                    row1, np.array([1, 2], dtype="int32"), self.flavor
                )
            )
            self.assertTrue(
                common.allequal(row2, np.array([], dtype="int32"), self.flavor)
            )
            self.assertTrue(
                common.allequal(
                    row3, np.array([7, 8, 9, 10], dtype="int32"), self.flavor
                )
            )
        elif self.flavor == "python":
            self.assertEqual(row1, [1, 2])
            self.assertEqual(row2, [])
            self.assertEqual(row3, [7, 8, 9, 10])
        self.assertEqual(len(row3), 4)

    def test04_get_row_size(self):
        """Checking get_row_size method."""

        self.h5file = tb.open_file(self.h5fname, "a")
        vlarray = self.h5file.get_node("/vlarray1")

        self.assertEqual(vlarray.get_row_size(0), 2 * vlarray.atom.size)
        self.assertEqual(vlarray.get_row_size(1), 3 * vlarray.atom.size)
        self.assertEqual(vlarray.get_row_size(2), 0 * vlarray.atom.size)
        self.assertEqual(vlarray.get_row_size(3), 4 * vlarray.atom.size)
        self.assertEqual(vlarray.get_row_size(4), 5 * vlarray.atom.size)


class BasicNumPyTestCase(BasicTestCase):
    flavor = "numpy"


class BasicPythonTestCase(BasicTestCase):
    flavor = "python"


class ZlibComprTestCase(BasicTestCase):
    compress = 1
    complib = "zlib"


@common.unittest.skipIf(
    not common.blosc_avail, "BLOSC compression library not available"
)
class BloscComprTestCase(BasicTestCase):
    compress = 9
    shuffle = 0
    complib = "blosc"


@common.unittest.skipIf(
    not common.blosc_avail, "BLOSC compression library not available"
)
class BloscShuffleComprTestCase(BasicTestCase):
    compress = 6
    shuffle = 1
    complib = "blosc"


@common.unittest.skipIf(
    not common.blosc_avail, "BLOSC compression library not available"
)
class BloscBitShuffleComprTestCase(BasicTestCase):
    compress = 9
    bitshuffle = 1
    complib = "blosc"


@common.unittest.skipIf(
    not common.blosc_avail, "BLOSC compression library not available"
)
class BloscBloscLZComprTestCase(BasicTestCase):
    compress = 9
    shuffle = 1
    complib = "blosc:blosclz"


@common.unittest.skipIf(
    not common.blosc_avail, "BLOSC compression library not available"
)
@common.unittest.skipIf(
    "lz4" not in tb.blosc_compressor_list(), "lz4 required"
)
class BloscLZ4ComprTestCase(BasicTestCase):
    compress = 9
    shuffle = 1
    complib = "blosc:lz4"


@common.unittest.skipIf(
    not common.blosc_avail, "BLOSC compression library not available"
)
@common.unittest.skipIf(
    "lz4" not in tb.blosc_compressor_list(), "lz4 required"
)
class BloscLZ4HCComprTestCase(BasicTestCase):
    compress = 9
    shuffle = 1
    complib = "blosc:lz4hc"


@common.unittest.skipIf(
    not common.blosc_avail, "BLOSC compression library not available"
)
@common.unittest.skipIf(
    "snappy" not in tb.blosc_compressor_list(), "snappy required"
)
class BloscSnappyComprTestCase(BasicTestCase):
    compress = 9
    shuffle = 1
    complib = "blosc:snappy"


@common.unittest.skipIf(
    not common.blosc_avail, "BLOSC compression library not available"
)
@common.unittest.skipIf(
    "zlib" not in tb.blosc_compressor_list(), "zlib required"
)
class BloscZlibComprTestCase(BasicTestCase):
    compress = 9
    shuffle = 1
    complib = "blosc:zlib"


@common.unittest.skipIf(
    not common.blosc_avail, "BLOSC compression library not available"
)
@common.unittest.skipIf(
    "zstd" not in tb.blosc_compressor_list(), "zstd required"
)
class BloscZstdComprTestCase(BasicTestCase):
    compress = 9
    shuffle = 1
    complib = "blosc:zstd"


@common.unittest.skipIf(
    not common.lzo_avail, "LZO compression library not available"
)
class LZOComprTestCase(BasicTestCase):
    compress = 1
    complib = "lzo"


@common.unittest.skipIf(
    not common.bzip2_avail, "BZIP2 compression library not available"
)
class Bzip2ComprTestCase(BasicTestCase):
    compress = 1
    complib = "bzip2"


class ShuffleComprTestCase(BasicTestCase):
    compress = 1
    shuffle = 1


class TypesTestCase(common.TempFileMixin, common.PyTablesTestCase):
    open_mode = "w"
    compress = 0
    complib = "zlib"  # Default compression library

    def test01_StringAtom(self):
        """Checking vlarray with NumPy string atoms ('numpy' flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01_StringAtom..." % self.__class__.__name__)

        vlarray = self.h5file.create_vlarray(
            "/",
            "stringAtom",
            atom=tb.StringAtom(itemsize=3),
            title="Ragged array of strings",
        )
        vlarray.flavor = "numpy"
        vlarray.append(np.array(["1", "12", "123", "1234", "12345"]))
        vlarray.append(np.array(["1", "12345"]))

        if self.reopen:
            name = vlarray._v_pathname
            self._reopen()
            vlarray = self.h5file.get_node(name)

        # Read all the rows:
        row = vlarray.read()
        if common.verbose:
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("First row in vlarray ==>", row[0])

        self.assertEqual(vlarray.nrows, 2)
        np.testing.assert_array_equal(
            row[0], np.array(["1", "12", "123", "123", "123"], "S")
        )
        np.testing.assert_array_equal(row[1], np.array(["1", "123"], "S"))
        self.assertEqual(len(row[0]), 5)
        self.assertEqual(len(row[1]), 2)

    def test01a_StringAtom(self):
        """Checking vlarray with NumPy string atoms ('numpy' flavor,
        strided)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01a_StringAtom..." % self.__class__.__name__)

        vlarray = self.h5file.create_vlarray(
            "/",
            "stringAtom",
            atom=tb.StringAtom(itemsize=3),
            title="Ragged array of strings",
        )
        vlarray.flavor = "numpy"
        vlarray.append(np.array(["1", "12", "123", "1234", "12345"][::2]))
        vlarray.append(np.array(["1", "12345", "2", "321"])[::3])

        if self.reopen:
            name = vlarray._v_pathname
            self._reopen()
            vlarray = self.h5file.get_node(name)

        # Read all the rows:
        row = vlarray.read()
        if common.verbose:
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("First row in vlarray ==>", row[0])

        self.assertEqual(vlarray.nrows, 2)
        np.testing.assert_array_equal(
            row[0], np.array(["1", "123", "123"], "S")
        )
        np.testing.assert_array_equal(row[1], np.array(["1", "321"], "S"))
        self.assertEqual(len(row[0]), 3)
        self.assertEqual(len(row[1]), 2)

    def test01a_2_StringAtom(self):
        """Checking vlarray with NumPy string atoms (NumPy flavor, no conv)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test01a_2_StringAtom..." % self.__class__.__name__
            )

        vlarray = self.h5file.create_vlarray(
            "/",
            "stringAtom",
            atom=tb.StringAtom(itemsize=3),
            title="Ragged array of strings",
        )
        vlarray.flavor = "numpy"
        vlarray.append(np.array(["1", "12", "123", "123"]))
        vlarray.append(np.array(["1", "2", "321"]))

        if self.reopen:
            name = vlarray._v_pathname
            self._reopen()
            vlarray = self.h5file.get_node(name)

        # Read all the rows:
        row = vlarray.read()
        if common.verbose:
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("First row in vlarray ==>", row[0])

        self.assertEqual(vlarray.nrows, 2)
        np.testing.assert_array_equal(
            row[0], np.array(["1", "12", "123", "123"], "S")
        )
        np.testing.assert_array_equal(row[1], np.array(["1", "2", "321"], "S"))
        self.assertEqual(len(row[0]), 4)
        self.assertEqual(len(row[1]), 3)

    def test01b_StringAtom(self):
        """Checking vlarray with NumPy string atoms (python flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01b_StringAtom..." % self.__class__.__name__)

        vlarray = self.h5file.create_vlarray(
            "/",
            "stringAtom2",
            atom=tb.StringAtom(itemsize=3),
            title="Ragged array of strings",
        )
        vlarray.flavor = "python"
        vlarray.append(["1", "12", "123", "1234", "12345"])
        vlarray.append(["1", "12345"])

        if self.reopen:
            name = vlarray._v_pathname
            self._reopen()
            vlarray = self.h5file.get_node(name)

        # Read all the rows:
        row = vlarray.read()
        if common.verbose:
            print("Testing String flavor")
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("First row in vlarray ==>", row[0])

        self.assertEqual(vlarray.nrows, 2)
        self.assertEqual(row[0], [b"1", b"12", b"123", b"123", b"123"])
        self.assertEqual(row[1], [b"1", b"123"])
        self.assertEqual(len(row[0]), 5)
        self.assertEqual(len(row[1]), 2)

    def test01c_StringAtom(self):
        """Checking updating vlarray with NumPy string atoms
        ('numpy' flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01c_StringAtom..." % self.__class__.__name__)

        vlarray = self.h5file.create_vlarray(
            "/",
            "stringAtom",
            atom=tb.StringAtom(itemsize=3),
            title="Ragged array of strings",
        )
        vlarray.flavor = "numpy"
        vlarray.append(np.array(["1", "12", "123", "1234", "12345"]))
        vlarray.append(np.array(["1", "12345"]))

        # Modify the rows
        vlarray[0] = np.array(["1", "123", "12", "", "12345"])
        vlarray[1] = np.array(["44", "4"])  # This should work as well

        if self.reopen:
            name = vlarray._v_pathname
            self._reopen()
            vlarray = self.h5file.get_node(name)

        # Read all the rows:
        row = vlarray.read()
        if common.verbose:
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("First row in vlarray ==>", row[0])

        self.assertEqual(vlarray.nrows, 2)
        self.assertTrue(
            common.allequal(
                row[0], np.array([b"1", b"123", b"12", b"", b"123"])
            )
        )
        self.assertTrue(
            common.allequal(row[1], np.array(["44", "4"], dtype="S3"))
        )
        self.assertEqual(len(row[0]), 5)
        self.assertEqual(len(row[1]), 2)

    def test01d_StringAtom(self):
        """Checking updating vlarray with string atoms (String flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01d_StringAtom..." % self.__class__.__name__)

        vlarray = self.h5file.create_vlarray(
            "/",
            "stringAtom2",
            atom=tb.StringAtom(itemsize=3),
            title="Ragged array of strings",
        )
        vlarray.flavor = "python"
        vlarray.append(["1", "12", "123", "1234", "12345"])
        vlarray.append(["1", "12345"])

        # Modify the rows
        vlarray[0] = ["1", "123", "12", "", "12345"]
        vlarray[1] = ["44", "4"]

        if self.reopen:
            name = vlarray._v_pathname
            self._reopen()
            vlarray = self.h5file.get_node(name)

        # Read all the rows:
        row = vlarray.read()
        if common.verbose:
            print("Testing String flavor")
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("First row in vlarray ==>", row[0])

        self.assertEqual(vlarray.nrows, 2)
        self.assertEqual(row[0], [b"1", b"123", b"12", b"", b"123"])
        self.assertEqual(row[1], [b"44", b"4"])
        self.assertEqual(len(row[0]), 5)
        self.assertEqual(len(row[1]), 2)

    def test02_BoolAtom(self):
        """Checking vlarray with boolean atoms."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02_BoolAtom..." % self.__class__.__name__)

        vlarray = self.h5file.create_vlarray(
            "/",
            "BoolAtom",
            atom=tb.BoolAtom(),
            title="Ragged array of Booleans",
        )
        vlarray.append([1, 0, 3])
        vlarray.append([1, 0])

        if self.reopen:
            name = vlarray._v_pathname
            self._reopen()
            vlarray = self.h5file.get_node(name)

        # Read all the rows:
        row = vlarray.read()
        if common.verbose:
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("First row in vlarray ==>", row[0])

        self.assertEqual(vlarray.nrows, 2)
        self.assertTrue(
            common.allequal(row[0], np.array([1, 0, 1], dtype="bool"))
        )
        self.assertTrue(
            common.allequal(row[1], np.array([1, 0], dtype="bool"))
        )
        self.assertEqual(len(row[0]), 3)
        self.assertEqual(len(row[1]), 2)

    def test02b_BoolAtom(self):
        """Checking setting vlarray with boolean atoms."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02b_BoolAtom..." % self.__class__.__name__)

        vlarray = self.h5file.create_vlarray(
            "/",
            "BoolAtom",
            atom=tb.BoolAtom(),
            title="Ragged array of Booleans",
        )
        vlarray.append([1, 0, 3])
        vlarray.append([1, 0])

        # Modify the rows
        vlarray[0] = (0, 1, 3)
        vlarray[1] = (0, 1)

        if self.reopen:
            name = vlarray._v_pathname
            self._reopen()
            vlarray = self.h5file.get_node(name)

        # Read all the rows:
        row = vlarray.read()
        if common.verbose:
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("First row in vlarray ==>", row[0])

        self.assertEqual(vlarray.nrows, 2)
        self.assertTrue(
            common.allequal(row[0], np.array([0, 1, 1], dtype="bool"))
        )
        self.assertTrue(
            common.allequal(row[1], np.array([0, 1], dtype="bool"))
        )
        self.assertEqual(len(row[0]), 3)
        self.assertEqual(len(row[1]), 2)

    def test03_IntAtom(self):
        """Checking vlarray with integer atoms."""

        ttypes = [
            "int8",
            "uint8",
            "int16",
            "uint16",
            "int32",
            "uint32",
            "int64",
            # "UInt64",  # Unavailable in some platforms
        ]
        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03_IntAtom..." % self.__class__.__name__)

        for atype in ttypes:
            vlarray = self.h5file.create_vlarray(
                "/", atype, atom=tb.Atom.from_sctype(atype)
            )
            vlarray.append([1, 2, 3])
            vlarray.append([1, 0])

            if self.reopen:
                name = vlarray._v_pathname
                self._reopen(mode="a")
                vlarray = self.h5file.get_node(name)

            # Read all the rows:
            row = vlarray.read()
            if common.verbose:
                print("Testing type:", atype)
                print("Object read:", row)
                print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
                print("First row in vlarray ==>", row[0])

            self.assertEqual(vlarray.nrows, 2)
            self.assertTrue(
                common.allequal(row[0], np.array([1, 2, 3], dtype=atype))
            )
            self.assertTrue(
                common.allequal(row[1], np.array([1, 0], dtype=atype))
            )
            self.assertEqual(len(row[0]), 3)
            self.assertEqual(len(row[1]), 2)

    def test03a_IntAtom(self):
        """Checking vlarray with integer atoms (byteorder swapped)"""

        ttypes = {
            "int8": np.int8,
            "uint8": np.uint8,
            "int16": np.int16,
            "uint16": np.uint16,
            "int32": np.int32,
            "uint32": np.uint32,
            "int64": np.int64,
            # "uint64": np.int64,  # Unavailable in some platforms
        }
        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03a_IntAtom..." % self.__class__.__name__)

        for atype in ttypes:
            vlarray = self.h5file.create_vlarray(
                "/", atype, atom=tb.Atom.from_sctype(ttypes[atype])
            )
            a0 = np.array([1, 2, 3], dtype=atype)
            a0 = a0.byteswap()
            a0 = a0.view(a0.dtype.newbyteorder())
            vlarray.append(a0)
            a1 = np.array([1, 0], dtype=atype)
            a1 = a1.byteswap()
            a1 = a1.view(a1.dtype.newbyteorder())
            vlarray.append(a1)

            if self.reopen:
                name = vlarray._v_pathname
                self._reopen(mode="a")
                vlarray = self.h5file.get_node(name)

            # Read all the rows:
            row = vlarray.read()
            if common.verbose:
                print("Testing type:", atype)
                print("Object read:", row)
                print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
                print("First row in vlarray ==>", row[0])

            self.assertEqual(vlarray.nrows, 2)
            self.assertTrue(
                common.allequal(
                    row[0], np.array([1, 2, 3], dtype=ttypes[atype])
                )
            )
            self.assertTrue(
                common.allequal(row[1], np.array([1, 0], dtype=ttypes[atype]))
            )
            self.assertEqual(len(row[0]), 3)
            self.assertEqual(len(row[1]), 2)

    def test03b_IntAtom(self):
        """Checking updating vlarray with integer atoms."""

        ttypes = [
            "int8",
            "uint8",
            "int16",
            "uint16",
            "int32",
            "uint32",
            "int64",
            # "UInt64",  # Unavailable in some platforms
        ]
        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03_IntAtom..." % self.__class__.__name__)

        for atype in ttypes:
            vlarray = self.h5file.create_vlarray(
                "/", atype, atom=tb.Atom.from_sctype(atype)
            )
            vlarray.append([1, 2, 3])
            vlarray.append([1, 0])

            # Modify rows
            vlarray[0] = (3, 2, 1)
            vlarray[1] = (0, 1)

            if self.reopen:
                name = vlarray._v_pathname
                self._reopen(mode="a")
                vlarray = self.h5file.get_node(name)

            # Read all the rows:
            row = vlarray.read()
            if common.verbose:
                print("Testing type:", atype)
                print("Object read:", row)
                print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
                print("First row in vlarray ==>", row[0])

            self.assertEqual(vlarray.nrows, 2)
            self.assertTrue(
                common.allequal(row[0], np.array([3, 2, 1], dtype=atype))
            )
            self.assertTrue(
                common.allequal(row[1], np.array([0, 1], dtype=atype))
            )
            self.assertEqual(len(row[0]), 3)
            self.assertEqual(len(row[1]), 2)

    def test03c_IntAtom(self):
        """Checking updating vlarray with integer atoms (byteorder swapped)"""

        ttypes = {
            "int8": np.int8,
            "uint8": np.uint8,
            "int16": np.int16,
            "uint16": np.uint16,
            "int32": np.int32,
            "uint32": np.uint32,
            "int64": np.int64,
            # "uint64": np.int64,  # Unavailable in some platforms
        }
        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03c_IntAtom..." % self.__class__.__name__)

        for atype in ttypes:
            vlarray = self.h5file.create_vlarray(
                "/", atype, atom=tb.Atom.from_sctype(ttypes[atype])
            )
            a0 = np.array([1, 2, 3], dtype=atype)
            vlarray.append(a0)
            a1 = np.array([1, 0], dtype=atype)
            vlarray.append(a1)

            # Modify rows
            a0 = np.array([3, 2, 1], dtype=atype)
            a0 = a0.byteswap()
            a0 = a0.view(a0.dtype.newbyteorder())
            vlarray[0] = a0
            a1 = np.array([0, 1], dtype=atype)
            a1 = a1.byteswap()
            a1 = a1.view(a1.dtype.newbyteorder())
            vlarray[1] = a1

            if self.reopen:
                name = vlarray._v_pathname
                self._reopen(mode="a")
                vlarray = self.h5file.get_node(name)

            # Read all the rows:
            row = vlarray.read()
            if common.verbose:
                print("Testing type:", atype)
                print("Object read:", row)
                print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
                print("First row in vlarray ==>", row[0])

            self.assertEqual(vlarray.nrows, 2)
            self.assertTrue(
                common.allequal(
                    row[0], np.array([3, 2, 1], dtype=ttypes[atype])
                )
            )
            self.assertTrue(
                common.allequal(row[1], np.array([0, 1], dtype=ttypes[atype]))
            )
            self.assertEqual(len(row[0]), 3)
            self.assertEqual(len(row[1]), 2)

    def test03d_IntAtom(self):
        """Checking updating vlarray with integer atoms (another byteorder)"""

        ttypes = {
            "int8": np.int8,
            "uint8": np.uint8,
            "int16": np.int16,
            "uint16": np.uint16,
            "int32": np.int32,
            "uint32": np.uint32,
            "int64": np.int64,
            # "uint64": np.int64,  # Unavailable in some platforms
        }
        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03d_IntAtom..." % self.__class__.__name__)

        byteorder = {"little": "big", "big": "little"}[sys.byteorder]
        for atype in ttypes:
            vlarray = self.h5file.create_vlarray(
                "/",
                atype,
                atom=tb.Atom.from_sctype(ttypes[atype]),
                byteorder=byteorder,
            )
            a0 = np.array([1, 2, 3], dtype=atype)
            vlarray.append(a0)
            a1 = np.array([1, 0], dtype=atype)
            vlarray.append(a1)

            # Modify rows
            a0 = np.array([3, 2, 1], dtype=atype)
            a0 = a0.byteswap()
            a0 = a0.view(a0.dtype.newbyteorder())
            vlarray[0] = a0
            a1 = np.array([0, 1], dtype=atype)
            a1 = a1.byteswap()
            a1 = a1.view(a1.dtype.newbyteorder())
            vlarray[1] = a1

            if self.reopen:
                name = vlarray._v_pathname
                self._reopen(mode="a")
                vlarray = self.h5file.get_node(name)

            # Read all the rows:
            row = vlarray.read()
            if common.verbose:
                print("Testing type:", atype)
                print("Object read:", row)
                print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
                print("First row in vlarray ==>", row[0])

            byteorder2 = tb.utils.byteorders[row[0].dtype.byteorder]
            if byteorder2 != "irrelevant":
                self.assertEqual(
                    tb.utils.byteorders[row[0].dtype.byteorder], sys.byteorder
                )
                self.assertEqual(vlarray.byteorder, byteorder)
            self.assertEqual(vlarray.nrows, 2)
            self.assertTrue(
                common.allequal(
                    row[0], np.array([3, 2, 1], dtype=ttypes[atype])
                )
            )
            self.assertTrue(
                common.allequal(row[1], np.array([0, 1], dtype=ttypes[atype]))
            )
            self.assertEqual(len(row[0]), 3)
            self.assertEqual(len(row[1]), 2)

    def test04_FloatAtom(self):
        """Checking vlarray with floating point atoms."""

        ttypes = [
            "float32",
            "float64",
        ]
        for name in ("float16", "float96", "float128"):
            atomname = name.capitalize() + "Atom"
            if hasattr(tb, atomname):
                ttypes.append(name)

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test04_FloatAtom..." % self.__class__.__name__)

        for atype in ttypes:
            vlarray = self.h5file.create_vlarray(
                "/", atype, atom=tb.Atom.from_sctype(atype)
            )
            vlarray.append([1.3, 2.2, 3.3])
            vlarray.append([5.96, 0.597])

            if self.reopen:
                name = vlarray._v_pathname
                self._reopen(mode="a")
                vlarray = self.h5file.get_node(name)

            # Read all the rows:
            row = vlarray.read()
            if common.verbose:
                print("Testing type:", atype)
                print("Object read:", row)
                print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
                print("First row in vlarray ==>", row[0])

            self.assertEqual(vlarray.nrows, 2)
            self.assertTrue(
                common.allequal(row[0], np.array([1.3, 2.2, 3.3], atype))
            )
            self.assertTrue(
                common.allequal(row[1], np.array([5.96, 0.597], atype))
            )
            self.assertEqual(len(row[0]), 3)
            self.assertEqual(len(row[1]), 2)

    def test04a_FloatAtom(self):
        """Checking vlarray with float atoms (byteorder swapped)"""

        ttypes = {
            "float32": np.float32,
            "float64": np.float64,
        }
        if hasattr(tb, "Float16Atom"):
            ttypes["float16"] = np.float16
        if hasattr(tb, "Float96Atom"):
            ttypes["float96"] = np.float96
        if hasattr(tb, "Float128Atom"):
            ttypes["float128"] = np.float128

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test04a_FloatAtom..." % self.__class__.__name__)

        for atype in ttypes:
            vlarray = self.h5file.create_vlarray(
                "/", atype, atom=tb.Atom.from_sctype(ttypes[atype])
            )
            a0 = np.array([1.3, 2.2, 3.3], dtype=atype)
            a0 = a0.byteswap()
            a0 = a0.view(a0.dtype.newbyteorder())
            vlarray.append(a0)
            a1 = np.array([5.96, 0.597], dtype=atype)
            a1 = a1.byteswap()
            a1 = a1.view(a1.dtype.newbyteorder())
            vlarray.append(a1)

            if self.reopen:
                name = vlarray._v_pathname
                self._reopen(mode="a")
                vlarray = self.h5file.get_node(name)

            # Read all the rows:
            row = vlarray.read()
            if common.verbose:
                print("Testing type:", atype)
                print("Object read:", row)
                print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
                print("First row in vlarray ==>", row[0])

            self.assertEqual(vlarray.nrows, 2)
            self.assertTrue(
                common.allequal(
                    row[0], np.array([1.3, 2.2, 3.3], dtype=ttypes[atype])
                )
            )
            self.assertTrue(
                common.allequal(
                    row[1], np.array([5.96, 0.597], dtype=ttypes[atype])
                )
            )
            self.assertEqual(len(row[0]), 3)
            self.assertEqual(len(row[1]), 2)

    def test04b_FloatAtom(self):
        """Checking updating vlarray with floating point atoms."""

        ttypes = [
            "float32",
            "float64",
        ]
        for name in ("float16", "float96", "float128"):
            atomname = name.capitalize() + "Atom"
            if hasattr(tb, atomname):
                ttypes.append(name)

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test04b_FloatAtom..." % self.__class__.__name__)

        for atype in ttypes:
            vlarray = self.h5file.create_vlarray(
                "/", atype, atom=tb.Atom.from_sctype(atype)
            )
            vlarray.append([1.3, 2.2, 3.3])
            vlarray.append([5.96, 0.597])

            # Modifiy some rows
            vlarray[0] = (4.3, 2.2, 4.3)
            vlarray[1] = (1.123, 1.1e-3)

            if self.reopen:
                name = vlarray._v_pathname
                self._reopen(mode="a")
                vlarray = self.h5file.get_node(name)

            # Read all the rows:
            row = vlarray.read()
            if common.verbose:
                print("Testing type:", atype)
                print("Object read:", row)
                print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
                print("First row in vlarray ==>", row[0])

            self.assertEqual(vlarray.nrows, 2)
            self.assertTrue(
                common.allequal(row[0], np.array([4.3, 2.2, 4.3], atype))
            )
            self.assertTrue(
                common.allequal(row[1], np.array([1.123, 1.1e-3], atype))
            )
            self.assertEqual(len(row[0]), 3)
            self.assertEqual(len(row[1]), 2)

    def test04c_FloatAtom(self):
        """Checking updating vlarray with float atoms (byteorder swapped)"""

        ttypes = {
            "float32": np.float32,
            "float64": np.float64,
        }
        if hasattr(tb, "Float16Atom"):
            ttypes["float16"] = np.float16
        if hasattr(tb, "Float96Atom"):
            ttypes["float96"] = np.float96
        if hasattr(tb, "Float128Atom"):
            ttypes["float128"] = np.float128

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test04c_FloatAtom..." % self.__class__.__name__)

        for atype in ttypes:
            vlarray = self.h5file.create_vlarray(
                "/", atype, atom=tb.Atom.from_sctype(ttypes[atype])
            )
            a0 = np.array([1.3, 2.2, 3.3], dtype=atype)
            vlarray.append(a0)
            a1 = np.array([1, 0], dtype=atype)
            vlarray.append(a1)

            # Modify rows
            a0 = np.array([4.3, 2.2, 4.3], dtype=atype)
            a0 = a0.byteswap()
            a0 = a0.view(a0.dtype.newbyteorder())
            vlarray[0] = a0
            a1 = np.array([1.123, 1.1e-3], dtype=atype)
            a1 = a1.byteswap()
            a1 = a1.view(a1.dtype.newbyteorder())
            vlarray[1] = a1

            if self.reopen:
                name = vlarray._v_pathname
                self._reopen(mode="a")
                vlarray = self.h5file.get_node(name)

            # Read all the rows:
            row = vlarray.read()
            if common.verbose:
                print("Testing type:", atype)
                print("Object read:", row)
                print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
                print("First row in vlarray ==>", row[0])

            self.assertEqual(vlarray.nrows, 2)
            self.assertTrue(
                common.allequal(
                    row[0], np.array([4.3, 2.2, 4.3], dtype=ttypes[atype])
                )
            )
            self.assertTrue(
                common.allequal(
                    row[1], np.array([1.123, 1.1e-3], dtype=ttypes[atype])
                )
            )
            self.assertEqual(len(row[0]), 3)
            self.assertEqual(len(row[1]), 2)

    def test04d_FloatAtom(self):
        """Checking updating vlarray with float atoms (another byteorder)"""

        ttypes = {
            "float32": np.float32,
            "float64": np.float64,
        }
        if hasattr(tb, "Float16Atom"):
            ttypes["float16"] = np.float16
        if hasattr(tb, "Float96Atom"):
            ttypes["float96"] = np.float96
        if hasattr(tb, "Float128Atom"):
            ttypes["float128"] = np.float128

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test04d_FloatAtom..." % self.__class__.__name__)

        byteorder = {"little": "big", "big": "little"}[sys.byteorder]
        for atype in ttypes:
            vlarray = self.h5file.create_vlarray(
                "/",
                atype,
                atom=tb.Atom.from_sctype(ttypes[atype]),
                byteorder=byteorder,
            )
            a0 = np.array([1.3, 2.2, 3.3], dtype=atype)
            vlarray.append(a0)
            a1 = np.array([1, 0], dtype=atype)
            vlarray.append(a1)

            # Modify rows
            a0 = np.array([4.3, 2.2, 4.3], dtype=atype)
            a0 = a0.byteswap()
            a0 = a0.view(a0.dtype.newbyteorder())
            vlarray[0] = a0
            a1 = np.array([1.123, 1.1e-3], dtype=atype)
            a1 = a1.byteswap()
            a1 = a1.view(a1.dtype.newbyteorder())
            vlarray[1] = a1

            if self.reopen:
                name = vlarray._v_pathname
                self._reopen(mode="a")
                vlarray = self.h5file.get_node(name)

            # Read all the rows:
            row = vlarray.read()
            if common.verbose:
                print("Testing type:", atype)
                print("Object read:", row)
                print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
                print("First row in vlarray ==>", row[0])

            self.assertEqual(vlarray.byteorder, byteorder)
            self.assertEqual(
                tb.utils.byteorders[row[0].dtype.byteorder], sys.byteorder
            )
            self.assertEqual(vlarray.nrows, 2)
            self.assertTrue(
                common.allequal(
                    row[0], np.array([4.3, 2.2, 4.3], dtype=ttypes[atype])
                )
            )
            self.assertTrue(
                common.allequal(
                    row[1], np.array([1.123, 1.1e-3], dtype=ttypes[atype])
                )
            )
            self.assertEqual(len(row[0]), 3)
            self.assertEqual(len(row[1]), 2)

    def test04_ComplexAtom(self):
        """Checking vlarray with numerical complex atoms."""

        ttypes = [
            "complex64",
            "complex128",
        ]

        if hasattr(tb, "Complex192Atom"):
            ttypes.append("complex192")
        if hasattr(tb, "Complex256Atom"):
            ttypes.append("complex256")

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test04_ComplexAtom..." % self.__class__.__name__)

        for atype in ttypes:
            vlarray = self.h5file.create_vlarray(
                "/", atype, atom=tb.Atom.from_sctype(atype)
            )
            vlarray.append([(1.3 + 0j), (0 + 2.2j), (3.3 + 3.3j)])
            vlarray.append([(0 - 5.96j), (0.597 + 0j)])

            if self.reopen:
                name = vlarray._v_pathname
                self._reopen(mode="a")
                vlarray = self.h5file.get_node(name)

            # Read all the rows:
            row = vlarray.read()
            if common.verbose:
                print("Testing type:", atype)
                print("Object read:", row)
                print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
                print("First row in vlarray ==>", row[0])

            self.assertEqual(vlarray.nrows, 2)
            self.assertTrue(
                common.allequal(
                    row[0],
                    np.array([(1.3 + 0j), (0 + 2.2j), (3.3 + 3.3j)], atype),
                )
            )
            self.assertTrue(
                common.allequal(
                    row[1], np.array([(0 - 5.96j), (0.597 + 0j)], atype)
                )
            )
            self.assertEqual(len(row[0]), 3)
            self.assertEqual(len(row[1]), 2)

    def test04b_ComplexAtom(self):
        """Checking modifying vlarray with numerical complex atoms."""

        ttypes = [
            "complex64",
            "complex128",
        ]

        if hasattr(tb, "Complex192Atom"):
            ttypes.append("complex192")
        if hasattr(tb, "Complex256Atom"):
            ttypes.append("complex256")

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test04b_ComplexAtom..." % self.__class__.__name__
            )

        for atype in ttypes:
            vlarray = self.h5file.create_vlarray(
                "/", atype, atom=tb.Atom.from_sctype(atype)
            )
            vlarray.append([(1.3 + 0j), (0 + 2.2j), (3.3 + 3.3j)])
            vlarray.append([(0 - 5.96j), (0.597 + 0j)])

            # Modify the rows
            vlarray[0] = ((1.4 + 0j), (0 + 4.2j), (3.3 + 4.3j))
            vlarray[1] = ((4 - 5.96j), (0.597 + 4j))

            if self.reopen:
                name = vlarray._v_pathname
                self._reopen(mode="a")
                vlarray = self.h5file.get_node(name)

            # Read all the rows:
            row = vlarray.read()
            if common.verbose:
                print("Testing type:", atype)
                print("Object read:", row)
                print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
                print("First row in vlarray ==>", row[0])

            self.assertEqual(vlarray.nrows, 2)
            self.assertTrue(
                common.allequal(
                    row[0],
                    np.array([(1.4 + 0j), (0 + 4.2j), (3.3 + 4.3j)], atype),
                )
            )
            self.assertTrue(
                common.allequal(
                    row[1], np.array([(4 - 5.96j), (0.597 + 4j)], atype)
                )
            )
            self.assertEqual(len(row[0]), 3)
            self.assertEqual(len(row[1]), 2)

    def test05_VLStringAtom(self):
        """Checking vlarray with variable length strings."""

        # Skip the test if the default encoding has been mangled.
        if sys.getdefaultencoding() != "ascii":
            return

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test05_VLStringAtom..." % self.__class__.__name__
            )

        vlarray = self.h5file.create_vlarray(
            "/", "VLStringAtom", atom=tb.VLStringAtom()
        )
        vlarray.append(b"asd")
        vlarray.append(b"asd\xe4")
        vlarray.append(b"aaana")
        vlarray.append(b"")
        # Check for ticket #62.
        self.assertRaises(TypeError, vlarray.append, [b"foo", b"bar"])
        # `VLStringAtom` makes no encoding assumptions.  See ticket #51.
        self.assertRaises(UnicodeEncodeError, vlarray.append, "asd\xe4")

        if self.reopen:
            name = vlarray._v_pathname
            self._reopen()
            vlarray = self.h5file.get_node(name)

        # Read all the rows:
        row = vlarray.read()
        if common.verbose:
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("First row in vlarray ==>", row[0])

        self.assertEqual(vlarray.nrows, 4)
        self.assertEqual(row[0], b"asd")
        self.assertEqual(row[1], b"asd\xe4")
        self.assertEqual(row[2], b"aaana")
        self.assertEqual(row[3], b"")
        self.assertEqual(len(row[0]), 3)
        self.assertEqual(len(row[1]), 4)
        self.assertEqual(len(row[2]), 5)
        self.assertEqual(len(row[3]), 0)

    def test05b_VLStringAtom(self):
        """Checking updating vlarray with variable length strings."""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test05b_VLStringAtom..." % self.__class__.__name__
            )

        vlarray = self.h5file.create_vlarray(
            "/", "VLStringAtom", atom=tb.VLStringAtom()
        )
        vlarray.append(b"asd")
        vlarray.append(b"aaana")

        # Modify values
        vlarray[0] = b"as4"
        vlarray[1] = b"aaanc"
        self.assertRaises(ValueError, vlarray.__setitem__, 1, b"shrt")
        self.assertRaises(ValueError, vlarray.__setitem__, 1, b"toolong")

        if self.reopen:
            name = vlarray._v_pathname
            self._reopen()
            vlarray = self.h5file.get_node(name)

        # Read all the rows:
        row = vlarray.read()
        if common.verbose:
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("First row in vlarray ==>", repr(row[0]))
            print("Second row in vlarray ==>", repr(row[1]))

        self.assertEqual(vlarray.nrows, 2)
        self.assertEqual(row[0], b"as4")
        self.assertEqual(row[1], b"aaanc")
        self.assertEqual(len(row[0]), 3)
        self.assertEqual(len(row[1]), 5)

    def test06a_Object(self):
        """Checking vlarray with object atoms."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test06a_Object..." % self.__class__.__name__)

        vlarray = self.h5file.create_vlarray(
            "/", "Object", atom=tb.ObjectAtom()
        )
        vlarray.append(
            [[1, 2, 3], "aaa", "aaa\xef\xbf\xbd\xef\xbf\xbd\xef\xbf\xbd"]
        )
        vlarray.append([3, 4, C()])
        vlarray.append(42)

        if self.reopen:
            name = vlarray._v_pathname
            self._reopen()
            vlarray = self.h5file.get_node(name)

        # Read all the rows:
        row = vlarray.read()
        if common.verbose:
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("First row in vlarray ==>", row[0])

        self.assertEqual(vlarray.nrows, 3)
        self.assertEqual(
            row[0],
            [[1, 2, 3], "aaa", "aaa\xef\xbf\xbd\xef\xbf\xbd\xef\xbf\xbd"],
        )
        list1 = list(row[1])
        obj = list1.pop()
        self.assertEqual(list1, [3, 4])
        self.assertEqual(obj.c, C().c)
        self.assertEqual(row[2], 42)
        self.assertEqual(len(row[0]), 3)
        self.assertEqual(len(row[1]), 3)
        self.assertRaises(TypeError, len, row[2])

    def test06b_Object(self):
        """Checking updating vlarray with object atoms."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test06b_Object..." % self.__class__.__name__)

        vlarray = self.h5file.create_vlarray(
            "/", "Object", atom=tb.ObjectAtom()
        )
        # When updating an object, this seems to change the number
        # of bytes that pickle.dumps generates
        # vlarray.append(
        #     ([1,2,3], "aaa", "aaa\xef\xbf\xbd\xef\xbf\xbd\xef\xbf\xbd"))
        vlarray.append(([1, 2, 3], "aaa", "\xef\xbf\xbd\xef\xbf\xbd4"))
        # vlarray.append([3,4, C()])
        vlarray.append([3, 4, [24]])

        # Modify the rows
        # vlarray[0] = ([1,2,4], "aa4", "aaa\xef\xbf\xbd\xef\xbf\xbd4")
        vlarray[0] = ([1, 2, 4], "aa4", "\xef\xbf\xbd\xef\xbf\xbd5")
        # vlarray[1] = (3,4, C())
        vlarray[1] = [4, 4, [24]]

        if self.reopen:
            name = vlarray._v_pathname
            self._reopen()
            vlarray = self.h5file.get_node(name)

        # Read all the rows:
        row = vlarray.read()
        if common.verbose:
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("First row in vlarray ==>", row[0])

        self.assertEqual(vlarray.nrows, 2)
        self.assertEqual(
            row[0], ([1, 2, 4], "aa4", "\xef\xbf\xbd\xef\xbf\xbd5")
        )
        list1 = list(row[1])
        obj = list1.pop()
        self.assertEqual(list1, [4, 4])

        # self.assertEqual(obj.c, C().c)
        self.assertEqual(obj, [24])
        self.assertEqual(len(row[0]), 3)
        self.assertEqual(len(row[1]), 3)

    def test06c_Object(self):
        """Checking vlarray with object atoms (numpy arrays as values)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test06c_Object..." % self.__class__.__name__)

        vlarray = self.h5file.create_vlarray(
            "/", "Object", atom=tb.ObjectAtom()
        )
        vlarray.append(np.array([[1, 2], [0, 4]], "i4"))
        vlarray.append(np.array([0, 1, 2, 3], "i8"))
        vlarray.append(np.array(42, "i1"))

        if self.reopen:
            name = vlarray._v_pathname
            self._reopen()
            vlarray = self.h5file.get_node(name)

        # Read all the rows:
        row = vlarray.read()
        if common.verbose:
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("First row in vlarray ==>", row[0])

        self.assertEqual(vlarray.nrows, 3)
        self.assertTrue(
            common.allequal(row[0], np.array([[1, 2], [0, 4]], "i4"))
        )
        self.assertTrue(common.allequal(row[1], np.array([0, 1, 2, 3], "i8")))
        self.assertTrue(common.allequal(row[2], np.array(42, "i1")))

    def test06d_Object(self):
        """Checking updating vlarray with object atoms (numpy arrays)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test06d_Object..." % self.__class__.__name__)

        vlarray = self.h5file.create_vlarray(
            "/", "Object", atom=tb.ObjectAtom()
        )
        vlarray.append(np.array([[1, 2], [0, 4]], "i4"))
        vlarray.append(np.array([0, 1, 2, 3], "i8"))
        vlarray.append(np.array(42, "i1"))

        # Modify the rows.  Since PyTables 2.2.1 we use a binary
        # pickle for arrays and ObjectAtoms, so the next should take
        # the same space than the above.
        vlarray[0] = np.array([[1, 0], [0, 4]], "i4")
        vlarray[1] = np.array([0, 1, 0, 3], "i8")
        vlarray[2] = np.array(22, "i1")

        if self.reopen:
            name = vlarray._v_pathname
            self._reopen()
            vlarray = self.h5file.get_node(name)

        # Read all the rows:
        row = vlarray.read()
        if common.verbose:
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("First row in vlarray ==>", row[0])

        self.assertEqual(vlarray.nrows, 3)
        self.assertTrue(
            common.allequal(row[0], np.array([[1, 0], [0, 4]], "i4"))
        )
        self.assertTrue(common.allequal(row[1], np.array([0, 1, 0, 3], "i8")))
        self.assertTrue(common.allequal(row[2], np.array(22, "i1")))

    def test07_VLUnicodeAtom(self):
        """Checking vlarray with variable length Unicode strings."""

        # Skip the test if the default encoding has been mangled.
        if sys.getdefaultencoding() != "ascii":
            return

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test07_VLUnicodeAtom..." % self.__class__.__name__
            )

        vlarray = self.h5file.create_vlarray(
            "/", "VLUnicodeAtom", atom=tb.VLUnicodeAtom()
        )
        vlarray.append("asd")
        vlarray.append("asd\u0140")
        vlarray.append("aaana")
        vlarray.append("")
        # Check for ticket #62.
        self.assertRaises(TypeError, vlarray.append, ["foo", "bar"])
        # `VLUnicodeAtom` makes no encoding assumptions.
        self.assertRaises(UnicodeDecodeError, vlarray.append, "asd\xe4")

        if self.reopen:
            name = vlarray._v_pathname
            self._reopen()
            vlarray = self.h5file.get_node(name)

        # Read all the rows:
        row = vlarray.read()
        if common.verbose:
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("First row in vlarray ==>", row[0])

        self.assertEqual(vlarray.nrows, 4)
        self.assertEqual(row[0], "asd")
        self.assertEqual(row[1], "asd\u0140")
        self.assertEqual(row[2], "aaana")
        self.assertEqual(row[3], "")
        self.assertEqual(len(row[0]), 3)
        self.assertEqual(len(row[1]), 4)
        self.assertEqual(len(row[2]), 5)
        self.assertEqual(len(row[3]), 0)

    def test07b_VLUnicodeAtom(self):
        """Checking updating vlarray with variable length Unicode strings."""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test07b_VLUnicodeAtom..." % self.__class__.__name__
            )

        vlarray = self.h5file.create_vlarray(
            "/", "VLUnicodeAtom", atom=tb.VLUnicodeAtom()
        )
        vlarray.append("asd")
        vlarray.append("aaan\xe4")

        # Modify values
        vlarray[0] = "as\xe4"
        vlarray[1] = "aaan\u0140"
        self.assertRaises(ValueError, vlarray.__setitem__, 1, "shrt")
        self.assertRaises(ValueError, vlarray.__setitem__, 1, "toolong")

        if self.reopen:
            name = vlarray._v_pathname
            self._reopen()
            vlarray = self.h5file.get_node(name)

        # Read all the rows:
        row = vlarray.read()
        if common.verbose:
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("First row in vlarray ==>", repr(row[0]))
            print("Second row in vlarray ==>", repr(row[1]))

        self.assertEqual(vlarray.nrows, 2)
        self.assertEqual(row[0], "as\xe4")
        self.assertEqual(row[1], "aaan\u0140")
        self.assertEqual(len(row[0]), 3)
        self.assertEqual(len(row[1]), 5)


class TypesReopenTestCase(TypesTestCase):
    title = "Reopen"
    reopen = True


class TypesNoReopenTestCase(TypesTestCase):
    title = "No reopen"
    reopen = False


class MDTypesTestCase(common.TempFileMixin, common.PyTablesTestCase):
    open_mode = "w"
    compress = 0
    complib = "zlib"  # Default compression library

    def setUp(self):
        super().setUp()
        self.rootgroup = self.h5file.root

    def test01_StringAtom(self):
        """Checking vlarray with MD NumPy string atoms."""

        root = self.rootgroup
        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01_StringAtom..." % self.__class__.__name__)

        # Create an string atom
        vlarray = self.h5file.create_vlarray(
            root,
            "stringAtom",
            tb.StringAtom(itemsize=3, shape=(2,)),
            "Ragged array of strings",
        )
        vlarray.append([["123", "45"], ["45", "123"]])
        vlarray.append([["s", "abc"], ["abc", "f"], ["s", "ab"], ["ab", "f"]])

        # Read all the rows:
        row = vlarray.read()
        if common.verbose:
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, 2)
        np.testing.assert_array_equal(
            row[0], np.array([["123", "45"], ["45", "123"]], "S")
        )
        np.testing.assert_array_equal(
            row[1],
            np.array(
                [["s", "abc"], ["abc", "f"], ["s", "ab"], ["ab", "f"]], "S"
            ),
        )
        self.assertEqual(len(row[0]), 2)
        self.assertEqual(len(row[1]), 4)

    def test01b_StringAtom(self):
        """Checking vlarray with MD NumPy string atoms ('python' flavor)"""

        root = self.rootgroup
        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01b_StringAtom..." % self.__class__.__name__)

        # Create an string atom
        vlarray = self.h5file.create_vlarray(
            root,
            "stringAtom",
            tb.StringAtom(itemsize=3, shape=(2,)),
            "Ragged array of strings",
        )
        vlarray.flavor = "python"
        vlarray.append([["123", "45"], ["45", "123"]])
        vlarray.append([["s", "abc"], ["abc", "f"], ["s", "ab"], ["ab", "f"]])

        # Read all the rows:
        row = vlarray.read()
        if common.verbose:
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, 2)
        self.assertEqual(row[0], [[b"123", b"45"], [b"45", b"123"]])
        self.assertEqual(
            row[1],
            [[b"s", b"abc"], [b"abc", b"f"], [b"s", b"ab"], [b"ab", b"f"]],
        )
        self.assertEqual(len(row[0]), 2)
        self.assertEqual(len(row[1]), 4)

    def test01c_StringAtom(self):
        """Checking vlarray with MD NumPy string atoms (with offset)"""

        root = self.rootgroup
        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01c_StringAtom..." % self.__class__.__name__)

        # Create an string atom
        vlarray = self.h5file.create_vlarray(
            root,
            "stringAtom",
            tb.StringAtom(itemsize=3, shape=(2,)),
            "Ragged array of strings",
        )
        vlarray.flavor = "python"
        a = np.array([["a", "b"], ["123", "45"], ["45", "123"]], dtype="S3")
        vlarray.append(a[1:])
        a = np.array(
            [
                ["s", "a"],
                ["ab", "f"],
                ["s", "abc"],
                ["abc", "f"],
                ["s", "ab"],
                ["ab", "f"],
            ]
        )
        vlarray.append(a[2:])

        # Read all the rows:
        row = vlarray.read()
        if common.verbose:
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, 2)
        self.assertEqual(row[0], [[b"123", b"45"], [b"45", b"123"]])
        self.assertEqual(
            row[1],
            [[b"s", b"abc"], [b"abc", b"f"], [b"s", b"ab"], [b"ab", b"f"]],
        )
        self.assertEqual(len(row[0]), 2)
        self.assertEqual(len(row[1]), 4)

    def test01d_StringAtom(self):
        """Checking vlarray with MD NumPy string atoms (with stride)"""

        root = self.rootgroup
        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01d_StringAtom..." % self.__class__.__name__)

        # Create an string atom
        vlarray = self.h5file.create_vlarray(
            root,
            "stringAtom",
            tb.StringAtom(itemsize=3, shape=(2,)),
            "Ragged array of strings",
        )
        vlarray.flavor = "python"
        a = np.array([["a", "b"], ["123", "45"], ["45", "123"]], dtype="S3")
        vlarray.append(a[1::2])
        a = np.array(
            [
                ["s", "a"],
                ["ab", "f"],
                ["s", "abc"],
                ["abc", "f"],
                ["s", "ab"],
                ["ab", "f"],
            ]
        )
        vlarray.append(a[::3])

        # Read all the rows:
        row = vlarray.read()
        if common.verbose:
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, 2)
        self.assertEqual(row[0], [[b"123", b"45"]])
        self.assertEqual(row[1], [[b"s", b"a"], [b"abc", b"f"]])
        self.assertEqual(len(row[0]), 1)
        self.assertEqual(len(row[1]), 2)

    def test02_BoolAtom(self):
        """Checking vlarray with MD boolean atoms."""

        root = self.rootgroup
        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02_BoolAtom..." % self.__class__.__name__)

        # Create an string atom
        vlarray = self.h5file.create_vlarray(
            root,
            "BoolAtom",
            tb.BoolAtom(shape=(3,)),
            "Ragged array of Booleans",
        )
        vlarray.append([(1, 0, 3), (1, 1, 1), (0, 0, 0)])
        vlarray.append([(1, 0, 0)])

        # Read all the rows:
        row = vlarray.read()
        if common.verbose:
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, 2)
        self.assertTrue(
            common.allequal(
                row[0],
                np.array([[1, 0, 1], [1, 1, 1], [0, 0, 0]], dtype="bool"),
            )
        )
        self.assertTrue(
            common.allequal(row[1], np.array([[1, 0, 0]], dtype="bool"))
        )
        self.assertEqual(len(row[0]), 3)
        self.assertEqual(len(row[1]), 1)

    def test02b_BoolAtom(self):
        """Checking vlarray with MD boolean atoms (with offset)"""

        root = self.rootgroup
        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02b_BoolAtom..." % self.__class__.__name__)

        # Create an string atom
        vlarray = self.h5file.create_vlarray(
            root,
            "BoolAtom",
            tb.BoolAtom(shape=(3,)),
            "Ragged array of Booleans",
        )
        a = np.array(
            [(0, 0, 0), (1, 0, 3), (1, 1, 1), (0, 0, 0)], dtype="bool"
        )
        vlarray.append(a[1:])  # Create an offset
        a = np.array([(1, 1, 1), (1, 0, 0)], dtype="bool")
        vlarray.append(a[1:])  # Create an offset

        # Read all the rows:
        row = vlarray.read()
        if common.verbose:
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, 2)
        self.assertTrue(
            common.allequal(
                row[0],
                np.array([[1, 0, 1], [1, 1, 1], [0, 0, 0]], dtype="bool"),
            )
        )
        self.assertTrue(
            common.allequal(row[1], np.array([[1, 0, 0]], dtype="bool"))
        )
        self.assertEqual(len(row[0]), 3)
        self.assertEqual(len(row[1]), 1)

    def test02c_BoolAtom(self):
        """Checking vlarray with MD boolean atoms (with strides)"""

        root = self.rootgroup
        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02c_BoolAtom..." % self.__class__.__name__)

        # Create an string atom
        vlarray = self.h5file.create_vlarray(
            root,
            "BoolAtom",
            tb.BoolAtom(shape=(3,)),
            "Ragged array of Booleans",
        )
        a = np.array(
            [(0, 0, 0), (1, 0, 3), (1, 1, 1), (0, 0, 0)], dtype="bool"
        )
        vlarray.append(a[1::2])  # Create an strided array
        a = np.array([(1, 1, 1), (1, 0, 0), (0, 0, 0)], dtype="bool")
        vlarray.append(a[::2])  # Create an strided array

        # Read all the rows:
        row = vlarray.read()
        if common.verbose:
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, 2)
        self.assertTrue(
            common.allequal(
                row[0], np.array([[1, 0, 1], [0, 0, 0]], dtype="bool")
            )
        )
        self.assertTrue(
            common.allequal(
                row[1], np.array([[1, 1, 1], [0, 0, 0]], dtype="bool")
            )
        )
        self.assertEqual(len(row[0]), 2)
        self.assertEqual(len(row[1]), 2)

    def test03_IntAtom(self):
        """Checking vlarray with MD integer atoms."""

        ttypes = [
            "int8",
            "uint8",
            "int16",
            "uint16",
            "int32",
            "uint32",
            "int64",
            # "UInt64",  # Unavailable in some platforms
        ]
        root = self.rootgroup
        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03_IntAtom..." % self.__class__.__name__)

        # Create an string atom
        for atype in ttypes:
            vlarray = self.h5file.create_vlarray(
                root, atype, atom=tb.Atom.from_sctype(atype, (2, 3))
            )
            vlarray.append([np.ones((2, 3), atype), np.zeros((2, 3), atype)])
            vlarray.append([np.ones((2, 3), atype) * 100])

            # Read all the rows:
            row = vlarray.read()
            if common.verbose:
                print("Testing type:", atype)
                print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
                print("Second row in vlarray ==>", repr(row[1]))

            self.assertEqual(vlarray.nrows, 2)
            self.assertTrue(
                common.allequal(
                    row[0],
                    np.array([np.ones((2, 3)), np.zeros((2, 3))], atype),
                )
            )
            self.assertTrue(
                common.allequal(
                    row[1], np.array([np.ones((2, 3)) * 100], atype)
                )
            )
            self.assertEqual(len(row[0]), 2)
            self.assertEqual(len(row[1]), 1)

    def test04_FloatAtom(self):
        """Checking vlarray with MD floating point atoms."""

        ttypes = [
            "float32",
            "float64",
            "complex64",
            "complex128",
        ]

        for name in ("float16", "float96", "float128"):
            atomname = name.capitalize() + "Atom"
            if hasattr(tb, atomname):
                ttypes.append(name)
        for itemsize in (192, 256):
            atomname = "Complex%dAtom" % itemsize
            if hasattr(tb, atomname):
                ttypes.append("complex%d" % (itemsize))

        root = self.rootgroup
        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test04_FloatAtom..." % self.__class__.__name__)

        # Create an string atom
        for atype in ttypes:
            vlarray = self.h5file.create_vlarray(
                root, atype, atom=tb.Atom.from_sctype(atype, (5, 2, 6))
            )
            vlarray.append(
                [np.ones((5, 2, 6), atype) * 1.3, np.zeros((5, 2, 6), atype)]
            )
            vlarray.append([np.ones((5, 2, 6), atype) * 2.0e4])

            # Read all the rows:
            row = vlarray.read()
            if common.verbose:
                print("Testing type:", atype)
                print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
                print("Second row in vlarray ==>", row[1])

            self.assertEqual(vlarray.nrows, 2)
            self.assertTrue(
                common.allequal(
                    row[0],
                    np.array(
                        [np.ones((5, 2, 6)) * 1.3, np.zeros((5, 2, 6))], atype
                    ),
                )
            )
            self.assertTrue(
                common.allequal(
                    row[1], np.array([np.ones((5, 2, 6)) * 2.0e4], atype)
                )
            )
            self.assertEqual(len(row[0]), 2)
            self.assertEqual(len(row[1]), 1)


class MDTypesNumPyTestCase(MDTypesTestCase):
    title = "MDTypes"


class AppendShapeTestCase(common.TempFileMixin, common.PyTablesTestCase):
    open_mode = "w"

    def setUp(self):
        super().setUp()
        self.rootgroup = self.h5file.root

    def test00_difinputs(self):
        """Checking vlarray.append() with different inputs."""

        root = self.rootgroup
        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test00_difinputs..." % self.__class__.__name__)

        # Create an string atom
        vlarray = self.h5file.create_vlarray(
            root, "vlarray", tb.Int32Atom(), "Ragged array of ints"
        )
        vlarray.flavor = "python"

        # Check different ways to input
        # All of the next should lead to the same rows
        vlarray.append((1, 2, 3))  # a tuple
        vlarray.append([1, 2, 3])  # a unique list
        vlarray.append(np.array([1, 2, 3], dtype="int32"))  # and array

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen()
            vlarray = self.h5file.root.vlarray

        # Read all the vlarray
        row = vlarray.read()
        if common.verbose:
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("First row in vlarray ==>", row[0])

        self.assertEqual(vlarray.nrows, 3)
        self.assertEqual(row[0], [1, 2, 3])
        self.assertEqual(row[1], [1, 2, 3])
        self.assertEqual(row[2], [1, 2, 3])

    def test01_toomanydims(self):
        """Checking vlarray.append() with too many dimensions."""

        root = self.rootgroup
        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01_toomanydims..." % self.__class__.__name__)

        # Create an string atom
        vlarray = self.h5file.create_vlarray(
            root,
            "vlarray",
            tb.StringAtom(itemsize=3),
            "Ragged array of strings",
        )
        # Adding an array with one dimensionality more than allowed
        with self.assertRaises(ValueError):
            vlarray.append([["123", "456", "3"]])

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen()
            vlarray = self.h5file.root.vlarray

        # Read all the rows (there should be none)
        row = vlarray.read()
        if common.verbose:
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)

        self.assertEqual(vlarray.nrows, 0)

    def test02_zerodims(self):
        """Checking vlarray.append() with a zero-dimensional array"""

        root = self.rootgroup
        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02_zerodims..." % self.__class__.__name__)

        # Create an string atom
        vlarray = self.h5file.create_vlarray(
            root, "vlarray", tb.Int32Atom(), "Ragged array of ints"
        )
        vlarray.append(np.zeros(dtype="int32", shape=(6, 0)))

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen()
            vlarray = self.h5file.root.vlarray

        # Read the only row in vlarray
        row = vlarray.read(0)[0]
        if common.verbose:
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("First row in vlarray ==>", repr(row))

        self.assertEqual(vlarray.nrows, 1)
        self.assertTrue(
            common.allequal(row, np.zeros(dtype="int32", shape=(0,)))
        )
        self.assertEqual(len(row), 0)

    def test03a_cast(self):
        """Checking vlarray.append() with a casted array (upgrading case)"""

        root = self.rootgroup
        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03a_cast..." % self.__class__.__name__)

        # Create an string atom
        vlarray = self.h5file.create_vlarray(
            root, "vlarray", tb.Int32Atom(), "Ragged array of ints"
        )
        # This type has to be upgraded
        vlarray.append(np.array([1, 2], dtype="int16"))

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen()
            vlarray = self.h5file.root.vlarray

        # Read the only row in vlarray
        row = vlarray.read(0)[0]
        if common.verbose:
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("First row in vlarray ==>", repr(row))

        self.assertEqual(vlarray.nrows, 1)
        self.assertTrue(common.allequal(row, np.array([1, 2], dtype="int32")))
        self.assertEqual(len(row), 2)

    def test03b_cast(self):
        """Checking vlarray.append() with a casted array (downgrading case)"""

        root = self.rootgroup
        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03b_cast..." % self.__class__.__name__)

        # Create an string atom
        vlarray = self.h5file.create_vlarray(
            root, "vlarray", tb.Int32Atom(), "Ragged array of ints"
        )
        # This type has to be downcasted
        vlarray.append(np.array([1, 2], dtype="float64"))

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen()
            vlarray = self.h5file.root.vlarray

        # Read the only row in vlarray
        row = vlarray.read(0)[0]
        if common.verbose:
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("First row in vlarray ==>", repr(row))

        self.assertEqual(vlarray.nrows, 1)
        self.assertTrue(common.allequal(row, np.array([1, 2], dtype="int32")))
        self.assertEqual(len(row), 2)


class OpenAppendShapeTestCase(AppendShapeTestCase):
    close = 0


class CloseAppendShapeTestCase(AppendShapeTestCase):
    close = 1


class FlavorTestCase(common.TempFileMixin, common.PyTablesTestCase):
    open_mode = "w"
    compress = 0
    complib = "zlib"  # Default compression library

    def setUp(self):
        super().setUp()
        self.rootgroup = self.h5file.root

    def test01a_EmptyVLArray(self):
        """Checking empty vlarrays with different flavors (closing the file)"""

        root = self.rootgroup
        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test01_EmptyVLArray..." % self.__class__.__name__
            )

        # Create an string atom
        vlarray = self.h5file.create_vlarray(
            root, "vlarray", tb.Atom.from_kind("int", itemsize=4)
        )
        vlarray.flavor = self.flavor
        self.h5file.close()
        self.h5file = tb.open_file(self.h5fname, "r")

        # Read all the rows (it should be empty):
        vlarray = self.h5file.root.vlarray
        row = vlarray.read()
        if common.verbose:
            print("Testing flavor:", self.flavor)
            print("Object read:", row, repr(row))
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)

        # Check that the object read is effectively empty
        self.assertEqual(vlarray.nrows, 0)
        self.assertEqual(row, [])

    def test01b_EmptyVLArray(self):
        """Checking empty vlarrays with different flavors (no closing file)"""

        root = self.rootgroup
        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test01_EmptyVLArray..." % self.__class__.__name__
            )

        # Create an string atom
        vlarray = self.h5file.create_vlarray(
            root, "vlarray", tb.Atom.from_kind("int", itemsize=4)
        )
        vlarray.flavor = self.flavor

        # Read all the rows (it should be empty):
        row = vlarray.read()
        if common.verbose:
            print("Testing flavor:", self.flavor)
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)

        # Check that the object read is effectively empty
        self.assertEqual(vlarray.nrows, 0)
        self.assertEqual(row, [])

    def test02_BooleanAtom(self):
        """Checking vlarray with different flavors (boolean versions)"""

        root = self.rootgroup
        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02_BoolAtom..." % self.__class__.__name__)

        # Create an string atom
        vlarray = self.h5file.create_vlarray(root, "Bool", tb.BoolAtom())
        vlarray.flavor = self.flavor
        vlarray.append([1, 2, 3])
        vlarray.append(())  # Empty row
        vlarray.append([100, 0])

        # Read all the rows:
        row = vlarray.read()
        if common.verbose:
            print("Testing flavor:", self.flavor)
            print("Object read:", row)
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("First row in vlarray ==>", row[0])

        self.assertEqual(vlarray.nrows, 3)
        self.assertEqual(len(row[0]), 3)
        self.assertEqual(len(row[1]), 0)
        self.assertEqual(len(row[2]), 2)
        if self.flavor == "python":
            arr1 = [1, 1, 1]
            arr2 = []
            arr3 = [1, 0]
        elif self.flavor == "numpy":
            arr1 = np.array([1, 1, 1], dtype="bool")
            arr2 = np.array([], dtype="bool")
            arr3 = np.array([1, 0], dtype="bool")

        if self.flavor == "numpy":
            self.assertTrue(common.allequal(row[0], arr1, self.flavor))
            self.assertTrue(common.allequal(row[1], arr2, self.flavor))
            self.assertTrue(common.allequal(row[1], arr2, self.flavor))
        else:
            # 'python' flavor
            self.assertEqual(row[0], arr1)
            self.assertEqual(row[1], arr2)
            self.assertEqual(row[2], arr3)

    def test03_IntAtom(self):
        """Checking vlarray with different flavors (integer versions)"""

        ttypes = [
            "int8",
            "uint8",
            "int16",
            "uint16",
            "int32",
            "uint32",
            "int64",
            # Not checked because some platforms does not support it
            # "UInt64",
        ]

        root = self.rootgroup
        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03_IntAtom..." % self.__class__.__name__)

        # Create an string atom
        for atype in ttypes:
            vlarray = self.h5file.create_vlarray(
                root, atype, tb.Atom.from_sctype(atype)
            )
            vlarray.flavor = self.flavor
            vlarray.append([1, 2, 3])
            vlarray.append(())
            vlarray.append([100, 0])

            # Read all the rows:
            row = vlarray.read()
            if common.verbose:
                print("Testing flavor:", self.flavor)
                print("Object read:", row)
                print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
                print("First row in vlarray ==>", row[0])

            self.assertEqual(vlarray.nrows, 3)
            self.assertEqual(len(row[0]), 3)
            self.assertEqual(len(row[1]), 0)
            self.assertEqual(len(row[2]), 2)
            if self.flavor == "python":
                arr1 = [1, 2, 3]
                arr2 = []
                arr3 = [100, 0]
            elif self.flavor == "numpy":
                arr1 = np.array([1, 2, 3], dtype=atype)
                arr2 = np.array([], dtype=atype)
                arr3 = np.array([100, 0], dtype=atype)

            if self.flavor == "numpy":
                self.assertTrue(common.allequal(row[0], arr1, self.flavor))
                self.assertTrue(common.allequal(row[1], arr2, self.flavor))
                self.assertTrue(common.allequal(row[2], arr3, self.flavor))
            else:
                # "python" flavor
                self.assertEqual(row[0], arr1)
                self.assertEqual(row[1], arr2)
                self.assertEqual(row[2], arr3)

    def test03b_IntAtom(self):
        """Checking vlarray flavors (integer versions and closed file)"""

        ttypes = [
            "int8",
            "uint8",
            "int16",
            "uint16",
            "int32",
            "uint32",
            "int64",
            # Not checked because some platforms does not support it
            # "UInt64",
        ]

        root = self.rootgroup
        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03_IntAtom..." % self.__class__.__name__)

        # Create an string atom
        for atype in ttypes:
            vlarray = self.h5file.create_vlarray(
                root, atype, tb.Atom.from_sctype(atype)
            )
            vlarray.flavor = self.flavor
            vlarray.append([1, 2, 3])
            vlarray.append(())
            vlarray.append([100, 0])
            self._reopen(mode="a")  # open in "a"ppend mode
            root = self.h5file.root  # Very important!
            vlarray = self.h5file.get_node(root, str(atype))

            # Read all the rows:
            row = vlarray.read()
            if common.verbose:
                print("Testing flavor:", self.flavor)
                print("Object read:", row)
                print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
                print("First row in vlarray ==>", row[0])

            self.assertEqual(vlarray.nrows, 3)
            self.assertEqual(len(row[0]), 3)
            self.assertEqual(len(row[1]), 0)
            self.assertEqual(len(row[2]), 2)
            if self.flavor == "python":
                arr1 = [1, 2, 3]
                arr2 = []
                arr3 = [100, 0]
            elif self.flavor == "numpy":
                arr1 = np.array([1, 2, 3], dtype=atype)
                arr2 = np.array([], dtype=atype)
                arr3 = np.array([100, 0], dtype=atype)

            if self.flavor == "numpy":
                self.assertTrue(common.allequal(row[0], arr1, self.flavor))
                self.assertTrue(common.allequal(row[1], arr2, self.flavor))
                self.assertTrue(common.allequal(row[2], arr3, self.flavor))
            else:
                # Tuple or List flavors
                self.assertEqual(row[0], arr1)
                self.assertEqual(row[1], arr2)
                self.assertEqual(row[2], arr3)

    def test04_FloatAtom(self):
        """Checking vlarray with different flavors (floating point versions)"""

        ttypes = [
            "float32",
            "float64",
            "complex64",
            "complex128",
        ]

        for name in ("float16", "float96", "float128"):
            atomname = name.capitalize() + "Atom"
            if hasattr(tb, atomname):
                ttypes.append(name)

        for itemsize in (192, 256):
            atomname = "Complex%dAtom" % itemsize
            if hasattr(tb, atomname):
                ttypes.append("complex%d" % (itemsize))

        root = self.rootgroup
        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test04_FloatAtom..." % self.__class__.__name__)

        # Create an string atom
        for atype in ttypes:
            vlarray = self.h5file.create_vlarray(
                root, atype, tb.Atom.from_sctype(atype)
            )
            vlarray.flavor = self.flavor
            vlarray.append([1.3, 2.2, 3.3])
            vlarray.append(())
            vlarray.append([5.96, 0.597])

            # Read all the rows:
            row = vlarray.read()
            if common.verbose:
                print("Testing flavor:", self.flavor)
                print("Object read:", row)
                print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
                print("First row in vlarray ==>", row[0])

            self.assertEqual(vlarray.nrows, 3)
            self.assertEqual(len(row[0]), 3)
            self.assertEqual(len(row[1]), 0)
            self.assertEqual(len(row[2]), 2)
            if self.flavor == "python":
                arr1 = list(np.array([1.3, 2.2, 3.3], atype))
                arr2 = list(np.array([], atype))
                arr3 = list(np.array([5.96, 0.597], atype))
            elif self.flavor == "numpy":
                arr1 = np.array([1.3, 2.2, 3.3], dtype=atype)
                arr2 = np.array([], dtype=atype)
                arr3 = np.array([5.96, 0.597], dtype=atype)

            if self.flavor == "numpy":
                self.assertTrue(common.allequal(row[0], arr1, self.flavor))
                self.assertTrue(common.allequal(row[1], arr2, self.flavor))
                self.assertTrue(common.allequal(row[2], arr3, self.flavor))
            else:
                # Tuple or List flavors
                self.assertEqual(row[0], arr1)
                self.assertEqual(row[1], arr2)
                self.assertEqual(row[2], arr3)


class NumPyFlavorTestCase(FlavorTestCase):
    flavor = "numpy"


class PythonFlavorTestCase(FlavorTestCase):
    flavor = "python"


class ReadRangeTestCase(common.TempFileMixin, common.PyTablesTestCase):
    nrows = 100
    mode = "w"
    compress = 0
    complib = "zlib"  # Default compression library

    def setUp(self):
        super().setUp()
        self.rootgroup = self.h5file.root
        self.populateFile()
        self._reopen()

    def populateFile(self):
        group = self.rootgroup
        filters = tb.Filters(complevel=self.compress, complib=self.complib)
        vlarray = self.h5file.create_vlarray(
            group,
            "vlarray",
            tb.Int32Atom(),
            "ragged array if ints",
            filters=filters,
            expectedrows=1000,
        )

        # Fill it with 100 rows with variable length
        for i in range(self.nrows):
            vlarray.append(list(range(i)))

    def test01_start(self):
        """Checking reads with only a start value"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01_start..." % self.__class__.__name__)

        vlarray = self.h5file.root.vlarray

        # Read some rows:
        row = []
        row.append(vlarray.read(0)[0])
        row.append(vlarray.read(10)[0])
        row.append(vlarray.read(99)[0])
        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, self.nrows)
        self.assertEqual(len(row[0]), 0)
        self.assertEqual(len(row[1]), 10)
        self.assertEqual(len(row[2]), 99)
        self.assertTrue(common.allequal(row[0], np.arange(0, dtype="int32")))
        self.assertTrue(common.allequal(row[1], np.arange(10, dtype="int32")))
        self.assertTrue(common.allequal(row[2], np.arange(99, dtype="int32")))

    def test01b_start(self):
        """Checking reads with only a start value in a slice"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01b_start..." % self.__class__.__name__)

        vlarray = self.h5file.root.vlarray

        # Read some rows:
        row = []
        row.append(vlarray[0])
        row.append(vlarray[10])
        row.append(vlarray[99])
        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, self.nrows)
        self.assertEqual(len(row[0]), 0)
        self.assertEqual(len(row[1]), 10)
        self.assertEqual(len(row[2]), 99)
        self.assertTrue(common.allequal(row[0], np.arange(0, dtype="int32")))
        self.assertTrue(common.allequal(row[1], np.arange(10, dtype="int32")))
        self.assertTrue(common.allequal(row[2], np.arange(99, dtype="int32")))

    def test01np_start(self):
        """Checking reads with only a start value in a slice (numpy indexes)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01np_start..." % self.__class__.__name__)

        vlarray = self.h5file.root.vlarray

        # Read some rows:
        row = []
        row.append(vlarray[np.int8(0)])
        row.append(vlarray[np.int32(10)])
        row.append(vlarray[np.int64(99)])
        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, self.nrows)
        self.assertEqual(len(row[0]), 0)
        self.assertEqual(len(row[1]), 10)
        self.assertEqual(len(row[2]), 99)
        self.assertTrue(common.allequal(row[0], np.arange(0, dtype="int32")))
        self.assertTrue(common.allequal(row[1], np.arange(10, dtype="int32")))
        self.assertTrue(common.allequal(row[2], np.arange(99, dtype="int32")))

    def test02_stop(self):
        """Checking reads with only a stop value"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02_stop..." % self.__class__.__name__)

        vlarray = self.h5file.root.vlarray
        # Choose a small value for buffer size
        vlarray._nrowsinbuf = 3

        # Read some rows:
        row = []
        row.append(vlarray.read(stop=1))
        row.append(vlarray.read(stop=10))
        row.append(vlarray.read(stop=99))
        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("First row in vlarray ==>", row[0])
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, self.nrows)
        self.assertEqual(len(row[0]), 1)
        self.assertEqual(len(row[1]), 10)
        self.assertEqual(len(row[2]), 99)
        self.assertTrue(
            common.allequal(row[0][0], np.arange(0, dtype="int32"))
        )
        for x in range(10):
            self.assertTrue(
                common.allequal(row[1][x], np.arange(x, dtype="int32"))
            )
        for x in range(99):
            self.assertTrue(
                common.allequal(row[2][x], np.arange(x, dtype="int32"))
            )

    def test02b_stop(self):
        """Checking reads with only a stop value in a slice"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02b_stop..." % self.__class__.__name__)

        vlarray = self.h5file.root.vlarray

        # Choose a small value for buffer size
        vlarray._nrowsinbuf = 3

        # Read some rows:
        row = []
        row.append(vlarray[:1])
        row.append(vlarray[:10])
        row.append(vlarray[:99])
        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, self.nrows)
        self.assertEqual(len(row[0]), 1)
        self.assertEqual(len(row[1]), 10)
        self.assertEqual(len(row[2]), 99)
        for x in range(1):
            self.assertTrue(
                common.allequal(row[0][x], np.arange(0, dtype="int32"))
            )
        for x in range(10):
            self.assertTrue(
                common.allequal(row[1][x], np.arange(x, dtype="int32"))
            )
        for x in range(99):
            self.assertTrue(
                common.allequal(row[2][x], np.arange(x, dtype="int32"))
            )

    def test03_startstop(self):
        """Checking reads with a start and stop values"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03_startstop..." % self.__class__.__name__)

        vlarray = self.h5file.root.vlarray

        # Choose a small value for buffer size
        vlarray._nrowsinbuf = 3

        # Read some rows:
        row = []
        row.append(vlarray.read(0, 10))
        row.append(vlarray.read(5, 15))
        row.append(vlarray.read(0, 100))  # read all the array
        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, self.nrows)
        self.assertEqual(len(row[0]), 10)
        self.assertEqual(len(row[1]), 10)
        self.assertEqual(len(row[2]), 100)
        for x in range(0, 10):
            self.assertTrue(
                common.allequal(row[0][x], np.arange(x, dtype="int32"))
            )
        for x in range(5, 15):
            self.assertTrue(
                common.allequal(row[1][x - 5], np.arange(x, dtype="int32"))
            )
        for x in range(0, 100):
            self.assertTrue(
                common.allequal(row[2][x], np.arange(x, dtype="int32"))
            )

    def test03b_startstop(self):
        """Checking reads with a start and stop values in slices"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03b_startstop..." % self.__class__.__name__)

        vlarray = self.h5file.root.vlarray

        # Choose a small value for buffer size
        vlarray._nrowsinbuf = 3

        # Read some rows:
        row = []
        row.append(vlarray[0:10])
        row.append(vlarray[5:15])
        row.append(vlarray[:])  # read all the array
        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, self.nrows)
        self.assertEqual(len(row[0]), 10)
        self.assertEqual(len(row[1]), 10)
        self.assertEqual(len(row[2]), 100)
        for x in range(0, 10):
            self.assertTrue(
                common.allequal(row[0][x], np.arange(x, dtype="int32"))
            )
        for x in range(5, 15):
            self.assertTrue(
                common.allequal(row[1][x - 5], np.arange(x, dtype="int32"))
            )
        for x in range(0, 100):
            self.assertTrue(
                common.allequal(row[2][x], np.arange(x, dtype="int32"))
            )

    def test04_startstopstep(self):
        """Checking reads with a start, stop & step values"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test04_startstopstep..." % self.__class__.__name__
            )

        vlarray = self.h5file.root.vlarray

        # Choose a small value for buffer size
        vlarray._nrowsinbuf = 3

        # Read some rows:
        row = []
        row.append(vlarray.read(0, 10, 2))
        row.append(vlarray.read(5, 15, 3))
        row.append(vlarray.read(0, 100, 20))
        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, self.nrows)
        self.assertEqual(len(row[0]), 5)
        self.assertEqual(len(row[1]), 4)
        self.assertEqual(len(row[2]), 5)
        for x in range(0, 10, 2):
            self.assertTrue(
                common.allequal(row[0][x // 2], np.arange(x, dtype="int32"))
            )
        for x in range(5, 15, 3):
            self.assertTrue(
                common.allequal(
                    row[1][(x - 5) // 3], np.arange(x, dtype="int32")
                )
            )
        for x in range(0, 100, 20):
            self.assertTrue(
                common.allequal(row[2][x // 20], np.arange(x, dtype="int32"))
            )

    def test04np_startstopstep(self):
        """Checking reads with a start, stop & step values (numpy indices)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test04np_startstopstep..."
                % self.__class__.__name__
            )

        vlarray = self.h5file.root.vlarray

        # Choose a small value for buffer size
        vlarray._nrowsinbuf = 3

        # Read some rows:
        row = []
        row.append(vlarray.read(np.int8(0), np.int8(10), np.int8(2)))
        row.append(vlarray.read(np.int8(5), np.int8(15), np.int8(3)))
        row.append(vlarray.read(np.int8(0), np.int8(100), np.int8(20)))
        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, self.nrows)
        self.assertEqual(len(row[0]), 5)
        self.assertEqual(len(row[1]), 4)
        self.assertEqual(len(row[2]), 5)
        for x in range(0, 10, 2):
            self.assertTrue(
                common.allequal(row[0][x // 2], np.arange(x, dtype="int32"))
            )
        for x in range(5, 15, 3):
            self.assertTrue(
                common.allequal(
                    row[1][(x - 5) // 3], np.arange(x, dtype="int32")
                )
            )
        for x in range(0, 100, 20):
            self.assertTrue(
                common.allequal(row[2][x // 20], np.arange(x, dtype="int32"))
            )

    def test04b_slices(self):
        """Checking reads with start, stop & step values in slices"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test04b_slices..." % self.__class__.__name__)

        vlarray = self.h5file.root.vlarray

        # Choose a small value for buffer size
        vlarray._nrowsinbuf = 3

        # Read some rows:
        row = []
        row.append(vlarray[0:10:2])
        row.append(vlarray[5:15:3])
        row.append(vlarray[0:100:20])
        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, self.nrows)
        self.assertEqual(len(row[0]), 5)
        self.assertEqual(len(row[1]), 4)
        self.assertEqual(len(row[2]), 5)
        for x in range(0, 10, 2):
            self.assertTrue(
                common.allequal(row[0][x // 2], np.arange(x, dtype="int32"))
            )
        for x in range(5, 15, 3):
            self.assertTrue(
                common.allequal(
                    row[1][(x - 5) // 3], np.arange(x, dtype="int32")
                )
            )
        for x in range(0, 100, 20):
            self.assertTrue(
                common.allequal(row[2][x // 20], np.arange(x, dtype="int32"))
            )

    def test04bnp_slices(self):
        """Checking reads with start, stop & step values in slices.

        (numpy indices)

        """

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test04bnp_slices..." % self.__class__.__name__)

        vlarray = self.h5file.root.vlarray

        # Choose a small value for buffer size
        vlarray._nrowsinbuf = 3

        # Read some rows:
        row = []
        row.append(vlarray[np.int16(0) : np.int16(10) : np.int32(2)])
        row.append(vlarray[np.int16(5) : np.int16(15) : np.int64(3)])
        row.append(vlarray[np.uint16(0) : np.int32(100) : np.int8(20)])
        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, self.nrows)
        self.assertEqual(len(row[0]), 5)
        self.assertEqual(len(row[1]), 4)
        self.assertEqual(len(row[2]), 5)
        for x in range(0, 10, 2):
            self.assertTrue(
                common.allequal(row[0][x // 2], np.arange(x, dtype="int32"))
            )
        for x in range(5, 15, 3):
            self.assertTrue(
                common.allequal(
                    row[1][(x - 5) // 3], np.arange(x, dtype="int32")
                )
            )
        for x in range(0, 100, 20):
            self.assertTrue(
                common.allequal(row[2][x // 20], np.arange(x, dtype="int32"))
            )

    def test05_out_of_range(self):
        """Checking out of range reads"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test05_out_of_range..." % self.__class__.__name__
            )

        vlarray = self.h5file.root.vlarray

        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)

        with self.assertRaises(IndexError):
            row = vlarray.read(1000)[0]
            print("row-->", row)


class GetItemRangeTestCase(common.TempFileMixin, common.PyTablesTestCase):
    nrows = 100
    open_mode = "w"
    compress = 0
    complib = "zlib"  # Default compression library

    def setUp(self):
        super().setUp()

        self.rootgroup = self.h5file.root
        self.populateFile()
        self._reopen()

    def populateFile(self):
        group = self.rootgroup
        filters = tb.Filters(complevel=self.compress, complib=self.complib)
        vlarray = self.h5file.create_vlarray(
            group,
            "vlarray",
            tb.Int32Atom(),
            "ragged array if ints",
            filters=filters,
            expectedrows=1000,
        )

        # Fill it with 100 rows with variable length
        for i in range(self.nrows):
            vlarray.append(list(range(i)))

    def test01_start(self):
        """Checking reads with only a start value"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01_start..." % self.__class__.__name__)

        vlarray = self.h5file.root.vlarray

        # Read some rows:
        row = []
        row.append(vlarray[0])

        # rank-0 array should work as a regular index (see #303)
        row.append(vlarray[np.array(10)])
        row.append(vlarray[99])
        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, self.nrows)
        self.assertEqual(len(row[0]), 0)
        self.assertEqual(len(row[1]), 10)
        self.assertEqual(len(row[2]), 99)
        self.assertTrue(common.allequal(row[0], np.arange(0, dtype="int32")))
        self.assertTrue(common.allequal(row[1], np.arange(10, dtype="int32")))
        self.assertTrue(common.allequal(row[2], np.arange(99, dtype="int32")))

    def test01b_start(self):
        """Checking reads with only a start value in a slice"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01b_start..." % self.__class__.__name__)

        vlarray = self.h5file.root.vlarray

        # Read some rows:
        row = []
        row.append(vlarray[0])
        row.append(vlarray[10])
        row.append(vlarray[99])
        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, self.nrows)
        self.assertEqual(len(row[0]), 0)
        self.assertEqual(len(row[1]), 10)
        self.assertEqual(len(row[2]), 99)
        self.assertTrue(common.allequal(row[0], np.arange(0, dtype="int32")))
        self.assertTrue(common.allequal(row[1], np.arange(10, dtype="int32")))
        self.assertTrue(common.allequal(row[2], np.arange(99, dtype="int32")))

    def test02_stop(self):
        """Checking reads with only a stop value"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02_stop..." % self.__class__.__name__)

        vlarray = self.h5file.root.vlarray

        # Choose a small value for buffer size
        vlarray._nrowsinbuf = 3

        # Read some rows:
        row = []
        row.append(vlarray[:1])
        row.append(vlarray[:10])
        row.append(vlarray[:99])
        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("First row in vlarray ==>", row[0])
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, self.nrows)
        self.assertEqual(len(row[0]), 1)
        self.assertEqual(len(row[1]), 10)
        self.assertEqual(len(row[2]), 99)
        self.assertTrue(
            common.allequal(row[0][0], np.arange(0, dtype="int32"))
        )
        for x in range(10):
            self.assertTrue(
                common.allequal(row[1][x], np.arange(x, dtype="int32"))
            )
        for x in range(99):
            self.assertTrue(
                common.allequal(row[2][x], np.arange(x, dtype="int32"))
            )

    def test02b_stop(self):
        """Checking reads with only a stop value in a slice"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02b_stop..." % self.__class__.__name__)

        vlarray = self.h5file.root.vlarray

        # Choose a small value for buffer size
        vlarray._nrowsinbuf = 3

        # Read some rows:
        row = []
        row.append(vlarray[:1])
        row.append(vlarray[:10])
        row.append(vlarray[:99])
        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, self.nrows)
        self.assertEqual(len(row[0]), 1)
        self.assertEqual(len(row[1]), 10)
        self.assertEqual(len(row[2]), 99)
        for x in range(1):
            self.assertTrue(
                common.allequal(row[0][x], np.arange(0, dtype="int32"))
            )
        for x in range(10):
            self.assertTrue(
                common.allequal(row[1][x], np.arange(x, dtype="int32"))
            )
        for x in range(99):
            self.assertTrue(
                common.allequal(row[2][x], np.arange(x, dtype="int32"))
            )

    def test03_startstop(self):
        """Checking reads with a start and stop values"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03_startstop..." % self.__class__.__name__)

        vlarray = self.h5file.root.vlarray

        # Choose a small value for buffer size
        vlarray._nrowsinbuf = 3

        # Read some rows:
        row = []
        row.append(vlarray[0:10])
        row.append(vlarray[5:15])
        row.append(vlarray[0:100])  # read all the array
        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, self.nrows)
        self.assertEqual(len(row[0]), 10)
        self.assertEqual(len(row[1]), 10)
        self.assertEqual(len(row[2]), 100)
        for x in range(0, 10):
            self.assertTrue(
                common.allequal(row[0][x], np.arange(x, dtype="int32"))
            )
        for x in range(5, 15):
            self.assertTrue(
                common.allequal(row[1][x - 5], np.arange(x, dtype="int32"))
            )
        for x in range(0, 100):
            self.assertTrue(
                common.allequal(row[2][x], np.arange(x, dtype="int32"))
            )

    def test03b_startstop(self):
        """Checking reads with a start and stop values in slices"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03b_startstop..." % self.__class__.__name__)

        vlarray = self.h5file.root.vlarray

        # Choose a small value for buffer size
        vlarray._nrowsinbuf = 3

        # Read some rows:
        row = []
        row.append(vlarray[0:10])
        row.append(vlarray[5:15])
        row.append(vlarray[:])  # read all the array
        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, self.nrows)
        self.assertEqual(len(row[0]), 10)
        self.assertEqual(len(row[1]), 10)
        self.assertEqual(len(row[2]), 100)
        for x in range(0, 10):
            self.assertTrue(
                common.allequal(row[0][x], np.arange(x, dtype="int32"))
            )
        for x in range(5, 15):
            self.assertTrue(
                common.allequal(row[1][x - 5], np.arange(x, dtype="int32"))
            )
        for x in range(0, 100):
            self.assertTrue(
                common.allequal(row[2][x], np.arange(x, dtype="int32"))
            )

    def test04_slices(self):
        """Checking reads with a start, stop & step values"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test04_slices..." % self.__class__.__name__)

        vlarray = self.h5file.root.vlarray

        # Choose a small value for buffer size
        vlarray._nrowsinbuf = 3

        # Read some rows:
        row = []
        row.append(vlarray[0:10:2])
        row.append(vlarray[5:15:3])
        row.append(vlarray[0:100:20])
        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, self.nrows)
        self.assertEqual(len(row[0]), 5)
        self.assertEqual(len(row[1]), 4)
        self.assertEqual(len(row[2]), 5)
        for x in range(0, 10, 2):
            self.assertTrue(
                common.allequal(row[0][x // 2], np.arange(x, dtype="int32"))
            )
        for x in range(5, 15, 3):
            self.assertTrue(
                common.allequal(
                    row[1][(x - 5) // 3], np.arange(x, dtype="int32")
                )
            )
        for x in range(0, 100, 20):
            self.assertTrue(
                common.allequal(row[2][x // 20], np.arange(x, dtype="int32"))
            )

    def test04bnp_slices(self):
        """Checking reads with start, stop & step values (numpy indices)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test04np_slices..." % self.__class__.__name__)

        vlarray = self.h5file.root.vlarray

        # Choose a small value for buffer size
        vlarray._nrowsinbuf = 3

        # Read some rows:
        row = []
        row.append(vlarray[np.int8(0) : np.int8(10) : np.int8(2)])
        row.append(vlarray[np.int8(5) : np.int8(15) : np.int8(3)])
        row.append(vlarray[np.int8(0) : np.int8(100) : np.int8(20)])
        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, self.nrows)
        self.assertEqual(len(row[0]), 5)
        self.assertEqual(len(row[1]), 4)
        self.assertEqual(len(row[2]), 5)
        for x in range(0, 10, 2):
            self.assertTrue(
                common.allequal(row[0][x // 2], np.arange(x, dtype="int32"))
            )
        for x in range(5, 15, 3):
            self.assertTrue(
                common.allequal(
                    row[1][(x - 5) // 3], np.arange(x, dtype="int32")
                )
            )
        for x in range(0, 100, 20):
            self.assertTrue(
                common.allequal(row[2][x // 20], np.arange(x, dtype="int32"))
            )

    def test05_out_of_range(self):
        """Checking out of range reads"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test05_out_of_range..." % self.__class__.__name__
            )

        vlarray = self.h5file.root.vlarray

        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)

        with self.assertRaises(IndexError):
            row = vlarray[1000]
            print("row-->", row)

    def test05np_out_of_range(self):
        """Checking out of range reads (numpy indexes)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test05np_out_of_range..." % self.__class__.__name__
            )

        vlarray = self.h5file.root.vlarray

        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)

        with self.assertRaises(IndexError):
            row = vlarray[np.int32(1000)]
            print("row-->", row)


class SetRangeTestCase(common.TempFileMixin, common.PyTablesTestCase):
    nrows = 100
    open_mode = "w"
    compress = 0
    complib = "zlib"  # Default compression library

    def setUp(self):
        super().setUp()
        self.rootgroup = self.h5file.root
        self.populateFile()
        self._reopen(mode="a")

    def populateFile(self):
        group = self.rootgroup
        filters = tb.Filters(complevel=self.compress, complib=self.complib)
        vlarray = self.h5file.create_vlarray(
            group,
            "vlarray",
            tb.Int32Atom(),
            "ragged array if ints",
            filters=filters,
            expectedrows=1000,
        )

        # Fill it with 100 rows with variable length
        for i in range(self.nrows):
            vlarray.append(list(range(i)))

    def test01_start(self):
        """Checking updates that modifies a complete row"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01_start..." % self.__class__.__name__)

        vlarray = self.h5file.root.vlarray

        # Modify some rows:
        vlarray[0] = vlarray[0] * 2 + 3
        vlarray[10] = vlarray[10] * 2 + 3
        vlarray[99] = vlarray[99] * 2 + 3

        # Read some rows:
        row = []
        row.append(vlarray.read(0)[0])
        row.append(vlarray.read(10)[0])
        row.append(vlarray.read(99)[0])
        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, self.nrows)
        self.assertEqual(len(row[0]), 0)
        self.assertEqual(len(row[1]), 10)
        self.assertEqual(len(row[2]), 99)
        self.assertTrue(
            common.allequal(row[0], np.arange(0, dtype="int32") * 2 + 3)
        )
        self.assertTrue(
            common.allequal(row[1], np.arange(10, dtype="int32") * 2 + 3)
        )
        self.assertTrue(
            common.allequal(row[2], np.arange(99, dtype="int32") * 2 + 3)
        )

    def test01np_start(self):
        """Checking updates that modifies a complete row"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01np_start..." % self.__class__.__name__)

        vlarray = self.h5file.root.vlarray

        # Modify some rows:
        vlarray[np.int8(0)] = vlarray[np.int16(0)] * 2 + 3
        vlarray[np.int8(10)] = vlarray[np.int8(10)] * 2 + 3
        vlarray[np.int32(99)] = vlarray[np.int64(99)] * 2 + 3

        # Read some rows:
        row = []
        row.append(vlarray.read(np.int8(0))[0])
        row.append(vlarray.read(np.int8(10))[0])
        row.append(vlarray.read(np.int8(99))[0])
        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, self.nrows)
        self.assertEqual(len(row[0]), 0)
        self.assertEqual(len(row[1]), 10)
        self.assertEqual(len(row[2]), 99)
        self.assertTrue(
            common.allequal(row[0], np.arange(0, dtype="int32") * 2 + 3)
        )
        self.assertTrue(
            common.allequal(row[1], np.arange(10, dtype="int32") * 2 + 3)
        )
        self.assertTrue(
            common.allequal(row[2], np.arange(99, dtype="int32") * 2 + 3)
        )

    def test02_partial(self):
        """Checking updates with only a part of a row"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02_partial..." % self.__class__.__name__)

        vlarray = self.h5file.root.vlarray

        # Modify some rows:
        vlarray[0] = vlarray[0] * 2 + 3
        vlarray[10] = vlarray[10] * 2 + 3
        vlarray[96] = vlarray[99][3:] * 2 + 3

        # Read some rows:
        row = []
        row.append(vlarray.read(0)[0])
        row.append(vlarray.read(10)[0])
        row.append(vlarray.read(96)[0])
        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, self.nrows)
        self.assertEqual(len(row[0]), 0)
        self.assertEqual(len(row[1]), 10)
        self.assertEqual(len(row[2]), 96)
        self.assertTrue(
            common.allequal(row[0], np.arange(0, dtype="int32") * 2 + 3)
        )
        self.assertTrue(
            common.allequal(row[1], np.arange(10, dtype="int32") * 2 + 3)
        )
        a = np.arange(3, 99, dtype="int32")
        a = a * 2 + 3
        self.assertTrue(common.allequal(row[2], a))

    def test03a_several_rows(self):
        """Checking updating several rows at once (slice style)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test03a_several_rows..." % self.__class__.__name__
            )

        vlarray = self.h5file.root.vlarray

        # Modify some rows:
        vlarray[3:6] = (
            vlarray[3] * 2 + 3,
            vlarray[4] * 2 + 3,
            vlarray[5] * 2 + 3,
        )

        # Read some rows:
        row = []
        row.append(vlarray.read(3)[0])
        row.append(vlarray.read(4)[0])
        row.append(vlarray.read(5)[0])
        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, self.nrows)
        self.assertEqual(len(row[0]), 3)
        self.assertEqual(len(row[1]), 4)
        self.assertEqual(len(row[2]), 5)
        self.assertTrue(
            common.allequal(row[0], np.arange(3, dtype="int32") * 2 + 3)
        )
        self.assertTrue(
            common.allequal(row[1], np.arange(4, dtype="int32") * 2 + 3)
        )
        self.assertTrue(
            common.allequal(row[2], np.arange(5, dtype="int32") * 2 + 3)
        )

    def test03b_several_rows(self):
        """Checking updating several rows at once (list style)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test03b_several_rows..." % self.__class__.__name__
            )

        vlarray = self.h5file.root.vlarray

        # Modify some rows:
        vlarray[[0, 10, 96]] = (
            vlarray[0] * 2 + 3,
            vlarray[10] * 2 + 3,
            vlarray[96] * 2 + 3,
        )

        # Read some rows:
        row = []
        row.append(vlarray.read(0)[0])
        row.append(vlarray.read(10)[0])
        row.append(vlarray.read(96)[0])
        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, self.nrows)
        self.assertEqual(len(row[0]), 0)
        self.assertEqual(len(row[1]), 10)
        self.assertEqual(len(row[2]), 96)
        self.assertTrue(
            common.allequal(row[0], np.arange(0, dtype="int32") * 2 + 3)
        )
        self.assertTrue(
            common.allequal(row[1], np.arange(10, dtype="int32") * 2 + 3)
        )
        self.assertTrue(
            common.allequal(row[2], np.arange(96, dtype="int32") * 2 + 3)
        )

    def test03c_several_rows(self):
        """Checking updating several rows at once (NumPy's where style)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test03c_several_rows..." % self.__class__.__name__
            )

        vlarray = self.h5file.root.vlarray

        # Modify some rows:
        vlarray[(np.array([0, 10, 96]),)] = (
            vlarray[0] * 2 + 3,
            vlarray[10] * 2 + 3,
            vlarray[96] * 2 + 3,
        )

        # Read some rows:
        row = []
        row.append(vlarray.read(0)[0])
        row.append(vlarray.read(10)[0])
        row.append(vlarray.read(96)[0])
        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)
            print("Second row in vlarray ==>", row[1])

        self.assertEqual(vlarray.nrows, self.nrows)
        self.assertEqual(len(row[0]), 0)
        self.assertEqual(len(row[1]), 10)
        self.assertEqual(len(row[2]), 96)
        self.assertTrue(
            common.allequal(row[0], np.arange(0, dtype="int32") * 2 + 3)
        )
        self.assertTrue(
            common.allequal(row[1], np.arange(10, dtype="int32") * 2 + 3)
        )
        self.assertTrue(
            common.allequal(row[2], np.arange(96, dtype="int32") * 2 + 3)
        )

    def test04_out_of_range(self):
        """Checking out of range updates (first index)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test04_out_of_range..." % self.__class__.__name__
            )

        vlarray = self.h5file.root.vlarray

        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)

        with self.assertRaises(IndexError):
            vlarray[1000] = [1]

    def test05_value_error(self):
        """Checking out value errors"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test05_value_error..." % self.__class__.__name__)

        vlarray = self.h5file.root.vlarray

        if common.verbose:
            print("Nrows in", vlarray._v_pathname, ":", vlarray.nrows)

        with self.assertRaises(ValueError):
            vlarray[10] = [1] * 100


class CopyTestCase(common.TempFileMixin, common.PyTablesTestCase):
    close = True

    def test01a_copy(self):
        """Checking VLArray.copy() method."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01a_copy..." % self.__class__.__name__)

        # Create an Vlarray
        arr = tb.Int16Atom(shape=2)
        array1 = self.h5file.create_vlarray(
            self.h5file.root, "array1", arr, "title array1"
        )
        array1.flavor = "python"
        array1.append([[2, 3]])
        array1.append(())  # an empty row
        array1.append([[3, 457], [2, 4]])

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
            print("array1-->", repr(array1))
            print("array2-->", repr(array2))
            print("array1[:]-->", repr(array1.read()))
            print("array2[:]-->", repr(array2.read()))
            print("attrs array1-->", repr(array1.attrs))
            print("attrs array2-->", repr(array2.attrs))

        # Check that all the elements are equal
        self.assertEqual(array1.read(), array2.read())

        # Assert other properties in array
        self.assertEqual(array1.nrows, array2.nrows)
        self.assertEqual(array1.shape, array2.shape)
        self.assertEqual(array1.flavor, array2.flavor)
        self.assertEqual(array1.atom.dtype, array2.atom.dtype)
        self.assertEqual(repr(array1.atom), repr(array2.atom))

        self.assertEqual(array1.title, array2.title)

    def test01b_copy(self):
        """Checking VLArray.copy() method (Pseudo-atom case)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01b_copy..." % self.__class__.__name__)

        # Create an Vlarray
        arr = tb.VLStringAtom()
        array1 = self.h5file.create_vlarray(
            self.h5file.root, "array1", arr, "title array1"
        )
        array1.flavor = "python"
        array1.append(b"a string")
        array1.append(b"")  # an empty row
        array1.append(b"another string")

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
            print("array1-->", repr(array1))
            print("array2-->", repr(array2))
            print("array1[:]-->", repr(array1.read()))
            print("array2[:]-->", repr(array2.read()))
            print("attrs array1-->", repr(array1.attrs))
            print("attrs array2-->", repr(array2.attrs))

        # Check that all the elements are equal
        self.assertEqual(array1.read(), array2.read())

        # Assert other properties in array
        self.assertEqual(array1.nrows, array2.nrows)
        self.assertEqual(array1.shape, array2.shape)
        self.assertEqual(array1.flavor, array2.flavor)
        self.assertEqual(array1.atom.type, array2.atom.type)
        self.assertEqual(repr(array1.atom), repr(array2.atom))

        self.assertEqual(array1.title, array2.title)

    def test02_copy(self):
        """Checking VLArray.copy() method (where specified)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02_copy..." % self.__class__.__name__)

        # Create an VLArray
        arr = tb.Int16Atom(shape=2)
        array1 = self.h5file.create_vlarray(
            self.h5file.root, "array1", arr, "title array1"
        )
        array1.flavor = "python"
        array1.append([[2, 3]])
        array1.append(())  # an empty row
        array1.append([[3, 457], [2, 4]])

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
            print("array1-->", repr(array1))
            print("array2-->", repr(array2))
            print("array1-->", array1.read())
            print("array2-->", array2.read())
            print("attrs array1-->", repr(array1.attrs))
            print("attrs array2-->", repr(array2.attrs))

        # Check that all the elements are equal
        self.assertEqual(array1.read(), array2.read())

        # Assert other properties in array
        self.assertEqual(array1.nrows, array2.nrows)
        self.assertEqual(array1.shape, array2.shape)
        self.assertEqual(array1.flavor, array2.flavor)
        self.assertEqual(array1.atom.dtype, array2.atom.dtype)
        self.assertEqual(repr(array1.atom), repr(array1.atom))
        self.assertEqual(array1.title, array2.title)

    def test03_copy(self):
        """Checking VLArray.copy() method ('python' flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03_copy..." % self.__class__.__name__)

        # Create an VLArray
        atom = tb.Int16Atom(shape=2)
        array1 = self.h5file.create_vlarray(
            self.h5file.root, "array1", atom, title="title array1"
        )
        array1.flavor = "python"
        array1.append(((2, 3),))
        array1.append(())  # an empty row
        array1.append(((3, 457), (2, 4)))

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

        # Assert other properties in array
        self.assertEqual(array1.nrows, array2.nrows)
        self.assertEqual(array1.shape, array2.shape)
        self.assertEqual(array1.flavor, array2.flavor)  # Very important here
        self.assertEqual(array1.atom.dtype, array2.atom.dtype)
        self.assertEqual(repr(array1.atom), repr(array1.atom))
        self.assertEqual(array1.title, array2.title)

    def test04_copy(self):
        """Checking VLArray.copy() method (checking title copying)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test04_copy..." % self.__class__.__name__)

        # Create an VLArray
        atom = tb.Int16Atom(shape=2)
        array1 = self.h5file.create_vlarray(
            self.h5file.root, "array1", atom=atom, title="title array1"
        )
        array1.append(((2, 3),))
        array1.append(())  # an empty row
        array1.append(((3, 457), (2, 4)))

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
        """Checking VLArray.copy() method (user attributes copied)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test05_copy..." % self.__class__.__name__)

        # Create an Array
        atom = tb.Int16Atom(shape=2)
        array1 = self.h5file.create_vlarray(
            self.h5file.root, "array1", atom=atom, title="title array1"
        )
        array1.append(((2, 3),))
        array1.append(())  # an empty row
        array1.append(((3, 457), (2, 4)))

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

    def notest05b_copy(self):
        """Checking VLArray.copy() method (user attributes not copied)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test05b_copy..." % self.__class__.__name__)

        # Create an VLArray
        atom = tb.Int16Atom(shape=2)
        array1 = self.h5file.create_vlarray(
            self.h5file.root, "array1", atom=atom, title="title array1"
        )
        array1.append(((2, 3),))
        array1.append(())  # an empty row
        array1.append(((3, 457), (2, 4)))

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
        self.assertEqual(array2.attrs.attr1, None)
        self.assertEqual(array2.attrs.attr2, None)


class CloseCopyTestCase(CopyTestCase):
    close = 1


class OpenCopyTestCase(CopyTestCase):
    close = 0


class CopyIndexTestCase(common.TempFileMixin, common.PyTablesTestCase):

    def test01_index(self):
        """Checking VLArray.copy() method with indexes."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01_index..." % self.__class__.__name__)

        # Create an VLArray
        atom = tb.Int32Atom(shape=(2,))
        array1 = self.h5file.create_vlarray(
            self.h5file.root, "array1", atom, "t array1"
        )
        array1.flavor = "python"

        # The next creates 20 rows of variable length
        r = []
        for row in range(20):
            r.append([[row, row + 1]])
            array1.append([row, row + 1])

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="a")
            array1 = self.h5file.root.array1

        # Copy to another array
        array2 = array1.copy(
            "/", "array2", start=self.start, stop=self.stop, step=self.step
        )

        r2 = r[self.start : self.stop : self.step]
        if common.verbose:
            print("r2-->", r2)
            print("array2-->", array2[:])
            print("attrs array1-->", repr(array1.attrs))
            print("attrs array2-->", repr(array2.attrs))
            print("nrows in array2-->", array2.nrows)
            print("and it should be-->", len(r2))

        # Check that all the elements are equal
        self.assertEqual(r2, array2[:])

        # Assert the number of rows in array
        self.assertEqual(len(r2), array2.nrows)


class CopyIndex1TestCase(CopyIndexTestCase):
    close = 0
    start = 0
    stop = 7
    step = 1


class CopyIndex2TestCase(CopyIndexTestCase):
    close = 1
    start = 0
    stop = -1
    step = 1


class CopyIndex3TestCase(CopyIndexTestCase):
    close = 0
    start = 1
    stop = 7
    step = 1


class CopyIndex4TestCase(CopyIndexTestCase):
    close = 1
    start = 0
    stop = 6
    step = 1


class CopyIndex5TestCase(CopyIndexTestCase):
    close = 0
    start = 3
    stop = 7
    step = 1


class CopyIndex6TestCase(CopyIndexTestCase):
    close = 1
    start = 3
    stop = 6
    step = 2


class CopyIndex7TestCase(CopyIndexTestCase):
    close = 0
    start = 0
    stop = 7
    step = 10


class CopyIndex8TestCase(CopyIndexTestCase):
    close = 1
    start = 6
    stop = -1  # Negative values means starting from the end
    step = 1


class CopyIndex9TestCase(CopyIndexTestCase):
    close = 0
    start = 3
    stop = 4
    step = 1


class CopyIndex10TestCase(CopyIndexTestCase):
    close = 1
    start = 3
    stop = 4
    step = 2


class CopyIndex11TestCase(CopyIndexTestCase):
    close = 0
    start = -3
    stop = -1
    step = 2


class CopyIndex12TestCase(CopyIndexTestCase):
    close = 1
    start = -1  # Should point to the last element
    stop = None  # None should mean the last element (including it)
    step = 1


class ChunkshapeTestCase(common.TempFileMixin, common.PyTablesTestCase):

    def setUp(self):
        super().setUp()
        atom = tb.Int32Atom(shape=(2,))
        self.h5file.create_vlarray(
            "/", "vlarray", atom=atom, title="t array1", chunkshape=13
        )

    def test00(self):
        """Test setting the chunkshape in a table (no reopen)."""

        vla = self.h5file.root.vlarray
        if common.verbose:
            print("chunkshape-->", vla.chunkshape)
        self.assertEqual(vla.chunkshape, (13,))

    def test01(self):
        """Test setting the chunkshape in a table (reopen)."""

        self.h5file.close()
        self.h5file = tb.open_file(self.h5fname, "r")
        vla = self.h5file.root.vlarray
        if common.verbose:
            print("chunkshape-->", vla.chunkshape)
        self.assertEqual(vla.chunkshape, (13,))


class VLUEndianTestCase(common.PyTablesTestCase):
    def setUp(self):
        super().setUp()
        self.h5fname = common.test_filename("vlunicode_endian.h5")
        self.h5file = tb.open_file(self.h5fname)

    def tearDown(self):
        self.h5file.close()
        super().tearDown()

    def test(self):
        """Accessing ``vlunicode`` data of a different endianness."""

        bedata = self.h5file.root.vlunicode_big[0]
        ledata = self.h5file.root.vlunicode_little[0]
        self.assertEqual(bedata, "para\u0140lel")
        self.assertEqual(ledata, "para\u0140lel")


class TruncateTestCase(common.TempFileMixin, common.PyTablesTestCase):

    def setUp(self):
        super().setUp()

        # Create an VLArray
        arr = tb.Int16Atom(dflt=3)
        array1 = self.h5file.create_vlarray(
            self.h5file.root, "array1", arr, "title array1"
        )

        # Add a couple of rows
        array1.append(np.array([456, 2], dtype="int16"))
        array1.append(np.array([3], dtype="int16"))

    def test00_truncate(self):
        """Checking VLArray.truncate() method (truncating to 0 rows)"""

        array1 = self.h5file.root.array1
        # Truncate to 0 elements
        array1.truncate(0)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen()
            array1 = self.h5file.root.array1

        if common.verbose:
            print("array1-->", array1.read())

        self.assertEqual(array1.nrows, 0)
        self.assertEqual(array1[:], [])

    def test01_truncate(self):
        """Checking VLArray.truncate() method (truncating to 1 rows)"""

        array1 = self.h5file.root.array1
        # Truncate to 1 element
        array1.truncate(1)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen()
            array1 = self.h5file.root.array1

        if common.verbose:
            print("array1-->", array1.read())

        self.assertEqual(array1.nrows, 1)
        self.assertTrue(
            common.allequal(array1[0], np.array([456, 2], dtype="int16"))
        )

    def test02_truncate(self):
        """Checking VLArray.truncate() method (truncating to == self.nrows)"""

        array1 = self.h5file.root.array1
        # Truncate to 2 elements
        array1.truncate(2)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen()
            array1 = self.h5file.root.array1

        if common.verbose:
            print("array1-->", array1.read())

        self.assertEqual(array1.nrows, 2)
        self.assertTrue(
            common.allequal(array1[0], np.array([456, 2], dtype="int16"))
        )
        self.assertTrue(
            common.allequal(array1[1], np.array([3], dtype="int16"))
        )

    def test03_truncate(self):
        """Checking VLArray.truncate() method (truncating to > self.nrows)"""

        array1 = self.h5file.root.array1
        # Truncate to 4 elements
        array1.truncate(4)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen()
            array1 = self.h5file.root.array1

        if common.verbose:
            print("array1-->", array1.read())

        self.assertEqual(array1.nrows, 4)

        # Check the original values
        self.assertTrue(
            common.allequal(array1[0], np.array([456, 2], dtype="int16"))
        )
        self.assertTrue(
            common.allequal(array1[1], np.array([3], dtype="int16"))
        )

        # Check that the added rows are empty
        self.assertTrue(
            common.allequal(array1[2], np.array([], dtype="int16"))
        )
        self.assertTrue(
            common.allequal(array1[3], np.array([], dtype="int16"))
        )


class TruncateOpenTestCase(TruncateTestCase):
    close = 0


class TruncateCloseTestCase(TruncateTestCase):
    close = 1


class PointSelectionTestCase(common.TempFileMixin, common.PyTablesTestCase):

    def setUp(self):
        super().setUp()

        # The next are valid selections for both NumPy and PyTables
        self.working_keyset = [
            [],  # empty list
            [2],  # single-entry list
            [0, 2],  # list
            [0, -2],  # negative values
            ([0, 2],),  # tuple of list
            np.array([], dtype="i4"),  # empty array
            np.array([1], dtype="i4"),  # single-entry array
            np.array([True, False, True]),  # array of bools
        ]

        # The next are invalid selections for VLArrays
        self.not_working_keyset = [
            [1, 2, 100],  # coordinate 100 > len(vlarray)
            ([True, False, True],),  # tuple of bools
        ]

        # Create a sample array
        arr1 = np.array([5, 6], dtype="i4")
        arr2 = np.array([5, 6, 7], dtype="i4")
        arr3 = np.array([5, 6, 9, 8], dtype="i4")
        self.nparr = np.array([arr1, arr2, arr3], dtype="object")

        # Create the VLArray
        self.vlarr = self.h5file.create_vlarray(
            self.h5file.root, "vlarray", tb.Int32Atom()
        )
        self.vlarr.append(arr1)
        self.vlarr.append(arr2)
        self.vlarr.append(arr3)

    def test01a_read(self):
        """Test for point-selections (read, boolean keys)."""

        nparr = self.nparr
        vlarr = self.vlarr
        for key in self.working_keyset:
            if common.verbose:
                print("Selection to test:", repr(key))
            a = nparr[key].tolist()
            b = vlarr[key]
            # if common.verbose:
            #     print "NumPy selection:", a, type(a)
            #     print "PyTables selection:", b, type(b)
            self.assertEqual(
                repr(a),
                repr(b),
                "NumPy array and PyTables selections does not match.",
            )

    def test01b_read(self):
        """Test for point-selections (not working selections, read)."""

        vlarr = self.vlarr
        for key in self.not_working_keyset:
            if common.verbose:
                print("Selection to test:", key)
            self.assertRaises(IndexError, vlarr.__getitem__, key)


class SizeInMemoryPropertyTestCase(
    common.TempFileMixin, common.PyTablesTestCase
):
    def create_array(self, atom, complevel):
        filters = tb.Filters(complevel=complevel, complib="blosc")
        self.array = self.h5file.create_vlarray(
            "/", "vlarray", atom=atom, filters=filters
        )

    def test_zero_length(self):
        atom = tb.Int32Atom()
        complevel = 0
        self.create_array(atom, complevel)
        self.assertEqual(self.array.size_in_memory, 0)

    def int_tests(self, complevel, flavor):
        atom = tb.Int32Atom()
        self.create_array(atom, complevel)
        self.array.flavor = flavor
        expected_size = 0
        for i in range(10):
            row = np.arange((i + 1) * 10, dtype="i4")
            self.array.append(row)
            expected_size += row.nbytes
        return expected_size

    def test_numpy_int_numpy_flavor(self):
        complevel = 0
        flavor = "numpy"
        expected_size = self.int_tests(complevel, flavor)
        self.assertEqual(self.array.size_in_memory, expected_size)

    # compression will have no effect, since this is uncompressed size
    def test_numpy_int_numpy_flavor_compressed(self):
        complevel = 1
        flavor = "numpy"
        expected_size = self.int_tests(complevel, flavor)
        self.assertEqual(self.array.size_in_memory, expected_size)

    # flavor will have no effect on what's stored in HDF5 file
    def test_numpy_int_python_flavor(self):
        complevel = 0
        flavor = "python"
        expected_size = self.int_tests(complevel, flavor)
        self.assertEqual(self.array.size_in_memory, expected_size)

    # this relies on knowledge of the implementation, so it's not
    # a great test
    def test_object_atom(self):
        atom = tb.ObjectAtom()
        complevel = 0
        self.create_array(atom, complevel)
        obj = [1, 2, 3]
        for i in range(10):
            self.array.append(obj)
        pickle_array = atom.toarray(obj)
        expected_size = 10 * pickle_array.nbytes
        self.assertEqual(self.array.size_in_memory, expected_size)


class SizeOnDiskPropertyTestCase(
    common.TempFileMixin, common.PyTablesTestCase
):
    def create_array(self, atom, complevel):
        filters = tb.Filters(complevel=complevel, complib="blosc")
        self.h5file.create_vlarray("/", "vlarray", atom, filters=filters)
        self.array = self.h5file.get_node("/", "vlarray")

    def test_not_implemented(self):
        atom = tb.IntAtom()
        complevel = 0
        self.create_array(atom, complevel)
        self.assertRaises(
            NotImplementedError, getattr, self.array, "size_on_disk"
        )


class AccessClosedTestCase(common.TempFileMixin, common.PyTablesTestCase):
    def setUp(self):
        super().setUp()
        self.array = self.h5file.create_vlarray(
            self.h5file.root, "array", atom=tb.StringAtom(8)
        )
        self.array.append([str(i) for i in range(5, 5005, 100)])

    def test_read(self):
        self.h5file.close()
        self.assertRaises(tb.ClosedNodeError, self.array.read)

    def test_getitem(self):
        self.h5file.close()
        self.assertRaises(tb.ClosedNodeError, self.array.__getitem__, 0)

    def test_setitem(self):
        self.h5file.close()
        self.assertRaises(tb.ClosedNodeError, self.array.__setitem__, 0, "0")

    def test_append(self):
        self.h5file.close()
        self.assertRaises(tb.ClosedNodeError, self.array.append, "xxxxxxxxx")


class TestCreateVLArrayArgs(common.TempFileMixin, common.PyTablesTestCase):
    obj = np.array([1, 2, 3])
    where = "/"
    name = "vlarray"
    atom = tb.Atom.from_dtype(obj.dtype)
    title = "title"
    filters = None
    expectedrows = None
    chunkshape = None
    byteorder = None
    createparents = False

    def test_positional_args_01(self):
        self.h5file.create_vlarray(
            self.where,
            self.name,
            self.atom,
            self.title,
            self.filters,
            self.expectedrows,
        )
        self.h5file.close()

        self.h5file = tb.open_file(self.h5fname)
        ptarr = self.h5file.get_node(self.where, self.name)

        self.assertEqual(ptarr.title, self.title)
        self.assertEqual(ptarr.shape, (0,))
        self.assertEqual(ptarr.nrows, 0)
        self.assertEqual(ptarr.atom, self.atom)
        self.assertEqual(ptarr.atom.dtype, self.atom.dtype)

    def test_positional_args_02(self):
        ptarr = self.h5file.create_vlarray(
            self.where,
            self.name,
            self.atom,
            self.title,
            self.filters,
            self.expectedrows,
        )
        ptarr.append(self.obj)
        self.h5file.close()

        self.h5file = tb.open_file(self.h5fname)
        ptarr = self.h5file.get_node(self.where, self.name)
        nparr = ptarr.read()[0]

        self.assertEqual(ptarr.title, self.title)
        self.assertEqual(ptarr.shape, (1,))
        self.assertEqual(ptarr[0].shape, self.obj.shape)
        self.assertEqual(ptarr.nrows, 1)
        self.assertEqual(ptarr.atom, self.atom)
        self.assertEqual(ptarr.atom.dtype, self.atom.dtype)
        self.assertTrue(common.allequal(self.obj, nparr))

    def test_positional_args_obj(self):
        self.h5file.create_vlarray(
            self.where,
            self.name,
            None,
            self.title,
            self.filters,
            self.expectedrows,
            self.chunkshape,
            self.byteorder,
            self.createparents,
            self.obj,
        )
        self.h5file.close()

        self.h5file = tb.open_file(self.h5fname)
        ptarr = self.h5file.get_node(self.where, self.name)
        nparr = ptarr.read()[0]

        self.assertEqual(ptarr.title, self.title)
        self.assertEqual(ptarr.shape, (1,))
        self.assertEqual(ptarr[0].shape, self.obj.shape)
        self.assertEqual(ptarr.nrows, 1)
        self.assertEqual(ptarr.atom, self.atom)
        self.assertEqual(ptarr.atom.dtype, self.atom.dtype)
        self.assertTrue(common.allequal(self.obj, nparr))

    def test_kwargs_obj(self):
        self.h5file.create_vlarray(
            self.where, self.name, title=self.title, obj=self.obj
        )
        self.h5file.close()

        self.h5file = tb.open_file(self.h5fname)
        ptarr = self.h5file.get_node(self.where, self.name)
        nparr = ptarr.read()[0]

        self.assertEqual(ptarr.title, self.title)
        self.assertEqual(ptarr.shape, (1,))
        self.assertEqual(ptarr[0].shape, self.obj.shape)
        self.assertEqual(ptarr.nrows, 1)
        self.assertEqual(ptarr.atom, self.atom)
        self.assertEqual(ptarr.atom.dtype, self.atom.dtype)
        self.assertTrue(common.allequal(self.obj, nparr))

    def test_kwargs_atom_01(self):
        ptarr = self.h5file.create_vlarray(
            self.where, self.name, title=self.title, atom=self.atom
        )
        ptarr.append(self.obj)
        self.h5file.close()

        self.h5file = tb.open_file(self.h5fname)
        ptarr = self.h5file.get_node(self.where, self.name)
        nparr = ptarr.read()[0]

        self.assertEqual(ptarr.title, self.title)
        self.assertEqual(ptarr.shape, (1,))
        self.assertEqual(ptarr[0].shape, self.obj.shape)
        self.assertEqual(ptarr.nrows, 1)
        self.assertEqual(ptarr.atom, self.atom)
        self.assertEqual(ptarr.atom.dtype, self.atom.dtype)
        self.assertTrue(common.allequal(self.obj, nparr))

    def test_kwargs_atom_02(self):
        ptarr = self.h5file.create_vlarray(
            self.where, self.name, title=self.title, atom=self.atom
        )
        # ptarr.append(self.obj)
        self.h5file.close()

        self.h5file = tb.open_file(self.h5fname)
        ptarr = self.h5file.get_node(self.where, self.name)

        self.assertEqual(ptarr.title, self.title)
        self.assertEqual(ptarr.shape, (0,))
        self.assertEqual(ptarr.nrows, 0)
        self.assertEqual(ptarr.atom, self.atom)
        self.assertEqual(ptarr.atom.dtype, self.atom.dtype)

    def test_kwargs_obj_atom(self):
        ptarr = self.h5file.create_vlarray(
            self.where,
            self.name,
            title=self.title,
            obj=self.obj,
            atom=self.atom,
        )
        self.h5file.close()

        self.h5file = tb.open_file(self.h5fname)
        ptarr = self.h5file.get_node(self.where, self.name)
        nparr = ptarr.read()[0]

        self.assertEqual(ptarr.title, self.title)
        self.assertEqual(ptarr.shape, (1,))
        self.assertEqual(ptarr[0].shape, self.obj.shape)
        self.assertEqual(ptarr.nrows, 1)
        self.assertEqual(ptarr.atom, self.atom)
        self.assertEqual(ptarr.atom.dtype, self.atom.dtype)
        self.assertTrue(common.allequal(self.obj, nparr))

    def test_kwargs_obj_atom_error(self):
        atom = tb.Atom.from_dtype(np.dtype("complex"))
        # shape = self.shape + self.shape
        self.assertRaises(
            TypeError,
            self.h5file.create_vlarray,
            self.where,
            self.name,
            title=self.title,
            obj=self.obj,
            atom=atom,
        )


def suite():
    theSuite = common.unittest.TestSuite()
    niter = 1

    for n in range(niter):
        theSuite.addTest(common.make_suite(BasicNumPyTestCase))
        theSuite.addTest(common.make_suite(BasicPythonTestCase))
        theSuite.addTest(common.make_suite(ZlibComprTestCase))
        theSuite.addTest(common.make_suite(BloscComprTestCase))
        theSuite.addTest(common.make_suite(BloscShuffleComprTestCase))
        theSuite.addTest(common.make_suite(BloscBitShuffleComprTestCase))
        theSuite.addTest(common.make_suite(BloscBloscLZComprTestCase))
        theSuite.addTest(common.make_suite(BloscLZ4ComprTestCase))
        theSuite.addTest(common.make_suite(BloscLZ4HCComprTestCase))
        theSuite.addTest(common.make_suite(BloscSnappyComprTestCase))
        theSuite.addTest(common.make_suite(BloscZlibComprTestCase))
        theSuite.addTest(common.make_suite(BloscZstdComprTestCase))
        theSuite.addTest(common.make_suite(LZOComprTestCase))
        theSuite.addTest(common.make_suite(Bzip2ComprTestCase))
        theSuite.addTest(common.make_suite(TypesReopenTestCase))
        theSuite.addTest(common.make_suite(TypesNoReopenTestCase))
        theSuite.addTest(common.make_suite(MDTypesNumPyTestCase))
        theSuite.addTest(common.make_suite(OpenAppendShapeTestCase))
        theSuite.addTest(common.make_suite(CloseAppendShapeTestCase))
        theSuite.addTest(common.make_suite(PythonFlavorTestCase))
        theSuite.addTest(common.make_suite(NumPyFlavorTestCase))
        theSuite.addTest(common.make_suite(ReadRangeTestCase))
        theSuite.addTest(common.make_suite(GetItemRangeTestCase))
        theSuite.addTest(common.make_suite(SetRangeTestCase))
        theSuite.addTest(common.make_suite(ShuffleComprTestCase))
        theSuite.addTest(common.make_suite(CloseCopyTestCase))
        theSuite.addTest(common.make_suite(OpenCopyTestCase))
        theSuite.addTest(common.make_suite(CopyIndex1TestCase))
        theSuite.addTest(common.make_suite(CopyIndex2TestCase))
        theSuite.addTest(common.make_suite(CopyIndex3TestCase))
        theSuite.addTest(common.make_suite(CopyIndex4TestCase))
        theSuite.addTest(common.make_suite(CopyIndex5TestCase))
        theSuite.addTest(common.make_suite(CopyIndex6TestCase))
        theSuite.addTest(common.make_suite(CopyIndex7TestCase))
        theSuite.addTest(common.make_suite(CopyIndex8TestCase))
        theSuite.addTest(common.make_suite(CopyIndex9TestCase))
        theSuite.addTest(common.make_suite(CopyIndex10TestCase))
        theSuite.addTest(common.make_suite(CopyIndex11TestCase))
        theSuite.addTest(common.make_suite(CopyIndex12TestCase))
        theSuite.addTest(common.make_suite(ChunkshapeTestCase))
        theSuite.addTest(common.make_suite(VLUEndianTestCase))
        theSuite.addTest(common.make_suite(TruncateOpenTestCase))
        theSuite.addTest(common.make_suite(TruncateCloseTestCase))
        theSuite.addTest(common.make_suite(PointSelectionTestCase))
        theSuite.addTest(common.make_suite(SizeInMemoryPropertyTestCase))
        theSuite.addTest(common.make_suite(SizeOnDiskPropertyTestCase))
        theSuite.addTest(common.make_suite(AccessClosedTestCase))
        theSuite.addTest(common.make_suite(TestCreateVLArrayArgs))

    return theSuite


if __name__ == "__main__":
    common.parse_argv(sys.argv)
    common.print_versions()
    common.unittest.main(defaultTest="suite")
