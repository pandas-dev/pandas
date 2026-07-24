import sys
import tempfile
from pathlib import Path

import numpy as np

import tables as tb
from tables.tests import common

typecodes = ["b", "h", "i", "l", "q", "f", "d"]
# UInt64 checking disabled on win platforms
# because this type is not supported
if sys.platform != "win32":
    typecodes += ["B", "H", "I", "L", "Q", "F", "D"]
else:
    typecodes += ["B", "H", "I", "L", "F", "D"]
typecodes += ["b1"]  # boolean

if hasattr(tb, "Float16Atom"):
    typecodes.append("e")
if hasattr(tb, "Float96Atom") or hasattr(tb, "Float128Atom"):
    typecodes.append("g")
if hasattr(tb, "Complex192Atom") or hasattr(tb, "Conplex256Atom"):
    typecodes.append("G")

byteorder = {"little": "<", "big": ">"}[sys.byteorder]


class BasicTestCase(common.PyTablesTestCase):
    """Basic test for all the supported typecodes present in NumPy.

    All of them are included on PyTables.

    """

    endiancheck = 0

    def WriteRead(self, testArray):
        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running test for array with typecode '%s'"
                % testArray.dtype.char,
                end=" ",
            )
            print("for class check:", self.title)

        # Create an instance of HDF5 Table
        self.h5fname = tempfile.mktemp(".h5")
        try:
            with tb.open_file(self.h5fname, mode="w") as self.h5file:
                self.root = self.h5file.root

                # Create the array under root and name 'somearray'
                a = testArray
                self.h5file.create_array(
                    self.root, "somearray", a, "Some array"
                )

            # Re-open the file in read-only mode
            with tb.open_file(self.h5fname, mode="r") as self.h5file:
                self.root = self.h5file.root

                # Read the saved array
                b = self.root.somearray.read()

                # For cases that read returns a python type instead of a
                # numpy type
                if not hasattr(b, "shape"):
                    b = np.np.array(b, dtype=a.dtype.str)

                # Compare them. They should be equal.
                # if not allequal(a,b, "numpy") and common.verbose:
                if common.verbose:
                    print("Array written:", a)
                    print("Array written shape:", a.shape)
                    print("Array written itemsize:", a.itemsize)
                    print("Array written type:", a.dtype.char)
                    print("Array read:", b)
                    print("Array read shape:", b.shape)
                    print("Array read itemsize:", b.itemsize)
                    print("Array read type:", b.dtype.char)

                type_ = self.root.somearray.atom.type

                # Check strictly the array equality
                self.assertEqual(type(a), type(b))
                self.assertEqual(a.shape, b.shape)
                self.assertEqual(a.shape, self.root.somearray.shape)
                self.assertEqual(a.dtype, b.dtype)
                if a.dtype.char[0] == "S":
                    self.assertEqual(type_, "string")
                else:
                    self.assertEqual(a.dtype.base.name, type_)

                self.assertTrue(common.allequal(a, b, "numpy"))
        finally:
            # Then, delete the file
            if Path(self.h5fname).is_file():
                Path(self.h5fname).unlink()

    def test00_char(self):
        """Data integrity during recovery (character objects)"""

        a = np.array(self.tupleChar, "S" + str(len(self.tupleChar)))
        self.WriteRead(a)

    def test01_char_nc(self):
        """Data integrity during recovery (non-contiguous character objects)"""

        a = np.array(self.tupleChar, "S" + str(len(self.tupleChar)))
        if a.shape == ():
            b = a  # We cannot use the indexing notation
        else:
            b = a[::2]
            # Ensure that this numpy string is non-contiguous
            if a.shape[0] > 2:
                self.assertEqual(b.flags["CONTIGUOUS"], False)
        self.WriteRead(b)

    def test02_types(self):
        """Data integrity during recovery (numerical types)"""

        for typecode in typecodes:
            if self.tupleInt.shape:
                a = self.tupleInt.astype(typecode)
            else:
                # shape is the empty tuple ()
                a = np.array(self.tupleInt, dtype=typecode)
            self.WriteRead(a)

    def test03_types_nc(self):
        """Data integrity during recovery (non-contiguous numerical types)"""

        for typecode in typecodes:
            if self.tupleInt.shape:
                a = self.tupleInt.astype(typecode)
            else:
                # shape is the empty tuple ()
                a = np.array(self.tupleInt, dtype=typecode)

            # This should not be tested for the rank-0 case
            if len(a.shape) == 0:
                raise common.unittest.SkipTest
            b = a[::2]

            # Ensure that this array is non-contiguous (for non-trivial case)
            if a.shape[0] > 2:
                self.assertEqual(b.flags["CONTIGUOUS"], False)
            self.WriteRead(b)


class Basic0DOneTestCase(BasicTestCase):
    # Rank-0 case
    title = "Rank-0 case 1"
    tupleInt = np.array(3)
    tupleChar = "4"


class Basic0DTwoTestCase(BasicTestCase):
    # Rank-0 case
    title = "Rank-0 case 2"
    tupleInt = np.array(33)
    tupleChar = "44"


class Basic1DOneTestCase(BasicTestCase):
    # 1D case
    title = "Rank-1 case 1"
    tupleInt = np.array((3,))
    tupleChar = ("a",)


class Basic1DTwoTestCase(BasicTestCase):
    # 1D case
    title = "Rank-1 case 2"
    tupleInt = np.array((0, 4))
    tupleChar = ("aaa",)


class Basic1DThreeTestCase(BasicTestCase):
    # 1D case
    title = "Rank-1 case 3"
    tupleInt = np.array((3, 4, 5))
    tupleChar = (
        "aaaa",
        "bbb",
    )


class Basic2DTestCase(BasicTestCase):
    # 2D case
    title = "Rank-2 case 1"
    # tupleInt = reshape(np.array(np.arange((4)**2)), (4,)*2)
    tupleInt = np.ones((4,) * 2)
    tupleChar = [["aaa", "ddddd"], ["d", "ss"], ["s", "tt"]]


class Basic10DTestCase(BasicTestCase):
    # 10D case
    title = "Rank-10 case 1"
    # tupleInt = reshape(np.array(np.arange((2)**10)), (2,)*10)
    tupleInt = np.ones((2,) * 10)
    # tupleChar = reshape(np.array([1],dtype="S1"),(1,)*10)
    # The next tuple consumes far more time, so this
    # test should be run in common.heavy mode.
    tupleChar = np.array(tupleInt, dtype="S1")


# class Basic32DTestCase(BasicTestCase):
#     # 32D case (maximum)
#     tupleInt = reshape(np.array((22,)), (1,)*32)
#     # Strings seems to be very slow with somewhat large dimensions
#     # This should not be run unless the numarray people address this problem
#     # F. Alted 2006-01-04
#     tupleChar = np.array(tupleInt, dtype="S1")


class GroupsArrayTestCase(common.TempFileMixin, common.PyTablesTestCase):
    """This test class checks combinations of arrays with groups.

    It also uses arrays ranks which ranges until 10.

    """

    def test00_iterativeGroups(self):
        """Checking combinations of arrays with groups

        It also uses arrays ranks which ranges until 10.

        """

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test00_iterativeGroups..."
                % self.__class__.__name__
            )

        # Get the root group
        group = self.h5file.root

        i = 1
        for typecode in typecodes:
            # Create an array of typecode, with incrementally bigger ranges
            a = np.ones((2,) * i, typecode)
            # Save it on the HDF5 file
            dsetname = "array_" + typecode
            if common.verbose:
                print("Creating dataset:", group._g_join(dsetname))
            self.h5file.create_array(group, dsetname, a, "Large array")
            # Create a new group
            group = self.h5file.create_group(group, "group" + str(i))
            # increment the range for next iteration
            i += 1

        self._reopen()

        # Get the root group
        group = self.h5file.root

        # Get the metadata on the previosly saved arrays
        for i in range(1, len(typecodes)):
            # Create an array for later comparison
            a = np.ones((2,) * i, typecodes[i - 1])
            # Get the dset object hanging from group
            dset = getattr(group, "array_" + typecodes[i - 1])
            # Get the actual array
            b = dset.read()
            if not common.allequal(a, b, "numpy") and common.verbose:
                print("Array a original. Shape: ==>", a.shape)
                print("Array a original. Data: ==>", a)
                print("Info from dataset:", dset._v_pathname)
                print("  shape ==>", dset.shape, end=" ")
                print("  dtype ==> %s" % dset.dtype)
                print("Array b read from file. Shape: ==>", b.shape, end=" ")
                print(". Type ==> %s" % b.dtype.char)

            self.assertEqual(a.shape, b.shape)
            if np.dtype("l").itemsize == 4:
                if a.dtype.char == "i" or a.dtype.char == "l":
                    # Special expection. We have no way to distinguish between
                    # "l" and "i" typecode, and we can consider them the same
                    # to all practical effects
                    self.assertIn(b.dtype.char, ("l", "i"))
                elif a.dtype.char == "I" or a.dtype.char == "L":
                    # Special expection. We have no way to distinguish between
                    # "L" and "I" typecode, and we can consider them the same
                    # to all practical effects
                    self.assertIn(b.dtype.char, ("L", "I"))
                else:
                    self.assertTrue(common.allequal(a, b, "numpy"))
            elif np.dtype("l").itemsize == 8:
                if a.dtype.char == "q" or a.dtype.char == "l":
                    # Special expection. We have no way to distinguish between
                    # "q" and "l" typecode in 64-bit platforms, and we can
                    # consider them the same to all practical effects
                    self.assertIn(b.dtype.char, ("l", "q"))
                elif a.dtype.char == "Q" or a.dtype.char == "L":
                    # Special expection. We have no way to distinguish between
                    # "Q" and "L" typecode in 64-bit platforms, and we can
                    # consider them the same to all practical effects
                    self.assertIn(b.dtype.char, ("L", "Q"))
                else:
                    self.assertTrue(common.allequal(a, b, "numpy"))

            # Iterate over the next group
            group = getattr(group, "group" + str(i))

    def test01_largeRankArrays(self):
        """Checking creation of large rank arrays (0 < rank <= 32)

        It also uses arrays ranks which ranges until maxrank.

        """

        # maximum level of recursivity (deepest group level) achieved:
        # maxrank = 32 (for an effective maximum rank of 32)
        # This limit is due to a limit in the HDF5 library.
        minrank = 1
        maxrank = 32

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test01_largeRankArrays..."
                % self.__class__.__name__
            )
            print("Maximum rank for tested arrays:", maxrank)

        group = self.h5file.root
        if common.verbose:
            print("Rank array writing progress: ", end=" ")
        for rank in range(minrank, maxrank + 1):
            # Create an array of integers, with incrementally bigger ranges
            a = np.ones((1,) * rank, "i")
            if common.verbose:
                print("%3d," % (rank), end=" ")
            self.h5file.create_array(group, "array", a, "Rank: %s" % rank)
            group = self.h5file.create_group(group, "group" + str(rank))

        # Flush the buffers
        self.h5file.flush()

        self._reopen()

        group = self.h5file.root
        if common.verbose:
            print()
            print("Rank array reading progress: ")
        # Get the metadata on the previously saved arrays
        for rank in range(minrank, maxrank + 1):
            # Create an array for later comparison
            a = np.ones((1,) * rank, "i")
            # Get the actual array
            b = group.array.read()
            if common.verbose:
                print("%3d," % (rank), end=" ")
            if not a.tolist() == b.tolist() and common.verbose:
                dset = group.array
                print("Info from dataset:", dset._v_pathname)
                print("  Shape: ==>", dset.shape, end=" ")
                print("  typecode ==> %c" % dset.typecode)
                print("Array b read from file. Shape: ==>", b.shape, end=" ")
                print(". Type ==> %c" % b.dtype.char)
            self.assertEqual(a.shape, b.shape)
            if a.dtype.char == "i":
                # Special expection. We have no way to distinguish between
                # "l" and "i" typecode, and we can consider them the same
                # to all practical effects
                self.assertIn(b.dtype.char, ("l", "i"))
            else:
                self.assertEqual(a.dtype.char, b.dtype.char)

            self.assertEqual(a, b)

            # Iterate over the next group
            group = self.h5file.get_node(group, "group" + str(rank))

        if common.verbose:
            print()  # This flush the stdout buffer


# Test Record class
class Record(tb.IsDescription):
    var1 = tb.StringCol(itemsize=4, dflt=b"abcd", pos=0)
    var2 = tb.StringCol(itemsize=1, dflt=b"a", pos=1)
    var3 = tb.BoolCol(dflt=1)
    var4 = tb.Int8Col(dflt=1)
    var5 = tb.UInt8Col(dflt=1)
    var6 = tb.Int16Col(dflt=1)
    var7 = tb.UInt16Col(dflt=1)
    var8 = tb.Int32Col(dflt=1)
    var9 = tb.UInt32Col(dflt=1)
    var10 = tb.Int64Col(dflt=1)
    var11 = tb.Float32Col(dflt=1.0)
    var12 = tb.Float64Col(dflt=1.0)
    var13 = tb.ComplexCol(itemsize=8, dflt=(1.0 + 0.0j))
    var14 = tb.ComplexCol(itemsize=16, dflt=(1.0 + 0.0j))
    if hasattr(tb, "Float16Col"):
        var15 = tb.Float16Col(dflt=1.0)
    if hasattr(tb, "Float96Col"):
        var16 = tb.Float96Col(dflt=1.0)
    if hasattr(tb, "Float128Col"):
        var17 = tb.Float128Col(dflt=1.0)
    if hasattr(tb, "Complex196Col"):
        var18 = tb.ComplexCol(itemsize=24, dflt=(1.0 + 0.0j))
    if hasattr(tb, "Complex256Col"):
        var19 = tb.ComplexCol(itemsize=32, dflt=(1.0 + 0.0j))


class TableReadTestCase(common.TempFileMixin, common.PyTablesTestCase):
    nrows = 100

    def setUp(self):
        super().setUp()

        # Create an instance of an HDF5 Table
        table = self.h5file.create_table(self.h5file.root, "table", Record)
        for i in range(self.nrows):
            table.row.append()  # Fill 100 rows with default values

        self._reopen(mode="a")

    def test01_readTableChar(self):
        """Checking column conversion into NumPy in read().

        Char flavor

        """

        table = self.h5file.root.table
        table.flavor = "numpy"
        for colname in table.colnames:
            numcol = table.read(field=colname)
            typecol = table.coltypes[colname]
            itemsizecol = table.description._v_dtypes[colname].base.itemsize
            nctypecode = numcol.dtype.char
            if typecol == "string":
                if itemsizecol > 1:
                    orignumcol = np.array(["abcd"] * self.nrows, dtype="S4")
                else:
                    orignumcol = np.array(["a"] * self.nrows, dtype="S1")
                if common.verbose:
                    print("Typecode of NumPy column read:", nctypecode)
                    print("Should look like:", "c")
                    print("Itemsize of column:", itemsizecol)
                    print("Shape of NumPy column read:", numcol.shape)
                    print("Should look like:", orignumcol.shape)
                    print("First 3 elements of read col:", numcol[:3])
                # Check that both NumPy objects are equal
                self.assertTrue(common.allequal(numcol, orignumcol, "numpy"))

    def test01_readTableNum(self):
        """Checking column conversion into NumPy in read().

        NumPy flavor

        """

        table = self.h5file.root.table
        table.flavor = "numpy"
        for colname in table.colnames:
            numcol = table.read(field=colname)
            typecol = table.coltypes[colname]
            nctypecode = np.dtype(numcol.dtype.char[0]).type
            if typecol != "string":
                if common.verbose:
                    print("Typecode of NumPy column read:", nctypecode)
                    print("Should look like:", typecol)
                orignumcol = np.ones(shape=self.nrows, dtype=numcol.dtype.char)
                # Check that both NumPy objects are equal
                self.assertTrue(common.allequal(numcol, orignumcol, "numpy"))

    def test02_readCoordsChar(self):
        """Column conversion into NumPy in readCoords().

        Chars

        """

        table = self.h5file.root.table
        table.flavor = "numpy"
        coords = [1, 2, 3]
        self.nrows = len(coords)
        for colname in table.colnames:
            numcol = table.read_coordinates(coords, field=colname)
            typecol = table.coltypes[colname]
            itemsizecol = table.description._v_dtypes[colname].base.itemsize
            nctypecode = numcol.dtype.char
            if typecol == "string":
                if itemsizecol > 1:
                    orignumcol = np.array(["abcd"] * self.nrows, dtype="S4")
                else:
                    orignumcol = np.array(["a"] * self.nrows, dtype="S1")
                if common.verbose:
                    print("Typecode of NumPy column read:", nctypecode)
                    print("Should look like:", "c")
                    print("Itemsize of column:", itemsizecol)
                    print("Shape of NumPy column read:", numcol.shape)
                    print("Should look like:", orignumcol.shape)
                    print("First 3 elements of read col:", numcol[:3])
                # Check that both NumPy objects are equal
                self.assertTrue(common.allequal(numcol, orignumcol, "numpy"))

    def test02_readCoordsNum(self):
        """Column conversion into NumPy in read_coordinates().

        NumPy.

        """

        table = self.h5file.root.table
        table.flavor = "numpy"
        coords = [1, 2, 3]
        self.nrows = len(coords)
        for colname in table.colnames:
            numcol = table.read_coordinates(coords, field=colname)
            typecol = table.coltypes[colname]
            type_ = numcol.dtype.type
            if typecol != "string":
                if typecol == "int64":
                    return
                if common.verbose:
                    print("Type of read NumPy column:", type_)
                    print("Should look like:", typecol)
                orignumcol = np.ones(shape=self.nrows, dtype=numcol.dtype.char)
                # Check that both NumPy objects are equal
                self.assertTrue(common.allequal(numcol, orignumcol, "numpy"))

    def test03_getIndexNumPy(self):
        """Getting table rows specified as NumPy scalar integers."""

        table = self.h5file.root.table
        coords = np.array([1, 2, 3], dtype="int8")
        for colname in table.colnames:
            numcol = [table[coord][colname] for coord in coords]
            typecol = table.coltypes[colname]
            if typecol != "string":
                if typecol == "int64":
                    return
                numcol = np.array(numcol, typecol)
                if common.verbose:
                    type_ = numcol.dtype.type
                    print("Type of read NumPy column:", type_)
                    print("Should look like:", typecol)
                orignumcol = np.ones(
                    shape=len(numcol), dtype=numcol.dtype.char
                )
                # Check that both NumPy objects are equal
                self.assertTrue(common.allequal(numcol, orignumcol, "numpy"))

    def test04_setIndexNumPy(self):
        """Setting table rows specified as NumPy integers."""

        self._reopen(mode="a")
        table = self.h5file.root.table
        table.flavor = "numpy"
        coords = np.array([1, 2, 3], dtype="int8")
        # Modify row 1
        # From PyTables 2.0 on, assignments to records can be done
        # only as tuples (see http://projects.scipy.org/scipy/numpy/ticket/315)
        # table[coords[0]] = ["aasa","x"]+[123]*12

        n = len(Record.columns) - 2

        table[coords[0]] = tuple(["aasa", "x"] + [123] * n)  # XXX
        # record = list(table[coords[0]])
        record = table.read(coords[0], coords[0] + 1)
        if common.verbose:
            print(
                "Original row:\n"
                "['aasa', 'x', True, 123, 123, 123, 123, 123, 123L, "
                "123, 123.0, 123.0, (123 + 0j), (123+0j), 123.0, "
                "(123+0j)]\n"
            )
            print("Read row:\n", record)
        self.assertEqual(record["var1"], b"aasa")
        self.assertEqual(record["var2"], b"x")
        self.assertEqual(record["var3"], True)
        self.assertEqual(record["var4"], 123)
        self.assertEqual(record["var7"], 123)


# The declaration of the nested table:
class Info(tb.IsDescription):
    _v_pos = 3
    Name = tb.StringCol(itemsize=2)
    Value = tb.ComplexCol(itemsize=16)


class TestTDescr(tb.IsDescription):
    """A description that has several nested columns."""

    x = tb.Int32Col(dflt=0, shape=2, pos=0)  # 0
    y = tb.FloatCol(dflt=1, shape=(2, 2))
    z = tb.UInt8Col(dflt=1)
    z3 = tb.EnumCol({"r": 4, "g": 2, "b": 1}, "r", "int32", shape=2)
    color = tb.StringCol(itemsize=4, dflt=b"ab", pos=2)
    info = Info()

    class Info(tb.IsDescription):  # 1
        _v_pos = 1
        name = tb.StringCol(itemsize=2)
        value = tb.ComplexCol(itemsize=16, pos=0)  # 0
        y2 = tb.FloatCol(pos=1)  # 1
        z2 = tb.UInt8Col()

        class Info2(tb.IsDescription):
            y3 = tb.Time64Col(shape=2)
            name = tb.StringCol(itemsize=2)
            value = tb.ComplexCol(itemsize=16, shape=2)


class TableNativeFlavorTestCase(common.TempFileMixin, common.PyTablesTestCase):
    nrows = 100

    def setUp(self):
        super().setUp()

        # Create an instance of an HDF5 Table
        table = self.h5file.create_table(
            self.h5file.root, "table", TestTDescr, expectedrows=self.nrows
        )
        table.flavor = "numpy"
        for i in range(self.nrows):
            table.row.append()  # Fill 100 rows with default values
        table.flush()

    def test01a_basicTableRead(self):
        """Checking the return of a NumPy in read()."""

        if self.close:
            self._reopen(mode="a")
        table = self.h5file.root.table
        data = table[:]
        if common.verbose:
            print("Type of read:", type(data))
            print("Description of the record:", data.dtype.descr)
            print("First 3 elements of read:", data[:3])

        # Check that both NumPy objects are equal
        self.assertIsInstance(data, np.ndarray)

        # Check the value of some columns
        # A flat column
        col = table.cols.x[:3]
        self.assertIsInstance(col, np.ndarray)
        npcol = np.zeros((3, 2), dtype="int32")
        self.assertTrue(common.allequal(col, npcol, "numpy"))

        # A nested column
        col = table.cols.Info[:3]
        self.assertIsInstance(col, np.ndarray)
        dtype = [
            ("value", "c16"),
            ("y2", "f8"),
            (
                "Info2",
                [("name", "S2"), ("value", "c16", (2,)), ("y3", "f8", (2,))],
            ),
            ("name", "S2"),
            ("z2", "u1"),
        ]
        npcol = np.zeros((3,), dtype=dtype)
        self.assertEqual(col.dtype.descr, npcol.dtype.descr)
        if common.verbose:
            print("col-->", col)
            print("npcol-->", npcol)

        # A copy() is needed in case the buffer can be in different segments
        self.assertEqual(bytes(col.copy().data), bytes(npcol.data))

    def test01b_basicTableRead(self):
        """Checking the return of a NumPy in read() (strided version)."""

        if self.close:
            self._reopen(mode="a")
        table = self.h5file.root.table
        data = table[::3]
        if common.verbose:
            print("Type of read:", type(data))
            print("Description of the record:", data.dtype.descr)
            print("First 3 elements of read:", data[:3])

        # Check that both NumPy objects are equal
        self.assertIsInstance(data, np.ndarray)

        # Check the value of some columns
        # A flat column
        col = table.cols.x[:9:3]
        self.assertIsInstance(col, np.ndarray)
        npcol = np.zeros((3, 2), dtype="int32")
        self.assertTrue(common.allequal(col, npcol, "numpy"))

        # A nested column
        col = table.cols.Info[:9:3]
        self.assertIsInstance(col, np.ndarray)
        dtype = [
            ("value", "%sc16" % byteorder),
            ("y2", "%sf8" % byteorder),
            (
                "Info2",
                [
                    ("name", "|S2"),
                    ("value", "%sc16" % byteorder, (2,)),
                    ("y3", "%sf8" % byteorder, (2,)),
                ],
            ),
            ("name", "|S2"),
            ("z2", "|u1"),
        ]
        npcol = np.zeros((3,), dtype=dtype)
        self.assertEqual(col.dtype.descr, npcol.dtype.descr)
        if common.verbose:
            print("col-->", col)
            print("npcol-->", npcol)

        # A copy() is needed in case the buffer can be in different segments
        self.assertEqual(bytes(col.copy().data), bytes(npcol.data))

    def test02_getWhereList(self):
        """Checking the return of NumPy in get_where_list method."""

        if self.close:
            self._reopen(mode="a")
        table = self.h5file.root.table
        data = table.get_where_list("z == 1")
        if common.verbose:
            print("Type of read:", type(data))
            print("Description of the record:", data.dtype.descr)
            print("First 3 elements of read:", data[:3])

        # Check that both NumPy objects are equal
        self.assertIsInstance(data, np.ndarray)

        # Check that all columns have been selected
        self.assertEqual(len(data), 100)

        # Finally, check that the contents are ok
        self.assertTrue(
            common.allequal(data, np.arange(100, dtype="i8"), "numpy")
        )

    def test03a_readWhere(self):
        """Checking the return of NumPy in read_where method (strings)."""

        table = self.h5file.root.table
        table.cols.color.create_index()
        if self.close:
            self._reopen(mode="a")
            table = self.h5file.root.table
        data = table.read_where('color == b"ab"')
        if common.verbose:
            print("Type of read:", type(data))
            print("Length of the data read:", len(data))

        # Check that both NumPy objects are equal
        self.assertIsInstance(data, np.ndarray)

        # Check that all columns have been selected
        self.assertEqual(len(data), self.nrows)

    def test03b_readWhere(self):
        """Checking the return of NumPy in read_where method (numeric)."""

        table = self.h5file.root.table
        table.cols.z.create_index()
        if self.close:
            self._reopen(mode="a")
            table = self.h5file.root.table
        data = table.read_where("z == 0")
        if common.verbose:
            print("Type of read:", type(data))
            print("Length of the data read:", len(data))

        # Check that both NumPy objects are equal
        self.assertIsInstance(data, np.ndarray)

        # Check that all columns have been selected
        self.assertEqual(len(data), 0)

    def test04a_createTable(self):
        """Checking the Table creation from a numpy recarray."""

        dtype = [
            ("value", "%sc16" % byteorder),
            ("y2", "%sf8" % byteorder),
            (
                "Info2",
                [
                    ("name", "|S2"),
                    ("value", "%sc16" % byteorder, (2,)),
                    ("y3", "%sf8" % byteorder, (2,)),
                ],
            ),
            ("name", "|S2"),
            ("z2", "|u1"),
        ]
        npdata = np.zeros((3,), dtype=dtype)
        table = self.h5file.create_table(self.h5file.root, "table2", npdata)
        if self.close:
            self._reopen(mode="a")
            table = self.h5file.root.table2
        data = table[:]
        if common.verbose:
            print("Type of read:", type(data))
            print("Description of the record:", data.dtype.descr)
            print("First 3 elements of read:", data[:3])
            print("Length of the data read:", len(data))

        # Check that both NumPy objects are equal
        self.assertIsInstance(data, np.ndarray)

        # Check the type
        self.assertEqual(data.dtype.descr, npdata.dtype.descr)
        if common.verbose:
            print("npdata-->", npdata)
            print("data-->", data)

        # A copy() is needed in case the buffer would be in different segments
        self.assertEqual(bytes(data.copy().data), bytes(npdata.data))

    def test04b_appendTable(self):
        """Checking appending a numpy recarray."""

        table = self.h5file.root.table
        npdata = table[3:6]
        table.append(npdata)
        if self.close:
            self._reopen(mode="a")
            table = self.h5file.root.table
        data = table[-3:]
        if common.verbose:
            print("Type of read:", type(data))
            print("Description of the record:", data.dtype.descr)
            print("Last 3 elements of read:", data[-3:])
            print("Length of the data read:", len(data))

        # Check that both NumPy objects are equal
        self.assertIsInstance(data, np.ndarray)

        # Check the type
        self.assertEqual(data.dtype.descr, npdata.dtype.descr)
        if common.verbose:
            print("npdata-->", npdata)
            print("data-->", data)

        # A copy() is needed in case the buffer would be in different segments
        self.assertEqual(bytes(data.copy().data), bytes(npdata.data))

    def test05a_assignColumn(self):
        """Checking assigning to a column."""

        table = self.h5file.root.table
        table.cols.z[:] = np.zeros((100,), dtype="u1")
        if self.close:
            self._reopen(mode="a")
            table = self.h5file.root.table
        data = table.cols.z[:]
        if common.verbose:
            print("Type of read:", type(data))
            print("Description of the record:", data.dtype.descr)
            print("First 3 elements of read:", data[:3])
            print("Length of the data read:", len(data))

        # Check that both NumPy objects are equal
        self.assertIsInstance(data, np.ndarray)

        # Check that all columns have been selected
        self.assertEqual(len(data), 100)

        # Finally, check that the contents are ok
        self.assertTrue(
            common.allequal(data, np.zeros((100,), dtype="u1"), "numpy")
        )

    def test05b_modifyingColumns(self):
        """Checking modifying several columns at once."""

        table = self.h5file.root.table
        xcol = np.ones((3, 2), "int32")
        ycol = np.zeros((3, 2, 2), "float64")
        zcol = np.zeros((3,), "uint8")
        table.modify_columns(3, 6, 1, [xcol, ycol, zcol], ["x", "y", "z"])
        if self.close:
            self._reopen(mode="a")
            table = self.h5file.root.table
        data = table.cols.y[3:6]
        if common.verbose:
            print("Type of read:", type(data))
            print("Description of the record:", data.dtype.descr)
            print("First 3 elements of read:", data[:3])
            print("Length of the data read:", len(data))

        # Check that both NumPy objects are equal
        self.assertIsInstance(data, np.ndarray)

        # Check the type
        self.assertEqual(data.dtype.descr, ycol.dtype.descr)
        if common.verbose:
            print("ycol-->", ycol)
            print("data-->", data)

        # A copy() is needed in case the buffer would be in different segments
        self.assertEqual(data.copy().data, ycol.data)

    def test05c_modifyingColumns(self):
        """Checking modifying several columns using a single numpy buffer."""

        table = self.h5file.root.table
        dtype = [("x", "i4", (2,)), ("y", "f8", (2, 2)), ("z", "u1")]
        nparray = np.zeros((3,), dtype=dtype)
        table.modify_columns(3, 6, 1, nparray, ["x", "y", "z"])
        if self.close:
            self._reopen(mode="a")
            table = self.h5file.root.table
        ycol = np.zeros((3, 2, 2), "float64")
        data = table.cols.y[3:6]
        if common.verbose:
            print("Type of read:", type(data))
            print("Description of the record:", data.dtype.descr)
            print("First 3 elements of read:", data[:3])
            print("Length of the data read:", len(data))

        # Check that both NumPy objects are equal
        self.assertIsInstance(data, np.ndarray)

        # Check the type
        self.assertEqual(data.dtype.descr, ycol.dtype.descr)
        if common.verbose:
            print("ycol-->", ycol)
            print("data-->", data)

        # A copy() is needed in case the buffer would be in different segments
        self.assertEqual(data.copy().data, ycol.data)

    def test06a_assignNestedColumn(self):
        """Checking assigning a nested column (using modify_column)."""

        table = self.h5file.root.table
        dtype = [
            ("value", "%sc16" % byteorder),
            ("y2", "%sf8" % byteorder),
            (
                "Info2",
                [
                    ("name", "|S2"),
                    ("value", "%sc16" % byteorder, (2,)),
                    ("y3", "%sf8" % byteorder, (2,)),
                ],
            ),
            ("name", "|S2"),
            ("z2", "|u1"),
        ]
        npdata = np.zeros((3,), dtype=dtype)
        data = table.cols.Info[3:6]
        table.modify_column(3, 6, 1, column=npdata, colname="Info")
        if self.close:
            self._reopen(mode="a")
            table = self.h5file.root.table
        data = table.cols.Info[3:6]
        if common.verbose:
            print("Type of read:", type(data))
            print("Description of the record:", data.dtype.descr)
            print("First 3 elements of read:", data[:3])
            print("Length of the data read:", len(data))

        # Check that both NumPy objects are equal
        self.assertIsInstance(data, np.ndarray)

        # Check the type
        self.assertEqual(data.dtype.descr, npdata.dtype.descr)
        if common.verbose:
            print("npdata-->", npdata)
            print("data-->", data)

        # A copy() is needed in case the buffer would be in different segments
        self.assertEqual(bytes(data.copy().data), bytes(npdata.data))

    def test06b_assignNestedColumn(self):
        """Checking assigning a nested column (using the .cols accessor)."""

        table = self.h5file.root.table
        dtype = [
            ("value", "%sc16" % byteorder),
            ("y2", "%sf8" % byteorder),
            (
                "Info2",
                [
                    ("name", "|S2"),
                    ("value", "%sc16" % byteorder, (2,)),
                    ("y3", "%sf8" % byteorder, (2,)),
                ],
            ),
            ("name", "|S2"),
            ("z2", "|u1"),
        ]
        npdata = np.zeros((3,), dtype=dtype)
        # self.assertRaises(NotImplementedError,
        #                   table.cols.Info.__setitem__, slice(3,6,1),  npdata)
        table.cols.Info[3:6] = npdata
        if self.close:
            self._reopen(mode="a")
            table = self.h5file.root.table
        data = table.cols.Info[3:6]
        if common.verbose:
            print("Type of read:", type(data))
            print("Description of the record:", data.dtype.descr)
            print("First 3 elements of read:", data[:3])
            print("Length of the data read:", len(data))

        # Check that both NumPy objects are equal
        self.assertIsInstance(data, np.ndarray)

        # Check the type
        self.assertEqual(data.dtype.descr, npdata.dtype.descr)
        if common.verbose:
            print("npdata-->", npdata)
            print("data-->", data)

        # A copy() is needed in case the buffer would be in different segments
        self.assertEqual(bytes(data.copy().data), bytes(npdata.data))

    def test07a_modifyingRows(self):
        """Checking modifying several rows at once (using modify_rows)."""

        table = self.h5file.root.table

        # Read a chunk of the table
        chunk = table[0:3]

        # Modify it somewhat
        chunk["y"][:] = -1
        table.modify_rows(3, 6, 1, rows=chunk)
        if self.close:
            self._reopen(mode="a")
            table = self.h5file.root.table
        ycol = np.zeros((3, 2, 2), "float64") - 1
        data = table.cols.y[3:6]
        if common.verbose:
            print("Type of read:", type(data))
            print("Description of the record:", data.dtype.descr)
            print("First 3 elements of read:", data[:3])
            print("Length of the data read:", len(data))

        # Check that both NumPy objects are equal
        self.assertIsInstance(data, np.ndarray)

        # Check the type
        self.assertEqual(data.dtype.descr, ycol.dtype.descr)
        if common.verbose:
            print("ycol-->", ycol)
            print("data-->", data)
        self.assertTrue(common.allequal(ycol, data, "numpy"))

    def test07b_modifyingRows(self):
        """Checking modifying several rows at once (using cols accessor)."""

        table = self.h5file.root.table

        # Read a chunk of the table
        chunk = table[0:3]

        # Modify it somewhat
        chunk["y"][:] = -1
        table.cols[3:6] = chunk
        if self.close:
            self._reopen(mode="a")
            table = self.h5file.root.table

        # Check that some column has been actually modified
        ycol = np.zeros((3, 2, 2), "float64") - 1
        data = table.cols.y[3:6]
        if common.verbose:
            print("Type of read:", type(data))
            print("Description of the record:", data.dtype.descr)
            print("First 3 elements of read:", data[:3])
            print("Length of the data read:", len(data))

        # Check that both NumPy objects are equal
        self.assertIsInstance(data, np.ndarray)

        # Check the type
        self.assertEqual(data.dtype.descr, ycol.dtype.descr)
        if common.verbose:
            print("ycol-->", ycol)
            print("data-->", data)
        self.assertTrue(common.allequal(ycol, data, "numpy"))

    def test08a_modifyingRows(self):
        """Checking modifying just one row at once (using modify_rows)."""

        table = self.h5file.root.table

        # Read a chunk of the table
        chunk = table[3:4]

        # Modify it somewhat
        chunk["y"][:] = -1
        table.modify_rows(6, 7, 1, chunk)
        if self.close:
            self._reopen(mode="a")
            table = self.h5file.root.table

        # Check that some column has been actually modified
        ycol = np.zeros((2, 2), "float64") - 1
        data = table.cols.y[6]
        if common.verbose:
            print("Type of read:", type(data))
            print("Description of the record:", data.dtype.descr)
            print("First 3 elements of read:", data[:3])
            print("Length of the data read:", len(data))

        # Check that both NumPy objects are equal
        self.assertIsInstance(data, np.ndarray)

        # Check the type
        self.assertEqual(data.dtype.descr, ycol.dtype.descr)
        if common.verbose:
            print("ycol-->", ycol)
            print("data-->", data)
        self.assertTrue(common.allequal(ycol, data, "numpy"))

    def test08b_modifyingRows(self):
        """Checking modifying just one row at once (using cols accessor)."""

        table = self.h5file.root.table

        # Read a chunk of the table
        chunk = table[3:4]

        # Modify it somewhat
        chunk["y"][:] = -1
        table.cols[6] = chunk
        if self.close:
            self._reopen(mode="a")
            table = self.h5file.root.table

        # Check that some column has been actually modified
        ycol = np.zeros((2, 2), "float64") - 1
        data = table.cols.y[6]
        if common.verbose:
            print("Type of read:", type(data))
            print("Description of the record:", data.dtype.descr)
            print("First 3 elements of read:", data[:3])
            print("Length of the data read:", len(data))

        # Check that both NumPy objects are equal
        self.assertIsInstance(data, np.ndarray)

        # Check the type
        self.assertEqual(data.dtype.descr, ycol.dtype.descr)
        if common.verbose:
            print("ycol-->", ycol)
            print("data-->", data)
        self.assertTrue(common.allequal(ycol, data, "numpy"))

    def test09a_getStrings(self):
        """Checking the return of string columns with spaces."""

        if self.close:
            self._reopen(mode="a")
        table = self.h5file.root.table
        rdata = table.get_where_list('color == b"ab"')
        data = table.read_coordinates(rdata)
        if common.verbose:
            print("Type of read:", type(data))
            print("Description of the record:", data.dtype.descr)
            print("First 3 elements of read:", data[:3])

        # Check that both NumPy objects are equal
        self.assertIsInstance(data, np.ndarray)

        # Check that all columns have been selected
        self.assertEqual(len(data), 100)

        # Finally, check that the contents are ok
        for idata in data["color"]:
            self.assertEqual(idata, np.array("ab", dtype="|S4"))

    def test09b_getStrings(self):
        """Checking the return of string columns with spaces.

        (modify)

        """

        if self.close:
            self._reopen(mode="a")
        table = self.h5file.root.table
        for i in range(50):
            table.cols.color[i] = "a  "
        table.flush()
        data = table[:]
        if common.verbose:
            print("Type of read:", type(data))
            print("Description of the record:", data.dtype.descr)
            print("First 3 elements of read:", data[:3])

        # Check that both NumPy objects are equal
        self.assertIsInstance(data, np.ndarray)

        # Check that all columns have been selected
        self.assertEqual(len(data), 100)

        # Finally, check that the contents are ok
        for i in range(100):
            idata = data["color"][i]
            if i >= 50:
                self.assertEqual(idata, np.array("ab", dtype="|S4"))
            else:
                self.assertEqual(idata, np.array("a  ", dtype="|S4"))

    def test09c_getStrings(self):
        """Checking the return of string columns with spaces.

        (append)

        """

        if self.close:
            self._reopen(mode="a")
        table = self.h5file.root.table
        row = table.row
        for i in range(50):
            row["color"] = "a  "  # note the trailing spaces
            row.append()
        table.flush()
        if self.close:
            self.h5file.close()
            self.h5file = tb.open_file(self.h5fname, "a")
        data = self.h5file.root.table[:]
        if common.verbose:
            print("Type of read:", type(data))
            print("Description of the record:", data.dtype.descr)
            print("First 3 elements of read:", data[:3])

        # Check that both NumPy objects are equal
        self.assertIsInstance(data, np.ndarray)

        # Check that all columns have been selected
        self.assertEqual(len(data), 150)

        # Finally, check that the contents are ok
        for i in range(150):
            idata = data["color"][i]
            if i < 100:
                self.assertEqual(idata, np.array("ab", dtype="|S4"))
            else:
                self.assertEqual(idata, np.array("a  ", dtype="|S4"))


class TableNativeFlavorOpenTestCase(TableNativeFlavorTestCase):
    close = False


class TableNativeFlavorCloseTestCase(TableNativeFlavorTestCase):
    close = True


class AttributesTestCase(common.TempFileMixin, common.PyTablesTestCase):
    def setUp(self):
        super().setUp()

        # Create an instance of an HDF5 Table
        self.h5file.create_group(self.h5file.root, "group")

    def test01_writeAttribute(self):
        """Checking the creation of a numpy attribute."""

        group = self.h5file.root.group
        g_attrs = group._v_attrs
        g_attrs.numpy1 = np.zeros((1, 1), dtype="int16")
        if self.close:
            self._reopen(mode="a")
            group = self.h5file.root.group
            g_attrs = group._v_attrs

        # Check that we can retrieve a numpy object
        data = g_attrs.numpy1
        npcomp = np.zeros((1, 1), dtype="int16")

        # Check that both NumPy objects are equal
        self.assertIsInstance(data, np.ndarray)

        # Check the type
        self.assertEqual(data.dtype.descr, npcomp.dtype.descr)
        if common.verbose:
            print("npcomp-->", npcomp)
            print("data-->", data)
        self.assertTrue(common.allequal(npcomp, data, "numpy"))

    def test02_updateAttribute(self):
        """Checking the modification of a numpy attribute."""

        group = self.h5file.root.group
        g_attrs = group._v_attrs
        g_attrs.numpy1 = np.zeros((1, 2), dtype="int16")
        if self.close:
            self._reopen(mode="a")
            group = self.h5file.root.group
            g_attrs = group._v_attrs

        # Update this attribute
        g_attrs.numpy1 = np.ones((1, 2), dtype="int16")

        # Check that we can retrieve a numpy object
        data = g_attrs.numpy1
        npcomp = np.ones((1, 2), dtype="int16")

        # Check that both NumPy objects are equal
        self.assertIsInstance(data, np.ndarray)

        # Check the type
        self.assertEqual(data.dtype.descr, npcomp.dtype.descr)
        if common.verbose:
            print("npcomp-->", npcomp)
            print("data-->", data)
        self.assertTrue(common.allequal(npcomp, data, "numpy"))


class AttributesOpenTestCase(AttributesTestCase):
    close = 0


class AttributesCloseTestCase(AttributesTestCase):
    close = 1


class StrlenTestCase(common.TempFileMixin, common.PyTablesTestCase):
    def setUp(self):
        super().setUp()

        # Create an instance of an HDF5 Table
        group = self.h5file.create_group(self.h5file.root, "group")
        tablelayout = {
            "Text": tb.StringCol(itemsize=1000),
        }
        self.table = self.h5file.create_table(group, "table", tablelayout)
        self.table.flavor = "numpy"
        row = self.table.row
        row["Text"] = "Hello Francesc!"  # XXX: check unicode --> bytes
        row.append()
        row["Text"] = "Hola Francesc!"  # XXX: check unicode --> bytes
        row.append()
        self.table.flush()

    def test01(self):
        """Checking the lengths of strings (read field)."""

        if self.close:
            self._reopen(mode="a")
            self.table = self.h5file.root.group.table

        # Get both strings
        str1 = self.table.col("Text")[0]
        str2 = self.table.col("Text")[1]
        if common.verbose:
            print("string1-->", str1)
            print("string2-->", str2)

        # Check that both NumPy objects are equal
        self.assertEqual(len(str1), len(b"Hello Francesc!"))
        self.assertEqual(len(str2), len(b"Hola Francesc!"))
        self.assertEqual(str1, b"Hello Francesc!")
        self.assertEqual(str2, b"Hola Francesc!")

    def test02(self):
        """Checking the lengths of strings (read recarray)."""

        if self.close:
            self._reopen(mode="a")
            self.table = self.h5file.root.group.table

        # Get both strings
        str1 = self.table[:]["Text"][0]
        str2 = self.table[:]["Text"][1]

        # Check that both NumPy objects are equal
        self.assertEqual(len(str1), len(b"Hello Francesc!"))
        self.assertEqual(len(str2), len(b"Hola Francesc!"))
        self.assertEqual(str1, b"Hello Francesc!")
        self.assertEqual(str2, b"Hola Francesc!")

    def test03(self):
        """Checking the lengths of strings (read recarray, row by row)."""

        if self.close:
            self._reopen(mode="a")
            self.table = self.h5file.root.group.table

        # Get both strings
        str1 = self.table[0]["Text"]
        str2 = self.table[1]["Text"]

        # Check that both NumPy objects are equal
        self.assertEqual(len(str1), len(b"Hello Francesc!"))
        self.assertEqual(len(str2), len(b"Hola Francesc!"))
        self.assertEqual(str1, b"Hello Francesc!")
        self.assertEqual(str2, b"Hola Francesc!")


class StrlenOpenTestCase(StrlenTestCase):
    close = 0


class StrlenCloseTestCase(StrlenTestCase):
    close = 1


def suite():
    theSuite = common.unittest.TestSuite()
    niter = 1

    # theSuite.addTest(make_suite(StrlenOpenTestCase))
    # theSuite.addTest(make_suite(Basic0DOneTestCase))
    # theSuite.addTest(make_suite(GroupsArrayTestCase))
    for i in range(niter):
        theSuite.addTest(common.make_suite(Basic0DOneTestCase))
        theSuite.addTest(common.make_suite(Basic0DTwoTestCase))
        theSuite.addTest(common.make_suite(Basic1DOneTestCase))
        theSuite.addTest(common.make_suite(Basic1DTwoTestCase))
        theSuite.addTest(common.make_suite(Basic1DThreeTestCase))
        theSuite.addTest(common.make_suite(Basic2DTestCase))
        theSuite.addTest(common.make_suite(GroupsArrayTestCase))
        theSuite.addTest(common.make_suite(TableReadTestCase))
        theSuite.addTest(common.make_suite(TableNativeFlavorOpenTestCase))
        theSuite.addTest(common.make_suite(TableNativeFlavorCloseTestCase))
        theSuite.addTest(common.make_suite(AttributesOpenTestCase))
        theSuite.addTest(common.make_suite(AttributesCloseTestCase))
        theSuite.addTest(common.make_suite(StrlenOpenTestCase))
        theSuite.addTest(common.make_suite(StrlenCloseTestCase))
        if common.heavy:
            theSuite.addTest(common.make_suite(Basic10DTestCase))
            # The 32 dimensions case takes forever to run!!
            # theSuite.addTest(make_suite(Basic32DTestCase))
    return theSuite


if __name__ == "__main__":
    common.parse_argv(sys.argv)
    common.print_versions()
    common.unittest.main(defaultTest="suite")
