import sys

import numpy as np

import tables as tb
from tables.tests import common

# It is important that columns are ordered according to their names
# to ease the comparison with structured arrays.


# Test Record class
class Record(tb.IsDescription):
    var0 = tb.StringCol(itemsize=4, dflt=b"", shape=2)  # 4-char str array
    var1 = tb.StringCol(itemsize=4, dflt=[b"abcd", b"efgh"], shape=(2, 2))
    var1_ = tb.IntCol(dflt=((1, 1),), shape=2)  # integer array
    var2 = tb.IntCol(dflt=((1, 1), (1, 1)), shape=(2, 2))  # integer array
    var3 = tb.Int16Col(dflt=2)  # short integer
    var4 = tb.FloatCol(dflt=3.1)  # double (double-precision)
    var5 = tb.Float32Col(dflt=4.2)  # float  (single-precision)
    var6 = tb.UInt16Col(dflt=5)  # unsigned short integer
    var7 = tb.StringCol(itemsize=1, dflt=b"e")  # 1-character String


#  Dictionary definition
RecordDescriptionDict = {
    "var0": tb.StringCol(itemsize=4, dflt=b"", shape=2),  # 4-char str array
    "var1": tb.StringCol(itemsize=4, dflt=[b"abcd", b"efgh"], shape=(2, 2)),
    # 'var0': StringCol(itemsize=4, shape=2),  # 4-character String
    # 'var1': StringCol(itemsize=4, shape=(2,2)),  # 4-character String
    "var1_": tb.IntCol(shape=2),  # integer array
    "var2": tb.IntCol(shape=(2, 2)),  # integer array
    "var3": tb.Int16Col(),  # short integer
    "var4": tb.FloatCol(),  # double (double-precision)
    "var5": tb.Float32Col(),  # float  (single-precision)
    "var6": tb.Int16Col(),  # unsigned short integer
    "var7": tb.StringCol(itemsize=1),  # 1-character String
}

# Record class with numpy dtypes (mixed shapes is checked here)


class RecordDT(tb.IsDescription):
    var0 = tb.Col.from_dtype(np.dtype("2S4"), dflt=b"")  # shape in dtype
    var1 = tb.Col.from_dtype(
        np.dtype(("S4", (2, 2))), dflt=[b"abcd", b"efgh"]
    )  # shape is a mix
    var1_ = tb.Col.from_dtype(
        np.dtype("2i4"), dflt=((1, 1),)
    )  # shape in dtype
    var2 = tb.Col.from_sctype(
        "i4", shape=(2, 2), dflt=((1, 1), (1, 1))
    )  # shape is a mix
    var3 = tb.Col.from_dtype(np.dtype("i2"), dflt=2)
    var4 = tb.Col.from_dtype(np.dtype("2f8"), dflt=3.1)
    var5 = tb.Col.from_dtype(np.dtype("f4"), dflt=4.2)
    var6 = tb.Col.from_dtype(np.dtype("()u2"), dflt=5)
    var7 = tb.Col.from_dtype(np.dtype("S1"), dflt=b"e")  # no shape


class BasicTestCase(common.TempFileMixin, common.PyTablesTestCase):
    # file  = "test.h5"
    open_mode = "w"
    title = "This is the table title"
    expectedrows = 100
    appendrows = 20
    compress = 0
    complib = "zlib"  # Default compression library
    record = Record
    recarrayinit = 0
    maxshort = 1 << 15

    def setUp(self):
        super().setUp()

        # Create an instance of an HDF5 Table
        self.rootgroup = self.h5file.root
        self.populateFile()
        self.h5file.close()

    def initRecArray(self):
        record = self.recordtemplate
        row = record[0]
        buflist = []
        # Fill the recarray
        for i in range(self.expectedrows):
            tmplist = []
            # Both forms (list or chararray) work
            var0 = ["%04d" % (self.expectedrows - i)] * 2
            tmplist.append(var0)
            var1 = [["%04d" % (self.expectedrows - i)] * 2] * 2
            tmplist.append(var1)
            var1_ = (i, 1)
            tmplist.append(var1_)
            var2 = ((i, 1), (1, 1))  # *-*
            tmplist.append(var2)
            var3 = i % self.maxshort
            tmplist.append(var3)
            if isinstance(row["var4"], np.ndarray):
                tmplist.append([float(i), float(i * i)])
            else:
                tmplist.append(float(i))
            if isinstance(row["var5"], np.ndarray):
                tmplist.append(np.array((float(i),) * 4))
            else:
                tmplist.append(float(i))
            # var6 will be like var3 but byteswaped
            tmplist.append(((var3 >> 8) & 0xFF) + ((var3 << 8) & 0xFF00))
            var7 = var1[0][0][-1]
            tmplist.append(var7)
            buflist.append(tuple(tmplist))

        self.record = np.rec.array(
            buflist, dtype=record.dtype, shape=self.expectedrows
        )

    def populateFile(self):
        group = self.rootgroup
        if self.recarrayinit:
            # Initialize a starting buffer, if any
            self.initRecArray()
        for j in range(3):
            # Create a table
            filters = tb.Filters(complevel=self.compress, complib=self.complib)
            if j < 2:
                byteorder = sys.byteorder
            else:
                # table2 will be byteswapped
                byteorder = {"little": "big", "big": "little"}[sys.byteorder]
            table = self.h5file.create_table(
                group,
                "table" + str(j),
                self.record,
                title=self.title,
                filters=filters,
                expectedrows=self.expectedrows,
                byteorder=byteorder,
            )
            if not self.recarrayinit:
                # Get the row object associated with the new table
                row = table.row

                # Fill the table
                for i in range(self.expectedrows):
                    s = "%04d" % (self.expectedrows - i)
                    row["var0"] = s.encode("ascii")
                    row["var1"] = s.encode("ascii")
                    row["var7"] = s[-1].encode("ascii")
                    row["var1_"] = (i, 1)
                    row["var2"] = ((i, 1), (1, 1))  # *-*
                    row["var3"] = i % self.maxshort
                    if isinstance(row["var4"], np.ndarray):
                        row["var4"] = [float(i), float(i * i)]
                    else:
                        row["var4"] = float(i)
                    if isinstance(row["var5"], np.ndarray):
                        row["var5"] = np.array((float(i),) * 4)
                    else:
                        row["var5"] = float(i)
                    # var6 will be like var3 but byteswaped
                    row["var6"] = ((row["var3"] >> 8) & 0xFF) + (
                        (row["var3"] << 8) & 0xFF00
                    )
                    row.append()

            # Flush the buffer for this table
            table.flush()
            # Create a new group (descendant of group)
            group2 = self.h5file.create_group(group, "group" + str(j))
            # Iterate over this new group (group2)
            group = group2

    def test00_description(self):
        """Checking table description and descriptive fields."""

        self.h5file = tb.open_file(self.h5fname)

        tbl = self.h5file.get_node("/table0")
        desc = tbl.description

        if isinstance(self.record, dict):
            columns = self.record
        elif isinstance(self.record, np.ndarray):
            descr, _ = tb.description.descr_from_dtype(self.record.dtype)
            columns = descr._v_colobjects
        elif isinstance(self.record, np.dtype):
            descr, _ = tb.description.descr_from_dtype(self.record)
            columns = descr._v_colobjects
        else:
            # This is an ordinary description.
            columns = self.record.columns

        # Check table and description attributes at the same time.
        # These checks are only valid for non-nested tb.

        # Column names.
        expectedNames = [
            "var0",
            "var1",
            "var1_",
            "var2",
            "var3",
            "var4",
            "var5",
            "var6",
            "var7",
        ]
        self.assertEqual(expectedNames, list(tbl.colnames))
        self.assertEqual(expectedNames, list(desc._v_names))

        # Column types.
        expectedTypes = [columns[colname].dtype for colname in expectedNames]
        self.assertEqual(
            expectedTypes, [tbl.coldtypes[v] for v in expectedNames]
        )
        self.assertEqual(
            expectedTypes, [desc._v_dtypes[v] for v in expectedNames]
        )

        # Column string types.
        expectedTypes = [columns[colname].type for colname in expectedNames]
        self.assertEqual(
            expectedTypes, [tbl.coltypes[v] for v in expectedNames]
        )
        self.assertEqual(
            expectedTypes, [desc._v_types[v] for v in expectedNames]
        )

        # Column defaults.
        for v in expectedNames:
            if common.verbose:
                print("dflt-->", columns[v].dflt)
                print("coldflts-->", tbl.coldflts[v])
                print("desc.dflts-->", desc._v_dflts[v])
            self.assertTrue(
                common.areArraysEqual(tbl.coldflts[v], columns[v].dflt)
            )
            self.assertTrue(
                common.areArraysEqual(desc._v_dflts[v], columns[v].dflt)
            )

        # Column path names.
        self.assertEqual(expectedNames, list(desc._v_pathnames))

        # Column objects.
        for colName in expectedNames:
            expectedCol = columns[colName]
            col = desc._v_colobjects[colName]
            self.assertEqual(expectedCol.dtype, col.dtype)
            self.assertEqual(expectedCol.type, col.type)

    def test01_readTable(self):
        """Checking table read and cuts."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01_readTable..." % self.__class__.__name__)

        # Create an instance of an HDF5 Table
        self.h5file = tb.open_file(self.h5fname, "r")
        table = self.h5file.get_node("/table0")

        # Choose a small value for buffer size
        table.nrowsinbuf = 3
        # Read the records and select those with "var2" file less than 20
        result = [
            r["var2"][0][0] for r in table.iterrows() if r["var2"][0][0] < 20
        ]

        if common.verbose:
            print("Table:", repr(table))
            print("Nrows in", table._v_pathname, ":", table.nrows)
            print("Last record in table ==>", table[-1])
            print("Total selected records in table ==> ", len(result))
        nrows = self.expectedrows - 1
        r = [r for r in table.iterrows() if r["var2"][0][0] < 20][-1]
        self.assertEqual(
            (
                r["var0"][0],
                r["var1"][0][0],
                r["var1_"][0],
                r["var2"][0][0],
                r["var7"],
            ),
            (b"0001", b"0001", nrows, nrows, b"1"),
        )
        if isinstance(r["var5"], np.ndarray):
            self.assertTrue(
                common.allequal(r["var5"], np.array((nrows,) * 4, np.float32))
            )
        else:
            self.assertEqual(r["var5"], float(nrows))
        self.assertEqual(len(result), 20)

    def test01b_readTable(self):
        """Checking table read and cuts (multidimensional columns case)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01b_readTable..." % self.__class__.__name__)

        # Create an instance of an HDF5 Table
        self.h5file = tb.open_file(self.h5fname, "r")
        table = self.h5file.get_node("/table0")

        # Choose a small value for buffer size
        table.nrowsinbuf = 3
        # Read the records and select those with "var2" file less than 20
        result1 = [r["var5"] for r in table.iterrows() if r["var2"][0][0] < 20]
        if common.verbose:
            print("Nrows in", table._v_pathname, ":", table.nrows)
            print("Last record in table ==>", table[-1])
            print("Total selected records in table ==> ", len(result1))
        nrows = table.nrows
        result2 = [r for r in table.iterrows() if r["var2"][0][0] < 20][-1]
        if isinstance(result2["var5"], np.ndarray):
            self.assertTrue(
                common.allequal(
                    result1[0], np.array((float(0),) * 4, np.float32)
                )
            )
            self.assertTrue(
                common.allequal(
                    result1[1], np.array((float(1),) * 4, np.float32)
                )
            )
            self.assertTrue(
                common.allequal(
                    result1[2], np.array((float(2),) * 4, np.float32)
                )
            )
            self.assertTrue(
                common.allequal(
                    result1[3], np.array((float(3),) * 4, np.float32)
                )
            )
            self.assertTrue(
                common.allequal(
                    result1[10], np.array((float(10),) * 4, np.float32)
                )
            )
            self.assertTrue(
                common.allequal(
                    result2["var5"],
                    np.array((float(nrows - 1),) * 4, np.float32),
                )
            )
        else:
            self.assertEqual(result2["var5"], float(nrows - 1))
        self.assertEqual(len(result1), 20)

        # Read the records and select those with "var2" file less than 20
        result1 = [r["var1"] for r in table.iterrows() if r["var2"][0][0] < 20]
        result2 = [r for r in table.iterrows() if r["var2"][0][0] < 20][-1]

        if result2["var1"].dtype.char == "S":
            a = np.array([["%04d" % (self.expectedrows - 0)] * 2] * 2, "S")
            self.assertTrue(common.allequal(result1[0], a))
            a = np.array([["%04d" % (self.expectedrows - 1)] * 2] * 2, "S")
            self.assertTrue(common.allequal(result1[1], a))
            a = np.array([["%04d" % (self.expectedrows - 2)] * 2] * 2, "S")
            self.assertTrue(common.allequal(result1[2], a))
            a = np.array([["%04d" % (self.expectedrows - 3)] * 2] * 2, "S")
            self.assertTrue(common.allequal(result1[3], a))
            a = np.array([["%04d" % (self.expectedrows - 10)] * 2] * 2, "S")
            self.assertTrue(common.allequal(result1[10], a))
            a = np.array([["%04d" % (1)] * 2] * 2, "S")
            self.assertTrue(common.allequal(result2["var1"], a))
        else:
            self.assertEqual(result1["var1"], "0001")
        self.assertEqual(len(result1), 20)

    def test01c_readTable(self):
        """Checking shape of multidimensional columns."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01c_readTable..." % self.__class__.__name__)

        # Create an instance of an HDF5 Table
        self.h5file = tb.open_file(self.h5fname, "r")
        table = self.h5file.get_node("/table0")

        if common.verbose:
            print("var2 col shape:", table.cols.var2.shape)
            print("Should be:", table.cols.var2[:].shape)
        self.assertEqual(table.cols.var2.shape, table.cols.var2[:].shape)

    def test02_AppendRows(self):
        """Checking whether appending record rows works or not."""

        # Now, open it, but in "append" mode
        self.h5file = tb.open_file(self.h5fname, mode="a")
        self.rootgroup = self.h5file.root
        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02_AppendRows..." % self.__class__.__name__)

        # Get a table
        table = self.h5file.get_node("/group0/table1")
        # Get their row object
        row = table.row
        if common.verbose:
            print("Nrows in old", table._v_pathname, ":", table.nrows)
            print("Record Format ==>", table.description._v_nested_formats)
            print("Record Size ==>", table.rowsize)
        # Append some rows
        for i in range(self.appendrows):
            s = "%04d" % (self.appendrows - i)
            row["var0"] = s.encode("ascii")
            row["var1"] = s.encode("ascii")
            row["var7"] = s[-1].encode("ascii")
            row["var1_"] = (i, 1)
            row["var2"] = ((i, 1), (1, 1))  # *-*
            row["var3"] = i % self.maxshort
            if isinstance(row["var4"], np.ndarray):
                row["var4"] = [float(i), float(i * i)]
            else:
                row["var4"] = float(i)
            if isinstance(row["var5"], np.ndarray):
                row["var5"] = np.array((float(i),) * 4)
            else:
                row["var5"] = float(i)
            row.append()

        # Flush the buffer for this table and read it
        table.flush()
        result = [
            r["var2"][0][0] for r in table.iterrows() if r["var2"][0][0] < 20
        ]
        row = [r for r in table.iterrows() if r["var2"][0][0] < 20][-1]

        nrows = self.appendrows - 1
        self.assertEqual(
            (
                row["var0"][0],
                row["var1"][0][0],
                row["var1_"][0],
                row["var2"][0][0],
                row["var7"],
            ),
            (b"0001", b"0001", nrows, nrows, b"1"),
        )
        if isinstance(row["var5"], np.ndarray):
            self.assertTrue(
                common.allequal(
                    row["var5"], np.array((float(nrows),) * 4, np.float32)
                )
            )
        else:
            self.assertEqual(row["var5"], float(nrows))
        if self.appendrows <= 20:
            add = self.appendrows
        else:
            add = 20
        self.assertEqual(len(result), 20 + add)  # because we appended new rows
        # del table

    # CAVEAT: The next test only works for tables with rows < 2**15
    def test03_endianess(self):
        """Checking if table is endianess aware."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03_endianess..." % self.__class__.__name__)

        # Create an instance of an HDF5 Table
        self.h5file = tb.open_file(self.h5fname, "r")
        table = self.h5file.get_node("/group0/group1/table2")

        # Read the records and select the ones with "var3" column less than 20
        result = [r["var2"] for r in table.iterrows() if r["var3"] < 20]
        if common.verbose:
            print("Nrows in", table._v_pathname, ":", table.nrows)
            print("On-disk byteorder ==>", table.byteorder)
            print("Last record in table ==>", table[-1])
            print("Total selected records in table ==>", len(result))
        nrows = self.expectedrows - 1
        r = list(table.iterrows())[-1]
        self.assertEqual((r["var1"][0][0], r["var3"]), (b"0001", nrows))
        self.assertEqual(len(result), 20)


class BasicWriteTestCase(BasicTestCase):
    title = "BasicWrite"


class DictWriteTestCase(BasicTestCase):
    # This checks also unidimensional arrays as columns
    title = "DictWrite"
    record = RecordDescriptionDict
    nrows = 21
    nrowsinbuf = 3  # Choose a small value for the buffer size
    start = 0
    stop = 10
    step = 3


class RecordDTWriteTestCase(BasicTestCase):
    title = "RecordDTWriteTestCase"
    record = RecordDT


# Pure NumPy dtype
class NumPyDTWriteTestCase(BasicTestCase):
    title = "NumPyDTWriteTestCase"
    record = np.dtype("(2,)S4,(2,2)S4,(2,)i4,(2,2)i4,i2,2f8,f4,i2,S1")
    record.names = "var0,var1,var1_,var2,var3,var4,var5,var6,var7".split(",")


class RecArrayOneWriteTestCase(BasicTestCase):
    title = "RecArrayOneWrite"
    record = np.rec.array(
        None,
        formats="(2,)S4,(2,2)S4,(2,)i4,(2,2)i4,i2,2f8,f4,i2,S1",
        names="var0,var1,var1_,var2,var3,var4,var5,var6,var7",
        shape=0,
    )


class RecArrayTwoWriteTestCase(BasicTestCase):
    title = "RecArrayTwoWrite"
    expectedrows = 100
    recarrayinit = 1
    recordtemplate = np.rec.array(
        None,
        formats="(2,)S4,(2,2)S4,(2,)i4,(2,2)i4,i2,f8,f4,i2,S1",
        names="var0,var1,var1_,var2,var3,var4,var5,var6,var7",
        shape=1,
    )


class RecArrayThreeWriteTestCase(BasicTestCase):
    title = "RecArrayThreeWrite"
    expectedrows = 100
    recarrayinit = 1
    recordtemplate = np.rec.array(
        None,
        formats="(2,)S4,(2,2)S4,(2,)i4,(2,2)i4,i2,2f8,4f4,i2,S1",
        names="var0,var1,var1_,var2,var3,var4,var5,var6,var7",
        shape=1,
    )


class RecArrayAlignedWriteTestCase(BasicTestCase):
    title = "RecArrayThreeWrite"
    expectedrows = 100
    recarrayinit = 1
    recordtemplate = np.rec.array(
        None,
        formats="(2,)S4,(2,2)S4,(2,)i4,(2,2)i4,i2,2f8,4f4,i2,S1",
        names="var0,var1,var1_,var2,var3,var4,var5,var6,var7",
        shape=1,
        aligned=True,
    )


@common.unittest.skipIf(
    not common.blosc_avail, "BLOSC compression library not available"
)
class CompressBloscTablesTestCase(BasicTestCase):
    title = "CompressBloscTables"
    compress = 1
    complib = "blosc"


@common.unittest.skipIf(
    not common.lzo_avail, "LZO compression library not available"
)
class CompressLZOTablesTestCase(BasicTestCase):
    title = "CompressLZOTables"
    compress = 1
    complib = "lzo"


@common.unittest.skipIf(
    not common.bzip2_avail, "BZIP2 compression library not available"
)
class CompressBzip2TablesTestCase(BasicTestCase):
    title = "CompressBzip2Tables"
    compress = 1
    complib = "bzip2"


class CompressZLIBTablesTestCase(BasicTestCase):
    title = "CompressOneTables"
    compress = 1
    complib = "zlib"


class CompressTwoTablesTestCase(BasicTestCase):
    title = "CompressTwoTables"
    compress = 1
    # This checks also unidimensional arrays as columns
    record = RecordDescriptionDict


class BigTablesTestCase(BasicTestCase):
    title = "BigTables"
    # 10000 rows takes much more time than we can afford for tests
    # reducing to 1000 would be more than enough
    # F. Alted 2004-01-19

    # expectedrows = 10000
    # appendrows = 1000
    expectedrows = 1000
    appendrows = 100


class BasicRangeTestCase(common.TempFileMixin, common.PyTablesTestCase):
    # file  = "test.h5"
    open_mode = "w"
    title = "This is the table title"
    record = Record
    maxshort = 1 << 15
    expectedrows = 100
    compress = 0
    # Default values
    nrows = 20
    nrowsinbuf = 3  # Choose a small value for the buffer size
    start = 1
    stop = nrows
    checkrecarray = 0
    checkgetCol = 0

    def setUp(self):
        super().setUp()

        # Create an instance of an HDF5 Table
        self.rootgroup = self.h5file.root
        self.populateFile()
        self.h5file.close()

    def populateFile(self):
        group = self.rootgroup
        for j in range(3):
            # Create a table
            table = self.h5file.create_table(
                group,
                "table" + str(j),
                self.record,
                title=self.title,
                filters=tb.Filters(self.compress),
                expectedrows=self.expectedrows,
            )
            # Get the row object associated with the new table
            row = table.row

            # Fill the table
            for i in range(self.expectedrows):
                row["var1"] = "%04d" % (self.expectedrows - i)
                row["var7"] = row["var1"][0][0][-1]
                row["var2"] = i
                row["var3"] = i % self.maxshort
                if isinstance(row["var4"], np.ndarray):
                    row["var4"] = [float(i), float(i * i)]
                else:
                    row["var4"] = float(i)
                if isinstance(row["var5"], np.ndarray):
                    row["var5"] = np.array((float(i),) * 4)
                else:
                    row["var5"] = float(i)
                # var6 will be like var3 but byteswaped
                row["var6"] = ((row["var3"] >> 8) & 0xFF) + (
                    (row["var3"] << 8) & 0xFF00
                )
                row.append()

            # Flush the buffer for this table
            table.flush()
            # Create a new group (descendant of group)
            group2 = self.h5file.create_group(group, "group" + str(j))
            # Iterate over this new group (group2)
            group = group2

    def check_range(self):
        # Create an instance of an HDF5 Table
        self.h5file = tb.open_file(self.h5fname, "r")
        table = self.h5file.get_node("/table0")

        table.nrowsinbuf = self.nrowsinbuf
        resrange = slice(self.start, self.stop, self.step).indices(table.nrows)
        reslength = len(list(range(*resrange)))
        if self.checkrecarray:
            recarray = table.read(self.start, self.stop, self.step)
            result = []
            for nrec in range(len(recarray)):
                if recarray["var2"][nrec][0][0] < self.nrows and 0 < self.step:
                    result.append(recarray["var2"][nrec][0][0])
                elif (
                    recarray["var2"][nrec][0][0] > self.nrows and 0 > self.step
                ):
                    result.append(recarray["var2"][nrec][0][0])
        elif self.checkgetCol:
            column = table.read(self.start, self.stop, self.step, "var2")
            result = []
            for nrec in range(len(column)):
                if column[nrec][0][0] < self.nrows and 0 < self.step:  # *-*
                    result.append(column[nrec][0][0])  # *-*
                elif column[nrec][0][0] > self.nrows and 0 > self.step:  # *-*
                    result.append(column[nrec][0][0])  # *-*
        else:
            if 0 < self.step:
                result = [
                    r["var2"][0][0]
                    for r in table.iterrows(self.start, self.stop, self.step)
                    if r["var2"][0][0] < self.nrows
                ]
            elif 0 > self.step:
                result = [
                    r["var2"][0][0]
                    for r in table.iterrows(self.start, self.stop, self.step)
                    if r["var2"][0][0] > self.nrows
                ]

        if self.start < 0:
            startr = self.expectedrows + self.start
        else:
            startr = self.start

        if self.stop is None:
            if self.checkrecarray or self.checkgetCol:
                # data read using the read method
                stopr = startr + 1
            else:
                # data read using the  iterrows method
                stopr = self.nrows
        elif self.stop < 0:
            stopr = self.expectedrows + self.stop
        else:
            stopr = self.stop

        if self.nrows < stopr:
            stopr = self.nrows

        if common.verbose:
            print("Nrows in", table._v_pathname, ":", table.nrows)
            if reslength:
                if self.checkrecarray:
                    print("Last record *read* in recarray ==>", recarray[-1])
                elif self.checkgetCol:
                    print("Last value *read* in getCol ==>", column[-1])
                else:
                    r = list(table.iterrows(self.start, self.stop, self.step))[
                        -1
                    ]
                    print("Last record *read* in table range ==>", r)
            print("Total number of selected records ==>", len(result))
            print("Selected records:\n", result)
            print(
                "Selected records should look like:\n",
                list(range(startr, stopr, self.step)),
            )
            print("start, stop, step ==>", startr, stopr, self.step)

        self.assertEqual(result, list(range(startr, stopr, self.step)))
        if not (self.checkrecarray or self.checkgetCol):
            if startr < stopr and 0 < self.step:
                r = [
                    r["var2"]
                    for r in table.iterrows(self.start, self.stop, self.step)
                    if r["var2"][0][0] < self.nrows
                ][-1]
                if self.nrows > self.expectedrows:
                    self.assertEqual(
                        r[0][0],
                        list(range(self.start, self.stop, self.step))[-1],
                    )
                else:
                    self.assertEqual(
                        r[0][0], list(range(startr, stopr, self.step))[-1]
                    )
            elif startr > stopr and 0 > self.step:
                r = [
                    r["var2"]
                    for r in table.iterrows(self.start, self.stop, self.step)
                    if r["var2"][0][0] > self.nrows
                ][0]
                if self.nrows < self.expectedrows:
                    self.assertEqual(
                        r[0][0],
                        list(range(self.start, self.stop or -1, self.step))[0],
                    )
                else:
                    self.assertEqual(
                        r[0][0], list(range(startr, stopr or -1, self.step))[0]
                    )

        # Close the file
        self.h5file.close()

    def test01_range(self):
        """Checking ranges in table iterators (case1)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01_range..." % self.__class__.__name__)

        # Case where step < nrowsinbuf < 2 * step
        self.nrows = 21
        self.nrowsinbuf = 3
        self.start = 0
        self.stop = self.expectedrows
        self.step = 2

        self.check_range()

    def test01a_range(self):
        """Checking ranges in table iterators (case1)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01_range..." % self.__class__.__name__)

        # Case where step < nrowsinbuf < 2 * step
        self.nrows = 21
        self.nrowsinbuf = 3
        self.start = self.expectedrows - 1
        self.stop = None
        self.step = -2

        self.check_range()

    def test02_range(self):
        """Checking ranges in table iterators (case2)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02_range..." % self.__class__.__name__)

        # Case where step < nrowsinbuf < 10 * step
        self.nrows = 21
        self.nrowsinbuf = 31
        self.start = 11
        self.stop = self.expectedrows
        self.step = 3

        self.check_range()

    def test03_range(self):
        """Checking ranges in table iterators (case3)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03_range..." % self.__class__.__name__)

        # Case where step < nrowsinbuf < 1.1 * step
        self.nrows = self.expectedrows
        self.nrowsinbuf = 11  # Choose a small value for the buffer size
        self.start = 0
        self.stop = self.expectedrows
        self.step = 10

        self.check_range()

    def test04_range(self):
        """Checking ranges in table iterators (case4)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test04_range..." % self.__class__.__name__)

        # Case where step == nrowsinbuf
        self.nrows = self.expectedrows
        self.nrowsinbuf = 11  # Choose a small value for the buffer size
        self.start = 1
        self.stop = self.expectedrows
        self.step = 11

        self.check_range()

    def test05_range(self):
        """Checking ranges in table iterators (case5)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test05_range..." % self.__class__.__name__)

        # Case where step > 1.1 * nrowsinbuf
        self.nrows = 21
        self.nrowsinbuf = 10  # Choose a small value for the buffer size
        self.start = 1
        self.stop = self.expectedrows
        self.step = 11

        self.check_range()

    def test06_range(self):
        """Checking ranges in table iterators (case6)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test06_range..." % self.__class__.__name__)

        # Case where step > 3 * nrowsinbuf
        self.nrows = 3
        self.nrowsinbuf = 3  # Choose a small value for the buffer size
        self.start = 2
        self.stop = self.expectedrows
        self.step = 10

        self.check_range()

    def test07_range(self):
        """Checking ranges in table iterators (case7)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test07_range..." % self.__class__.__name__)

        # Case where start == stop
        self.nrows = 2
        self.nrowsinbuf = 3  # Choose a small value for the buffer size
        self.start = self.nrows
        self.stop = self.nrows
        self.step = 10

        self.check_range()

    def test08_range(self):
        """Checking ranges in table iterators (case8)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test08_range..." % self.__class__.__name__)

        # Case where start > stop
        self.nrows = 2
        self.nrowsinbuf = 3  # Choose a small value for the buffer size
        self.start = self.nrows + 1
        self.stop = self.nrows
        self.step = 1

        self.check_range()

    def test09_range(self):
        """Checking ranges in table iterators (case9)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test09_range..." % self.__class__.__name__)

        # Case where stop = None
        self.nrows = 100
        self.nrowsinbuf = 3  # Choose a small value for the buffer size
        self.start = 1
        self.stop = 2
        self.step = 1

        self.check_range()

    def test10_range(self):
        """Checking ranges in table iterators (case10)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test10_range..." % self.__class__.__name__)

        # Case where start < 0 and stop = 0
        self.nrows = self.expectedrows
        self.nrowsinbuf = 5  # Choose a small value for the buffer size
        self.start = -6
        self.startr = self.expectedrows + self.start
        self.stop = 0
        self.stopr = self.expectedrows + self.stop
        self.step = 2

        self.check_range()

    def test11_range(self):
        """Checking ranges in table iterators (case11)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test11_range..." % self.__class__.__name__)

        # Case where start < 0 and stop < 0
        self.nrows = self.expectedrows
        self.nrowsinbuf = 5  # Choose a small value for the buffer size
        self.start = -6
        self.startr = self.expectedrows + self.start
        self.stop = -2
        self.stopr = self.expectedrows + self.stop
        self.step = 1

        self.check_range()

    def test12_range(self):
        """Checking ranges in table iterators (case12)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test12_range..." % self.__class__.__name__)

        # Case where start < 0 and stop < 0 and start > stop
        self.nrows = self.expectedrows
        self.nrowsinbuf = 5  # Choose a small value for the buffer size
        self.start = -1
        self.startr = self.expectedrows + self.start
        self.stop = -2
        self.stopr = self.expectedrows + self.stop
        self.step = 1

        self.check_range()

    def test13_range(self):
        """Checking ranges in table iterators (case13)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test13_range..." % self.__class__.__name__)

        # Case where step < 0
        self.step = -11
        try:
            self.check_range()
        except ValueError:
            if common.verbose:
                (type, value, traceback) = sys.exc_info()
                print("\nGreat!, the next ValueError was catched!")
            self.h5file.close()
        # else:
        #     self.fail("expected a ValueError")

        # Case where step == 0
        self.step = 0
        try:
            self.check_range()
        except ValueError:
            if common.verbose:
                (type, value, traceback) = sys.exc_info()
                print("\nGreat!, the next ValueError was catched!")
            self.h5file.close()
        # else:
        #     self.fail("expected a ValueError")


class IterRangeTestCase(BasicRangeTestCase):
    pass


class RecArrayRangeTestCase(BasicRangeTestCase):
    checkrecarray = 1


class GetColRangeTestCase(BasicRangeTestCase):
    checkgetCol = 1

    def test01_nonexistentField(self):
        """Checking non-existing Field in getCol method"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test01_nonexistentField..."
                % self.__class__.__name__
            )

        # Create an instance of an HDF5 Table
        self.h5file = tb.open_file(self.h5fname, "r")
        self.root = self.h5file.root
        table = self.h5file.get_node("/table0")

        with self.assertRaises(KeyError):
            table.read(field="non-existent-column")


class Rec(tb.IsDescription):
    col1 = tb.IntCol(pos=1, shape=(2,))
    col2 = tb.StringCol(itemsize=3, pos=2, shape=(3,))
    col3 = tb.FloatCol(pos=3, shape=(3, 2))


class RecArrayIO(common.TempFileMixin, common.PyTablesTestCase):

    def test00(self):
        """Checking saving a normal recarray"""

        # Create a recarray
        intlist1 = [[456, 23] * 3] * 2
        intlist2 = np.array([[2, 2] * 3] * 2, dtype=int)
        arrlist1 = [["dbe"] * 2] * 3
        arrlist2 = [["de"] * 2] * 3
        floatlist1 = [[1.2, 2.3] * 3] * 4
        floatlist2 = np.array([[4.5, 2.4] * 3] * 4)
        b = [
            (intlist1, arrlist1, floatlist1),
            (intlist2, arrlist2, floatlist2),
        ]
        r = np.rec.array(
            b, formats="(2,6)i4,(3,2)S3,(4,6)f8", names="col1,col2,col3"
        )

        # Save it in a table:
        self.h5file.create_table(self.h5file.root, "recarray", r)

        # Read it again
        r2 = self.h5file.root.recarray.read()

        self.assertEqual(r.tobytes(), r2.tobytes())

    def test01(self):
        """Checking saving a recarray with an offset in its buffer"""

        # Create a recarray
        intlist1 = [[456, 23] * 3] * 2
        intlist2 = np.array([[2, 2] * 3] * 2, dtype=int)
        arrlist1 = [["dbe"] * 2] * 3
        arrlist2 = [["de"] * 2] * 3
        floatlist1 = [[1.2, 2.3] * 3] * 4
        floatlist2 = np.array([[4.5, 2.4] * 3] * 4)
        b = [
            (intlist1, arrlist1, floatlist1),
            (intlist2, arrlist2, floatlist2),
        ]
        r = np.rec.array(
            b, formats="(2,6)i4,(3,2)S3,(4,6)f8", names="col1,col2,col3"
        )

        # Get a view of the recarray
        r1 = r[1:]

        # Save it in a table:
        self.h5file.create_table(self.h5file.root, "recarray", r1)

        # Read it again
        r2 = self.h5file.root.recarray.read()

        self.assertEqual(r1.tobytes(), r2.tobytes())

    def test02(self):
        """Checking saving a slice of a large recarray"""

        # Create a recarray
        intlist1 = [[[23, 24, 35] * 6] * 6]
        intlist2 = np.array([[[2, 3, 4] * 6] * 6], dtype=int)
        arrlist1 = [["dbe"] * 2] * 3
        arrlist2 = [["de"] * 2] * 3
        floatlist1 = [[1.2, 2.3] * 3] * 4
        floatlist2 = np.array([[4.5, 2.4] * 3] * 4)
        b = [
            (intlist1, arrlist1, floatlist1),
            (intlist2, arrlist2, floatlist2),
        ]
        r = np.rec.array(
            b * 300,
            formats="(1,6,18)i4,(3,2)S3,(4,6)f8",
            names="col1,col2,col3",
        )

        # Get a slice of recarray
        r1 = r[290:292]

        # Save it in a table:
        self.h5file.create_table(self.h5file.root, "recarray", r1)

        # Read it again
        r2 = self.h5file.root.recarray.read()

        self.assertEqual(r1.tobytes(), r2.tobytes())

    def test03(self):
        """Checking saving a slice of a strided recarray"""

        # Create a recarray
        intlist1 = [[[23, 24, 35] * 6] * 6]
        intlist2 = np.array([[[2, 3, 4] * 6] * 6], dtype=int)
        arrlist1 = [["dbe"] * 2] * 3
        arrlist2 = [["de"] * 2] * 3
        floatlist1 = [[1.2, 2.3] * 3] * 4
        floatlist2 = np.array([[4.5, 2.4] * 3] * 4)
        b = [
            (intlist1, arrlist1, floatlist1),
            (intlist2, arrlist2, floatlist2),
        ]
        r = np.rec.array(
            b * 300,
            formats="(1,6,18)i4,(3,2)S3,(4,6)f8",
            names="col1,col2,col3",
            shape=600,
        )

        # Get a strided recarray
        r2 = r[::2]

        # Get a slice
        r1 = r2[148:]

        # Save it in a table:
        self.h5file.create_table(self.h5file.root, "recarray", r1)

        # Read it again
        r2 = self.h5file.root.recarray.read()

        self.assertEqual(r1.tobytes(), r2.tobytes())

    def test08a(self):
        """Checking modifying one column (single column version, list)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test08a..." % self.__class__.__name__)

        # Create a new table:
        table = self.h5file.create_table(self.h5file.root, "recarray", Rec)

        # Append new rows
        s0, s1, s2, s3 = ["dbe"] * 3, ["ded"] * 3, ["db1"] * 3, ["de1"] * 3
        f0, f1, f2, f3 = (
            [[1.2] * 2] * 3,
            [[1.3] * 2] * 3,
            [[1.4] * 2] * 3,
            [[1.5] * 2] * 3,
        )
        r = np.rec.array(
            [([456, 457], s0, f0), ([2, 3], s1, f1)],
            formats="(2,)i4,(3,)S3,(3,2)f8",
        )
        table.append(r)
        table.append([([457, 458], s2, f2), ([5, 6], s3, f3)])

        # Modify just one existing column
        table.cols.col1[1:] = [[[2, 3], [3, 4], [4, 5]]]

        # Create the modified recarray
        r1 = np.rec.array(
            [
                ([456, 457], s0, f0),
                ([2, 3], s1, f1),
                ([3, 4], s2, f2),
                ([4, 5], s3, f3),
            ],
            formats="(2,)i4,(3,)S3,(3,2)f8",
            names="col1,col2,col3",
        )

        # Read the modified table
        r2 = table.read()
        if common.verbose:
            print("Original table-->", repr(r2))
            print("Should look like-->", repr(r1))
        self.assertEqual(r1.tobytes(), r2.tobytes())
        self.assertEqual(table.nrows, 4)

    def test08b(self):
        """Checking modifying one column (single column version, recarray)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test08b..." % self.__class__.__name__)

        # Create a new table:
        table = self.h5file.create_table(self.h5file.root, "recarray", Rec)

        # Append new rows
        s0, s1, s2, s3 = ["dbe"] * 3, ["ded"] * 3, ["db1"] * 3, ["de1"] * 3
        f0, f1, f2, f3 = (
            [[1.2] * 2] * 3,
            [[1.3] * 2] * 3,
            [[1.4] * 2] * 3,
            [[1.5] * 2] * 3,
        )
        r = np.rec.array(
            [([456, 457], s0, f0), ([2, 3], s1, f1)],
            formats="(2,)i4,(3,)S3,(3,2)f8",
        )
        table.append(r)
        table.append([([457, 458], s2, f2), ([5, 6], s3, f3)])

        # Modify just one existing column
        columns = np.rec.fromarrays(
            np.array([[[2, 3], [3, 4], [4, 5]]]), formats="i4"
        )
        table.modify_columns(start=1, columns=columns, names=["col1"])

        # Create the modified recarray
        r1 = np.rec.array(
            [
                ([456, 457], s0, f0),
                ([2, 3], s1, f1),
                ([3, 4], s2, f2),
                ([4, 5], s3, f3),
            ],
            formats="(2,)i4,(3,)S3,(3,2)f8",
            names="col1,col2,col3",
        )

        # Read the modified table
        r2 = table.read()
        if common.verbose:
            print("Original table-->", repr(r2))
            print("Should look like-->", repr(r1))
        self.assertEqual(r1.tobytes(), r2.tobytes())
        self.assertEqual(table.nrows, 4)

    def test08b2(self):
        """Checking modifying one column (single column version, recarray,
        modify_column)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test08b2..." % self.__class__.__name__)

        # Create a new table:
        table = self.h5file.create_table(self.h5file.root, "recarray", Rec)

        # Append new rows
        s0, s1, s2, s3 = ["dbe"] * 3, ["ded"] * 3, ["db1"] * 3, ["de1"] * 3
        f0, f1, f2, f3 = (
            [[1.2] * 2] * 3,
            [[1.3] * 2] * 3,
            [[1.4] * 2] * 3,
            [[1.5] * 2] * 3,
        )
        r = np.rec.array(
            [([456, 457], s0, f0), ([2, 3], s1, f1)],
            formats="(2,)i4,(3,)S3,(3,2)f8",
        )
        table.append(r)
        table.append([([457, 458], s2, f2), ([5, 6], s3, f3)])

        # Modify just one existing column
        columns = np.rec.fromarrays(
            np.array([[[2, 3], [3, 4], [4, 5]]]), formats="i4"
        )
        table.modify_column(start=1, column=columns, colname="col1")

        # Create the modified recarray
        r1 = np.rec.array(
            [
                ([456, 457], s0, f0),
                ([2, 3], s1, f1),
                ([3, 4], s2, f2),
                ([4, 5], s3, f3),
            ],
            formats="(2,)i4,(3,)S3,(3,2)f8",
            names="col1,col2,col3",
        )

        # Read the modified table
        r2 = table.read()
        if common.verbose:
            print("Original table-->", repr(r2))
            print("Should look like-->", repr(r1))
        self.assertEqual(r1.tobytes(), r2.tobytes())
        self.assertEqual(table.nrows, 4)


class DefaultValues(common.TempFileMixin, common.PyTablesTestCase):

    def test00(self):
        """Checking saving a Table MD with default values"""

        # Create a table
        table = self.h5file.create_table(self.h5file.root, "table", Record)

        # Take a number of records a bit large
        # nrows = int(table.nrowsinbuf * 1.1)
        nrows = 5  # for test
        # Fill the table with nrows records
        for i in range(nrows):
            if i == 3 or i == 4:
                table.row["var2"] = ((2, 2), (2, 2))  # *-*
            # This injects the row values.
            table.row.append()

        # We need to flush the buffers in table in order to get an
        # accurate number of records on it.
        table.flush()

        # Create a recarray with the same default values
        buffer = [
            (
                ["\x00"]
                * 2,  # just "" does not initialize the buffer properly
                [["abcd", "efgh"]] * 2,
                (1, 1),
                ((1, 1), (1, 1)),
                2,
                3.1,
                4.2,
                5,
                "e",
            )
        ]
        r = np.rec.array(
            buffer * nrows,
            formats="(2,)S4,(2,2)S4,(2,)i4,(2,2)i4,i2,f8,f4,u2,S1",
            names=[
                "var0",
                "var1",
                "var1_",
                "var2",
                "var3",
                "var4",
                "var5",
                "var6",
                "var7",
            ],
        )  # *-*

        # Assign the value exceptions
        r["var2"][3] = ((2, 2), (2, 2))  # *-*
        r["var2"][4] = ((2, 2), (2, 2))  # *-*

        # Read the table in another recarray
        r2 = table.read()

        # This generates too much output. Activate only when
        # self.nrowsinbuf is very small (<10)
        if common.verbose and 1:
            print("Table values:")
            print(r2)
            print("Record values:")
            print(r)

        # Both checks do work, however, tobytes() seems more stringent.
        self.assertEqual(r.tobytes(), r2.tobytes())
        # self.assertTrue(common.areArraysEqual(r,r2))


class RecordT(tb.IsDescription):
    var0 = tb.IntCol(dflt=1, shape=())  # native int
    var1 = tb.IntCol(dflt=[1], shape=(1,))  # 1-D int (one element)
    var2_s = tb.IntCol(dflt=[1, 1], shape=2)  # 1-D int (two elements)
    var2 = tb.IntCol(dflt=[1, 1], shape=(2,))  # 1-D int (two elements)
    var3 = tb.IntCol(dflt=[[0, 0], [1, 1]], shape=(2, 2))  # 2-D int


class ShapeTestCase(common.TempFileMixin, common.PyTablesTestCase):

    def setUp(self):
        super().setUp()
        self.populateFile()

    def populateFile(self):
        table = self.h5file.create_table(self.h5file.root, "table", RecordT)
        row = table.row
        # Fill the table with some rows with default values
        for i in range(1):
            row.append()

        # Flush the buffer for this table
        table.flush()

    def test00(self):
        """Checking scalar shapes"""

        if self.reopen:
            self._reopen()
        table = self.h5file.root.table

        if common.verbose:
            print("The values look like:", table.cols.var0[:])
            print("They should look like:", [1])

        # The real check
        self.assertEqual(table.cols.var0[:].tolist(), [1])

    def test01(self):
        """Checking undimensional (one element) shapes"""

        if self.reopen:
            self._reopen()
        table = self.h5file.root.table

        if common.verbose:
            print("The values look like:", table.cols.var1[:])
            print("They should look like:", [[1]])

        # The real check
        self.assertEqual(table.cols.var1[:].tolist(), [[1]])

    def test02(self):
        """Checking undimensional (two elements) shapes"""

        if self.reopen:
            self._reopen()
        table = self.h5file.root.table

        if common.verbose:
            print("The values look like:", table.cols.var2[:])
            print("They should look like:", [[1, 1]])

        # The real check
        self.assertEqual(table.cols.var2[:].tolist(), [[1, 1]])
        self.assertEqual(table.cols.var2_s[:].tolist(), [[1, 1]])

    def test03(self):
        """Checking bidimensional shapes"""

        if self.reopen:
            self._reopen()
        table = self.h5file.root.table

        if common.verbose:
            print("The values look like:", table.cols.var3[:])
            print("They should look like:", [[[0, 0], [1, 1]]])

        # The real check
        self.assertEqual(table.cols.var3[:].tolist(), [[[0, 0], [1, 1]]])


class ShapeTestCase1(ShapeTestCase):
    reopen = 0


class ShapeTestCase2(ShapeTestCase):
    reopen = 1


class SetItemTestCase(common.TempFileMixin, common.PyTablesTestCase):

    def setUp(self):
        super().setUp()

        # Create a new table:
        self.table = self.h5file.create_table(
            self.h5file.root, "recarray", Rec
        )
        self.table.nrowsinbuf = self.buffersize  # set buffer value

    def test01(self):
        """Checking modifying one table row with __setitem__"""

        table = self.table
        formats = table.description._v_nested_formats

        # append new rows
        r = np.rec.array([(456, "dbe", 1.2), (2, "ded", 1.3)], formats=formats)
        table.append(r)
        table.append([(457, "db1", 1.2), (5, "de1", 1.3)])

        # Modify just one existing row
        table[2] = (456, "db2", 1.2)

        # Create the modified recarray
        r1 = np.rec.array(
            [
                (456, "dbe", 1.2),
                (2, "ded", 1.3),
                (456, "db2", 1.2),
                (5, "de1", 1.3),
            ],
            formats=formats,
            names="col1,col2,col3",
        )
        # Read the modified table
        if self.reopen:
            self._reopen()
            table = self.h5file.root.recarray
            table.nrowsinbuf = self.buffersize  # set buffer value
        r2 = table.read()
        if common.verbose:
            print("Original table-->", repr(r2))
            print("Should look like-->", repr(r1))
        self.assertEqual(r1.tobytes(), r2.tobytes())
        self.assertEqual(table.nrows, 4)

    def test01b(self):
        """Checking modifying one table row with __setitem__ (long index)"""

        table = self.table
        formats = table.description._v_nested_formats

        # append new rows
        r = np.rec.array([(456, "dbe", 1.2), (2, "ded", 1.3)], formats=formats)
        table.append(r)
        table.append([(457, "db1", 1.2), (5, "de1", 1.3)])

        # Modify just one existing row
        table[2] = (456, "db2", 1.2)

        # Create the modified recarray
        r1 = np.rec.array(
            [
                (456, "dbe", 1.2),
                (2, "ded", 1.3),
                (456, "db2", 1.2),
                (5, "de1", 1.3),
            ],
            formats=formats,
            names="col1,col2,col3",
        )

        # Read the modified table
        if self.reopen:
            self._reopen()
            table = self.h5file.root.recarray
            table.nrowsinbuf = self.buffersize  # set buffer value
        r2 = table.read()
        if common.verbose:
            print("Original table-->", repr(r2))
            print("Should look like-->", repr(r1))
        self.assertEqual(r1.tobytes(), r2.tobytes())
        self.assertEqual(table.nrows, 4)

    def test02(self):
        """Modifying one row, with a step (__setitem__)"""

        table = self.table
        formats = table.description._v_nested_formats

        # append new rows
        r = np.rec.array([(456, "dbe", 1.2), (2, "ded", 1.3)], formats=formats)
        table.append(r)
        table.append([(457, "db1", 1.2), (5, "de1", 1.3)])

        # Modify two existing rows
        rows = np.rec.array([(457, "db1", 1.2)], formats=formats)
        table[1:3:2] = rows

        # Create the modified recarray
        r1 = np.rec.array(
            [
                (456, "dbe", 1.2),
                (457, "db1", 1.2),
                (457, "db1", 1.2),
                (5, "de1", 1.3),
            ],
            formats=formats,
            names="col1,col2,col3",
        )

        # Read the modified table
        if self.reopen:
            self._reopen()
            table = self.h5file.root.recarray
            table.nrowsinbuf = self.buffersize  # set buffer value
        r2 = table.read()
        if common.verbose:
            print("Original table-->", repr(r2))
            print("Should look like-->", repr(r1))
        self.assertEqual(r1.tobytes(), r2.tobytes())
        self.assertEqual(table.nrows, 4)

    def test03(self):
        """Checking modifying several rows at once (__setitem__)"""

        table = self.table
        formats = table.description._v_nested_formats

        # append new rows
        r = np.rec.array([(456, "dbe", 1.2), (2, "ded", 1.3)], formats=formats)
        table.append(r)
        table.append([(457, "db1", 1.2), (5, "de1", 1.3)])

        # Modify two existing rows
        rows = np.rec.array(
            [(457, "db1", 1.2), (5, "de1", 1.3)], formats=formats
        )

        # table.modify_rows(start=1, rows=rows)
        table[1:3] = rows

        # Create the modified recarray
        r1 = np.rec.array(
            [
                (456, "dbe", 1.2),
                (457, "db1", 1.2),
                (5, "de1", 1.3),
                (5, "de1", 1.3),
            ],
            formats=formats,
            names="col1,col2,col3",
        )

        # Read the modified table
        if self.reopen:
            self._reopen()
            table = self.h5file.root.recarray
            table.nrowsinbuf = self.buffersize  # set buffer value
        r2 = table.read()
        if common.verbose:
            print("Original table-->", repr(r2))
            print("Should look like-->", repr(r1))
        self.assertEqual(r1.tobytes(), r2.tobytes())
        self.assertEqual(table.nrows, 4)

    def test04(self):
        """Modifying several rows at once, with a step (__setitem__)"""

        table = self.table
        formats = table.description._v_nested_formats

        # append new rows
        r = np.rec.array([(456, "dbe", 1.2), (2, "ded", 1.3)], formats=formats)
        table.append(r)
        table.append([(457, "db1", 1.2), (5, "de1", 1.3)])

        # Modify two existing rows
        rows = np.rec.array(
            [(457, "db1", 1.2), (6, "de2", 1.3)], formats=formats
        )
        # table[1:4:2] = rows
        table[1::2] = rows

        # Create the modified recarray
        r1 = np.rec.array(
            [
                (456, "dbe", 1.2),
                (457, "db1", 1.2),
                (457, "db1", 1.2),
                (6, "de2", 1.3),
            ],
            formats=formats,
            names="col1,col2,col3",
        )

        # Read the modified table
        if self.reopen:
            self._reopen()
            table = self.h5file.root.recarray
            table.nrowsinbuf = self.buffersize  # set buffer value
        r2 = table.read()
        if common.verbose:
            print("Original table-->", repr(r2))
            print("Should look like-->", repr(r1))
        self.assertEqual(r1.tobytes(), r2.tobytes())
        self.assertEqual(table.nrows, 4)

    def test05(self):
        """Checking modifying one column (single element, __setitem__)"""

        table = self.table
        formats = table.description._v_nested_formats

        # append new rows
        r = np.rec.array([(456, "dbe", 1.2), (2, "ded", 1.3)], formats=formats)
        table.append(r)
        table.append([(457, "db1", 1.2), (5, "de1", 1.3)])

        # Modify just one existing column
        table.cols.col1[1] = -1

        # Create the modified recarray
        r1 = np.rec.array(
            [
                (456, "dbe", 1.2),
                (-1, "ded", 1.3),
                (457, "db1", 1.2),
                (5, "de1", 1.3),
            ],
            formats=formats,
            names="col1,col2,col3",
        )

        # Read the modified table
        if self.reopen:
            self._reopen()
            table = self.h5file.root.recarray
            table.nrowsinbuf = self.buffersize  # set buffer value
        r2 = table.read()
        if common.verbose:
            print("Original table-->", repr(r2))
            print("Should look like-->", repr(r1))
        self.assertEqual(r1.tobytes(), r2.tobytes())
        self.assertEqual(table.nrows, 4)

    def test06a(self):
        """Checking modifying one column (several elements, __setitem__)"""

        table = self.table
        formats = table.description._v_nested_formats

        # append new rows
        r = np.rec.array([(456, "dbe", 1.2), (2, "ded", 1.3)], formats=formats)
        table.append(r)
        table.append([(457, "db1", 1.2), (5, "de1", 1.3)])

        # Modify just one existing column
        table.cols.col1[1:4] = [(2, 2), (3, 3), (4, 4)]

        # Create the modified recarray
        r1 = np.rec.array(
            [
                (456, "dbe", 1.2),
                (2, "ded", 1.3),
                (3, "db1", 1.2),
                (4, "de1", 1.3),
            ],
            formats=formats,
            names="col1,col2,col3",
        )
        # Read the modified table
        if self.reopen:
            self._reopen()
            table = self.h5file.root.recarray
            table.nrowsinbuf = self.buffersize  # set buffer value
        r2 = table.read()
        if common.verbose:
            print("Original table-->", repr(r2))
            print("Should look like-->", repr(r1))
        self.assertEqual(r1.tobytes(), r2.tobytes())
        self.assertEqual(table.nrows, 4)

    def test06b(self):
        """Checking modifying one column (iterator, __setitem__)"""

        table = self.table
        formats = table.description._v_nested_formats

        # append new rows
        r = np.rec.array([(456, "dbe", 1.2), (2, "ded", 1.3)], formats=formats)
        table.append(r)
        table.append([(457, "db1", 1.2), (5, "de1", 1.3)])

        # Modify just one existing column
        with self.assertRaises(NotImplementedError):
            for row in table.iterrows():
                row["col1"] = row.nrow + 1
                row.append()
            table.flush()

    def test07(self):
        """Modifying one column (several elements, __setitem__, step)"""

        table = self.table
        formats = table.description._v_nested_formats

        # append new rows
        r = np.rec.array([(456, "dbe", 1.2), (1, "ded", 1.3)], formats=formats)
        table.append(r)
        table.append([(457, "db1", 1.2), (5, "de1", 1.3)])

        # Modify just one existing column
        table.cols.col1[1:4:2] = [(2, 2), (3, 3)]

        # Create the modified recarray
        r1 = np.rec.array(
            [
                (456, "dbe", 1.2),
                (2, "ded", 1.3),
                (457, "db1", 1.2),
                (3, "de1", 1.3),
            ],
            formats=formats,
            names="col1,col2,col3",
        )

        # Read the modified table
        if self.reopen:
            self._reopen()
            table = self.h5file.root.recarray
            table.nrowsinbuf = self.buffersize  # set buffer value
        r2 = table.read()
        if common.verbose:
            print("Original table-->", repr(r2))
            print("Should look like-->", repr(r1))
        self.assertEqual(r1.tobytes(), r2.tobytes())
        self.assertEqual(table.nrows, 4)

    def test08(self):
        """Modifying one column (one element, __setitem__, step)"""

        table = self.table
        formats = table.description._v_nested_formats

        # append new rows
        r = np.rec.array([(456, "dbe", 1.2), (2, "ded", 1.3)], formats=formats)
        table.append(r)
        table.append([(457, "db1", 1.2), (5, "de1", 1.3)])

        # Modify just one existing column
        table.cols.col1[1:4:3] = [(2, 2)]

        # Create the modified recarray
        r1 = np.rec.array(
            [
                (456, "dbe", 1.2),
                (2, "ded", 1.3),
                (457, "db1", 1.2),
                (5, "de1", 1.3),
            ],
            formats=formats,
            names="col1,col2,col3",
        )

        # Read the modified table
        if self.reopen:
            self._reopen()
            table = self.h5file.root.recarray
            table.nrowsinbuf = self.buffersize  # set buffer value
        r2 = table.read()
        if common.verbose:
            print("Original table-->", repr(r2))
            print("Should look like-->", repr(r1))
        self.assertEqual(r1.tobytes(), r2.tobytes())
        self.assertEqual(table.nrows, 4)

    def test09(self):
        """Modifying beyond the table extend (__setitem__, step)"""

        table = self.table
        formats = table.description._v_nested_formats

        # append new rows
        r = np.rec.array([(456, "dbe", 1.2), (2, "ded", 1.3)], formats=formats)
        table.append(r)
        table.append([(457, "db1", 1.2), (5, "de1", 1.3)])

        # Try to modify beyond the extent
        # This will silently exclude the non-fitting rows
        rows = np.rec.array(
            [(457, "db1", 1.2), (6, "de2", 1.3)], formats=formats
        )
        table[1::2] = rows

        # How it should look like
        r1 = np.rec.array(
            [
                (456, "dbe", 1.2),
                (457, "db1", 1.2),
                (457, "db1", 1.2),
                (6, "de2", 1.3),
            ],
            formats=formats,
        )

        # Read the modified table
        if self.reopen:
            self._reopen()
            table = self.h5file.root.recarray
            table.nrowsinbuf = self.buffersize  # set buffer value
        r2 = table.read()
        if common.verbose:
            print("Original table-->", repr(r2))
            print("Should look like-->", repr(r1))
        self.assertEqual(r1.tobytes(), r2.tobytes())
        self.assertEqual(table.nrows, 4)


class SetItemTestCase1(SetItemTestCase):
    reopen = 0
    buffersize = 1


class SetItemTestCase2(SetItemTestCase):
    reopen = 1
    buffersize = 2


class SetItemTestCase3(SetItemTestCase):
    reopen = 0
    buffersize = 1000


class SetItemTestCase4(SetItemTestCase):
    reopen = 1
    buffersize = 1000


class UpdateRowTestCase(common.TempFileMixin, common.PyTablesTestCase):

    def setUp(self):
        super().setUp()

        # Create a new table:
        self.table = self.h5file.create_table(
            self.h5file.root, "recarray", Rec
        )
        self.table.nrowsinbuf = self.buffersize  # set buffer value

    def test01(self):
        """Checking modifying one table row with Row.update"""

        table = self.table
        formats = table.description._v_nested_formats

        # append new rows
        r = np.rec.array([(456, "dbe", 1.2), (2, "ded", 1.3)], formats=formats)
        table.append(r)
        table.append([(457, "db1", 1.2), (5, "de1", 1.3)])

        # Modify just one existing row
        for row in table.iterrows(2, 3):
            (row["col1"], row["col2"], row["col3"]) = [456, "db2", 1.2]
            row.update()

        # Create the modified recarray
        r1 = np.rec.array(
            [
                (456, "dbe", 1.2),
                (2, "ded", 1.3),
                (456, "db2", 1.2),
                (5, "de1", 1.3),
            ],
            formats=formats,
            names="col1,col2,col3",
        )

        # Read the modified table
        if self.reopen:
            self._reopen()
            table = self.h5file.root.recarray
            table.nrowsinbuf = self.buffersize  # set buffer value
        r2 = table.read()
        if common.verbose:
            print("Original table-->", repr(r2))
            print("Should look like-->", repr(r1))
        self.assertEqual(r1.tobytes(), r2.tobytes())
        self.assertEqual(table.nrows, 4)

    def test02(self):
        """Modifying one row, with a step (Row.update)"""

        table = self.table
        formats = table.description._v_nested_formats

        # append new rows
        r = np.rec.array([(456, "dbe", 1.2), (2, "ded", 1.3)], formats=formats)
        table.append(r)
        table.append([(457, "db1", 1.2), (5, "de1", 1.3)])

        # Modify two existing rows
        for row in table.iterrows(1, 3, 2):
            if row.nrow == 1:
                (row["col1"], row["col2"], row["col3"]) = (457, "db1", 1.2)
            elif row.nrow == 3:
                (row["col1"], row["col2"], row["col3"]) = (6, "de2", 1.3)
            row.update()

        # Create the modified recarray
        r1 = np.rec.array(
            [
                (456, "dbe", 1.2),
                (457, "db1", 1.2),
                (457, "db1", 1.2),
                (5, "de1", 1.3),
            ],
            formats=formats,
            names="col1,col2,col3",
        )

        # Read the modified table
        if self.reopen:
            self._reopen()
            table = self.h5file.root.recarray
            table.nrowsinbuf = self.buffersize  # set buffer value
        r2 = table.read()
        if common.verbose:
            print("Original table-->", repr(r2))
            print("Should look like-->", repr(r1))
        self.assertEqual(r1.tobytes(), r2.tobytes())
        self.assertEqual(table.nrows, 4)

    def test03(self):
        """Checking modifying several rows at once (Row.update)"""

        table = self.table
        formats = table.description._v_nested_formats

        # append new rows
        r = np.rec.array([(456, "dbe", 1.2), (2, "ded", 1.3)], formats=formats)
        table.append(r)
        table.append([(457, "db1", 1.2), (5, "de1", 1.3)])

        # Modify two existing rows
        for row in table.iterrows(1, 3):
            if row.nrow == 1:
                (row["col1"], row["col2"], row["col3"]) = (457, "db1", 1.2)
            elif row.nrow == 2:
                (row["col1"], row["col2"], row["col3"]) = (5, "de1", 1.3)
            row.update()

        # Create the modified recarray
        r1 = np.rec.array(
            [
                (456, "dbe", 1.2),
                (457, "db1", 1.2),
                (5, "de1", 1.3),
                (5, "de1", 1.3),
            ],
            formats=formats,
            names="col1,col2,col3",
        )

        # Read the modified table
        if self.reopen:
            self._reopen()
            table = self.h5file.root.recarray
            table.nrowsinbuf = self.buffersize  # set buffer value
        r2 = table.read()
        if common.verbose:
            print("Original table-->", repr(r2))
            print("Should look like-->", repr(r1))
        self.assertEqual(r1.tobytes(), r2.tobytes())
        self.assertEqual(table.nrows, 4)

    def test04(self):
        """Modifying several rows at once, with a step (Row.update)"""

        table = self.table
        formats = table.description._v_nested_formats

        # append new rows
        r = np.rec.array([(456, "dbe", 1.2), (2, "ded", 1.3)], formats=formats)
        table.append(r)
        table.append([(457, "db1", 1.2), (5, "de1", 1.3)])

        # Modify two existing rows
        for row in table.iterrows(1, stop=4, step=2):
            if row.nrow == 1:
                (row["col1"], row["col2"], row["col3"]) = (457, "db1", 1.2)
            elif row.nrow == 3:
                (row["col1"], row["col2"], row["col3"]) = (6, "de2", 1.3)
            row.update()

        # Create the modified recarray
        r1 = np.rec.array(
            [
                (456, "dbe", 1.2),
                (457, "db1", 1.2),
                (457, "db1", 1.2),
                (6, "de2", 1.3),
            ],
            formats=formats,
            names="col1,col2,col3",
        )

        # Read the modified table
        if self.reopen:
            self._reopen()
            table = self.h5file.root.recarray
            table.nrowsinbuf = self.buffersize  # set buffer value
        r2 = table.read()
        if common.verbose:
            print("Original table-->", repr(r2))
            print("Should look like-->", repr(r1))
        self.assertEqual(r1.tobytes(), r2.tobytes())
        self.assertEqual(table.nrows, 4)

    def test05(self):
        """Checking modifying one column (single element, Row.update)"""

        table = self.table
        formats = table.description._v_nested_formats

        # append new rows
        r = np.rec.array([(456, "dbe", 1.2), (2, "ded", 1.3)], formats=formats)
        table.append(r)
        table.append([(457, "db1", 1.2), (5, "de1", 1.3)])

        # Modify just one existing column
        for row in table.iterrows(1, 2):
            row["col1"] = -1
            row.update()

        # Create the modified recarray
        r1 = np.rec.array(
            [
                (456, "dbe", 1.2),
                (-1, "ded", 1.3),
                (457, "db1", 1.2),
                (5, "de1", 1.3),
            ],
            formats=formats,
            names="col1,col2,col3",
        )

        # Read the modified table
        if self.reopen:
            self._reopen()
            table = self.h5file.root.recarray
            table.nrowsinbuf = self.buffersize  # set buffer value
        r2 = table.read()
        if common.verbose:
            print("Original table-->", repr(r2))
            print("Should look like-->", repr(r1))
        self.assertEqual(r1.tobytes(), r2.tobytes())
        self.assertEqual(table.nrows, 4)

    def test06(self):
        """Checking modifying one column (several elements, Row.update)"""

        table = self.table
        formats = table.description._v_nested_formats

        # append new rows
        r = np.rec.array([(456, "dbe", 1.2), (2, "ded", 1.3)], formats=formats)
        table.append(r)
        table.append([(457, "db1", 1.2), (5, "de1", 1.3)])

        # Modify just one existing column
        for row in table.iterrows(1, 4):
            row["col1"] = row.nrow + 1
            row.update()

        # Create the modified recarray
        r1 = np.rec.array(
            [
                (456, "dbe", 1.2),
                (2, "ded", 1.3),
                (3, "db1", 1.2),
                (4, "de1", 1.3),
            ],
            formats=formats,
            names="col1,col2,col3",
        )

        # Read the modified table
        if self.reopen:
            self._reopen()
            table = self.h5file.root.recarray
            table.nrowsinbuf = self.buffersize  # set buffer value
        r2 = table.read()
        if common.verbose:
            print("Original table-->", repr(r2))
            print("Should look like-->", repr(r1))
        self.assertEqual(r1.tobytes(), r2.tobytes())
        self.assertEqual(table.nrows, 4)

    def test07(self):
        """Modifying values from a selection"""

        table = self.table
        formats = table.description._v_nested_formats

        # append new rows
        r = np.rec.array([(456, "dbe", 1.2), (1, "ded", 1.3)], formats=formats)
        table.append(r)
        table.append([(457, "db1", 1.2), (5, "de1", 1.3)])

        # Modify just rows with col1 < 456
        for row in table.iterrows():
            if row["col1"][0] < 456:
                row["col1"] = 2
                row["col2"] = "ada"
                row.update()

        # Create the modified recarray
        r1 = np.rec.array(
            [
                (456, "dbe", 1.2),
                (2, "ada", 1.3),
                (457, "db1", 1.2),
                (2, "ada", 1.3),
            ],
            formats=formats,
            names="col1,col2,col3",
        )

        # Read the modified table
        if self.reopen:
            self._reopen()
            table = self.h5file.root.recarray
            table.nrowsinbuf = self.buffersize  # set buffer value
        r2 = table.read()
        if common.verbose:
            print("Original table-->", repr(r2))
            print("Should look like-->", repr(r1))
        self.assertEqual(r1.tobytes(), r2.tobytes())
        self.assertEqual(table.nrows, 4)

    def test08(self):
        """Modifying a large table (Row.update)"""

        table = self.table
        formats = table.description._v_nested_formats

        nrows = 100

        # append new rows
        row = table.row
        for i in range(nrows):
            row["col1"] = i - 1
            row["col2"] = "a" + str(i - 1)
            row["col3"] = -1.0
            row.append()
        table.flush()

        # Modify all the rows
        for row in table.iterrows():
            row["col1"] = row.nrow
            row["col2"] = "b" + str(row.nrow)
            row["col3"] = 0.0
            row.update()

        # Create the modified recarray
        r1 = np.rec.array(
            None, shape=nrows, formats=formats, names="col1,col2,col3"
        )
        for i in range(nrows):
            r1["col1"][i] = i
            r1["col2"][i] = "b" + str(i)
            r1["col3"][i] = 0.0

        # Read the modified table
        if self.reopen:
            self._reopen()
            table = self.h5file.root.recarray
            table.nrowsinbuf = self.buffersize  # set buffer value
        r2 = table.read()
        if common.verbose:
            print("Original table-->", repr(r2))
            print("Should look like-->", repr(r1))
        self.assertEqual(r1.tobytes(), r2.tobytes())
        self.assertEqual(table.nrows, nrows)

    def test08b(self):
        """Setting values on a large table without calling Row.update"""

        table = self.table
        formats = table.description._v_nested_formats

        nrows = 100

        # append new rows
        row = table.row
        for i in range(nrows):
            row["col1"] = i - 1
            row["col2"] = "a" + str(i - 1)
            row["col3"] = -1.0
            row.append()
        table.flush()

        # Modify all the rows (actually don't)
        for row in table.iterrows():
            row["col1"] = row.nrow
            row["col2"] = "b" + str(row.nrow)
            row["col3"] = 0.0
            # row.update()

        # Create the modified recarray
        r1 = np.rec.array(
            None, shape=nrows, formats=formats, names="col1,col2,col3"
        )
        for i in range(nrows):
            r1["col1"][i] = i - 1
            r1["col2"][i] = "a" + str(i - 1)
            r1["col3"][i] = -1.0

        # Read the modified table
        if self.reopen:
            self._reopen()
            table = self.h5file.root.recarray
            table.nrowsinbuf = self.buffersize  # set buffer value
        r2 = table.read()
        if common.verbose:
            print("Original table-->", repr(r2))
            print("Should look like-->", repr(r1))
        self.assertEqual(r1.tobytes(), r2.tobytes())
        self.assertEqual(table.nrows, nrows)

    def test09(self):
        """Modifying selected values on a large table"""

        table = self.table
        formats = table.description._v_nested_formats

        nrows = 100

        # append new rows
        row = table.row
        for i in range(nrows):
            row["col1"] = i - 1
            row["col2"] = "a" + str(i - 1)
            row["col3"] = -1.0
            row.append()
        table.flush()

        # Modify selected rows
        for row in table.iterrows():
            if row["col1"][0] > nrows - 3:
                row["col1"] = row.nrow
                row["col2"] = "b" + str(row.nrow)
                row["col3"] = 0.0
                row.update()

        # Create the modified recarray
        r1 = np.rec.array(
            None, shape=nrows, formats=formats, names="col1,col2,col3"
        )
        for i in range(nrows):
            r1["col1"][i] = i - 1
            r1["col2"][i] = "a" + str(i - 1)
            r1["col3"][i] = -1.0

        # modify just the last line
        r1["col1"][i] = i
        r1["col2"][i] = "b" + str(i)
        r1["col3"][i] = 0.0

        # Read the modified table
        if self.reopen:
            self._reopen()
            table = self.h5file.root.recarray
            table.nrowsinbuf = self.buffersize  # set buffer value
        r2 = table.read()
        if common.verbose:
            print("Original table-->", repr(r2))
            print("Should look like-->", repr(r1))
        self.assertEqual(r1.tobytes(), r2.tobytes())
        self.assertEqual(table.nrows, nrows)

    def test09b(self):
        """Modifying selected values on a large table (alternate values)"""

        table = self.table
        formats = table.description._v_nested_formats

        nrows = 100

        # append new rows
        row = table.row
        for i in range(nrows):
            row["col1"] = i - 1
            row["col2"] = "a" + str(i - 1)
            row["col3"] = -1.0
            row.append()
        table.flush()

        # Modify selected rows
        for row in table.iterrows(step=10):
            row["col1"] = row.nrow
            row["col2"] = "b" + str(row.nrow)
            row["col3"] = 0.0
            row.update()

        # Create the modified recarray
        r1 = np.rec.array(
            None, shape=nrows, formats=formats, names="col1,col2,col3"
        )
        for i in range(nrows):
            if i % 10 > 0:
                r1["col1"][i] = i - 1
                r1["col2"][i] = "a" + str(i - 1)
                r1["col3"][i] = -1.0
            else:
                r1["col1"][i] = i
                r1["col2"][i] = "b" + str(i)
                r1["col3"][i] = 0.0

        # Read the modified table
        if self.reopen:
            self._reopen()
            table = self.h5file.root.recarray
            table.nrowsinbuf = self.buffersize  # set buffer value
        r2 = table.read()
        if common.verbose:
            print("Original table-->", repr(r2))
            print("Should look like-->", repr(r1))
        self.assertEqual(r1.tobytes(), r2.tobytes())
        self.assertEqual(table.nrows, nrows)


class UpdateRowTestCase1(UpdateRowTestCase):
    reopen = 0
    buffersize = 1


class UpdateRowTestCase2(UpdateRowTestCase):
    reopen = 1
    buffersize = 2


class UpdateRowTestCase3(UpdateRowTestCase):
    reopen = 0
    buffersize = 1000


class UpdateRowTestCase4(UpdateRowTestCase):
    reopen = 1
    buffersize = 1000


def suite():
    theSuite = common.unittest.TestSuite()
    niter = 1
    # common.heavy = 1  # Uncomment this only for testing purposes

    for n in range(niter):
        theSuite.addTest(common.make_suite(BasicWriteTestCase))
        theSuite.addTest(common.make_suite(DictWriteTestCase))
        theSuite.addTest(common.make_suite(RecordDTWriteTestCase))
        theSuite.addTest(common.make_suite(NumPyDTWriteTestCase))
        theSuite.addTest(common.make_suite(RecArrayOneWriteTestCase))
        theSuite.addTest(common.make_suite(RecArrayTwoWriteTestCase))
        theSuite.addTest(common.make_suite(RecArrayThreeWriteTestCase))
        theSuite.addTest(common.make_suite(RecArrayAlignedWriteTestCase))
        theSuite.addTest(common.make_suite(CompressZLIBTablesTestCase))
        theSuite.addTest(common.make_suite(CompressTwoTablesTestCase))
        theSuite.addTest(common.make_suite(IterRangeTestCase))
        theSuite.addTest(common.make_suite(RecArrayRangeTestCase))
        theSuite.addTest(common.make_suite(GetColRangeTestCase))
        theSuite.addTest(common.make_suite(DefaultValues))
        theSuite.addTest(common.make_suite(RecArrayIO))
        theSuite.addTest(common.make_suite(ShapeTestCase1))
        theSuite.addTest(common.make_suite(ShapeTestCase2))
        theSuite.addTest(common.make_suite(SetItemTestCase1))
        theSuite.addTest(common.make_suite(SetItemTestCase2))
        theSuite.addTest(common.make_suite(SetItemTestCase3))
        theSuite.addTest(common.make_suite(SetItemTestCase4))
        theSuite.addTest(common.make_suite(UpdateRowTestCase1))
        theSuite.addTest(common.make_suite(UpdateRowTestCase2))
        theSuite.addTest(common.make_suite(UpdateRowTestCase3))
        theSuite.addTest(common.make_suite(UpdateRowTestCase4))
        theSuite.addTest(common.make_suite(CompressBloscTablesTestCase))
        theSuite.addTest(common.make_suite(CompressLZOTablesTestCase))
    if common.heavy:
        theSuite.addTest(common.make_suite(CompressBzip2TablesTestCase))
        theSuite.addTest(common.make_suite(BigTablesTestCase))

    return theSuite


if __name__ == "__main__":
    common.parse_argv(sys.argv)
    common.print_versions()
    common.unittest.main(defaultTest="suite")
