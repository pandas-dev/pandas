"""Test module for nested types under PyTables."""

import sys
import itertools

import numpy as np

import tables as tb
from tables.tests import common

minRowIndex = 10


# This is the structure of the table used for testing (DON'T PANIC!):
#
# +-+---------------------------------+-----+----------+-+-+
# |x|Info                             |color|info      |y|z|
# | +-----+--+----------------+----+--+     +----+-----+ | |
# | |value|y2|Info2           |name|z2|     |Name|Value| | |
# | |     |  +----+-----+--+--+    |  |     |    |     | | |
# | |     |  |name|value|y3|z3|    |  |     |    |     | | |
# +-+-----+--+----+-----+--+--+----+--+-----+----+-----+-+-+
#
# Please note that some fields are explicitly ordered while others are
# ordered alphabetically by name.
# The declaration of the nested table:
class Info(tb.IsDescription):
    _v_pos = 3
    Name = tb.StringCol(itemsize=2)
    Value = tb.ComplexCol(itemsize=16)


class TestTDescr(tb.IsDescription):
    """A description that has several nested columns."""

    x = tb.Int32Col(dflt=0, shape=2, pos=0)  # 0
    y = tb.Float64Col(dflt=1, shape=(2, 2))
    z = tb.UInt8Col(dflt=1)
    color = tb.StringCol(itemsize=2, dflt=b" ", pos=2)
    info = Info()

    class Info(tb.IsDescription):  # 1
        _v_pos = 1
        name = tb.StringCol(itemsize=2)
        value = tb.ComplexCol(itemsize=16, pos=0)  # 0
        y2 = tb.Float64Col(dflt=1, pos=1)  # 1
        z2 = tb.UInt8Col(dflt=1)

        class Info2(tb.IsDescription):
            y3 = tb.Time64Col(dflt=1, shape=2)
            z3 = tb.EnumCol({"r": 4, "g": 2, "b": 1}, "r", "int32", shape=2)
            name = tb.StringCol(itemsize=2)
            value = tb.ComplexCol(itemsize=16, shape=2)


# The corresponding nested array description:
testADescr = [
    ("x", "(2,)int32"),
    (
        "Info",
        [
            ("value", "complex128"),
            ("y2", "float64"),
            (
                "Info2",
                [
                    ("name", "S2"),
                    ("value", "(2,)complex128"),
                    ("y3", "(2,)float64"),
                    ("z3", "(2,)int32"),
                ],
            ),
            ("name", "S2"),
            ("z2", "uint8"),
        ],
    ),
    ("color", "S2"),
    ("info", [("Name", "S2"), ("Value", "complex128")]),
    ("y", "(2,2)float64"),
    ("z", "uint8"),
]

# The corresponding nested array description (brief version):
testADescr2 = [
    ("x", "(2,)i4"),
    (
        "Info",
        [
            ("value", "()c16"),
            ("y2", "()f8"),
            (
                "Info2",
                [
                    ("name", "()S2"),
                    ("value", "(2,)c16"),
                    ("y3", "(2,)f8"),
                    ("z3", "(2,)i4"),
                ],
            ),
            ("name", "()S2"),
            ("z2", "()u1"),
        ],
    ),
    ("color", "()S2"),
    ("info", [("Name", "()S2"), ("Value", "()c16")]),
    ("y", "(2, 2)f8"),
    ("z", "()u1"),
]

# A nested array for testing:
testABuffer = [
    # x     Info    color info      y       z
    #       value y2 Info2      name z2         Name Value
    #                name   value    y3       z3
    (
        (3, 2),
        (6j, 6.0, ("nn", (6j, 4j), (6.0, 4.0), (1, 2)), "NN", 8),
        "cc",
        ("NN", 6j),
        ((6.0, 4.0), (6.0, 4.0)),
        8,
    ),
    (
        (4, 3),
        (7j, 7.0, ("oo", (7j, 5j), (7.0, 5.0), (2, 1)), "OO", 9),
        "dd",
        ("OO", 7j),
        ((7.0, 5.0), (7.0, 5.0)),
        9,
    ),
]
testAData = np.array(testABuffer, dtype=testADescr)
# The name of the column to be searched:
testCondCol = "Info/z2"
# The name of a nested column (it can not be searched):
testNestedCol = "Info"
# The condition to be applied on the column (all but the last row match it):
testCondition = "(2 < col) & (col < 9)"


def areDescriptionsEqual(desc1, desc2):
    """Are both `desc1` and `desc2` equivalent descriptions?

    The arguments may be description objects (``IsDescription``,
    ``Description``) or dictionaries.

    """

    if isinstance(desc1, tb.Col):
        # This is a rough comparison but it suffices here.
        return (
            desc1.type == desc2.type
            and desc2.dtype == desc2.dtype
            and desc1._v_pos == desc2._v_pos
            # and desc1.dflt == desc2.dflt)
            and common.areArraysEqual(desc1.dflt, desc2.dflt)
        )

    if hasattr(desc1, "_v_colobjects"):  # quacks like a Description
        cols1 = desc1._v_colobjects
    elif hasattr(desc1, "columns"):  # quacks like an IsDescription
        cols1 = desc1.columns
    else:  # hope it quacks like a dictionary
        cols1 = desc1

    if hasattr(desc2, "_v_colobjects"):  # quacks like a Description
        cols2 = desc2._v_colobjects
    elif hasattr(desc2, "columns"):  # quacks like an IsDescription
        cols2 = desc2.columns
    else:  # hope it quacks like a dictionary
        cols2 = desc2

    if len(cols1) != len(cols2):
        return False

    for colName, colobj1 in cols1.items():
        colobj2 = cols2[colName]
        if colName == "_v_pos":
            # The comparison may not be quite exhaustive!
            return colobj1 == colobj2
        if not areDescriptionsEqual(colobj1, colobj2):
            return False

    return True


# Test creating nested column descriptions
class DescriptionTestCase(common.PyTablesTestCase):
    _TestTDescr = TestTDescr
    _testADescr = testADescr
    _testADescr2 = testADescr2
    _testAData = testAData

    def test00_instance(self):
        """Creating an instance of a nested description."""

        self.assertTrue(
            areDescriptionsEqual(self._TestTDescr, self._TestTDescr()),
            "Table description does not match the given one.",
        )

    def test01_instance(self):
        """Checking attrs of an instance of a nested description."""

        descr = tb.description.Description(self._TestTDescr().columns)
        if common.verbose:
            print("Generated description:", descr._v_nested_descr)
            print("Should look like:", self._testADescr2)
        self.assertEqual(
            self._testADescr2,
            descr._v_nested_descr,
            "Description._v_nested_descr does not match.",
        )


# Test creating a nested table and opening it
class CreateTestCase(common.TempFileMixin, common.PyTablesTestCase):
    _TestTDescr = TestTDescr
    _testABuffer = testABuffer
    _testAData = testAData

    def _checkColumns(self, cols, desc):
        """Check that `cols` has all the accessors for `self._TestTDescr`."""

        # ``_desc`` is a leaf column and ``cols`` a ``Column``.
        if isinstance(desc, tb.Col):
            return isinstance(cols, tb.Column)

        # ``_desc`` is a description object and ``cols`` a ``Cols``.
        descColumns = desc._v_colobjects
        for colName in descColumns:
            if colName not in cols._v_colnames:
                return False
            if not self._checkColumns(
                cols._f_col(colName), descColumns[colName]
            ):
                return False

        return True

    def _checkDescription(self, table):
        """Check that description of `table` matches `self._TestTDescr`."""

        # Compare descriptions.
        self.assertTrue(
            areDescriptionsEqual(self._TestTDescr, table.description),
            "Table description does not match the given one.",
        )
        # Check access to columns.
        self._checkColumns(table.cols, table.description)

    def _checkColinstances(self, table):
        """Check that ``colinstances`` and ``cols`` of `table` match."""
        for colpathname in table.description._v_pathnames:
            self.assertTrue(
                table.colinstances[colpathname]
                is table.cols._f_col(colpathname)
            )

    def test00_create(self):
        """Creating a nested table."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )
        self._checkDescription(tbl)
        self._checkColinstances(tbl)

    def test01_open(self):
        """Opening a nested table."""

        self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )
        self._reopen()
        tbl = self.h5file.root.test
        self._checkDescription(tbl)
        self._checkColinstances(tbl)

    def test02_NestedRecArrayCompat(self):
        """Creating a compatible nested record array``."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )

        nrarr = np.array(testABuffer, dtype=tbl.description._v_nested_descr)
        self.assertTrue(
            common.areArraysEqual(nrarr, self._testAData),
            "Can not create a compatible structured array.",
        )

    def test03_NRA(self):
        """Creating a table from a nested record array object."""

        tbl = self.h5file.create_table(
            "/", "test", self._testAData, title=self._getMethodName()
        )
        tbl.flush()
        readAData = tbl.read()
        if common.verbose:
            print("Read data:", readAData)
            print("Should look like:", self._testAData)
        self.assertTrue(
            common.areArraysEqual(self._testAData, readAData),
            "Written and read values differ.",
        )

    def test04_NRA2(self):
        """Creating a table from a generated nested record array object."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )
        tbl.append(self._testAData)
        readAData = tbl.read()

        tbl2 = self.h5file.create_table(
            "/", "test2", readAData, title=self._getMethodName()
        )
        readAData2 = tbl2.read()

        self.assertTrue(
            common.areArraysEqual(self._testAData, readAData2),
            "Written and read values differ.",
        )


# Test writing data in a nested table
class WriteTestCase(common.TempFileMixin, common.PyTablesTestCase):
    _TestTDescr = TestTDescr
    _testAData = testAData
    _testCondition = testCondition
    _testCondCol = testCondCol
    _testNestedCol = testNestedCol

    def _testCondVars(self, table):
        """Get condition variables for the given `table`."""
        return {"col": table.cols._f_col(self._testCondCol)}

    def _testNestedCondVars(self, table):
        """Get condition variables for the given `table`."""
        return {"col": table.cols._f_col(self._testNestedCol)}

    def _appendRow(self, row, index):
        """
        Append the `index`-th row in `self._testAData` to `row`.

        Values are set field-by-field (be it nested or not).
        """

        record = self._testAData[index]
        for fieldName in self._testAData.dtype.names:
            row[fieldName] = record[fieldName]
        row.append()

    def test00_append(self):
        """Appending a set of rows."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )
        tbl.append(self._testAData)
        tbl.flush()

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        readAData = tbl.read()
        self.assertTrue(
            common.areArraysEqual(self._testAData, readAData),
            "Written and read values differ.",
        )

    def test01_row(self):
        """Appending individual rows."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )

        row = tbl.row
        # Add the first row
        self._appendRow(row, 0)
        # Add the rest of the rows field by field.
        for i in range(1, len(self._testAData)):
            self._appendRow(row, i)
        tbl.flush()

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        readAData = tbl.read()
        self.assertTrue(
            common.areArraysEqual(self._testAData, readAData),
            "Written and read values differ.",
        )

    def test02_where(self):
        """Searching nested data."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )
        tbl.append(self._testAData)
        tbl.flush()

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        searchedCoords = tbl.get_where_list(
            self._testCondition, self._testCondVars(tbl)
        )

        # All but the last row match the condition.
        searchedCoords.sort()
        self.assertEqual(
            searchedCoords.tolist(),
            list(range(len(self._testAData) - 1)),
            "Search returned incorrect results.",
        )

    def test02b_whereAppend(self):
        """Searching nested data and appending it to another table."""

        tbl1 = self.h5file.create_table(
            "/", "test1", self._TestTDescr, title=self._getMethodName()
        )
        tbl1.append(self._testAData)
        tbl1.flush()

        tbl2 = self.h5file.create_table(
            "/", "test2", self._TestTDescr, title=self._getMethodName()
        )
        tbl1.append_where(tbl2, self._testCondition, self._testCondVars(tbl1))

        if self.reopen:
            self._reopen()
            tbl1 = self.h5file.root.test1
            tbl2 = self.h5file.root.test2

        searchedCoords = tbl2.get_where_list(
            self._testCondition, self._testCondVars(tbl2)
        )

        # All but the last row match the condition.
        searchedCoords.sort()
        self.assertEqual(
            searchedCoords.tolist(),
            list(range(len(self._testAData) - 1)),
            "Search returned incorrect results.",
        )

    def test03_colscond(self):
        """Searching on a column with nested columns."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )
        tbl.append(self._testAData)
        tbl.flush()

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        self.assertRaises(
            TypeError,
            tbl.get_where_list,
            self._testCondition,
            self._testNestedCondVars(tbl),
        )

    def test04_modifyColumn(self):
        """Modifying one single nested column (modify_column)."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )
        tbl.append(self._testAData)
        tbl.flush()

        nColumn = self._testNestedCol
        # Get the nested column data and swap the first and last rows.
        raTable = self._testAData.copy()
        raColumn = raTable[nColumn]
        # The next will not work until NestedRecords supports copies
        (raColumn[0], raColumn[-1]) = (raColumn[-1], raColumn[0])

        # Write the resulting column and re-read the whole table.
        tbl.modify_column(colname=nColumn, column=raColumn)
        tbl.flush()

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        raReadTable = tbl.read()
        if common.verbose:
            print("Table read:", raReadTable)
            print("Should look like:", raTable)

        # Compare it to the written one.
        self.assertTrue(
            common.areArraysEqual(raTable, raReadTable),
            "Written and read values differ.",
        )

    def test05a_modifyColumns(self):
        """Modifying one nested column (modify_columns)."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )
        tbl.append(self._testAData)
        tbl.flush()

        nColumn = self._testNestedCol
        # Get the nested column data and swap the first and last rows.
        raTable = self._testAData.copy()
        raColumn = raTable[nColumn]
        (raColumn[0], raColumn[-1]) = (raColumn[-1].copy(), raColumn[0].copy())
        newdtype = np.dtype([(nColumn, raTable.dtype.fields[nColumn][0])])
        self.assertIsNotNone(newdtype)

        # Write the resulting column and re-read the whole table.
        tbl.modify_columns(names=[nColumn], columns=raColumn)
        tbl.flush()

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        raReadTable = tbl.read()
        if common.verbose:
            print("Table read:", raReadTable)
            print("Should look like:", raTable)

        # Compare it to the written one.
        self.assertTrue(
            common.areArraysEqual(raTable, raReadTable),
            "Written and read values differ.",
        )

    def test05b_modifyColumns(self):
        """Modifying two nested columns (modify_columns)."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )
        tbl.append(self._testAData)
        tbl.flush()

        # Get the nested column data and swap the first and last rows.
        colnames = ["x", "color"]  # Get the first two columns
        raCols = np.rec.fromarrays(
            [self._testAData["x"].copy(), self._testAData["color"].copy()],
            dtype=[("x", "(2,)i4"), ("color", "S2")],
        )
        # descr=tbl.description._v_nested_descr[0:2])
        # or...
        # names=tbl.description._v_nested_names[0:2],
        # formats=tbl.description._v_nested_formats[0:2])
        (raCols[0], raCols[-1]) = (raCols[-1].copy(), raCols[0].copy())

        # Write the resulting columns
        tbl.modify_columns(names=colnames, columns=raCols)
        tbl.flush()

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        # Re-read the appropriate columns
        raCols2 = np.rec.fromarrays(
            [tbl.cols._f_col("x"), tbl.cols._f_col("color")],
            dtype=raCols.dtype,
        )
        if common.verbose:
            print("Table read:", raCols2)
            print("Should look like:", raCols)

        # Compare it to the written one.
        self.assertTrue(
            common.areArraysEqual(raCols, raCols2),
            "Written and read values differ.",
        )

    def test06_modifyRows(self):
        """Checking modifying several rows at once (using nested rec array)"""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )
        tbl.append(self._testAData)
        tbl.flush()

        # Get the nested record and swap the first and last rows.
        raTable = self._testAData.copy()
        (raTable[0], raTable[-1]) = (raTable[-1].copy(), raTable[0].copy())

        # Write the resulting nested record and re-read the whole table.
        tbl.modify_rows(start=0, stop=2, rows=raTable)
        tbl.flush()

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        raReadTable = tbl.read()
        if common.verbose:
            print("Table read:", raReadTable)
            print("Should look like:", raTable)

        # Compare it to the written one.
        self.assertTrue(
            common.areArraysEqual(raTable, raReadTable),
            "Written and read values differ.",
        )

    def test07_index(self):
        """Checking indexes of nested columns."""

        tbl = self.h5file.create_table(
            "/",
            "test",
            self._TestTDescr,
            title=self._getMethodName(),
            expectedrows=minRowIndex * 2,
        )
        for i in range(minRowIndex):
            tbl.append(self._testAData)
        tbl.flush()
        coltoindex = tbl.cols._f_col(self._testCondCol)
        indexrows = coltoindex.create_index()
        self.assertIsNotNone(indexrows)

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test
            coltoindex = tbl.cols._f_col(self._testCondCol)

        if common.verbose:
            print("Number of written rows:", tbl.nrows)
            print("Number of indexed rows:", coltoindex.index.nelements)

        # Check indexing flags:
        self.assertEqual(tbl.indexed, True, "Table not indexed")
        self.assertNotEqual(coltoindex.index, None, "Column not indexed")
        self.assertTrue(
            tbl.colindexed[self._testCondCol], "Column not indexed"
        )
        # Do a look-up for values
        searchedCoords = tbl.get_where_list(
            self._testCondition, self._testCondVars(tbl)
        )
        searchedCoords.sort()

        expectedCoords = np.arange(0, minRowIndex * 2, 2, tb.utils.SizeType)
        if common.verbose:
            print("Searched coords:", searchedCoords)
            print("Expected coords:", expectedCoords)
        # All even rows match the condition.
        self.assertEqual(
            searchedCoords.tolist(),
            expectedCoords.tolist(),
            "Search returned incorrect results.",
        )

    def test08_setNestedField(self):
        """Checking modifying a nested field via natural naming."""
        # See ticket #93 (http://www.pytables.org/trac/ticket/93).

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )
        tbl.append(self._testAData)
        tbl.flush()

        oldvalue = tbl.cols.Info.z2[0]
        tbl.cols.Info.z2[0] = oldvalue + 1
        tbl.flush()

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        newvalue = tbl.cols.Info.z2[0]
        self.assertEqual(newvalue, oldvalue + 1)


class WriteNoReopen(WriteTestCase):
    reopen = 0


class WriteReopen(WriteTestCase):
    reopen = 1


class ReadTestCase(common.TempFileMixin, common.PyTablesTestCase):
    _TestTDescr = TestTDescr
    _testABuffer = testABuffer
    _testAData = testAData
    _testNestedCol = testNestedCol

    def test00a_repr(self):
        """Checking representation of a nested Table."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title="test00"
        )
        tbl.append(self._testAData)

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        if common.verbose:
            print("str(tbl)-->", str(tbl))
            print("repr(tbl)-->", repr(tbl))

        self.assertEqual(
            str(tbl), f"/test (Table({np.int64(2)!r},)) {np.str_('test00')!r}"
        )
        tblrepr = repr(tbl)
        # Remove the platform-dependent information (i.e. byteorder)
        tblrepr = "\n".join(tblrepr.split("\n")[:-2]) + "\n"
        template = f"""/test (Table({np.int64(2)!r},)) {np.str_('test00')!r}
  description := {{
  "x": Int32Col(shape=({np.int64(2)!r},), dflt={np.int32(0)!r}, pos=0),
  "Info": {{
    "value": ComplexCol(itemsize=16, shape=(), dflt={np.complex128(0j)!r}, pos=0),
    "y2": Float64Col(shape=(), dflt={np.float64(1.0)!r}, pos=1),
    "Info2": {{
      "name": StringCol(itemsize=2, shape=(), dflt={np.bytes_(b'')!r}, pos=0),
      "value": ComplexCol(itemsize=16, shape=({np.int64(2)!r},), dflt={np.complex128(0j)!r}, pos=1),
      "y3": Time64Col(shape=({np.int64(2)!r},), dflt={np.float64(1.0)!r}, pos=2),
      "z3": EnumCol(enum=Enum({{%(value)s}}), dflt='%(default)s', base=Int32Atom(shape=(), dflt={np.int32(0)!r}), shape=({np.int64(2)!r},), pos=3)}},
    "name": StringCol(itemsize=2, shape=(), dflt={np.bytes_(b'')!r}, pos=3),
    "z2": UInt8Col(shape=(), dflt={np.uint8(1)!r}, pos=4)}},
  "color": StringCol(itemsize=2, shape=(), dflt={np.bytes_(b' ')!r}, pos=2),
  "info": {{
    "Name": StringCol(itemsize=2, shape=(), dflt={np.bytes_(b'')!r}, pos=0),
    "Value": ComplexCol(itemsize=16, shape=(), dflt={np.complex128(0j)!r}, pos=1)}},
  "y": Float64Col(shape=({np.int64(2)!r}, {np.int64(2)!r}), dflt={np.float64(1.0)!r}, pos=4),
  "z": UInt8Col(shape=(), dflt={np.uint8(1)!r}, pos=5)}}
"""

        # The problem here is that the order in which items are stored in a
        # dict can't be assumed to be stable.
        # From python 3.3 on it is actually no more stable since the
        # "Hash randomization" feature is enable by default.
        #
        # For this reason we generate a representation string for each of the
        # permutations of the Enum items.
        #
        # Also the default value of enum types is not preserved in HDF5.
        # It is assumed that the default value is the first one in the array
        # of Enum names and hence it is also affected by the issue related to
        # the "Hash randomization" feature.
        #
        # Also in this case it is generated a representation string for each
        # of the possible default values.
        enums = [
            ", ".join(items)
            for items in itertools.permutations(
                (
                    f"'r': {np.int32(4)!r}",
                    f"'b': {np.int32(1)!r}",
                    f"'g': {np.int32(2)!r}",
                )
                if self.reopen
                else ("'r': 4", "'b': 1", "'g': 2")
            )
        ]
        defaults = ("r", "b", "g")
        values = [
            template % {"value": v, "default": d}
            for v, d in itertools.product(enums, defaults)
        ]
        self.assertIn(tblrepr, values)

    def test00b_repr(self):
        """Checking representation of a root Column."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title="test00"
        )
        tbl.append(self._testAData)

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        if common.verbose:
            print("str(tbl.cols.y)-->'%s'" % str(tbl.cols.y))
            print("repr(tbl.cols.y)-->'%s'" % repr(tbl.cols.y))

        self.assertEqual(
            str(tbl.cols.y),
            f"/test.cols.y (Column({np.int64(2)!r}, 2, 2), float64, idx=None)",
        )
        self.assertEqual(
            repr(tbl.cols.y),
            f"/test.cols.y (Column({np.int64(2)!r}, 2, 2), float64, idx=None)",
        )

    def test00c_repr(self):
        """Checking representation of a nested Column."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title="test00"
        )
        tbl.append(self._testAData)

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        if common.verbose:
            print("str(tbl.cols.Info.z2)-->'%s'" % str(tbl.cols.Info.z2))
            print("repr(tbl.cols.Info.z2)-->'%s'" % repr(tbl.cols.Info.z2))

        self.assertEqual(
            str(tbl.cols.Info.z2),
            f"/test.cols.Info.z2 (Column({np.int64(2)!r},), uint8, idx=None)",
        )
        self.assertEqual(
            repr(tbl.cols.Info.z2),
            f"/test.cols.Info.z2 (Column({np.int64(2)!r},), uint8, idx=None)",
        )

    def test01_read(self):
        """Checking Table.read with subgroups with a range index with step."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )
        tbl.append(self._testAData)

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        nrarr = np.rec.array(
            testABuffer, dtype=tbl.description._v_nested_descr
        )
        tblcols = tbl.read(start=0, step=2, field="Info")
        nrarrcols = nrarr["Info"][0::2]
        if common.verbose:
            print("Read cols:", tblcols)
            print("Should look like:", nrarrcols)
        self.assertTrue(
            common.areArraysEqual(nrarrcols, tblcols),
            "Original array are retrieved doesn't match.",
        )

    def test01_read_out_arg(self):
        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )
        tbl.append(self._testAData)

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        nrarr = np.rec.array(
            testABuffer, dtype=tbl.description._v_nested_descr
        )
        # When reading an entire nested column, the output array must contain
        # all fields in the table.  The output buffer will contain the contents
        # of all fields.  The selected column alone will be returned from the
        # method call.
        all_cols = np.empty(1, tbl.dtype)
        tblcols = tbl.read(start=0, step=2, field="Info", out=all_cols)
        nrarrcols = nrarr["Info"][0::2]
        if common.verbose:
            print("Read cols:", tblcols)
            print("Should look like:", nrarrcols)
        self.assertTrue(
            common.areArraysEqual(nrarrcols, tblcols),
            "Original array are retrieved doesn't match.",
        )
        self.assertTrue(
            common.areArraysEqual(nrarr[0::2], all_cols),
            "Output buffer does not match full table.",
        )

    def test02_read(self):
        """Checking Table.read with a nested Column."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )
        tbl.append(self._testAData)

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        tblcols = tbl.read(start=0, step=2, field="Info/value")
        nrarr = np.rec.array(
            testABuffer, dtype=tbl.description._v_nested_descr
        )
        nrarrcols = nrarr["Info"]["value"][0::2]
        self.assertTrue(
            common.areArraysEqual(nrarrcols, tblcols),
            "Original array are retrieved doesn't match.",
        )

    def test02_read_out_arg(self):
        """Checking Table.read with a nested Column."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )
        tbl.append(self._testAData)

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        tblcols = np.empty(1, dtype="c16")
        tbl.read(start=0, step=2, field="Info/value", out=tblcols)
        nrarr = np.rec.array(
            testABuffer, dtype=tbl.description._v_nested_descr
        )
        nrarrcols = nrarr["Info"]["value"][0::2]
        self.assertTrue(
            common.areArraysEqual(nrarrcols, tblcols),
            "Original array are retrieved doesn't match.",
        )


class ReadNoReopen(ReadTestCase):
    reopen = 0


class ReadReopen(ReadTestCase):
    reopen = 1


# Checking the Table.Cols accessor
class ColsTestCase(common.TempFileMixin, common.PyTablesTestCase):
    _TestTDescr = TestTDescr
    _testABuffer = testABuffer
    _testAData = testAData
    _testNestedCol = testNestedCol

    def test00a_repr(self):
        """Checking string representation of Cols."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title="test00"
        )

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        if common.verbose:
            print("str(tbl.cols)-->", str(tbl.cols))
            print("repr(tbl.cols)-->", repr(tbl.cols))

        self.assertEqual(str(tbl.cols), "/test.cols (Cols), 6 columns")
        try:
            self.assertEqual(
                repr(tbl.cols),
                """/test.cols (Cols), 6 columns
  x (Column(0, 2), ('int32',(2,)))
  Info (Cols(), Description)
  color (Column(0,), |S2)
  info (Cols(), Description)
  y (Column(0, 2, 2), ('float64',(2, 2)))
  z (Column(0,), uint8)
""",
            )
        except AssertionError:
            self.assertEqual(
                repr(tbl.cols),
                f"""/test.cols (Cols), 6 columns
  x (Column({np.int64(0)!r}, 2), ('{np.int32(0).dtype.str}', (2,)))
  Info (Cols(), Description)
  color (Column({np.int64(0)!r},), |S2)
  info (Cols(), Description)
  y (Column({np.int64(0)!r}, 2, 2), ('{np.float64(0).dtype.str}', (2, 2)))
  z (Column({np.int64(0)!r},), uint8)
""",
            )

    def test00b_repr(self):
        """Checking string representation of nested Cols."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        if common.verbose:
            print("str(tbl.cols.Info)-->", str(tbl.cols.Info))
            print("repr(tbl.cols.Info)-->", repr(tbl.cols.Info))

        self.assertEqual(
            str(tbl.cols.Info), "/test.cols.Info (Cols), 5 columns"
        )
        self.assertEqual(
            repr(tbl.cols.Info),
            f"""/test.cols.Info (Cols), 5 columns
  value (Column({np.int64(0)!r},), complex128)
  y2 (Column({np.int64(0)!r},), float64)
  Info2 (Cols(), Description)
  name (Column({np.int64(0)!r},), |S2)
  z2 (Column({np.int64(0)!r},), uint8)
""",
        )

    def test01a_f_col(self):
        """Checking cols._f_col() with a subgroup."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        tblcol = tbl.cols._f_col(self._testNestedCol)
        if common.verbose:
            print("Column group name:", tblcol._v_desc._v_pathname)
        self.assertEqual(
            tblcol._v_desc._v_pathname,
            self._testNestedCol,
            "Column group name doesn't match.",
        )

    def test01b_f_col(self):
        """Checking cols._f_col() with a column."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        tblcol = tbl.cols._f_col(self._testNestedCol + "/name")
        if common.verbose:
            print("Column name:", tblcol.name)
        self.assertEqual(tblcol.name, "name", "Column name doesn't match.")

    def test01c_f_col(self):
        """Checking cols._f_col() with a nested subgroup."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )

        tblcol = tbl.cols._f_col(self._testNestedCol + "/Info2")
        if common.verbose:
            print("Column group name:", tblcol._v_desc._v_pathname)
        self.assertEqual(
            tblcol._v_desc._v_pathname,
            self._testNestedCol + "/Info2",
            "Column group name doesn't match.",
        )

    def test02a__len__(self):
        """Checking cols.__len__() in root level."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        length = len(tbl.cols)
        if common.verbose:
            print("Column group length:", length)
        self.assertEqual(
            length, len(tbl.colnames), "Column group length doesn't match."
        )

    def test02b__len__(self):
        """Checking cols.__len__() in subgroup level."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        length = len(tbl.cols.Info)
        if common.verbose:
            print("Column group length:", length)
        self.assertEqual(
            length,
            len(tbl.cols.Info._v_colnames),
            "Column group length doesn't match.",
        )

    def test03a__getitem__(self):
        """Checking cols.__getitem__() with a single index."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )
        tbl.append(self._testAData)

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        nrarr = np.array(testABuffer, dtype=tbl.description._v_nested_descr)
        tblcols = tbl.cols[1]
        nrarrcols = nrarr[1]
        if common.verbose:
            print("Read cols:", tblcols)
            print("Should look like:", nrarrcols)
        self.assertTrue(
            common.areArraysEqual(nrarrcols, tblcols),
            "Original array are retrieved doesn't match.",
        )

    def test03b__getitem__(self):
        """Checking cols.__getitem__() with a range index."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )
        tbl.append(self._testAData)

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        nrarr = np.array(testABuffer, dtype=tbl.description._v_nested_descr)
        tblcols = tbl.cols[0:2]
        nrarrcols = nrarr[0:2]
        if common.verbose:
            print("Read cols:", tblcols)
            print("Should look like:", nrarrcols)
        self.assertTrue(
            common.areArraysEqual(nrarrcols, tblcols),
            "Original array are retrieved doesn't match.",
        )

    def test03c__getitem__(self):
        """Checking cols.__getitem__() with a range index with step."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )
        tbl.append(self._testAData)

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        nrarr = np.array(testABuffer, dtype=tbl.description._v_nested_descr)
        tblcols = tbl.cols[0::2]
        nrarrcols = nrarr[0::2]
        if common.verbose:
            print("Read cols:", tblcols)
            print("Should look like:", nrarrcols)
        self.assertTrue(
            common.areArraysEqual(nrarrcols, tblcols),
            "Original array are retrieved doesn't match.",
        )

    def test04a__getitem__(self):
        """Checking cols.__getitem__() with subgroups with a single index."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )
        tbl.append(self._testAData)

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        nrarr = np.array(testABuffer, dtype=tbl.description._v_nested_descr)
        tblcols = tbl.cols._f_col("Info")[1]
        nrarrcols = nrarr["Info"][1]
        if common.verbose:
            print("Read cols:", tblcols)
            print("Should look like:", nrarrcols)
        self.assertTrue(
            common.areArraysEqual(nrarrcols, tblcols),
            "Original array are retrieved doesn't match.",
        )

    def test04b__getitem__(self):
        """Checking cols.__getitem__() with subgroups with a range index."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )
        tbl.append(self._testAData)

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        nrarr = np.array(testABuffer, dtype=tbl.description._v_nested_descr)
        tblcols = tbl.cols._f_col("Info")[0:2]
        nrarrcols = nrarr["Info"][0:2]
        if common.verbose:
            print("Read cols:", tblcols)
            print("Should look like:", nrarrcols)
        self.assertTrue(
            common.areArraysEqual(nrarrcols, tblcols),
            "Original array are retrieved doesn't match.",
        )

    def test04c__getitem__(self):
        """Checking cols.__getitem__() with subgroups with a range index with
        step."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )
        tbl.append(self._testAData)

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        nrarr = np.array(testABuffer, dtype=tbl.description._v_nested_descr)
        tblcols = tbl.cols._f_col("Info")[0::2]
        nrarrcols = nrarr["Info"][0::2]
        if common.verbose:
            print("Read cols:", tblcols)
            print("Should look like:", nrarrcols)
        self.assertTrue(
            common.areArraysEqual(nrarrcols, tblcols),
            "Original array are retrieved doesn't match.",
        )

    def test05a__getitem__(self):
        """Checking cols.__getitem__() with a column with a single index."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )
        tbl.append(self._testAData)

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        nrarr = np.array(testABuffer, dtype=tbl.description._v_nested_descr)
        tblcols = tbl.cols._f_col("Info/value")[1]
        nrarrcols = nrarr["Info"]["value"][1]
        if common.verbose:
            print("Read cols:", tblcols)
            print("Should look like:", nrarrcols)
        self.assertEqual(
            nrarrcols, tblcols, "Original array are retrieved doesn't match."
        )

    def test05b__getitem__(self):
        """Checking cols.__getitem__() with a column with a range index."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )
        tbl.append(self._testAData)

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        nrarr = np.array(testABuffer, dtype=tbl.description._v_nested_descr)
        tblcols = tbl.cols._f_col("Info/value")[0:2]
        nrarrcols = nrarr["Info"]["value"][0:2]
        if common.verbose:
            print("Read cols:", tblcols)
            print("Should look like:", nrarrcols)
        self.assertTrue(
            common.areArraysEqual(nrarrcols, tblcols),
            "Original array are retrieved doesn't match.",
        )

    def test05c__getitem__(self):
        """Checking cols.__getitem__() with a column with a range index with
        step."""

        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )
        tbl.append(self._testAData)

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        nrarr = np.array(testABuffer, dtype=tbl.description._v_nested_descr)
        tblcols = tbl.cols._f_col("Info/value")[0::2]
        nrarrcols = nrarr["Info"]["value"][0::2]
        if common.verbose:
            print("Read cols:", tblcols)
            print("Should look like:", nrarrcols)
        self.assertTrue(
            common.areArraysEqual(nrarrcols, tblcols),
            "Original array are retrieved doesn't match.",
        )

    def test_01a__iter__(self):
        tbl = self.h5file.create_table(
            "/", "test", self._TestTDescr, title=self._getMethodName()
        )
        tbl.append(self._testAData)

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        nrarr = np.array(testABuffer, dtype=tbl.description._v_nested_descr)
        row_num = 0
        for item in tbl.cols.Info.value:
            self.assertEqual(item, nrarr["Info"]["value"][row_num])
            row_num += 1
        self.assertEqual(row_num, len(nrarr))


class ColsNoReopen(ColsTestCase):
    reopen = 0


class ColsReopen(ColsTestCase):
    reopen = 1


class Nested(tb.IsDescription):
    uid = tb.IntCol(pos=1)
    value = tb.FloatCol(pos=2)


class A_Candidate(tb.IsDescription):
    nested1 = Nested()
    nested2 = Nested()


class B_Candidate(tb.IsDescription):
    nested1 = Nested
    nested2 = Nested


class C_Candidate(tb.IsDescription):
    nested1 = Nested()
    nested2 = Nested


Dnested = {
    "uid": tb.IntCol(pos=1),
    "value": tb.FloatCol(pos=2),
}

D_Candidate = {
    "nested1": Dnested,
    "nested2": Dnested,
}

E_Candidate = {
    "nested1": Nested,
    "nested2": Dnested,
}

F_Candidate = {
    "nested1": Nested(),
    "nested2": Dnested,
}

# Checking several nested columns declared in the same way


class SameNestedTestCase(common.TempFileMixin, common.PyTablesTestCase):
    correct_names = [
        "",  # The root of columns
        "nested1",
        "nested1/uid",
        "nested1/value",
        "nested2",
        "nested2/uid",
        "nested2/value",
    ]

    def test01a(self):
        """Checking same nested columns (instance flavor)."""

        tbl = self.h5file.create_table(
            "/", "test", A_Candidate, title=self._getMethodName()
        )

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        names = [
            col._v_pathname for col in tbl.description._f_walk(type="All")
        ]
        if common.verbose:
            print("Pathnames of columns:", names)
            print("Should look like:", self.correct_names)
        self.assertEqual(
            names, self.correct_names, "Column nested names doesn't match."
        )

    def test01b(self):
        """Checking same nested columns (class flavor)."""

        tbl = self.h5file.create_table(
            "/", "test", B_Candidate, title=self._getMethodName()
        )

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        names = [
            col._v_pathname for col in tbl.description._f_walk(type="All")
        ]
        if common.verbose:
            print("Pathnames of columns:", names)
            print("Should look like:", self.correct_names)
        self.assertEqual(
            names, self.correct_names, "Column nested names doesn't match."
        )

    def test01c(self):
        """Checking same nested columns (mixed instance/class flavor)."""

        tbl = self.h5file.create_table(
            "/", "test", C_Candidate, title=self._getMethodName()
        )

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        names = [
            col._v_pathname for col in tbl.description._f_walk(type="All")
        ]
        if common.verbose:
            print("Pathnames of columns:", names)
            print("Should look like:", self.correct_names)
        self.assertEqual(
            names, self.correct_names, "Column nested names doesn't match."
        )

    def test01d(self):
        """Checking same nested columns (dictionary flavor)."""

        tbl = self.h5file.create_table(
            "/", "test", D_Candidate, title=self._getMethodName()
        )

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        names = [
            col._v_pathname for col in tbl.description._f_walk(type="All")
        ]
        if common.verbose:
            print("Pathnames of columns:", names)
            print("Should look like:", self.correct_names)
        self.assertEqual(
            names, self.correct_names, "Column nested names doesn't match."
        )

    def test01e(self):
        """Checking same nested columns (mixed dictionary/class flavor)."""

        tbl = self.h5file.create_table(
            "/", "test", E_Candidate, title=self._getMethodName()
        )

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        names = [
            col._v_pathname for col in tbl.description._f_walk(type="All")
        ]
        if common.verbose:
            print("Pathnames of columns:", names)
            print("Should look like:", self.correct_names)
        self.assertEqual(
            names, self.correct_names, "Column nested names doesn't match."
        )

    def test01f(self):
        """Checking same nested columns (mixed dictionary/instance flavor)."""

        tbl = self.h5file.create_table(
            "/", "test", F_Candidate, title=self._getMethodName()
        )

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test

        names = [
            col._v_pathname for col in tbl.description._f_walk(type="All")
        ]
        if common.verbose:
            print("Pathnames of columns:", names)
            print("Should look like:", self.correct_names)
        self.assertEqual(
            names, self.correct_names, "Column nested names doesn't match."
        )

    def test02a(self):
        """Indexing two simple columns under the same nested column."""

        desc = {"nested": {"i1": tb.Int32Col(), "i2": tb.Int32Col()}}

        i1 = "nested/i1"
        i2 = "nested/i2"
        tbl = self.h5file.create_table(
            "/", "test", desc, title=self._getMethodName()
        )

        row = tbl.row
        for i in range(1000):
            row[i1] = i
            row[i2] = i * 2
            row.append()
        tbl.flush()

        cols = {
            "i1": tbl.cols.nested.i1,
            "i2": tbl.cols.nested.i2,
        }
        cols["i1"].create_index()
        cols["i2"].create_index()

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test
            # Redefine the cols dictionary
            cols = {
                "i1": tbl.cols.nested.i1,
                "i2": tbl.cols.nested.i2,
            }

        i1res = [r[i1] for r in tbl.where("i1 < 10", cols)]
        i2res = [r[i2] for r in tbl.where("i2 < 10", cols)]

        if common.verbose:
            print("Retrieved values (i1):", i1res)
            print("Should look like:", list(range(10)))
            print("Retrieved values (i2):", i2res)
            print("Should look like:", list(range(0, 10, 2)))

        self.assertEqual(
            i1res,
            list(range(10)),
            "Select for nested column (i1) doesn't match.",
        )
        self.assertEqual(
            i2res,
            list(range(0, 10, 2)),
            "Select for nested column (i2) doesn't match.",
        )

    def test02b(self):
        """Indexing two simple columns under the same (very) nested column."""

        desc = {
            "nested1": {
                "nested2": {
                    "nested3": {"i1": tb.Int32Col(), "i2": tb.Int32Col()}
                }
            }
        }

        i1 = "nested1/nested2/nested3/i1"
        i2 = "nested1/nested2/nested3/i2"

        tbl = self.h5file.create_table(
            "/", "test", desc, title=self._getMethodName()
        )

        row = tbl.row
        for i in range(1000):
            row[i1] = i
            row[i2] = i * 2
            row.append()
        tbl.flush()

        cols = {
            "i1": tbl.cols.nested1.nested2.nested3.i1,
            "i2": tbl.cols.nested1.nested2.nested3.i2,
        }
        cols["i1"].create_index()
        cols["i2"].create_index()

        if self.reopen:
            self._reopen()
            tbl = self.h5file.root.test
            # Redefine the cols dictionary
            cols = {
                "i1": tbl.cols.nested1.nested2.nested3.i1,
                "i2": tbl.cols.nested1.nested2.nested3.i2,
            }

        i1res = [r[i1] for r in tbl.where("i1 < 10", cols)]
        i2res = [r[i2] for r in tbl.where("i2 < 10", cols)]

        if common.verbose:
            print("Retrieved values (i1):", i1res)
            print("Should look like:", list(range(10)))
            print("Retrieved values (i2):", i2res)
            print("Should look like:", list(range(0, 10, 2)))

        self.assertEqual(
            i1res,
            list(range(10)),
            "Select for nested column (i1) doesn't match.",
        )
        self.assertEqual(
            i2res,
            list(range(0, 10, 2)),
            "Select for nested column (i2) doesn't match.",
        )


class SameNestedNoReopen(SameNestedTestCase):
    reopen = 0


class SameNestedReopen(SameNestedTestCase):
    reopen = 1


class NestedTypesWithGaps(common.TestFileMixin, common.PyTablesTestCase):
    h5fname = common.test_filename("nested-type-with-gaps.h5")

    correct_descr = f"""{{
  "float": Float32Col(shape=(), dflt={np.float32(0.0)!r}, pos=0),
  "compound": {{
    "char": Int8Col(shape=(), dflt={np.int8(0)!r}, pos=0),
    "double": Float64Col(shape=(), dflt={np.float64(0.0)!r}, pos=1)}}}}"""

    def test01(self):
        """Opening a table with nested types with gaps."""

        tbl = self.h5file.get_node("/nestedtype")
        type_descr = repr(tbl.description)
        if common.verbose:
            print("Type size with gaps:", tbl.description._v_itemsize)
            print("And should be: 16")
            print("Representation of the nested type:\n", type_descr)
            print("And should be:\n", self.correct_descr)
            print("Here are the offsets: ", tbl.description._v_offsets)

        self.assertEqual(tbl.description._v_itemsize, 16)
        self.assertEqual(type_descr, self.correct_descr)

        if common.verbose:
            print("Great!  Nested types with gaps recognized correctly.")


def suite():
    """Return a test suite consisting of all the test cases in the module."""

    theSuite = common.unittest.TestSuite()
    niter = 1
    # common.heavy = 1  # uncomment this only for testing purposes

    for i in range(niter):
        theSuite.addTest(common.make_suite(DescriptionTestCase))
        theSuite.addTest(common.make_suite(CreateTestCase))
        theSuite.addTest(common.make_suite(WriteNoReopen))
        theSuite.addTest(common.make_suite(WriteReopen))
        theSuite.addTest(common.make_suite(ColsNoReopen))
        theSuite.addTest(common.make_suite(ColsReopen))
        theSuite.addTest(common.make_suite(ReadNoReopen))
        theSuite.addTest(common.make_suite(ReadReopen))
        theSuite.addTest(common.make_suite(SameNestedNoReopen))
        theSuite.addTest(common.make_suite(SameNestedReopen))
        theSuite.addTest(common.make_suite(NestedTypesWithGaps))

    return theSuite


if __name__ == "__main__":
    common.parse_argv(sys.argv)
    common.print_versions()
    common.unittest.main(defaultTest="suite")
