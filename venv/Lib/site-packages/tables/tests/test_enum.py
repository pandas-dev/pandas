"""Test module for enumerated types under PyTables."""

import operator
import itertools

import numpy as np

import tables as tb
from tables.tests import common


class CreateColTestCase(common.PyTablesTestCase):
    """Test creating enumerated column descriptions."""

    def _createCol(self, enum, dflt, base="uint32", shape=()):
        """Create and check an enumerated column description."""

        enumcol = tb.EnumCol(enum, dflt, base=base, shape=shape)
        sameEnum = tb.Enum(enum)
        self.assertEqual(enumcol.type, "enum")
        self.assertEqual(enumcol.dtype.base.name, enumcol.base.type)
        # To avoid 'LongInt' vs 'Int' issues
        # self.assertEqual(enumcol.dflt, sameEnum[dflt])
        self.assertEqual(int(enumcol.dflt), int(sameEnum[dflt]))
        self.assertEqual(enumcol.dtype.shape, shape)
        self.assertEqual(enumcol.enum, sameEnum)

    def test00a_validFromEnum(self):
        """Describing an enumerated column from an enumeration."""

        colors = tb.Enum(["red", "green", "blue"])
        self._createCol(colors, "red")

    def test00b_validFromDict(self):
        """Describing an enumerated column from a dictionary."""

        colors = {"red": 4, "green": 2, "blue": 1}
        self._createCol(colors, "red")

    def test00c_validFromList(self):
        """Describing an enumerated column from a list."""

        colors = ["red", "green", "blue"]
        self._createCol(colors, "red")

    def test00d_invalidFromType(self):
        """Describing an enumerated column from an invalid object."""

        colors = 123
        self.assertRaises(TypeError, self._createCol, colors, "red")

    def test01_invalidDflt(self):
        """Describing an enumerated column with an invalid default object."""

        colors = {"red": 4, "green": 2, "blue": 1}
        self.assertRaises(KeyError, self._createCol, colors, "black")

    def test02a_validDtypeBroader(self):
        """Describing an enumerated column with a broader type."""

        colors = {"red": 4, "green": 2, "blue": 1}
        self._createCol(colors, "red", "int64")

    def test02b_invalidDtypeTooNarrow(self):
        """Describing an enumerated column with a too narrow type."""

        colors = ["e%d" % i for i in range(300)]
        self.assertRaises(TypeError, self._createCol, colors, "e0", "uint8")

    def test03a_validShapeMD(self):
        """Describing an enumerated column with multidimensional shape."""

        colors = ["red", "green", "blue"]
        self._createCol(colors, "red", shape=(2,))

    def test04a_validReprEnum(self):
        """Checking the string representation of an enumeration."""

        colors = tb.Enum(["red", "green", "blue"])
        enumcol = tb.EnumCol(colors, "red", base="uint32", shape=())

        # needed due to "Hash randomization" (default on python 3.3)
        template = (
            "EnumCol(enum=Enum({%s}), dflt='red', base=UInt32Atom(shape=(), "
            f"dflt={np.uint32(0)!r}), shape=(), pos=None)"
        )
        permitations = [
            template % ", ".join(items)
            for items in itertools.permutations(
                ("'blue': 2", "'green': 1", "'red': 0")
            )
        ]
        self.assertIn(repr(enumcol), permitations)

    def test99a_nonIntEnum(self):
        """Describing an enumerated column of floats (not implemented)."""

        colors = {"red": 1.0}
        self.assertRaises(
            NotImplementedError,
            self._createCol,
            colors,
            "red",
            base=tb.FloatAtom(),
        )

    def test99b_nonIntDtype(self):
        """Describing an enumerated column encoded as floats.

        (not implemented).

        """

        colors = ["red", "green", "blue"]
        self.assertRaises(
            NotImplementedError, self._createCol, colors, "red", "float64"
        )

    def test99b_nonScalarEnum(self):
        """Describing an enumerated column of non-scalars (not implemented)."""

        colors = {"red": (1, 2, 3)}
        self.assertRaises(
            NotImplementedError,
            self._createCol,
            colors,
            "red",
            base=tb.IntAtom(shape=3),
        )


class CreateAtomTestCase(common.PyTablesTestCase):
    """Test creating enumerated atoms."""

    def _createAtom(self, enum, dflt, base="uint32", shape=()):
        """Create and check an enumerated atom."""

        enumatom = tb.EnumAtom(enum, dflt, base=base, shape=shape)
        sameEnum = tb.Enum(enum)
        self.assertEqual(enumatom.type, "enum")
        self.assertEqual(enumatom.dtype.base.name, enumatom.base.type)
        self.assertEqual(enumatom.shape, shape)
        self.assertEqual(enumatom.enum, sameEnum)

    def test00a_validFromEnum(self):
        """Describing an enumerated atom from an enumeration."""

        colors = tb.Enum(["red", "green", "blue"])
        self._createAtom(colors, "red")

    def test00b_validFromDict(self):
        """Describing an enumerated atom from a dictionary."""

        colors = {"red": 4, "green": 2, "blue": 1}
        self._createAtom(colors, "red")

    def test00c_validFromList(self):
        """Describing an enumerated atom from a list."""

        colors = ["red", "green", "blue"]
        self._createAtom(colors, "red")

    def test00d_invalidFromType(self):
        """Describing an enumerated atom from an invalid object."""

        colors = 123
        self.assertRaises(TypeError, self._createAtom, colors, "red")

    def test02a_validDtypeBroader(self):
        """Describing an enumerated atom with a broader type."""

        colors = {"red": 4, "green": 2, "blue": 1}
        self._createAtom(colors, "red", base="int64")

    def test02b_invalidDtypeTooNarrow(self):
        """Describing an enumerated atom with a too narrow type."""

        colors = ["e%d" % i for i in range(300)]
        self.assertRaises(TypeError, self._createAtom, colors, "red", "uint8")

    def test03a_validShapeMD(self):
        """Describing an enumerated atom with multidimensional shape."""

        colors = ["red", "green", "blue"]
        self._createAtom(colors, "red", shape=(2,))

    def test99a_nonIntEnum(self):
        """Describing an enumerated atom of floats (not implemented)."""

        colors = {"red": 1.0}
        self.assertRaises(
            NotImplementedError,
            self._createAtom,
            colors,
            "red",
            base=tb.FloatAtom(),
        )

    def test99b_nonIntDtype(self):
        """Describing an enumerated atom encoded as a float.

        (not implemented).

        """

        colors = ["red", "green", "blue"]
        self.assertRaises(
            NotImplementedError, self._createAtom, colors, "red", "float64"
        )

    def test99b_nonScalarEnum(self):
        """Describing an enumerated atom of non-scalars (not implemented)."""

        colors = {"red": (1, 2, 3)}
        self.assertRaises(
            NotImplementedError,
            self._createAtom,
            colors,
            "red",
            base=tb.IntAtom(shape=3),
        )


class EnumTableTestCase(common.TempFileMixin, common.PyTablesTestCase):
    """Test tables with enumerated columns."""

    enum = tb.Enum({"red": 4, "green": 2, "blue": 1, "black": 0})
    defaultName = "black"
    valueInEnum = enum.red
    valueOutOfEnum = 1234
    enumType = "uint16"

    def _description(self, shape=()):
        class TestDescription(tb.IsDescription):
            rid = tb.IntCol(pos=0)
            rcolor = tb.EnumCol(
                self.enum,
                self.defaultName,
                base=self.enumType,
                shape=shape,
                pos=1,
            )

        return TestDescription

    def test00a_reopen(self):
        """Reopening a file with tables using enumerated data."""

        self.h5file.create_table(
            "/", "test", self._description(), title=self._getMethodName()
        )

        self._reopen()

        self.assertEqual(
            self.h5file.root.test.get_enum("rcolor"),
            self.enum,
            "Enumerated type was not restored correctly from disk.",
        )

    def test00b_reopenMD(self):
        """Reopening a file with tables using enumerated multi-dimensional
        data."""

        self.h5file.create_table(
            "/", "test", self._description((2,)), title=self._getMethodName()
        )

        self._reopen()

        self.assertEqual(
            self.h5file.root.test.get_enum("rcolor"),
            self.enum,
            "Enumerated type was not restored correctly from disk.",
        )

    def test01_rowAppend(self):
        """Appending enumerated values using ``row.append()``."""

        tbl = self.h5file.create_table(
            "/", "test", self._description(), title=self._getMethodName()
        )

        appended = [(10, self.valueInEnum), (20, self.valueOutOfEnum)]

        row = tbl.row

        row["rid"] = appended[0][0]
        row["rcolor"] = appended[0][1]
        row.append()

        row["rid"] = appended[1][0]
        self.assertRaises(
            ValueError, operator.setitem, row, "rcolor", appended[1][1]
        )

        tbl.flush()
        tbl.flavor = "python"
        read = tbl.read()
        common.verbosePrint(
            "* appended value: %s\n"
            "* read value: %s\n" % (appended[:-1], read)
        )
        self.assertEqual(
            appended[:-1], read, "Written and read values differ."
        )

    def test02_append(self):
        """Appending enumerated values using ``table.append()``."""

        tbl = self.h5file.create_table(
            "/", "test", self._description(), title=self._getMethodName()
        )

        appended = [(10, self.valueInEnum), (20, self.valueOutOfEnum)]

        tbl.append(appended)
        tbl.flush()
        tbl.flavor = "python"
        read = tbl.read()
        common.verbosePrint(
            "* appended value: %s\n" "* read value: %s\n" % (appended, read)
        )
        self.assertEqual(appended, read, "Written and read values differ.")

    def test03_setitem(self):
        """Changing enumerated values using ``table.__setitem__()``."""

        tbl = self.h5file.create_table(
            "/", "test", self._description(), title=self._getMethodName()
        )

        appended = [(10, self.valueInEnum), (20, self.valueInEnum)]
        tbl.append(appended)

        written = [(10, self.valueInEnum), (20, self.valueOutOfEnum)]
        tbl[:] = written
        tbl.flavor = "python"
        read = tbl.read()
        common.verbosePrint(
            "* written value: %s\n" "* read value: %s\n" % (written, read)
        )
        self.assertEqual(written, read, "Written and read values differ.")

    def test04_multidim(self):
        """Appending multi-dimensional enumerated data."""

        tbl = self.h5file.create_table(
            "/", "test", self._description((2,)), title=self._getMethodName()
        )

        appended = [
            (10, (self.valueInEnum, self.valueOutOfEnum)),
            (20, (self.valueInEnum, self.valueOutOfEnum)),
        ]

        row = tbl.row
        row["rid"] = appended[0][0]
        self.assertRaises(
            ValueError, operator.setitem, row, "rcolor", appended[0][1]
        )

        tbl.append(appended)
        tbl.flush()
        tbl.flavor = "python"
        read = tbl.read()
        for x_appended, x_read in zip(appended, read):
            self.assertEqual(
                x_appended[0], x_read[0], "Written and read values differ."
            )
            self.assertEqual(
                x_appended[1][0],
                x_read[1][0],
                "Written and read values differ.",
            )
            self.assertEqual(
                x_appended[1][1],
                x_read[1][1],
                "Written and read values differ.",
            )

    def test05_where(self):
        """Searching enumerated data."""

        tbl = self.h5file.create_table(
            "/", "test", self._description(), title=self._getMethodName()
        )

        appended = [
            (10, self.valueInEnum),
            (20, self.valueInEnum),
            (30, self.valueOutOfEnum),
        ]
        tbl.append(appended)
        tbl.flush()

        searched = [
            (row["rid"], row["rcolor"])
            for row in tbl.where("rcolor == v", {"v": self.valueInEnum})
        ]
        common.verbosePrint(
            "* ``valueInEnum``: %s\n"
            "* ``rcolor`` column: ``%s``\n"
            "* ``searched``: %s\n"
            "* Should look like: %s\n"
            % (self.valueInEnum, tbl.cols.rcolor, searched, appended[:-1])
        )
        self.assertEqual(
            searched, appended[:-1], "Search returned incorrect results."
        )


class EnumEArrayTestCase(common.TempFileMixin, common.PyTablesTestCase):
    """Test extendable arrays of enumerated values."""

    enum = tb.Enum({"red": 4, "green": 2, "blue": 1, "black": 0})
    valueInEnum = enum.red
    valueOutOfEnum = 1234
    enumType = "uint16"

    def _atom(self, shape=()):
        return tb.EnumAtom(self.enum, "red", base=self.enumType, shape=shape)

    def test00a_reopen(self):
        """Reopening a file with extendable arrays using enumerated data."""

        self.h5file.create_earray(
            "/", "test", self._atom(), shape=(0,), title=self._getMethodName()
        )
        self.h5file.root.test.flavor = "python"

        self._reopen()

        self.assertEqual(
            self.h5file.root.test.get_enum(),
            self.enum,
            "Enumerated type was not restored correctly from disk.",
        )

    def test00b_reopenMD(self):
        """Reopening a file with extendable arrays using enumerated
        multi-dimensional data."""

        self.h5file.create_earray(
            "/",
            "test",
            self._atom(),
            shape=(0, 2),
            title=self._getMethodName(),
        )
        self.h5file.root.test.flavor = "python"

        self._reopen()

        self.assertEqual(
            self.h5file.root.test.get_enum(),
            self.enum,
            "Enumerated type was not restored correctly from disk.",
        )

    def test_enum_default_persistence_red(self):
        dflt = "red"
        atom = tb.EnumAtom(self.enum, dflt, base=self.enumType, shape=())

        self.h5file.create_earray(
            "/", "test", atom, shape=(0,), title=self._getMethodName()
        )
        self._reopen()

        self.assertEqual(
            self.h5file.root.test.get_enum(),
            self.enum,
            "Enumerated type was not restored correctly from disk.",
        )

        self.assertEqual(
            self.h5file.root.test.atom.dflt,
            self.enum[dflt],
            "The default value of enumerated type was not restored correctly "
            "from disk.",
        )

    def test_enum_default_persistence_green(self):
        dflt = "green"
        atom = tb.EnumAtom(self.enum, dflt, base=self.enumType, shape=())

        self.h5file.create_earray(
            "/", "test", atom, shape=(0,), title=self._getMethodName()
        )
        self._reopen()

        self.assertEqual(
            self.h5file.root.test.get_enum(),
            self.enum,
            "Enumerated type was not restored correctly from disk.",
        )

        self.assertEqual(
            self.h5file.root.test.atom.dflt,
            self.enum[dflt],
            "The default value of enumerated type was not restored correctly "
            "from disk.",
        )

    def test_enum_default_persistence_blue(self):
        dflt = "blue"
        atom = tb.EnumAtom(self.enum, dflt, base=self.enumType, shape=())

        self.h5file.create_earray(
            "/", "test", atom, shape=(0,), title=self._getMethodName()
        )
        self._reopen()

        self.assertEqual(
            self.h5file.root.test.get_enum(),
            self.enum,
            "Enumerated type was not restored correctly from disk.",
        )

        self.assertEqual(
            self.h5file.root.test.atom.dflt,
            self.enum[dflt],
            "The default value of enumerated type was not restored correctly "
            "from disk.",
        )

    def test_enum_default_persistence_black(self):
        dflt = "black"
        atom = tb.EnumAtom(self.enum, dflt, base=self.enumType, shape=())

        self.h5file.create_earray(
            "/", "test", atom, shape=(0,), title=self._getMethodName()
        )
        self._reopen()

        self.assertEqual(
            self.h5file.root.test.get_enum(),
            self.enum,
            "Enumerated type was not restored correctly from disk.",
        )

        self.assertEqual(
            self.h5file.root.test.atom.dflt,
            self.enum[dflt],
            "The default value of enumerated type was not restored correctly "
            "from disk.",
        )

    def test01_append(self):
        """Appending scalar elements of enumerated values."""

        earr = self.h5file.create_earray(
            "/", "test", self._atom(), shape=(0,), title=self._getMethodName()
        )
        earr.flavor = "python"

        appended = [self.valueInEnum, self.valueOutOfEnum]

        earr.append(appended)
        earr.flush()
        read = earr.read()
        self.assertEqual(appended, read, "Written and read values differ.")

    def test02_appendMD(self):
        """Appending multi-dimensional elements of enumerated values."""

        earr = self.h5file.create_earray(
            "/",
            "test",
            self._atom(),
            shape=(0, 2),
            title=self._getMethodName(),
        )
        earr.flavor = "python"

        appended = [
            [self.valueInEnum, self.valueOutOfEnum],
            [self.valueInEnum, self.valueOutOfEnum],
        ]

        earr.append(appended)
        earr.flush()
        read = earr.read()
        self.assertEqual(appended, read, "Written and read values differ.")

    def test03_setitem(self):
        """Changing enumerated values using ``earray.__setitem__()``."""

        earr = self.h5file.create_earray(
            "/", "test", self._atom(), shape=(0,), title=self._getMethodName()
        )
        earr.flavor = "python"

        appended = (self.valueInEnum, self.valueInEnum)
        earr.append(appended)

        written = [self.valueInEnum, self.valueOutOfEnum]
        earr[:] = written
        read = earr.read()
        self.assertEqual(written, read, "Written and read values differ.")


class EnumVLArrayTestCase(common.TempFileMixin, common.PyTablesTestCase):
    """Test variable-length arrays of enumerated values."""

    enum = tb.Enum({"red": 4, "green": 2, "blue": 1, "black": 0})
    valueInEnum = enum.red
    valueOutOfEnum = 1234
    enumType = "uint16"

    def _atom(self, shape=()):
        return tb.EnumAtom(self.enum, "red", base=self.enumType, shape=shape)

    def test00a_reopen(self):
        """Reopening a file with variable-length arrays using
        enumerated data."""

        self.h5file.create_vlarray(
            "/", "test", self._atom(), title=self._getMethodName()
        )
        self.h5file.root.test.flavor = "python"

        self._reopen()

        self.assertEqual(
            self.h5file.root.test.get_enum(),
            self.enum,
            "Enumerated type was not restored correctly from disk.",
        )

    def test00b_reopenMD(self):
        """Reopening a file with variable-length arrays using enumerated
        multi-dimensional data."""

        self.h5file.create_vlarray(
            "/", "test", self._atom((2,)), title=self._getMethodName()
        )
        self.h5file.root.test.flavor = "python"

        self._reopen()

        self.assertEqual(
            self.h5file.root.test.get_enum(),
            self.enum,
            "Enumerated type was not restored correctly from disk.",
        )

    def test01_append(self):
        """Appending scalar elements of enumerated values."""

        vlarr = self.h5file.create_vlarray(
            "/", "test", self._atom(), title=self._getMethodName()
        )
        vlarr.flavor = "python"

        appended = [
            [
                self.valueInEnum,
            ],
            [self.valueInEnum, self.valueOutOfEnum],
        ]

        vlarr.append(appended[0])
        vlarr.append(appended[1])
        vlarr.flush()
        read = vlarr.read()
        common.verbosePrint(
            "* appended value: %s\n" "* read value: %s\n" % (appended, read)
        )
        self.assertEqual(appended, read, "Written and read values differ.")

    def test02_appendMD(self):
        """Appending multi-dimensional elements of enumerated values."""

        vlarr = self.h5file.create_vlarray(
            "/", "test", self._atom((2,)), title=self._getMethodName()
        )
        vlarr.flavor = "python"

        appended = [
            [
                [self.valueInEnum, self.valueInEnum],
            ],
            [
                [self.valueInEnum, self.valueOutOfEnum],
                [self.valueInEnum, self.valueInEnum],
            ],
        ]

        vlarr.append(appended[0])
        vlarr.append(appended[1])
        vlarr.flush()
        read = vlarr.read()
        common.verbosePrint(
            "* appended value: %s\n" "* read value: %s\n" % (appended, read)
        )
        self.assertEqual(appended, read, "Written and read values differ.")

    def test03_setitem(self):
        """Changing enumerated values using ``vlarray.__setitem__()``."""

        vlarr = self.h5file.create_vlarray(
            "/", "test", self._atom(), title=self._getMethodName()
        )
        vlarr.flavor = "python"

        appended = (self.valueInEnum, self.valueInEnum)
        vlarr.append(appended)

        written = [self.valueInEnum, self.valueOutOfEnum]
        vlarr[0] = written
        read = vlarr.read()
        common.verbosePrint(
            "* written value: %s\n" "* read value: %s\n" % (written, read)
        )
        self.assertEqual(written, read[0], "Written and read values differ.")


def suite():
    """Return a test suite consisting of all the test cases in the module."""

    # These two are for including Enum's doctests here.
    import doctest

    theSuite = common.unittest.TestSuite()
    niter = 1

    # theSuite.addTest(make_suite(EnumTableTestCase))
    for i in range(niter):
        theSuite.addTest(doctest.DocTestSuite(tb.misc.enum))
        theSuite.addTest(common.make_suite(CreateColTestCase))
        theSuite.addTest(common.make_suite(CreateAtomTestCase))
        theSuite.addTest(common.make_suite(EnumTableTestCase))
        theSuite.addTest(common.make_suite(EnumEArrayTestCase))
        theSuite.addTest(common.make_suite(EnumVLArrayTestCase))

    return theSuite


if __name__ == "__main__":
    import sys

    common.parse_argv(sys.argv)
    common.print_versions()
    common.unittest.main(defaultTest="suite")
