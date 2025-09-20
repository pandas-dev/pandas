import shutil
import tempfile
import warnings
from pathlib import Path

import numpy as np

import tables as tb
from tables.tests import common


# Check read Tables from pytables version 0.8
class BackCompatTablesTestCase(common.PyTablesTestCase):
    def test01_readTable(self):
        """Checking backward compatibility of old formats of tables."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01_readTable..." % self.__class__.__name__)

        # Create an instance of an HDF5 Table
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            h5file = tb.open_file(common.test_filename(self.h5fname), "r")

        try:
            table = h5file.get_node("/tuple0")

            # Read the 100 records
            result = [rec["var2"] for rec in table]
            if common.verbose:
                print("Nrows in", table._v_pathname, ":", table.nrows)
                print("Last record in table ==>", table[-1])
                print("Total selected records in table ==> ", len(result))

            self.assertEqual(len(result), 100)
        finally:
            h5file.close()


@common.unittest.skipIf(not common.lzo_avail, "lzo not available")
class Table2_1LZO(BackCompatTablesTestCase):
    # pytables 0.8.x versions and after
    h5fname = "Table2_1_lzo_nrv2e_shuffle.h5"


@common.unittest.skipIf(not common.lzo_avail, "lzo not available")
class Tables_LZO1(BackCompatTablesTestCase):
    h5fname = "Tables_lzo1.h5"  # files compressed with LZO1


@common.unittest.skipIf(not common.lzo_avail, "lzo not available")
class Tables_LZO1_shuffle(BackCompatTablesTestCase):
    # files compressed with LZO1 and shuffle
    h5fname = "Tables_lzo1_shuffle.h5"


@common.unittest.skipIf(not common.lzo_avail, "lzo not available")
class Tables_LZO2(BackCompatTablesTestCase):
    h5fname = "Tables_lzo2.h5"  # files compressed with LZO2


@common.unittest.skipIf(not common.lzo_avail, "lzo not available")
class Tables_LZO2_shuffle(BackCompatTablesTestCase):
    # files compressed with LZO2 and shuffle
    h5fname = "Tables_lzo2_shuffle.h5"


# Check read attributes from PyTables >= 1.0 properly
class BackCompatAttrsTestCase(common.TestFileMixin, common.PyTablesTestCase):
    FILENAME = "zerodim-attrs-%s.h5"

    def setUp(self):
        self.h5fname = common.test_filename(self.FILENAME % self.format)
        super().setUp()

    def test01_readAttr(self):
        """Checking backward compatibility of old formats for attributes."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01_readAttr..." % self.__class__.__name__)

        # Read old formats
        a = self.h5file.get_node("/a")
        scalar = np.array(1, dtype="int32")
        vector = np.array([1], dtype="int32")
        if self.format == "1.3":
            self.assertTrue(common.allequal(a.attrs.arrdim1, vector))
            self.assertTrue(common.allequal(a.attrs.arrscalar, scalar))
            self.assertEqual(a.attrs.pythonscalar, 1)
        elif self.format == "1.4":
            self.assertTrue(common.allequal(a.attrs.arrdim1, vector))
            self.assertTrue(common.allequal(a.attrs.arrscalar, scalar))
            self.assertTrue(common.allequal(a.attrs.pythonscalar, scalar))


class Attrs_1_3(BackCompatAttrsTestCase):
    format = "1.3"  # pytables 1.0.x versions and earlier


class Attrs_1_4(BackCompatAttrsTestCase):
    format = "1.4"  # pytables 1.1.x versions and later


class VLArrayTestCase(common.TestFileMixin, common.PyTablesTestCase):
    h5fname = common.test_filename("flavored_vlarrays-format1.6.h5")

    def test01_backCompat(self):
        """Checking backward compatibility with old flavors of VLArray."""

        # Check that we can read the contents without problems (nor warnings!)
        vlarray1 = self.h5file.root.vlarray1
        self.assertEqual(vlarray1.flavor, "numeric")
        vlarray2 = self.h5file.root.vlarray2
        self.assertEqual(vlarray2.flavor, "python")
        self.assertEqual(vlarray2[1], [b"5", b"6", b"77"])


# Make sure that 1.x files with TimeXX types continue to be readable
# and that its byteorder is correctly retrieved.
class TimeTestCase(common.TestFileMixin, common.PyTablesTestCase):
    # Open a PYTABLES_FORMAT_VERSION=1.x file
    h5fname = common.test_filename("time-table-vlarray-1_x.h5")

    def test00_table(self):
        """Checking backward compatibility with old TimeXX types (tables)."""

        # Check that we can read the contents without problems (nor warnings!)
        table = self.h5file.root.table
        self.assertEqual(table.byteorder, "little")

    def test01_vlarray(self):
        """Checking backward compatibility with old TimeXX types (vlarrays)."""

        # Check that we can read the contents without problems (nor warnings!)
        vlarray4 = self.h5file.root.vlarray4
        self.assertEqual(vlarray4.byteorder, "little")
        vlarray8 = self.h5file.root.vlarray4
        self.assertEqual(vlarray8.byteorder, "little")


class OldFlavorsTestCase01(common.PyTablesTestCase):
    close = False

    # numeric
    def test01_open(self):
        """Checking opening of (X)Array (old 'numeric' flavor)"""

        # Open the HDF5 with old numeric flavor
        h5fname = common.test_filename("oldflavor_numeric.h5")
        with tb.open_file(h5fname) as h5file:

            # Assert other properties in array
            self.assertEqual(h5file.root.array1.flavor, "numeric")
            self.assertEqual(h5file.root.array2.flavor, "python")
            self.assertEqual(h5file.root.carray1.flavor, "numeric")
            self.assertEqual(h5file.root.carray2.flavor, "python")
            self.assertEqual(h5file.root.vlarray1.flavor, "numeric")
            self.assertEqual(h5file.root.vlarray2.flavor, "python")

    def test02_copy(self):
        """Checking (X)Array.copy() method ('numetic' flavor)"""

        srcfile = common.test_filename("oldflavor_numeric.h5")
        tmpfile = tempfile.mktemp(".h5")
        shutil.copy(srcfile, tmpfile)
        try:
            # Open the HDF5 with old numeric flavor
            with tb.open_file(tmpfile, "r+") as h5file:
                # Copy to another location
                self.assertWarns(
                    tb.exceptions.FlavorWarning,
                    h5file.root.array1.copy,
                    "/",
                    "array1copy",
                )
                h5file.root.array2.copy("/", "array2copy")
                h5file.root.carray1.copy("/", "carray1copy")
                h5file.root.carray2.copy("/", "carray2copy")
                h5file.root.vlarray1.copy("/", "vlarray1copy")
                h5file.root.vlarray2.copy("/", "vlarray2copy")

                if self.close:
                    h5file.close()
                    h5file = tb.open_file(tmpfile)
                else:
                    h5file.flush()

                # Assert other properties in array
                self.assertEqual(h5file.root.array1copy.flavor, "numeric")
                self.assertEqual(h5file.root.array2copy.flavor, "python")
                self.assertEqual(h5file.root.carray1copy.flavor, "numeric")
                self.assertEqual(h5file.root.carray2copy.flavor, "python")
                self.assertEqual(h5file.root.vlarray1copy.flavor, "numeric")
                self.assertEqual(h5file.root.vlarray2copy.flavor, "python")
        finally:
            Path(tmpfile).unlink()


class OldFlavorsTestCase02(common.PyTablesTestCase):
    close = True


def suite():
    theSuite = common.unittest.TestSuite()
    niter = 1

    for n in range(niter):
        theSuite.addTest(common.make_suite(VLArrayTestCase))
        theSuite.addTest(common.make_suite(TimeTestCase))
        theSuite.addTest(common.make_suite(OldFlavorsTestCase01))
        theSuite.addTest(common.make_suite(OldFlavorsTestCase02))
        theSuite.addTest(common.make_suite(Table2_1LZO))
        theSuite.addTest(common.make_suite(Tables_LZO1))
        theSuite.addTest(common.make_suite(Tables_LZO1_shuffle))
        theSuite.addTest(common.make_suite(Tables_LZO2))
        theSuite.addTest(common.make_suite(Tables_LZO2_shuffle))

    return theSuite


if __name__ == "__main__":
    import sys

    common.parse_argv(sys.argv)
    common.print_versions()
    common.unittest.main(defaultTest="suite")
