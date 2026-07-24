import sys

import numpy as np

import tables as tb
from tables.tests import common


class LargeTable(tb.IsDescription):
    time = tb.Int32Col()


class BasicTestCase(common.TempFileMixin, common.PyTablesTestCase):
    # file  = "test.h5"
    open_mode = "w"
    title = "This is the table title"
    dim1, dim2, dim3 = 24, 721, 1440
    nrows = dim1 * dim2 * dim3  # rows for a day
    chunkshape = nrows
    complib = "blosc2"  # default

    def setUp(self):
        super().setUp()

        # Create an instance of an HDF5 Table
        self.populateFile()
        self.h5file.close()

    def populateFile(self):
        group = self.h5file.root
        table = self.h5file.create_table(
            group,
            "table",
            LargeTable,
            "Large table",
            tb.Filters(complevel=1, complib=self.complib),
            chunkshape=self.chunkshape,
        )

        # Structured NumPy buffer for every day
        self.day_block = day_block = np.empty(self.nrows, dtype=table.dtype)
        day_block["time"] = np.arange(self.nrows)

        # Append groups of rows ("days") so that we have more than 2**31
        # (see https://github.com/PyTables/PyTables/issues/995)
        self.ndays = ndays = 90
        self.assertTrue(ndays * self.nrows > 2**31)
        if common.verbose:
            print(f"Writing {ndays} days...")
        for day in range(ndays):
            table.append(day_block)
        table.flush()

    def test00_values(self):
        """Check that written values are correct."""

        self.h5file = tb.open_file(self.h5fname)
        table = self.h5file.root.table
        nrows = self.nrows
        day_block = self.day_block
        if common.verbose:
            print(f"Checking {self.ndays} days...")
        for nday in range(self.ndays):
            day_block2 = table[nday * nrows : (nday + 1) * nrows]
            self.assertEqual(
                np.sum(day_block2["time"] == day_block["time"]),
                nrows,
                f"Values differ in day {nday}",
            )


@common.unittest.skipIf(
    not common.blosc_avail, "BLOSC compression library not available"
)
class BloscTestCase(BasicTestCase):
    title = "Blosc table"
    complib = "blosc"


@common.unittest.skipIf(
    not common.blosc2_avail, "BLOSC2 compression library not available"
)
class Blosc2TestCase(BasicTestCase):
    title = "Blosc2 table"
    complib = "blosc2"


class ZlibTestCase(BasicTestCase):
    title = "Zlib table"
    complib = "zlib"


def suite():
    theSuite = common.unittest.TestSuite()
    niter = 1
    # common.heavy = 1  # Uncomment this only for testing purposes

    for n in range(niter):
        theSuite.addTest(common.make_suite(BloscTestCase))
        theSuite.addTest(common.make_suite(Blosc2TestCase))
        if common.heavy:
            theSuite.addTest(common.make_suite(ZlibTestCase))

    return theSuite


if __name__ == "__main__":
    common.parse_argv(sys.argv)
    common.print_versions()
    common.unittest.main(defaultTest="suite")
