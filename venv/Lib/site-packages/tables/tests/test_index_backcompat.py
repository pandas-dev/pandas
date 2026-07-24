from tables.tests import common


# Check indexes from PyTables version 2.0
class IndexesTestCase(common.TestFileMixin, common.PyTablesTestCase):

    def setUp(self):
        super().setUp()
        self.table1 = self.h5file.root.table1
        self.table2 = self.h5file.root.table2
        self.il = 0
        self.sl = self.table1.cols.var1.index.slicesize

    def test00_version(self):
        """Checking index version."""

        t1var1 = self.table1.cols.var1
        if "2_0" in str(self.h5fname):
            self.assertEqual(t1var1.index._v_version, "2.0")
        elif "2_1" in str(self.h5fname):
            self.assertEqual(t1var1.index._v_version, "2.1")

    def test01_string(self):
        """Checking string indexes."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01_string..." % self.__class__.__name__)

        table1 = self.table1
        table2 = self.table2

        # Convert the limits to the appropriate type
        il = str(self.il).encode("ascii")
        sl = str(self.sl).encode("ascii")

        # Do some selections and check the results
        # First selection
        t1var1 = table1.cols.var1
        self.assertIsNotNone(t1var1)
        results1 = [
            p["var1"] for p in table1.where("(il<=t1var1)&(t1var1<=sl)")
        ]
        results2 = [p["var1"] for p in table2 if il <= p["var1"] <= sl]
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Should look like:", results2)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test02_bool(self):
        """Checking bool indexes."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02_bool..." % self.__class__.__name__)

        table1 = self.table1
        table2 = self.table2

        # Do some selections and check the results
        t1var2 = table1.cols.var2
        self.assertIsNotNone(t1var2)
        results1 = [p["var2"] for p in table1.where("t1var2 == True")]
        results2 = [p["var2"] for p in table2 if p["var2"] is True]
        if common.verbose:
            print("Selection results (index):", results1)
            print("Should look like:", results2)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test03_int(self):
        """Checking int indexes."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03_int..." % self.__class__.__name__)

        table1 = self.table1
        table2 = self.table2

        # Convert the limits to the appropriate type
        il = int(self.il)
        sl = int(self.sl)

        # Do some selections and check the results
        t1col = table1.cols.var3
        self.assertIsNotNone(t1col)

        # First selection
        results1 = [p["var3"] for p in table1.where("(il<=t1col)&(t1col<=sl)")]
        results2 = [p["var3"] for p in table2 if il <= p["var3"] <= sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test04_float(self):
        """Checking float indexes."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test04_float..." % self.__class__.__name__)

        table1 = self.table1
        table2 = self.table2

        # Convert the limits to the appropriate type
        il = float(self.il)
        sl = float(self.sl)

        # Do some selections and check the results
        t1col = table1.cols.var4
        self.assertIsNotNone(t1col)

        # First selection
        results1 = [p["var4"] for p in table1.where("(il<=t1col)&(t1col<=sl)")]
        results2 = [p["var4"] for p in table2 if il <= p["var4"] <= sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1.sort(), results2.sort())


# Check indexes from PyTables version 2.0
class Indexes2_0TestCase(IndexesTestCase):
    h5fname = common.test_filename("indexes_2_0.h5")


# Check indexes from PyTables version 2.1
class Indexes2_1TestCase(IndexesTestCase):
    h5fname = common.test_filename("indexes_2_1.h5")


def suite():
    theSuite = common.unittest.TestSuite()
    niter = 1

    for n in range(niter):
        theSuite.addTest(common.make_suite(Indexes2_0TestCase))
        theSuite.addTest(common.make_suite(Indexes2_1TestCase))

    return theSuite


if __name__ == "__main__":
    import sys

    common.parse_argv(sys.argv)
    common.print_versions()
    common.unittest.main(defaultTest="suite")
