import random
import tempfile
from pathlib import Path

import numpy as np

import tables as tb
from tables.tests import common

# An alias for frozenset
fzset = frozenset

# To make the values in the tests reproducible
random.seed(19)

# Sensible parameters for indexing with small blocksizes
small_blocksizes = (16, 8, 4, 2)  # The smaller set of parameters...
# The size for medium indexes
minRowIndex = 1000


class Small(tb.IsDescription):
    var1 = tb.StringCol(itemsize=4, dflt=b"")
    var2 = tb.BoolCol(dflt=0)
    var3 = tb.IntCol(dflt=0)
    var4 = tb.FloatCol(dflt=0)


class SelectValuesTestCase(common.TempFileMixin, common.PyTablesTestCase):
    compress = 1
    complib = "zlib"
    shuffle = 1
    fletcher32 = 0
    chunkshape = 10
    buffersize = 0
    random = 0
    values = None
    reopen = False

    def setUp(self):
        super().setUp()

        # Create an instance of an HDF5 Table
        if common.verbose:
            print("Checking index kind-->", self.kind)
        self.rootgroup = self.h5file.root
        self.populateFile()

    def populateFile(self):
        # Set a seed for the random generator if needed.
        # This is useful when one needs reproducible results.
        if self.random and hasattr(self, "seed"):
            random.seed(self.seed)
        group = self.rootgroup
        # Create a table
        title = "This is the IndexArray title"
        filters = tb.Filters(
            complevel=self.compress,
            complib=self.complib,
            shuffle=self.shuffle,
            fletcher32=self.fletcher32,
        )
        table1 = self.h5file.create_table(
            group,
            "table1",
            Small,
            title,
            filters,
            self.nrows,
            chunkshape=(self.chunkshape,),
        )
        table2 = self.h5file.create_table(
            group,
            "table2",
            Small,
            title,
            filters,
            self.nrows,
            chunkshape=(self.chunkshape,),
        )
        count = 0
        for i in range(0, self.nrows, self.nrep):
            for j in range(self.nrep):
                if self.random:
                    k = random.randrange(self.nrows)
                elif self.values is not None:
                    lenvalues = len(self.values)
                    if i >= lenvalues:
                        i %= lenvalues
                    k = self.values[i]
                else:
                    k = i
                bk = str(k).encode("ascii")
                table1.row["var1"] = bk
                table2.row["var1"] = bk
                table1.row["var2"] = k % 2
                table2.row["var2"] = k % 2
                table1.row["var3"] = k
                table2.row["var3"] = k
                table1.row["var4"] = float(self.nrows - k - 1)
                table2.row["var4"] = float(self.nrows - k - 1)
                table1.row.append()
                table2.row.append()
                count += 1
        table1.flush()
        table2.flush()
        if self.buffersize:
            # Change the buffersize by default
            table1.nrowsinbuf = self.buffersize
            # Make sure nrowsinbuf is a multiple of chunkshape
            table1.nrowsinbuf -= table1.nrowsinbuf % self.chunkshape
        # Index all entries:
        for col in table1.colinstances.values():
            indexrows = col.create_index(
                kind=self.kind, _blocksizes=self.blocksizes
            )
        if common.verbose:
            print("Number of written rows:", table1.nrows)
            print("Number of indexed rows:", indexrows)

        if self.reopen:
            self._reopen(mode="a")  # flavor changes
            self.table1 = self.h5file.root.table1
            self.table2 = self.h5file.root.table1

    def test01a(self):
        """Checking selecting values from an Index (string flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01a..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Convert the limits to the appropriate type
        il = str(self.il).encode("ascii")
        sl = str(self.sl).encode("ascii")

        # Do some selections and check the results
        # First selection
        t1var1 = table1.cols.var1
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

        # Second selection
        t1var1 = table1.cols.var1
        results1 = [
            p["var1"] for p in table1.where("(il<=t1var1)&(t1var1<sl)")
        ]
        results2 = [p["var1"] for p in table2 if il <= p["var1"] < sl]
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Third selection
        t1var1 = table1.cols.var1
        results1 = [
            p["var1"] for p in table1.where("(il<t1var1)&(t1var1<=sl)")
        ]
        results2 = [p["var1"] for p in table2 if il < p["var1"] <= sl]
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Forth selection
        t1var1 = table1.cols.var1
        self.assertIsNotNone(t1var1)
        results1 = [p["var1"] for p in table1.where("(il<t1var1)&(t1var1<sl)")]
        results2 = [p["var1"] for p in table2 if il < p["var1"] < sl]
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test01b(self):
        """Checking selecting values from an Index (string flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01b..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Convert the limits to the appropriate type
        # il = str(self.il).encode('ascii')
        sl = str(self.sl).encode("ascii")

        # Do some selections and check the results
        # First selection
        t1var1 = table1.cols.var1
        results1 = [p["var1"] for p in table1.where("t1var1 < sl")]
        results2 = [p["var1"] for p in table2 if p["var1"] < sl]
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Second selection
        t1var1 = table1.cols.var1
        results1 = [p["var1"] for p in table1.where("t1var1 <= sl")]
        results2 = [p["var1"] for p in table2 if p["var1"] <= sl]
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Third selection
        t1var1 = table1.cols.var1
        results1 = [p["var1"] for p in table1.where("t1var1 > sl")]
        results2 = [p["var1"] for p in table2 if p["var1"] > sl]
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Fourth selection
        t1var1 = table1.cols.var1
        self.assertIsNotNone(t1var1)
        results1 = [p["var1"] for p in table1.where("t1var1 >= sl")]
        results2 = [p["var1"] for p in table2 if p["var1"] >= sl]
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test02a(self):
        """Checking selecting values from an Index (bool flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02a..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Do some selections and check the results
        t1var2 = table1.cols.var2
        self.assertIsNotNone(t1var2)
        results1 = [p["var2"] for p in table1.where("t1var2 == True")]
        results2 = [p["var2"] for p in table2 if p["var2"] is True]
        if common.verbose:
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test02b(self):
        """Checking selecting values from an Index (bool flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02b..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Do some selections and check the results
        t1var2 = table1.cols.var2
        self.assertIsNotNone(t1var2)
        results1 = [p["var2"] for p in table1.where("t1var2 == False")]
        results2 = [p["var2"] for p in table2 if p["var2"] is False]
        if common.verbose:
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test03a(self):
        """Checking selecting values from an Index (int flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03a..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

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

        # Second selection
        results1 = [p["var3"] for p in table1.where("(il<=t1col)&(t1col<sl)")]
        results2 = [p["var3"] for p in table2 if il <= p["var3"] < sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Third selection
        results1 = [p["var3"] for p in table1.where("(il<t1col)&(t1col<=sl)")]
        results2 = [p["var3"] for p in table2 if il < p["var3"] <= sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Fourth selection
        results1 = [p["var3"] for p in table1.where("(il<t1col)&(t1col<sl)")]
        results2 = [p["var3"] for p in table2 if il < p["var3"] < sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test03b(self):
        """Checking selecting values from an Index (int flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03b..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Convert the limits to the appropriate type
        # il = int(self.il)
        sl = int(self.sl)

        # Do some selections and check the results
        t1col = table1.cols.var3
        self.assertIsNotNone(t1col)

        # First selection
        results1 = [p["var3"] for p in table1.where("t1col < sl")]
        results2 = [p["var3"] for p in table2 if p["var3"] < sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Second selection
        results1 = [p["var3"] for p in table1.where("t1col <= sl")]
        results2 = [p["var3"] for p in table2 if p["var3"] <= sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Third selection
        results1 = [p["var3"] for p in table1.where("t1col > sl")]
        results2 = [p["var3"] for p in table2 if p["var3"] > sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Fourth selection
        results1 = [p["var3"] for p in table1.where("t1col >= sl")]
        results2 = [p["var3"] for p in table2 if p["var3"] >= sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test03c(self):
        """Checking selecting values from an Index (long flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03c..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Convert the limits to the appropriate type
        # il = long(self.il)
        sl = int(self.sl)

        # Do some selections and check the results
        t1col = table1.cols.var3
        self.assertIsNotNone(t1col)

        # First selection
        results1 = [p["var3"] for p in table1.where("t1col < sl")]
        results2 = [p["var3"] for p in table2 if p["var3"] < sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Second selection
        results1 = [p["var3"] for p in table1.where("t1col <= sl")]
        results2 = [p["var3"] for p in table2 if p["var3"] <= sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Third selection
        results1 = [p["var3"] for p in table1.where("t1col > sl")]
        results2 = [p["var3"] for p in table2 if p["var3"] > sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Fourth selection
        results1 = [p["var3"] for p in table1.where("t1col >= sl")]
        results2 = [p["var3"] for p in table2 if p["var3"] >= sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test03d(self):
        """Checking selecting values from an Index (long and int flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03d..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Convert the limits to the appropriate type
        # il = int(self.il)
        sl = int(self.sl)

        # Do some selections and check the results
        t1col = table1.cols.var3
        self.assertIsNotNone(t1col)

        # First selection
        results1 = [p["var3"] for p in table1.where("t1col < sl")]
        results2 = [p["var3"] for p in table2 if p["var3"] < sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Second selection
        results1 = [p["var3"] for p in table1.where("t1col <= sl")]
        results2 = [p["var3"] for p in table2 if p["var3"] <= sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Third selection
        results1 = [p["var3"] for p in table1.where("t1col > sl")]
        results2 = [p["var3"] for p in table2 if p["var3"] > sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Fourth selection
        results1 = [p["var3"] for p in table1.where("t1col >= sl")]
        results2 = [p["var3"] for p in table2 if p["var3"] >= sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test04a(self):
        """Checking selecting values from an Index (float flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test04a..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

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

        # Second selection
        results1 = [p["var4"] for p in table1.where("(il<=t1col)&(t1col<sl)")]
        results2 = [p["var4"] for p in table2 if il <= p["var4"] < sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Third selection
        results1 = [p["var4"] for p in table1.where("(il<t1col)&(t1col<=sl)")]
        results2 = [p["var4"] for p in table2 if il < p["var4"] <= sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        self.assertEqual(results1, results2)

        # Fourth selection
        results1 = [p["var4"] for p in table1.where("(il<t1col)&(t1col<sl)")]
        results2 = [p["var4"] for p in table2 if il < p["var4"] < sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test04b(self):
        """Checking selecting values from an Index (float flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test04b..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Convert the limits to the appropriate type
        # il = float(self.il)
        sl = float(self.sl)

        # Do some selections and check the results
        t1col = table1.cols.var4
        self.assertIsNotNone(t1col)

        # First selection
        results1 = [p["var4"] for p in table1.where("t1col < sl")]
        results2 = [p["var4"] for p in table2 if p["var4"] < sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Second selection
        results1 = [p["var4"] for p in table1.where("t1col <= sl")]
        results2 = [p["var4"] for p in table2 if p["var4"] <= sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Third selection
        results1 = [p["var4"] for p in table1.where("t1col > sl")]
        results2 = [p["var4"] for p in table2 if p["var4"] > sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Fourth selection
        results1 = [p["var4"] for p in table1.where("t1col >= sl")]
        results2 = [p["var4"] for p in table2 if p["var4"] >= sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test05a(self):
        """Checking get_where_list & itersequence (string, python flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test05a..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Convert the limits to the appropriate type
        il = str(self.il).encode("ascii")
        sl = str(self.sl).encode("ascii")

        # Do some selections and check the results
        t1col = table1.cols.var1
        # First selection
        condition = "(il<=t1col)&(t1col<=sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        table1.flavor = "python"
        rowList1 = table1.get_where_list(condition)
        results1 = [p["var1"] for p in table1.itersequence(rowList1)]
        results2 = [p["var1"] for p in table2 if il <= p["var1"] <= sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1.sort(), results2.sort())

        # Second selection
        condition = "(il<=t1col)&(t1col<sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        table1.flavor = "python"
        rowList1 = table1.get_where_list(condition)
        results1 = [p["var1"] for p in table1.itersequence(rowList1)]
        results2 = [p["var1"] for p in table2 if il <= p["var1"] < sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Third selection
        condition = "(il<t1col)&(t1col<=sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        table1.flavor = "python"
        rowList1 = table1.get_where_list(condition)
        results1 = [p["var1"] for p in table1.itersequence(rowList1)]
        results2 = [p["var1"] for p in table2 if il < p["var1"] <= sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        self.assertEqual(results1, results2)

        # Fourth selection
        condition = "(il<t1col)&(t1col<sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        table1.flavor = "python"
        rowList1 = table1.get_where_list(condition)
        results1 = [p["var1"] for p in table1.itersequence(rowList1)]
        results2 = [p["var1"] for p in table2 if il < p["var1"] < sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test05b(self):
        """Checking get_where_list & itersequence (numpy string lims & python
        flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test05b..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Convert the limits to the appropriate type
        # il = np.bytes_(self.il)
        sl = np.bytes_(self.sl)

        # Do some selections and check the results
        t1col = table1.cols.var1

        # First selection
        condition = "t1col<sl"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        table1.flavor = "python"
        rowList1 = table1.get_where_list(condition)
        results1 = [p["var1"] for p in table1.itersequence(rowList1)]
        results2 = [p["var1"] for p in table2 if p["var1"] < sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Second selection
        condition = "t1col<=sl"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        rowList1 = table1.get_where_list(condition)
        results1 = [p["var1"] for p in table1.itersequence(rowList1)]
        results2 = [p["var1"] for p in table2 if p["var1"] <= sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Third selection
        condition = "t1col>sl"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        rowList1 = table1.get_where_list(condition)
        results1 = [p["var1"] for p in table1.itersequence(rowList1)]
        results2 = [p["var1"] for p in table2 if p["var1"] > sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Fourth selection
        condition = "t1col>=sl"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        rowList1 = table1.get_where_list(condition)
        results1 = [p["var1"] for p in table1.itersequence(rowList1)]
        results2 = [p["var1"] for p in table2 if p["var1"] >= sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test06a(self):
        """Checking get_where_list & itersequence (bool flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test06a..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Do some selections and check the results
        t1var2 = table1.cols.var2
        condition = "t1var2==True"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1var2.pathname])
        )
        table1.flavor = "python"
        rowList1 = table1.get_where_list(condition)
        results1 = [p["var2"] for p in table1.itersequence(rowList1)]
        results2 = [p["var2"] for p in table2 if p["var2"] is True]
        if common.verbose:
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test06b(self):
        """Checking get_where_list & itersequence (numpy bool limits &
        flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test06b..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Do some selections and check the results
        t1var2 = table1.cols.var2
        false = np.bool_(False)
        self.assertFalse(false)  # silence pyflakes
        condition = "t1var2==false"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1var2.pathname])
        )
        table1.flavor = "python"
        rowList1 = table1.get_where_list(condition)
        results1 = [p["var2"] for p in table1.itersequence(rowList1)]
        results2 = [p["var2"] for p in table2 if p["var2"] is False]
        if common.verbose:
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test07a(self):
        """Checking get_where_list & itersequence (int flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test07a..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Convert the limits to the appropriate type
        il = int(self.il)
        sl = int(self.sl)

        # Do some selections and check the results
        t1col = table1.cols.var3
        # First selection
        condition = "(il<=t1col)&(t1col<=sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        table1.flavor = "python"
        rowList1 = table1.get_where_list(condition)
        results1 = [p["var3"] for p in table1.itersequence(rowList1)]
        results2 = [p["var3"] for p in table2 if il <= p["var3"] <= sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1.sort(), results2.sort())

        # Second selection
        condition = "(il<=t1col)&(t1col<sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        table1.flavor = "python"
        rowList1 = table1.get_where_list(condition)
        results1 = [p["var3"] for p in table1.itersequence(rowList1)]
        results2 = [p["var3"] for p in table2 if il <= p["var3"] < sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Third selection
        condition = "(il<t1col)&(t1col<=sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        table1.flavor = "python"
        rowList1 = table1.get_where_list(condition)
        results1 = [p["var3"] for p in table1.itersequence(rowList1)]
        results2 = [p["var3"] for p in table2 if il < p["var3"] <= sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        self.assertEqual(results1, results2)

        # Fourth selection
        condition = "(il<t1col)&(t1col<sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        table1.flavor = "python"
        rowList1 = table1.get_where_list(condition)
        results1 = [p["var3"] for p in table1.itersequence(rowList1)]
        results2 = [p["var3"] for p in table2 if il < p["var3"] < sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test07b(self):
        """Checking get_where_list & itersequence (numpy int limits &
        flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test07b..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Convert the limits to the appropriate type
        # il = np.int32(self.il)
        sl = np.uint16(self.sl)

        # Do some selections and check the results
        t1col = table1.cols.var3

        # First selection
        condition = "t1col<sl"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        table1.flavor = "python"
        rowList1 = table1.get_where_list(condition)
        results1 = [p["var3"] for p in table1.itersequence(rowList1)]
        results2 = [p["var3"] for p in table2 if p["var3"] < sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Second selection
        condition = "t1col<=sl"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        rowList1 = table1.get_where_list(condition)
        results1 = [p["var3"] for p in table1.itersequence(rowList1)]
        results2 = [p["var3"] for p in table2 if p["var3"] <= sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Third selection
        condition = "t1col>sl"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        rowList1 = table1.get_where_list(condition)
        results1 = [p["var3"] for p in table1.itersequence(rowList1)]
        results2 = [p["var3"] for p in table2 if p["var3"] > sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Fourth selection
        condition = "t1col>=sl"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        rowList1 = table1.get_where_list(condition)
        results1 = [p["var3"] for p in table1.itersequence(rowList1)]
        results2 = [p["var3"] for p in table2 if p["var3"] >= sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test08a(self):
        """Checking get_where_list & itersequence (float flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test08a..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Convert the limits to the appropriate type
        il = float(self.il)
        sl = float(self.sl)

        # Do some selections and check the results
        t1col = table1.cols.var4
        # First selection
        condition = "(il<=t1col)&(t1col<=sl)"
        # results1 = [p["var4"] for p in table1.where(condition)]
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        table1.flavor = "python"
        rowList1 = table1.get_where_list(condition)
        results1 = [p["var4"] for p in table1.itersequence(rowList1)]
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

        # Second selection
        condition = "(il<=t1col)&(t1col<sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        table1.flavor = "python"
        rowList1 = table1.get_where_list(condition)
        results1 = [p["var4"] for p in table1.itersequence(rowList1)]
        results2 = [p["var4"] for p in table2 if il <= p["var4"] < sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Third selection
        condition = "(il<t1col)&(t1col<=sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        table1.flavor = "python"
        rowList1 = table1.get_where_list(condition)
        results1 = [p["var4"] for p in table1.itersequence(rowList1)]
        results2 = [p["var4"] for p in table2 if il < p["var4"] <= sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        self.assertEqual(results1, results2)

        # Fourth selection
        condition = "(il<t1col)&(t1col<sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        table1.flavor = "python"
        rowList1 = table1.get_where_list(condition)
        results1 = [p["var4"] for p in table1.itersequence(rowList1)]
        results2 = [p["var4"] for p in table2 if il < p["var4"] < sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test08b(self):
        """Checking get_where_list & itersequence (numpy float limits &
        flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test08b..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Convert the limits to the appropriate type
        # il = np.float32(self.il)
        sl = np.float64(self.sl)

        # Do some selections and check the results
        t1col = table1.cols.var4

        # First selection
        condition = "t1col<sl"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        rowList1 = table1.get_where_list(condition)
        results1 = [p["var4"] for p in table1.itersequence(rowList1)]
        results2 = [p["var4"] for p in table2 if p["var4"] < sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Second selection
        condition = "t1col<=sl"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        rowList1 = table1.get_where_list(condition)
        results1 = [p["var4"] for p in table1.itersequence(rowList1)]
        results2 = [p["var4"] for p in table2 if p["var4"] <= sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Third selection
        condition = "t1col>sl"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        rowList1 = table1.get_where_list(condition)
        results1 = [p["var4"] for p in table1.itersequence(rowList1)]
        results2 = [p["var4"] for p in table2 if p["var4"] > sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Fourth selection
        condition = "t1col>=sl"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        rowList1 = table1.get_where_list(condition)
        results1 = [p["var4"] for p in table1.itersequence(rowList1)]
        results2 = [p["var4"] for p in table2 if p["var4"] >= sl]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test09a(self):
        """Checking non-indexed where() (string flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test09a..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        table1._disable_indexing_in_queries()

        # Convert the limits to the appropriate type
        il = str(self.il).encode("ascii")
        sl = str(self.sl).encode("ascii")

        # Do some selections and check the results
        t1col = table1.cols.var1
        self.assertIsNotNone(t1col)

        # First selection
        condition = "t1col<=sl"
        self.assertTrue(not table1.will_query_use_indexing(condition))
        results1 = [
            p["var1"] for p in table1.where(condition, start=2, stop=10)
        ]
        results2 = [
            p["var1"] for p in table2.iterrows(2, 10) if p["var1"] <= sl
        ]
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Second selection
        condition = "(il<t1col)&(t1col<sl)"
        self.assertTrue(not table1.will_query_use_indexing(condition))
        results1 = [
            p["var1"]
            for p in table1.where(condition, start=2, stop=30, step=2)
        ]
        results2 = [
            p["var1"] for p in table2.iterrows(2, 30, 2) if il < p["var1"] < sl
        ]
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Third selection
        condition = "(il>t1col)&(t1col>sl)"
        self.assertTrue(not table1.will_query_use_indexing(condition))
        results1 = [
            p["var1"] for p in table1.where(condition, start=2, stop=-5)
        ]
        results2 = [
            p["var1"]
            for p in table2.iterrows(2, -5)  # Negative indices
            if (il > p["var1"] > sl)
        ]
        if common.verbose:
            print("Limits:", il, sl)
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # This selection to be commented out
        #         condition = 't1col>=sl'
        #         self.assertTrue(not table1.will_query_use_indexing(condition))
        #         results1 = [p['var1'] for p in table1.where(condition,start=2,
        #                                                     stop=-1,step=1)]
        #         results2 = [p["var1"] for p in table2.iterrows(2, -1, 1)
        #                     if p["var1"] >= sl]
        #         if verbose:
        #             print "Limit:", sl
        #             print "Selection results (in-kernel):", results1
        #             print "Should look like:", results2
        #             print "Length results:", len(results1)
        #             print "Should be:", len(results2)
        #         self.assertEqual(len(results1), len(results2))
        #         self.assertEqual(results1, results2)

        # Fourth selection
        # results1 = [p['var1'] for p in
        # table1.where(condition,start=2,stop=-1,step=3)]
        condition = "t1col>=sl"
        self.assertTrue(not table1.will_query_use_indexing(condition))
        results1 = [
            p["var1"]
            for p in table1.where(condition, start=2, stop=-1, step=3)
        ]
        results2 = [
            p["var1"] for p in table2.iterrows(2, -1, 3) if p["var1"] >= sl
        ]
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Re-enable the indexing in queries basically to unnail the
        # condition cache and not raising the performance warning
        # about some indexes being dirty
        table1._enable_indexing_in_queries()

    def test09b(self):
        """Checking non-indexed where() (float flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test09b..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        table1._disable_indexing_in_queries()

        # Convert the limits to the appropriate type
        il = float(self.il)
        sl = float(self.sl)

        # Do some selections and check the results
        t1col = table1.cols.var4
        self.assertIsNotNone(t1col)

        # First selection
        condition = "t1col<sl"
        self.assertTrue(not table1.will_query_use_indexing(condition))
        results1 = [
            p["var4"] for p in table1.where(condition, start=2, stop=5)
        ]
        results2 = [p["var4"] for p in table2.iterrows(2, 5) if p["var4"] < sl]
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Second selection
        condition = "(il<t1col)&(t1col<=sl)"
        self.assertTrue(not table1.will_query_use_indexing(condition))
        results1 = [
            p["var4"]
            for p in table1.where(condition, start=2, stop=-1, step=2)
        ]
        results2 = [
            p["var4"]
            for p in table2.iterrows(2, -1, 2)
            if il < p["var4"] <= sl
        ]
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Third selection
        condition = "(il<=t1col)&(t1col<=sl)"
        self.assertTrue(not table1.will_query_use_indexing(condition))
        results1 = [
            p["var4"] for p in table1.where(condition, start=2, stop=-5)
        ]
        results2 = [
            p["var4"]
            for p in table2.iterrows(2, -5)  # Negative indices
            if il <= p["var4"] <= sl
        ]
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Fourth selection
        condition = "t1col>=sl"
        self.assertTrue(not table1.will_query_use_indexing(condition))
        results1 = [
            p["var4"]
            for p in table1.where(condition, start=0, stop=-1, step=3)
        ]
        results2 = [
            p["var4"] for p in table2.iterrows(0, -1, 3) if p["var4"] >= sl
        ]
        if common.verbose:
            print("Limit:", sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Re-enable the indexing in queries basically to unnail the
        # condition cache and not raising the performance warning
        # about some indexes being dirty
        table1._enable_indexing_in_queries()

    def test09c(self):
        """Check non-indexed where() w/ ranges, changing step
        (string flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test09c..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        table1._disable_indexing_in_queries()

        # Convert the limits to the appropriate type
        il = str(self.il).encode("ascii")
        sl = str(self.sl).encode("ascii")

        # Do some selections and check the results
        t1col = table1.cols.var1
        self.assertIsNotNone(t1col)

        # First selection
        condition = "t1col>=sl"
        self.assertTrue(not table1.will_query_use_indexing(condition))
        results1 = [
            p["var1"]
            for p in table1.where(condition, start=2, stop=-1, step=3)
        ]
        results2 = [
            p["var1"] for p in table2.iterrows(2, -1, 3) if p["var1"] >= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Second selection
        condition = "t1col>=sl"
        self.assertTrue(not table1.will_query_use_indexing(condition))
        results1 = [
            p["var1"]
            for p in table1.where(condition, start=5, stop=-1, step=10)
        ]
        results2 = [
            p["var1"] for p in table2.iterrows(5, -1, 10) if p["var1"] >= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Third selection
        condition = "t1col>=sl"
        self.assertTrue(not table1.will_query_use_indexing(condition))
        results1 = [
            p["var1"]
            for p in table1.where(condition, start=5, stop=-3, step=11)
        ]
        results2 = [
            p["var1"] for p in table2.iterrows(5, -3, 11) if p["var1"] >= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Fourth selection
        condition = "t1col>=sl"
        self.assertTrue(not table1.will_query_use_indexing(condition))
        results1 = [
            p["var1"]
            for p in table1.where(condition, start=2, stop=-1, step=300)
        ]
        results2 = [
            p["var1"] for p in table2.iterrows(2, -1, 300) if p["var1"] >= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Re-enable the indexing in queries basically to unnail the
        # condition cache and not raising the performance warning
        # about some indexes being dirty
        table1._enable_indexing_in_queries()

    def test09d(self):
        """Checking non-indexed where() w/ ranges, changing step
        (int flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test09d..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        table1._disable_indexing_in_queries()

        # Convert the limits to the appropriate type
        il = int(self.il)
        sl = int(self.sl)

        # Do some selections and check the results
        t3col = table1.cols.var3
        self.assertIsNotNone(t3col)

        # First selection
        condition = "t3col>=sl"
        self.assertTrue(not table1.will_query_use_indexing(condition))
        results1 = [
            p["var3"]
            for p in table1.where(condition, start=2, stop=-1, step=3)
        ]
        results2 = [
            p["var3"] for p in table2.iterrows(2, -1, 3) if p["var3"] >= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Second selection
        condition = "t3col>=sl"
        self.assertTrue(not table1.will_query_use_indexing(condition))
        results1 = [
            p["var3"]
            for p in table1.where(condition, start=5, stop=-1, step=10)
        ]
        results2 = [
            p["var3"] for p in table2.iterrows(5, -1, 10) if p["var3"] >= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Third selection
        condition = "t3col>=sl"
        self.assertTrue(not table1.will_query_use_indexing(condition))
        results1 = [
            p["var3"]
            for p in table1.where(condition, start=5, stop=-3, step=11)
        ]
        results2 = [
            p["var3"] for p in table2.iterrows(5, -3, 11) if p["var3"] >= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Fourth selection
        condition = "t3col>=sl"
        self.assertTrue(not table1.will_query_use_indexing(condition))
        results1 = [
            p["var3"]
            for p in table1.where(condition, start=2, stop=-1, step=300)
        ]
        results2 = [
            p["var3"] for p in table2.iterrows(2, -1, 300) if p["var3"] >= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Re-enable the indexing in queries basically to unnail the
        # condition cache and not raising the performance warning
        # about some indexes being dirty
        table1._enable_indexing_in_queries()

    def test10a(self):
        """Checking indexed where() with ranges (string flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test10a..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Convert the limits to the appropriate type
        il = str(self.il).encode("ascii")
        sl = str(self.sl).encode("ascii")

        # Do some selections and check the results
        t1col = table1.cols.var1
        # First selection
        condition = "t1col<=sl"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        results1 = [
            p["var1"] for p in table1.where(condition, start=2, stop=10)
        ]
        results2 = [
            p["var1"] for p in table2.iterrows(2, 10) if p["var1"] <= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Second selection
        condition = "(il<=t1col)&(t1col<=sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        results1 = [
            p["var1"]
            for p in table1.where(condition, start=2, stop=30, step=1)
        ]
        results2 = [
            p["var1"]
            for p in table2.iterrows(2, 30, 1)
            if il <= p["var1"] <= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Repeat second selection (testing caches)
        condition = "(il<=t1col)&(t1col<=sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        results1 = [
            p["var1"]
            for p in table1.where(condition, start=2, stop=30, step=2)
        ]
        results2 = [
            p["var1"]
            for p in table2.iterrows(2, 30, 2)
            if il <= p["var1"] <= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Selection results (indexed):", results1)
            print("Should look like:", results2)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Third selection
        condition = "(il<t1col)&(t1col<sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        results1 = [
            p["var1"] for p in table1.where(condition, start=2, stop=-5)
        ]
        results2 = [
            p["var1"]
            for p in table2.iterrows(2, -5)  # Negative indices
            if (il < p["var1"] < sl)
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Fourth selection
        condition = "t1col>=sl"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        results1 = [
            p["var1"]
            for p in table1.where(condition, start=1, stop=-1, step=3)
        ]
        results2 = [
            p["var1"] for p in table2.iterrows(1, -1, 3) if p["var1"] >= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test10b(self):
        """Checking indexed where() with ranges (int flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test10b..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Convert the limits to the appropriate type
        il = int(self.il)
        sl = int(self.sl)

        # Do some selections and check the results
        t3col = table1.cols.var3
        # First selection
        condition = "t3col<=sl"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t3col.pathname])
        )
        results1 = [
            p["var3"] for p in table1.where(condition, start=2, stop=10)
        ]
        results2 = [
            p["var3"] for p in table2.iterrows(2, 10) if p["var3"] <= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Second selection
        condition = "(il<=t3col)&(t3col<=sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t3col.pathname])
        )
        results1 = [
            p["var3"]
            for p in table1.where(condition, start=2, stop=30, step=2)
        ]
        results2 = [
            p["var3"]
            for p in table2.iterrows(2, 30, 2)
            if il <= p["var3"] <= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Third selection
        condition = "(il<t3col)&(t3col<sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t3col.pathname])
        )
        results1 = [
            p["var3"] for p in table1.where(condition, start=2, stop=-5)
        ]
        results2 = [
            p["var3"]
            for p in table2.iterrows(2, -5)  # Negative indices
            if (il < p["var3"] < sl)
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Fourth selection
        condition = "t3col>=sl"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t3col.pathname])
        )
        results1 = [
            p["var3"]
            for p in table1.where(condition, start=1, stop=-1, step=3)
        ]
        results2 = [
            p["var3"] for p in table2.iterrows(1, -1, 3) if p["var3"] >= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test10c(self):
        """Checking indexed where() with ranges, changing step (string
        flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test10c..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Convert the limits to the appropriate type
        il = str(self.il).encode("ascii")
        sl = str(self.sl).encode("ascii")

        # Do some selections and check the results
        t1col = table1.cols.var1

        # First selection
        condition = "t1col>=sl"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        results1 = [
            p["var1"]
            for p in table1.where(condition, start=2, stop=-1, step=3)
        ]
        results2 = [
            p["var1"] for p in table2.iterrows(2, -1, 3) if p["var1"] >= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Second selection
        condition = "t1col>=sl"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        results1 = [
            p["var1"]
            for p in table1.where(condition, start=5, stop=-1, step=10)
        ]
        results2 = [
            p["var1"] for p in table2.iterrows(5, -1, 10) if p["var1"] >= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Third selection
        condition = "t1col>=sl"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        results1 = [
            p["var1"]
            for p in table1.where(condition, start=5, stop=-3, step=11)
        ]
        results2 = [
            p["var1"] for p in table2.iterrows(5, -3, 11) if p["var1"] >= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Fourth selection
        condition = "t1col>=sl"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        results1 = [
            p["var1"]
            for p in table1.where(condition, start=2, stop=-1, step=300)
        ]
        results2 = [
            p["var1"] for p in table2.iterrows(2, -1, 300) if p["var1"] >= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test10d(self):
        """Checking indexed where() with ranges, changing step (int flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test10d..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Convert the limits to the appropriate type
        il = int(self.il)
        sl = int(self.sl)

        # Do some selections and check the results
        t3col = table1.cols.var3

        # First selection
        condition = "t3col>=sl"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t3col.pathname])
        )
        results1 = [
            p["var3"]
            for p in table1.where(condition, start=2, stop=-1, step=3)
        ]
        results2 = [
            p["var3"] for p in table2.iterrows(2, -1, 3) if p["var3"] >= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Second selection
        condition = "t3col>=sl"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t3col.pathname])
        )
        results1 = [
            p["var3"]
            for p in table1.where(condition, start=5, stop=-1, step=10)
        ]
        results2 = [
            p["var3"] for p in table2.iterrows(5, -1, 10) if p["var3"] >= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Third selection
        condition = "t3col>=sl"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t3col.pathname])
        )
        results1 = [
            p["var3"]
            for p in table1.where(condition, start=5, stop=-3, step=11)
        ]
        results2 = [
            p["var3"] for p in table2.iterrows(5, -3, 11) if p["var3"] >= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Fourth selection
        condition = "t3col>=sl"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t3col.pathname])
        )
        results1 = [
            p["var3"]
            for p in table1.where(condition, start=2, stop=-1, step=300)
        ]
        results2 = [
            p["var3"] for p in table2.iterrows(2, -1, 300) if p["var3"] >= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test11a(self):
        """Checking selecting values from an Index via read_coordinates()"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test11a..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Convert the limits to the appropriate type
        il = str(self.il).encode("ascii")
        sl = str(self.sl).encode("ascii")

        # Do a selection and check the result
        t1var1 = table1.cols.var1
        condition = "(il<=t1var1)&(t1var1<=sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1var1.pathname])
        )
        coords1 = table1.get_where_list(condition)
        table1.flavor = "python"
        results1 = table1.read_coordinates(coords1, field="var1")
        results2 = [p["var1"] for p in table2 if il <= p["var1"] <= sl]
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test12a(self):
        """Checking selecting values after a Table.append() operation."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test12a..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Append more rows in already created indexes
        count = 0
        for i in range(0, self.nrows // 2, self.nrep):
            for j in range(self.nrep):
                if self.random:
                    k = random.randrange(self.nrows)
                elif self.values is not None:
                    lenvalues = len(self.values)
                    if i >= lenvalues:
                        i %= lenvalues
                    k = self.values[i]
                else:
                    k = i
                table1.row["var1"] = str(k)
                table2.row["var1"] = str(k)
                table1.row["var2"] = k % 2
                table2.row["var2"] = k % 2
                table1.row["var3"] = k
                table2.row["var3"] = k
                table1.row["var4"] = float(self.nrows - k - 1)
                table2.row["var4"] = float(self.nrows - k - 1)
                table1.row.append()
                table2.row.append()
                count += 1
        table1.flush()
        table2.flush()

        t1var1 = table1.cols.var1
        t1var2 = table1.cols.var2
        t1var3 = table1.cols.var3
        t1var4 = table1.cols.var4
        self.assertFalse(t1var1.index.dirty)
        self.assertFalse(t1var2.index.dirty)
        self.assertFalse(t1var3.index.dirty)
        self.assertFalse(t1var4.index.dirty)

        # Do some selections and check the results
        # First selection: string
        # Convert the limits to the appropriate type
        il = str(self.il).encode("ascii")
        sl = str(self.sl).encode("ascii")

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

        # Second selection: bool
        results1 = [p["var2"] for p in table1.where("t1var2 == True")]
        results2 = [p["var2"] for p in table2 if p["var2"] is True]
        t2var1_vals = [p["var1"] for p in table2]
        t2var2_vals = [p["var2"] for p in table2]
        msg = (
            f"Incorrect results for t1var2[n] == True where\n"
            f"t2var1_vals={repr(t2var1_vals)}\nt2var2_vals={repr(t2var2_vals)}\n"
            f"\n{results1=}\n{results2=}"
        )
        self.assertEqual(len(results1), len(results2), msg=msg)
        self.assertEqual(results1, results2, msg=msg)

        # Third selection: int
        # Convert the limits to the appropriate type
        il = int(self.il)
        sl = int(self.sl)

        t1var3 = table1.cols.var3
        results1 = [
            p["var3"] for p in table1.where("(il<=t1var3)&(t1var3<=sl)")
        ]
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

        # Fourth selection: float
        # Convert the limits to the appropriate type
        il = float(self.il)
        sl = float(self.sl)

        # Do some selections and check the results
        results1 = [
            p["var4"] for p in table1.where("(il<=t1var4)&(t1var4<=sl)")
        ]
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

    def test13a(self):
        """Checking repeated queries (checking caches)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test13a..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Convert the limits to the appropriate type
        il = str(self.il).encode("ascii")
        sl = str(self.sl).encode("ascii")

        # Do some selections and check the results
        t1col = table1.cols.var1
        condition = "(il<=t1col)&(t1col<=sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        results1 = [
            p["var1"]
            for p in table1.where(condition, start=2, stop=30, step=1)
        ]
        results2 = [
            p["var1"]
            for p in table2.iterrows(2, 30, 1)
            if il <= p["var1"] <= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Repeat the selection (testing caches)
        condition = "(il<=t1col)&(t1col<=sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        results1 = [
            p["var1"]
            for p in table1.where(condition, start=2, stop=30, step=2)
        ]
        results2 = [
            p["var1"]
            for p in table2.iterrows(2, 30, 2)
            if il <= p["var1"] <= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test13b(self):
        """Checking repeated queries, varying step (checking caches)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test13b..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Convert the limits to the appropriate type
        il = str(self.il).encode("ascii")
        sl = str(self.sl).encode("ascii")

        # Do some selections and check the results
        t1col = table1.cols.var1
        condition = "(il<=t1col)&(t1col<=sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        results1 = [
            p["var1"]
            for p in table1.where(condition, start=2, stop=30, step=1)
        ]
        results2 = [
            p["var1"]
            for p in table2.iterrows(2, 30, 1)
            if il <= p["var1"] <= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Repeat the selection (testing caches)
        condition = "(il<=t1col)&(t1col<=sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        results1 = [
            p["var1"]
            for p in table1.where(condition, start=2, stop=30, step=2)
        ]
        results2 = [
            p["var1"]
            for p in table2.iterrows(2, 30, 2)
            if il <= p["var1"] <= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test13c(self):
        """Checking repeated queries, varying start, stop, step."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test13c..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Convert the limits to the appropriate type
        il = str(self.il).encode("ascii")
        sl = str(self.sl).encode("ascii")

        # Do some selections and check the results
        t1col = table1.cols.var1
        condition = "(il<=t1col)&(t1col<=sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        results1 = [
            p["var1"] for p in table1.where(condition, start=0, stop=1, step=2)
        ]
        results2 = [
            p["var1"]
            for p in table2.iterrows(0, 1, 2)
            if il <= p["var1"] <= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Repeat the selection (testing caches)
        condition = "(il<=t1col)&(t1col<=sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        results1 = [
            p["var1"] for p in table1.where(condition, start=0, stop=5, step=1)
        ]
        results2 = [
            p["var1"]
            for p in table2.iterrows(0, 5, 1)
            if il <= p["var1"] <= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test13d(self):
        """Checking repeated queries, varying start, stop, step (another
        twist)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test13d..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Convert the limits to the appropriate type
        il = str(self.il).encode("ascii")
        sl = str(self.sl).encode("ascii")

        # Do some selections and check the results
        t1col = table1.cols.var1
        condition = "(il<=t1col)&(t1col<=sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        results1 = [
            p["var1"] for p in table1.where(condition, start=0, stop=1, step=1)
        ]
        results2 = [
            p["var1"]
            for p in table2.iterrows(0, 1, 1)
            if il <= p["var1"] <= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Repeat the selection (testing caches)
        condition = "(il<=t1col)&(t1col<=sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        results1 = [
            p["var1"] for p in table1.where(condition, start=0, stop=1, step=1)
        ]
        results2 = [
            p["var1"]
            for p in table2.iterrows(0, 1, 1)
            if il <= p["var1"] <= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test13e(self):
        """Checking repeated queries, with varying condition."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test13e..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Convert the limits to the appropriate type
        il = str(self.il).encode("ascii")
        sl = str(self.sl).encode("ascii")

        # Do some selections and check the results
        t1col = table1.cols.var1
        condition = "(il<=t1col)&(t1col<=sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        results1 = [
            p["var1"]
            for p in table1.where(condition, start=0, stop=10, step=1)
        ]
        results2 = [
            p["var1"]
            for p in table2.iterrows(0, 10, 1)
            if il <= p["var1"] <= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Repeat the selection with a more complex condition
        t2col = table1.cols.var2
        condition = "(il<=t1col)&(t1col<=sl)&(t2col==True)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname, t2col.pathname])
        )
        results1 = [
            p["var1"]
            for p in table1.where(condition, start=0, stop=10, step=1)
        ]
        results2 = [
            p["var1"]
            for p in table2.iterrows(0, 10, 1)
            if il <= p["var1"] <= sl and p["var2"] is True
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test13f(self):
        """Checking repeated queries, with varying condition."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test13f..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Remove indexes in var2 column
        table1.cols.var2.remove_index()
        table2.cols.var2.remove_index()

        # Convert the limits to the appropriate type
        il = str(self.il).encode("ascii")
        sl = str(self.sl).encode("ascii")

        # Do some selections and check the results
        t1col = table1.cols.var1
        t2col = table1.cols.var2
        self.assertIsNotNone(t2col)
        condition = "(il<=t1col)&(t1col<=sl)&(t2col==True)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        results1 = [
            p["var1"]
            for p in table1.where(condition, start=0, stop=10, step=1)
        ]
        results2 = [
            p["var1"]
            for p in table2.iterrows(0, 10, 1)
            if il <= p["var1"] <= sl and p["var2"] is True
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Repeat the selection with a simpler condition
        condition = "(il<=t1col)&(t1col<=sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        results1 = [
            p["var1"]
            for p in table1.where(condition, start=0, stop=10, step=1)
        ]
        results2 = [
            p["var1"]
            for p in table2.iterrows(0, 10, 1)
            if il <= p["var1"] <= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Repeat again with the original condition, but with a constant
        constant = True
        condition = "(il<=t1col)&(t1col<=sl)&(t2col==constant)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        results1 = [
            p["var1"]
            for p in table1.where(condition, start=0, stop=10, step=1)
        ]
        results2 = [
            p["var1"]
            for p in table2.iterrows(0, 10, 1)
            if il <= p["var1"] <= sl and p["var2"] == constant
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

    def test13g(self):
        """Checking repeated queries, with different limits."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test13g..." % self.__class__.__name__)

        table1 = self.h5file.root.table1
        table2 = self.h5file.root.table2

        # Convert the limits to the appropriate type
        il = str(self.il).encode("ascii")
        sl = str(self.sl).encode("ascii")

        # Do some selections and check the results
        t1col = table1.cols.var1
        condition = "(il<=t1col)&(t1col<=sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        results1 = [
            p["var1"]
            for p in table1.where(condition, start=0, stop=10, step=1)
        ]
        results2 = [
            p["var1"]
            for p in table2.iterrows(0, 10, 1)
            if il <= p["var1"] <= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)

        # Repeat the selection with different limits
        il, sl = (
            str(self.il + 1).encode("ascii"),
            str(self.sl - 2).encode("ascii"),
        )
        t2col = table1.cols.var2
        self.assertIsNotNone(t2col)
        condition = "(il<=t1col)&(t1col<=sl)"
        self.assertTrue(
            table1.will_query_use_indexing(condition)
            == fzset([t1col.pathname])
        )
        results1 = [
            p["var1"]
            for p in table1.where(condition, start=0, stop=10, step=1)
        ]
        results2 = [
            p["var1"]
            for p in table2.iterrows(0, 10, 1)
            if il <= p["var1"] <= sl
        ]
        # sort lists (indexing does not guarantee that rows are returned in
        # order)
        results1.sort()
        results2.sort()
        if common.verbose:
            print("Limits:", il, sl)
            print("Length results:", len(results1))
            print("Should be:", len(results2))
        self.assertEqual(len(results1), len(results2))
        self.assertEqual(results1, results2)


class SV1aTestCase(SelectValuesTestCase):
    blocksizes = small_blocksizes
    chunkshape = 1
    buffersize = 2
    ss = blocksizes[2]
    nrows = ss
    reopen = 0
    nrep = ss
    il = 0
    sl = ss


class SV1bTestCase(SV1aTestCase):
    blocksizes = tb.idxutils.calc_chunksize(minRowIndex, memlevel=1)
    chunkshape = blocksizes[2] // 2**9
    buffersize = chunkshape * 5


class SV2aTestCase(SelectValuesTestCase):
    blocksizes = small_blocksizes
    chunkshape = 2
    buffersize = 2
    ss = blocksizes[2]
    nrows = ss * 2 - 1
    reopen = 1
    nrep = 1
    il = 0
    sl = 2


class SV2bTestCase(SV2aTestCase):
    blocksizes = tb.idxutils.calc_chunksize(minRowIndex, memlevel=1)
    chunkshape = blocksizes[2] // 2**7
    buffersize = chunkshape * 20


class SV3aTestCase(SelectValuesTestCase):
    blocksizes = small_blocksizes
    chunkshape = 2
    buffersize = 3
    ss = blocksizes[2]
    nrows = ss * 5 - 1
    reopen = 1
    nrep = 3
    il = 0
    sl = 3


class SV3bTestCase(SV3aTestCase):
    blocksizes = tb.idxutils.calc_chunksize(minRowIndex, memlevel=1)
    #    chunkshape = 4
    #    buffersize = 16
    chunkshape = 3
    buffersize = 9


class SV4aTestCase(SelectValuesTestCase):
    blocksizes = small_blocksizes
    buffersize = 10
    ss = blocksizes[2]
    nrows = ss * 3
    reopen = 0
    nrep = 1
    # il = nrows-cs
    il = 0
    sl = nrows


class SV4bTestCase(SV4aTestCase):
    blocksizes = tb.idxutils.calc_chunksize(minRowIndex, memlevel=1)
    chunkshape = 500
    buffersize = 1000


class SV5aTestCase(SelectValuesTestCase):
    blocksizes = small_blocksizes
    ss = blocksizes[2]
    nrows = ss * 5
    reopen = 0
    nrep = 1
    il = 0
    sl = nrows


class SV5bTestCase(SV5aTestCase):
    blocksizes = tb.idxutils.calc_chunksize(minRowIndex, memlevel=1)


class SV6aTestCase(SelectValuesTestCase):
    blocksizes = small_blocksizes
    ss = blocksizes[2]
    nrows = ss * 5 + 1
    reopen = 0
    cs = blocksizes[3]
    nrep = cs + 1
    il = -1
    sl = nrows


class SV6bTestCase(SV6aTestCase):
    blocksizes = tb.idxutils.calc_chunksize(minRowIndex, memlevel=1)


class SV7aTestCase(SelectValuesTestCase):
    random = 1
    blocksizes = small_blocksizes
    ss = blocksizes[2]
    nrows = ss * 5 + 3
    reopen = 0
    cs = blocksizes[3]
    nrep = cs - 1
    il = -10
    sl = nrows


class SV7bTestCase(SV7aTestCase):
    blocksizes = tb.idxutils.calc_chunksize(minRowIndex, memlevel=1)


class SV8aTestCase(SelectValuesTestCase):
    random = 0
    chunkshape = 1
    blocksizes = small_blocksizes
    ss = blocksizes[2]
    nrows = ss * 5 - 3
    reopen = 0
    cs = blocksizes[3]
    nrep = cs - 1
    il = 10
    sl = nrows - 10


class SV8bTestCase(SV8aTestCase):
    random = 0
    blocksizes = tb.idxutils.calc_chunksize(minRowIndex, memlevel=1)


class SV9aTestCase(SelectValuesTestCase):
    random = 1
    blocksizes = small_blocksizes
    ss = blocksizes[2]
    nrows = ss * 5 + 11
    reopen = 0
    cs = blocksizes[3]
    nrep = cs - 1
    il = 10
    sl = nrows - 10


class SV9bTestCase(SV9aTestCase):
    blocksizes = tb.idxutils.calc_chunksize(minRowIndex, memlevel=1)


class SV10aTestCase(SelectValuesTestCase):
    random = 1
    blocksizes = small_blocksizes
    chunkshape = 1
    buffersize = 1
    ss = blocksizes[2]
    nrows = ss
    reopen = 0
    nrep = ss
    il = 0
    sl = ss


class SV10bTestCase(SV10aTestCase):
    blocksizes = tb.idxutils.calc_chunksize(minRowIndex, memlevel=1)
    chunkshape = 5
    buffersize = 6


class SV11aTestCase(SelectValuesTestCase):
    # This checks a special case that failed. It was discovered in a
    # random test above (SV10a). It is explicitely put here as a way
    # to always check that specific case.
    values = [1, 7, 6, 7, 0, 7, 4, 4, 9, 5]
    blocksizes = small_blocksizes
    chunkshape = 1
    buffersize = 1
    ss = blocksizes[2]
    nrows = ss
    reopen = 0
    nrep = ss
    il = 0
    sl = ss


class SV11bTestCase(SelectValuesTestCase):
    # This checks a special case that failed. It was discovered in a
    # random test above (SV10a). It is explicitely put here as a way
    # to always check that specific case.
    values = [1, 7, 6, 7, 0, 7, 4, 4, 9, 5]
    chunkshape = 2
    buffersize = 2
    blocksizes = tb.idxutils.calc_chunksize(minRowIndex, memlevel=1)
    ss = blocksizes[2]
    nrows = ss
    reopen = 0
    nrep = ss
    il = 0
    sl = ss


class SV12aTestCase(SelectValuesTestCase):
    # This checks a special case that failed. It was discovered in a
    # random test above (SV10b). It is explicitely put here as a way
    # to always check that specific case.
    # values = [0, 7, 0, 6, 5, 1, 6, 7, 0, 0]
    values = [4, 4, 1, 5, 2, 0, 1, 4, 3, 9]
    blocksizes = small_blocksizes
    chunkshape = 1
    buffersize = 1
    ss = blocksizes[2]
    nrows = ss
    reopen = 0
    nrep = ss
    il = 0
    sl = ss


class SV12bTestCase(SelectValuesTestCase):
    # This checks a special case that failed. It was discovered in a
    # random test above (SV10b). It is explicitely put here as a way
    # to always check that specific case.
    # values = [0, 7, 0, 6, 5, 1, 6, 7, 0, 0]
    values = [4, 4, 1, 5, 2, 0, 1, 4, 3, 9]
    blocksizes = tb.idxutils.calc_chunksize(minRowIndex, memlevel=1)
    chunkshape = 2
    buffersize = 2
    ss = blocksizes[2]
    nrows = ss
    reopen = 1
    nrep = ss
    il = 0
    sl = ss


class SV13aTestCase(SelectValuesTestCase):
    values = [0, 7, 0, 6, 5, 1, 6, 7, 0, 0]
    blocksizes = small_blocksizes
    chunkshape = 3
    buffersize = 5
    ss = blocksizes[2]
    nrows = ss
    reopen = 0
    nrep = ss
    il = 0
    sl = ss


class SV13bTestCase(SelectValuesTestCase):
    values = [0, 7, 0, 6, 5, 1, 6, 7, 0, 0]
    blocksizes = tb.idxutils.calc_chunksize(minRowIndex, memlevel=1)
    chunkshape = 5
    buffersize = 10
    ss = blocksizes[2]
    nrows = ss
    reopen = 1
    nrep = ss
    il = 0
    sl = ss


class SV14aTestCase(SelectValuesTestCase):
    values = [1, 7, 6, 7, 0, 7, 4, 4, 9, 5]
    blocksizes = small_blocksizes
    chunkshape = 2
    buffersize = 5
    ss = blocksizes[2]
    nrows = ss
    reopen = 0
    cs = blocksizes[3]
    nrep = cs
    il = -5
    sl = 500


class SV14bTestCase(SelectValuesTestCase):
    values = [1, 7, 6, 7, 0, 7, 4, 4, 9, 5]
    blocksizes = tb.idxutils.calc_chunksize(minRowIndex, memlevel=1)
    chunkshape = 9
    buffersize = 10
    ss = blocksizes[2]
    nrows = ss
    reopen = 1
    nrep = 9
    il = 0
    cs = blocksizes[3]
    sl = ss - cs + 1


class SV15aTestCase(SelectValuesTestCase):
    # Test that checks for case where there are not valid values in
    # the indexed part, but they exist in the non-indexed region.
    # At least, test01b takes account of that
    random = 1
    # Both values of seed below triggers a fail in indexing code
    # seed = 1885
    seed = 183
    blocksizes = small_blocksizes
    ss = blocksizes[2]
    nrows = ss * 5 + 1
    reopen = 0
    cs = blocksizes[3]
    nrep = cs - 1
    il = -10
    sl = nrows


class SV15bTestCase(SelectValuesTestCase):
    # Test that checks for case where there are not valid values in
    # the indexed part, but they exist in the non-indexed region.
    # At least, test01b takes account of that
    random = 1
    # Both values of seed below triggers a fail in indexing code
    seed = 1885
    # seed = 183
    blocksizes = tb.idxutils.calc_chunksize(minRowIndex, memlevel=1)
    ss = blocksizes[2]
    nrows = ss * 5 + 1
    reopen = 1
    cs = blocksizes[3]
    nrep = cs - 1
    il = -10
    sl = nrows


class LastRowReuseBuffers(common.PyTablesTestCase):
    # Test that checks for possible reuse of buffers coming
    # from last row in the sorted part of indexes
    nelem = 1221
    np.random.seed(1)
    random.seed(1)

    class Record(tb.IsDescription):
        id1 = tb.Int16Col()

    def setUp(self):
        super().setUp()
        self.h5fname = tempfile.mktemp(".h5")
        self.h5file = None

    def tearDown(self):
        if self.h5file is not None:
            self.h5file.close()
        if Path(self.h5fname).is_file():
            Path(self.h5fname).unlink()
        super().tearDown()

    def test00_lrucache(self):
        self.h5file = tb.open_file(self.h5fname, "w", node_cache_slots=64)
        ta = self.h5file.create_table(
            "/", "table", self.Record, filters=tb.Filters(1)
        )
        id1 = np.random.randint(0, 2**15, self.nelem)
        ta.append([id1])

        ta.cols.id1.create_index()

        for i in range(self.nelem):
            nrow = random.randrange(self.nelem)
            value = id1[nrow]
            idx = ta.get_where_list("id1 == %s" % value)
            self.assertGreater(len(idx), 0, f"idx--> {idx} {i} {nrow} {value}")
            self.assertTrue(
                nrow in idx, f"nrow not found: {idx} != {nrow}, {value}"
            )

    def test01_nocache(self):
        self.h5file = tb.open_file(self.h5fname, "w", node_cache_slots=0)
        ta = self.h5file.create_table(
            "/", "table", self.Record, filters=tb.Filters(1)
        )
        id1 = np.random.randint(0, 2**15, self.nelem)
        ta.append([id1])

        ta.cols.id1.create_index()

        for i in range(self.nelem):
            nrow = random.randrange(self.nelem)
            value = id1[nrow]
            idx = ta.get_where_list("id1 == %s" % value)
            self.assertGreater(len(idx), 0, f"idx--> {idx} {i} {nrow} {value}")
            self.assertTrue(
                nrow in idx, f"nrow not found: {idx} != {nrow}, {value}"
            )

    def test02_dictcache(self):
        self.h5file = tb.open_file(self.h5fname, "w", node_cache_slots=-64)
        ta = self.h5file.create_table(
            "/", "table", self.Record, filters=tb.Filters(1)
        )
        id1 = np.random.randint(0, 2**15, self.nelem)
        ta.append([id1])

        ta.cols.id1.create_index()

        for i in range(self.nelem):
            nrow = random.randrange(self.nelem)
            value = id1[nrow]
            idx = ta.get_where_list("id1 == %s" % value)
            self.assertGreater(len(idx), 0, f"idx--> {idx} {i} {nrow} {value}")
            self.assertTrue(
                nrow in idx, f"nrow not found: {idx} != {nrow}, {value}"
            )


normal_tests = (
    "SV1aTestCase",
    "SV2aTestCase",
    "SV3aTestCase",
)

heavy_tests = (
    # The next are too hard to be in the 'normal' suite
    "SV1bTestCase",
    "SV2bTestCase",
    "SV3bTestCase",
    "SV4aTestCase",
    "SV5aTestCase",
    "SV6aTestCase",
    "SV7aTestCase",
    "SV8aTestCase",
    "SV9aTestCase",
    "SV10aTestCase",
    "SV11aTestCase",
    "SV12aTestCase",
    "SV13aTestCase",
    "SV14aTestCase",
    "SV15aTestCase",
    # This are properly heavy
    "SV4bTestCase",
    "SV5bTestCase",
    "SV6bTestCase",
    "SV7bTestCase",
    "SV8bTestCase",
    "SV9bTestCase",
    "SV10bTestCase",
    "SV11bTestCase",
    "SV12bTestCase",
    "SV13bTestCase",
    "SV14bTestCase",
    "SV15bTestCase",
)


# Base classes for the different type indexes.
class UltraLightITableMixin:
    kind = "ultralight"


class LightITableMixin:
    kind = "light"


class MediumITableMixin:
    kind = "medium"


class FullITableMixin:
    kind = "full"


# Parameters for indexed queries.
ckinds = ["UltraLight", "Light", "Medium", "Full"]
testlevels = ["Normal", "Heavy"]

# Indexed queries: ``[ULMF]I[NH]SVXYTestCase``, where:
#
# 1. U is for 'UltraLight', L for 'Light', M for 'Medium', F for 'Full' indexes
# 2. N is for 'Normal', H for 'Heavy' tests


def iclassdata():
    for ckind in ckinds:
        for ctest in normal_tests + heavy_tests:
            classname = f"{ckind[0]}I{testlevels[common.heavy][0]}{ctest}"
            # Uncomment the next one and comment the past one if one
            # don't want to include the methods (testing purposes only)
            # cbasenames = ( '%sITableMixin' % ckind, "object")
            cbasenames = ("%sITableMixin" % ckind, ctest)
            classdict = dict(heavy=bool(ctest in heavy_tests))
            yield (classname, cbasenames, classdict)


# Create test classes.
for cname, cbasenames, cdict in iclassdata():
    cbases = tuple(eval(cbase) for cbase in cbasenames)
    class_ = type(cname, cbases, cdict)
    exec("%s = class_" % cname)


# Test case for issue #319
class BuffersizeMultipleChunksize(
    common.TempFileMixin, common.PyTablesTestCase
):
    open_mode = "w"

    def test01(self):
        np.random.seed(2)
        n = 700_000
        cs = 50_000
        nchunks = n // cs

        arr = np.zeros(
            (n,), dtype=[("index", "i8"), ("o", "i8"), ("value", "f8")]
        )
        arr["index"] = np.arange(n)
        arr["o"] = np.random.randint(-20_000, -15_000, size=n)
        arr["value"] = np.random.randn(n)

        node = self.h5file.create_group("/", "foo")
        table = self.h5file.create_table(
            node,
            "table",
            dict(
                index=tb.Int64Col(),
                o=tb.Int64Col(),
                value=tb.FloatCol(shape=()),
            ),
            expectedrows=10_000_000,
        )

        table.append(arr)

        self._reopen("a")

        v1 = np.unique(arr["o"])[0]
        v2 = np.unique(arr["o"])[1]
        res = np.array([v1, v2])
        selector = f"((o == {v1}) | (o == {v2}))"
        if common.verbose:
            print("selecting values: %s" % selector)

        table = self.h5file.root.foo.table

        result = np.unique(table.read_where(selector)["o"])
        np.testing.assert_almost_equal(result, res)
        if common.verbose:
            print("select entire table:")
            print(f"result: {result}\texpected: {res}")

        if common.verbose:
            print("index the column o")
        table.cols.o.create_index()  # this was triggering the issue

        if common.verbose:
            print("select via chunks")
        for i in range(nchunks):
            result = table.read_where(
                selector, start=i * cs, stop=(i + 1) * cs
            )
            result = np.unique(result["o"])
            np.testing.assert_almost_equal(np.unique(result), res)
            if common.verbose:
                print(f"result: {result}\texpected: {res}")


# Test case for issue #441
class SideEffectNumPyQuicksort(common.PyTablesTestCase):

    def test01(self):
        bug_file = common.test_filename("bug-idx.h5")
        tmp_file = tempfile.mktemp(".h5")
        tb.copy_file(bug_file, tmp_file)
        h5 = tb.open_file(tmp_file, "a")
        o = h5.root.table
        vals = o.cols.path[:]
        npvals = set(np.where(vals == 6)[0])

        # Setting the chunkshape is critical for reproducing the bug
        t = o.copy(newname="table2", chunkshape=2730)
        t.cols.path.create_index()
        indexed = {r.nrow for r in t.where("path == 6")}

        if common.verbose:
            diffs = sorted(npvals - indexed)
            print("ndiff:", len(diffs), diffs)
        self.assertEqual(len(npvals), len(indexed))

        h5.close()
        if Path(tmp_file).is_file():
            Path(tmp_file).unlink()


# -----------------------------


def suite():
    theSuite = common.unittest.TestSuite()

    niter = 1

    for n in range(niter):
        for cdata in iclassdata():
            class_ = eval(cdata[0])
            if not class_.heavy:
                suite_ = common.make_suite(class_)
                theSuite.addTest(suite_)
            elif common.heavy:
                suite_ = common.make_suite(class_)
                theSuite.addTest(suite_)
        theSuite.addTest(common.make_suite(LastRowReuseBuffers))
        theSuite.addTest(common.make_suite(BuffersizeMultipleChunksize))
        theSuite.addTest(common.make_suite(SideEffectNumPyQuicksort))
    return theSuite


if __name__ == "__main__":
    import sys

    common.parse_argv(sys.argv)
    common.print_versions()
    common.unittest.main(defaultTest="suite")
