import copy
import tempfile
from pathlib import Path

import numpy as np

import tables as tb
from tables.tests import common

# Sensible parameters for indexing with small blocksizes
minRowIndex = 10
small_blocksizes = (96, 24, 6, 3)


class TDescr(tb.IsDescription):
    var1 = tb.StringCol(itemsize=4, dflt=b"", pos=1)
    var2 = tb.BoolCol(dflt=0, pos=2)
    var3 = tb.IntCol(dflt=0, pos=3)
    var4 = tb.FloatCol(dflt=0, pos=4)


class BasicTestCase(common.TempFileMixin, common.PyTablesTestCase):
    compress = 0
    complib = "zlib"
    shuffle = 0
    fletcher32 = 0
    nrows = minRowIndex
    ss = small_blocksizes[2]

    def setUp(self):
        super().setUp()

        self.rootgroup = self.h5file.root
        self.populateFile()
        # Close the file
        self.h5file.close()

    def populateFile(self):
        group = self.rootgroup
        # Create a table
        title = "This is the IndexArray title"
        self.filters = tb.Filters(
            complevel=self.compress,
            complib=self.complib,
            shuffle=self.shuffle,
            fletcher32=self.fletcher32,
        )
        table = self.h5file.create_table(
            group, "table", TDescr, title, self.filters, self.nrows
        )
        for i in range(self.nrows):
            table.row["var1"] = str(i).encode("ascii")
            # table.row['var2'] = i > 2
            table.row["var2"] = i % 2
            table.row["var3"] = i
            table.row["var4"] = float(self.nrows - i - 1)
            table.row.append()
        table.flush()
        # Index all entries:
        for col in table.colinstances.values():
            indexrows = col.create_index(_blocksizes=small_blocksizes)
        if common.verbose:
            print("Number of written rows:", self.nrows)
            print("Number of indexed rows:", indexrows)

        return

    def test00_flushLastRow(self):
        """Checking flushing an Index incrementing only the last row."""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test00_flushLastRow..." % self.__class__.__name__
            )

        # Open the HDF5 file in append mode
        self.h5file = tb.open_file(self.h5fname, mode="a")
        table = self.h5file.root.table
        # Add just 3 rows more
        for i in range(3):
            table.row["var1"] = str(i).encode("ascii")
            table.row.append()
        table.flush()  # redo the indexes
        idxcol = table.cols.var1.index
        if common.verbose:
            print("Max rows in buf:", table.nrowsinbuf)
            print("Number of elements per slice:", idxcol.slicesize)
            print("Chunk size:", idxcol.sorted.chunksize)
            print("Elements in last row:", idxcol.indicesLR[-1])

        # Do a selection
        results = [p["var1"] for p in table.where('var1 == b"1"')]
        self.assertEqual(len(results), 2)
        self.assertEqual(results, [b"1"] * 2)

    def test00_update(self):
        """Checking automatic re-indexing after an update operation."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test00_update..." % self.__class__.__name__)

        # Open the HDF5 file in append mode
        self.h5file = tb.open_file(self.h5fname, mode="a")
        table = self.h5file.root.table
        # Modify a couple of columns
        for i, row in enumerate(table.where("(var3>1) & (var3<5)")):
            row["var1"] = str(i)
            row["var3"] = i
            row.update()
        table.flush()  # redo the indexes
        idxcol1 = table.cols.var1.index
        idxcol3 = table.cols.var3.index
        if common.verbose:
            print("Dirtyness of var1 col:", idxcol1.dirty)
            print("Dirtyness of var3 col:", idxcol3.dirty)
        self.assertEqual(idxcol1.dirty, False)
        self.assertEqual(idxcol3.dirty, False)

        # Do a couple of selections
        results = [p["var1"] for p in table.where('var1 == b"1"')]
        self.assertEqual(len(results), 2)
        self.assertEqual(results, [b"1"] * 2)
        results = [p["var3"] for p in table.where("var3 == 0")]
        self.assertEqual(len(results), 2)
        self.assertEqual(results, [0] * 2)

    def test01_readIndex(self):
        """Checking reading an Index (string flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01_readIndex..." % self.__class__.__name__)

        # Open the HDF5 file in read-only mode
        self.h5file = tb.open_file(self.h5fname, mode="r")
        table = self.h5file.root.table
        idxcol = table.cols.var1.index
        if common.verbose:
            print("Max rows in buf:", table.nrowsinbuf)
            print("Number of elements per slice:", idxcol.slicesize)
            print("Chunk size:", idxcol.sorted.chunksize)

        # Do a selection
        results = [p["var1"] for p in table.where('var1 == b"1"')]
        self.assertEqual(len(results), 1)
        self.assertEqual(results, [b"1"])

    def test02_readIndex(self):
        """Checking reading an Index (bool flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02_readIndex..." % self.__class__.__name__)

        # Open the HDF5 file in read-only mode
        self.h5file = tb.open_file(self.h5fname, mode="r")
        table = self.h5file.root.table
        idxcol = table.cols.var2.index
        if common.verbose:
            print("Rows in table:", table.nrows)
            print("Max rows in buf:", table.nrowsinbuf)
            print("Number of elements per slice:", idxcol.slicesize)
            print("Chunk size:", idxcol.sorted.chunksize)

        # Do a selection
        results = [p["var2"] for p in table.where("var2 == True")]
        if common.verbose:
            print("Selected values:", results)
        self.assertEqual(len(results), self.nrows // 2)
        self.assertEqual(results, [True] * (self.nrows // 2))

    def test03_readIndex(self):
        """Checking reading an Index (int flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03_readIndex..." % self.__class__.__name__)

        # Open the HDF5 file in read-only mode
        self.h5file = tb.open_file(self.h5fname, mode="r")
        table = self.h5file.root.table
        idxcol = table.cols.var3.index
        if common.verbose:
            print("Max rows in buf:", table.nrowsinbuf)
            print("Number of elements per slice:", idxcol.slicesize)
            print("Chunk size:", idxcol.sorted.chunksize)

        # Do a selection
        results = [p["var3"] for p in table.where("(1<var3)&(var3<10)")]
        if common.verbose:
            print("Selected values:", results)
        self.assertEqual(len(results), min(10, table.nrows) - 2)
        self.assertEqual(results, list(range(2, min(10, table.nrows))))

    def test04_readIndex(self):
        """Checking reading an Index (float flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test04_readIndex..." % self.__class__.__name__)

        # Open the HDF5 file in read-only mode
        self.h5file = tb.open_file(self.h5fname, mode="r")
        table = self.h5file.root.table
        idxcol = table.cols.var4.index
        if common.verbose:
            print("Max rows in buf:", table.nrowsinbuf)
            print("Number of rows in table:", table.nrows)
            print("Number of elements per slice:", idxcol.slicesize)
            print("Chunk size:", idxcol.sorted.chunksize)

        # Do a selection
        results = [p["var4"] for p in table.where("var4 < 10")]
        # results = [p["var4"] for p in table.where('(1<var4)&(var4<10)')]
        if common.verbose:
            print("Selected values:", results)
        self.assertEqual(len(results), min(10, table.nrows))
        self.assertEqual(
            results,
            [float(i) for i in reversed(list(range(min(10, table.nrows))))],
        )

    def test05_getWhereList(self):
        """Checking reading an Index with get_where_list (string flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test05_getWhereList..." % self.__class__.__name__
            )

        # Open the HDF5 file in read-write mode
        self.h5file = tb.open_file(self.h5fname, mode="a")
        table = self.h5file.root.table
        idxcol = table.cols.var4.index
        if common.verbose:
            print("Max rows in buf:", table.nrowsinbuf)
            print("Number of elements per slice:", idxcol.slicesize)
            print("Chunk size:", idxcol.sorted.chunksize)

        # Do a selection
        table.flavor = "python"
        rowList1 = table.get_where_list('var1 < b"10"')
        rowList2 = [p.nrow for p in table if p["var1"] < b"10"]
        if common.verbose:
            print("Selected values:", rowList1)
            print("Should look like:", rowList2)
        self.assertEqual(len(rowList1), len(rowList2))
        self.assertEqual(rowList1, rowList2)

    def test06_getWhereList(self):
        """Checking reading an Index with get_where_list (bool flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test06_getWhereList..." % self.__class__.__name__
            )

        # Open the HDF5 file in read-write mode
        self.h5file = tb.open_file(self.h5fname, mode="a")
        table = self.h5file.root.table
        idxcol = table.cols.var2.index
        if common.verbose:
            print("Max rows in buf:", table.nrowsinbuf)
            print("Rows in tables:", table.nrows)
            print("Number of elements per slice:", idxcol.slicesize)
            print("Chunk size:", idxcol.sorted.chunksize)

        # Do a selection
        table.flavor = "numpy"
        rowList1 = table.get_where_list("var2 == False", sort=True)
        rowList2 = [p.nrow for p in table if p["var2"] is False]
        # Convert to a NumPy object
        rowList2 = np.array(rowList2, np.int64)
        if common.verbose:
            print("Selected values:", rowList1)
            print("Should look like:", rowList2)
        self.assertEqual(len(rowList1), len(rowList2))
        self.assertTrue(common.allequal(rowList1, rowList2))

    def test07_getWhereList(self):
        """Checking reading an Index with get_where_list (int flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test07_getWhereList..." % self.__class__.__name__
            )

        # Open the HDF5 file in read-write mode
        self.h5file = tb.open_file(self.h5fname, mode="a")
        table = self.h5file.root.table
        idxcol = table.cols.var4.index
        if common.verbose:
            print("Max rows in buf:", table.nrowsinbuf)
            print("Number of elements per slice:", idxcol.slicesize)
            print("Chunk size:", idxcol.sorted.chunksize)

        # Do a selection
        table.flavor = "python"
        rowList1 = table.get_where_list("var3 < 15", sort=True)
        rowList2 = [p.nrow for p in table if p["var3"] < 15]
        if common.verbose:
            print("Selected values:", rowList1)
            print("Should look like:", rowList2)
        self.assertEqual(len(rowList1), len(rowList2))
        self.assertEqual(rowList1, rowList2)

    def test08_getWhereList(self):
        """Checking reading an Index with get_where_list (float flavor)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test08_getWhereList..." % self.__class__.__name__
            )

        # Open the HDF5 file in read-write mode
        self.h5file = tb.open_file(self.h5fname, mode="a")
        table = self.h5file.root.table
        idxcol = table.cols.var4.index
        if common.verbose:
            print("Max rows in buf:", table.nrowsinbuf)
            print("Number of elements per slice:", idxcol.slicesize)
            print("Chunk size:", idxcol.sorted.chunksize)

        # Do a selection
        table.flavor = "python"
        rowList1 = table.get_where_list("var4 < 10", sort=True)
        rowList2 = [p.nrow for p in table if p["var4"] < 10]
        if common.verbose:
            print("Selected values:", rowList1)
            print("Should look like:", rowList2)
        self.assertEqual(len(rowList1), len(rowList2))
        self.assertEqual(rowList1, rowList2)

    def test09a_removeIndex(self):
        """Checking removing an index."""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test09a_removeIndex..." % self.__class__.__name__
            )

        # Open the HDF5 file in read-write mode
        self.h5file = tb.open_file(self.h5fname, mode="a")
        table = self.h5file.root.table
        idxcol = table.cols.var1.index
        if common.verbose:
            print("Before deletion")
            print("var1 column:", table.cols.var1)
        self.assertEqual(table.colindexed["var1"], 1)
        self.assertIsNotNone(idxcol)

        # delete the index
        table.cols.var1.remove_index()
        if common.verbose:
            print("After deletion")
            print("var1 column:", table.cols.var1)
        self.assertIsNone(table.cols.var1.index)
        self.assertEqual(table.colindexed["var1"], 0)

        # re-create the index again
        indexrows = table.cols.var1.create_index(_blocksizes=small_blocksizes)
        self.assertIsNotNone(indexrows)
        idxcol = table.cols.var1.index
        if common.verbose:
            print("After re-creation")
            print("var1 column:", table.cols.var1)
        self.assertIsNotNone(idxcol)
        self.assertEqual(table.colindexed["var1"], 1)

    def test09b_removeIndex(self):
        """Checking removing an index (persistent version)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test09b_removeIndex..." % self.__class__.__name__
            )

        # Open the HDF5 file in read-write mode
        self.h5file = tb.open_file(self.h5fname, mode="a")
        table = self.h5file.root.table
        idxcol = table.cols.var1.index
        if common.verbose:
            print("Before deletion")
            print("var1 index column:", table.cols.var1)
        self.assertIsNotNone(idxcol)
        self.assertEqual(table.colindexed["var1"], 1)
        # delete the index
        table.cols.var1.remove_index()

        # close and reopen the file
        self._reopen(mode="a")
        table = self.h5file.root.table
        idxcol = table.cols.var1.index

        if common.verbose:
            print("After deletion")
            print("var1 column:", table.cols.var1)
        self.assertIsNone(table.cols.var1.index)
        self.assertEqual(table.colindexed["var1"], 0)

        # re-create the index again
        indexrows = table.cols.var1.create_index(_blocksizes=small_blocksizes)
        self.assertIsNotNone(indexrows)
        idxcol = table.cols.var1.index
        if common.verbose:
            print("After re-creation")
            print("var1 column:", table.cols.var1)
        self.assertIsNotNone(idxcol)
        self.assertEqual(table.colindexed["var1"], 1)

    def test10a_moveIndex(self):
        """Checking moving a table with an index."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test10a_moveIndex..." % self.__class__.__name__)

        # Open the HDF5 file in read-write mode
        self.h5file = tb.open_file(self.h5fname, mode="a")
        table = self.h5file.root.table
        idxcol = table.cols.var1.index
        if common.verbose:
            print("Before move")
            print("var1 column:", idxcol)
        self.assertEqual(table.colindexed["var1"], 1)
        self.assertIsNotNone(idxcol)

        # Create a new group called "agroup"
        agroup = self.h5file.create_group("/", "agroup")

        # move the table to "agroup"
        table.move(agroup, "table2")
        if common.verbose:
            print("After move")
            print("var1 column:", idxcol)
        self.assertIsNotNone(table.cols.var1.index)
        self.assertEqual(table.colindexed["var1"], 1)

        # Some sanity checks
        table.flavor = "python"
        rowList1 = table.get_where_list('var1 < b"10"')
        rowList2 = [p.nrow for p in table if p["var1"] < b"10"]
        if common.verbose:
            print("Selected values:", rowList1)
            print("Should look like:", rowList2)
        self.assertEqual(len(rowList1), len(rowList2))
        self.assertEqual(rowList1, rowList2)

    def test10b_moveIndex(self):
        """Checking moving a table with an index (persistent version)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test10b_moveIndex..." % self.__class__.__name__)

        # Open the HDF5 file in read-write mode
        self.h5file = tb.open_file(self.h5fname, mode="a")
        table = self.h5file.root.table
        idxcol = table.cols.var1.index
        if common.verbose:
            print("Before move")
            print("var1 index column:", idxcol)
        self.assertIsNotNone(idxcol)
        self.assertEqual(table.colindexed["var1"], 1)
        # Create a new group called "agroup"
        agroup = self.h5file.create_group("/", "agroup")

        # move the table to "agroup"
        table.move(agroup, "table2")

        # close and reopen the file
        self._reopen(mode="a")
        table = self.h5file.root.agroup.table2
        idxcol = table.cols.var1.index

        if common.verbose:
            print("After move")
            print("var1 column:", idxcol)
        self.assertIsNotNone(table.cols.var1.index)
        self.assertEqual(table.colindexed["var1"], 1)

        # Some sanity checks
        table.flavor = "python"
        rowList1 = table.get_where_list('var1 < b"10"')
        rowList2 = [p.nrow for p in table if p["var1"] < b"10"]
        if common.verbose:
            print("Selected values:", rowList1, type(rowList1))
            print("Should look like:", rowList2, type(rowList2))
        self.assertEqual(len(rowList1), len(rowList2))
        self.assertEqual(rowList1, rowList2)

    def test10c_moveIndex(self):
        """Checking moving a table with an index (small node cache)."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test10c_moveIndex..." % self.__class__.__name__)

        # Open the HDF5 file in read-write mode
        self.h5file = tb.open_file(self.h5fname, mode="a", node_cache_slots=10)
        table = self.h5file.root.table
        idxcol = table.cols.var1.index
        if common.verbose:
            print("Before move")
            print("var1 column:", idxcol)
        self.assertEqual(table.colindexed["var1"], 1)
        self.assertIsNotNone(idxcol)

        # Create a new group called "agroup"
        agroup = self.h5file.create_group("/", "agroup")

        # move the table to "agroup"
        table.move(agroup, "table2")
        if common.verbose:
            print("After move")
            print("var1 column:", idxcol)
        self.assertIsNotNone(table.cols.var1.index)
        self.assertEqual(table.colindexed["var1"], 1)

        # Some sanity checks
        table.flavor = "python"
        rowList1 = table.get_where_list('var1 < b"10"')
        rowList2 = [p.nrow for p in table if p["var1"] < b"10"]
        if common.verbose:
            print("Selected values:", rowList1)
            print("Should look like:", rowList2)
        self.assertEqual(len(rowList1), len(rowList2))
        self.assertEqual(rowList1, rowList2)

    def test10d_moveIndex(self):
        """Checking moving a table with an index (no node cache)."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test10d_moveIndex..." % self.__class__.__name__)

        # Open the HDF5 file in read-write mode
        self.h5file = tb.open_file(self.h5fname, mode="a", node_cache_slots=0)
        table = self.h5file.root.table
        idxcol = table.cols.var1.index
        if common.verbose:
            print("Before move")
            print("var1 column:", idxcol)
        self.assertEqual(table.colindexed["var1"], 1)
        self.assertIsNotNone(idxcol)

        # Create a new group called "agroup"
        agroup = self.h5file.create_group("/", "agroup")

        # move the table to "agroup"
        table.move(agroup, "table2")
        if common.verbose:
            print("After move")
            print("var1 column:", idxcol)
        self.assertIsNotNone(table.cols.var1.index)
        self.assertEqual(table.colindexed["var1"], 1)

        # Some sanity checks
        table.flavor = "python"
        rowList1 = table.get_where_list('var1 < b"10"')
        rowList2 = [p.nrow for p in table if p["var1"] < b"10"]
        if common.verbose:
            print("Selected values:", rowList1)
            print("Should look like:", rowList2)
        self.assertEqual(len(rowList1), len(rowList2))
        self.assertEqual(rowList1, rowList2)

    def test11a_removeTableWithIndex(self):
        """Checking removing a table with indexes."""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test11a_removeTableWithIndex..."
                % self.__class__.__name__
            )

        # Open the HDF5 file in read-write mode
        self.h5file = tb.open_file(self.h5fname, mode="a")
        table = self.h5file.root.table
        idxcol = table.cols.var1.index
        if common.verbose:
            print("Before deletion")
            print("var1 column:", table.cols.var1)
        self.assertEqual(table.colindexed["var1"], 1)
        self.assertIsNotNone(idxcol)

        # delete the table
        self.h5file.remove_node("/table")
        if common.verbose:
            print("After deletion")
        self.assertNotIn("table", self.h5file.root)

        # re-create the table and the index again
        table = self.h5file.create_table(
            "/", "table", TDescr, "New table", self.filters, self.nrows
        )
        for i in range(self.nrows):
            table.row["var1"] = str(i)
            table.row["var2"] = i % 2
            table.row["var3"] = i
            table.row["var4"] = float(self.nrows - i - 1)
            table.row.append()
        table.flush()
        # Index all entries:
        for col in table.colinstances.values():
            indexrows = col.create_index(_blocksizes=small_blocksizes)
            self.assertIsNotNone(indexrows)
        idxcol = table.cols.var1.index
        if common.verbose:
            print("After re-creation")
            print("var1 column:", table.cols.var1)
        self.assertIsNotNone(idxcol)
        self.assertEqual(table.colindexed["var1"], 1)

    def test11b_removeTableWithIndex(self):
        """Checking removing a table with indexes (persistent version 2)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test11b_removeTableWithIndex..."
                % self.__class__.__name__
            )

        self.h5file = tb.open_file(self.h5fname, mode="a")
        table = self.h5file.root.table
        idxcol = table.cols.var1.index
        if common.verbose:
            print("Before deletion")
            print("var1 column:", table.cols.var1)
        self.assertEqual(table.colindexed["var1"], 1)
        self.assertIsNotNone(idxcol)

        # delete the table
        self.h5file.remove_node("/table")
        if common.verbose:
            print("After deletion")
        self.assertNotIn("table", self.h5file.root)

        # close and reopen the file
        self._reopen(mode="r+")

        # re-create the table and the index again
        table = self.h5file.create_table(
            "/", "table", TDescr, "New table", self.filters, self.nrows
        )
        for i in range(self.nrows):
            table.row["var1"] = str(i)
            table.row["var2"] = i % 2
            table.row["var3"] = i
            table.row["var4"] = float(self.nrows - i - 1)
            table.row.append()
        table.flush()
        # Index all entries:
        for col in table.colinstances.values():
            indexrows = col.create_index(_blocksizes=small_blocksizes)
            self.assertIsNotNone(indexrows)
        idxcol = table.cols.var1.index
        if common.verbose:
            print("After re-creation")
            print("var1 column:", table.cols.var1)
        self.assertIsNotNone(idxcol)
        self.assertEqual(table.colindexed["var1"], 1)

    # Test provided by Andrew Straw
    def test11c_removeTableWithIndex(self):
        """Checking removing a table with indexes (persistent version 3)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test11c_removeTableWithIndex..."
                % self.__class__.__name__
            )

        class Distance(tb.IsDescription):
            frame = tb.Int32Col(pos=0)
            distance = tb.FloatCol(pos=1)

        # Delete the old temporal file
        Path(self.h5fname).unlink()

        self.h5fname = tempfile.mktemp(".h5")
        self.h5file = tb.open_file(self.h5fname, mode="w")
        table = self.h5file.create_table(
            self.h5file.root, "distance_table", Distance
        )
        table.cols.frame.create_index(_blocksizes=small_blocksizes)
        r = table.row
        for i in range(10):
            r["frame"] = i
            r["distance"] = float(i**2)
            r.append()
        table.flush()

        self._reopen(mode="r+")

        self.h5file.remove_node(self.h5file.root.distance_table)

    def test12_doubleIterate(self):
        self.h5file = tb.open_file(self.h5fname, mode="r")
        table = self.h5file.root.table
        tests = [1, 4, self.nrows]
        if self.nrows > 500:
            tests.append(self.nrows - 500)
        for limit in tests:
            handle_a = [0, table.where("(var3 < e)", dict(e=limit))]
            handle_b = [0, table.where("(var3 < e)", dict(e=limit))]

            try:
                while True:
                    next(handle_b[1])
                    handle_b[0] += 1
            except StopIteration:
                for _ in handle_a[1]:
                    handle_a[0] += 1
                for _ in handle_b[1]:
                    handle_b[0] += 1

            self.assertEqual(handle_a[0], limit)
            self.assertEqual(handle_b[0], limit)
            self.assertEqual(
                len(list(table.where("(var3 < e)", dict(e=limit)))), limit
            )


small_ss = small_blocksizes[2]


class BasicReadTestCase(BasicTestCase):
    compress = 0
    complib = "zlib"
    shuffle = 0
    fletcher32 = 0
    nrows = small_ss


class ZlibReadTestCase(BasicTestCase):
    compress = 1
    complib = "zlib"
    shuffle = 0
    fletcher32 = 0
    nrows = small_ss


@common.unittest.skipIf(
    not common.blosc_avail, "BLOSC compression library not available"
)
class BloscReadTestCase(BasicTestCase):
    compress = 1
    complib = "blosc"
    shuffle = 0
    fletcher32 = 0
    nrows = small_ss


@common.unittest.skipIf(
    not common.lzo_avail, "LZO compression library not available"
)
class LZOReadTestCase(BasicTestCase):
    compress = 1
    complib = "lzo"
    shuffle = 0
    fletcher32 = 0
    nrows = small_ss


@common.unittest.skipIf(
    not common.bzip2_avail, "BZIP2 compression library not available"
)
class Bzip2ReadTestCase(BasicTestCase):
    compress = 1
    complib = "bzip2"
    shuffle = 0
    fletcher32 = 0
    nrows = small_ss


class ShuffleReadTestCase(BasicTestCase):
    compress = 1
    complib = "zlib"
    shuffle = 1
    fletcher32 = 0
    nrows = small_ss


class Fletcher32ReadTestCase(BasicTestCase):
    compress = 1
    complib = "zlib"
    shuffle = 0
    fletcher32 = 1
    nrows = small_ss


class ShuffleFletcher32ReadTestCase(BasicTestCase):
    compress = 1
    complib = "zlib"
    shuffle = 1
    fletcher32 = 1
    nrows = small_ss


class OneHalfTestCase(BasicTestCase):
    nrows = small_ss + small_ss // 2


class UpperBoundTestCase(BasicTestCase):
    nrows = small_ss + 1


class LowerBoundTestCase(BasicTestCase):
    nrows = small_ss * 2 - 1


class DeepTableIndexTestCase(common.TempFileMixin, common.PyTablesTestCase):
    nrows = minRowIndex

    def test01(self):
        """Checking the indexing of a table in a 2nd level hierarchy"""

        # Create an instance of an HDF5 Table
        group = self.h5file.create_group(self.h5file.root, "agroup")
        # Create a table
        title = "This is the IndexArray title"
        table = self.h5file.create_table(
            group, "table", TDescr, title, None, self.nrows
        )
        for i in range(self.nrows):
            # Fill rows with defaults
            table.row.append()
        table.flush()
        # Index some column
        indexrows = table.cols.var1.create_index()
        self.assertIsNotNone(indexrows)
        idxcol = table.cols.var1.index
        # Some sanity checks
        self.assertEqual(table.colindexed["var1"], 1)
        self.assertIsNotNone(idxcol)
        self.assertEqual(idxcol.nelements, self.nrows)

    def test01b(self):
        """Checking the indexing of a table in 2nd level
        (persistent version)"""

        # Create an instance of an HDF5 Table
        group = self.h5file.create_group(self.h5file.root, "agroup")

        # Create a table
        title = "This is the IndexArray title"
        table = self.h5file.create_table(
            group, "table", TDescr, title, None, self.nrows
        )
        for i in range(self.nrows):
            # Fill rows with defaults
            table.row.append()
        table.flush()

        # Index some column
        indexrows = table.cols.var1.create_index()
        self.assertIsNotNone(indexrows)
        idxcol = table.cols.var1.index

        # Close and re-open this file
        self._reopen(mode="a")

        table = self.h5file.root.agroup.table
        idxcol = table.cols.var1.index
        # Some sanity checks
        self.assertEqual(table.colindexed["var1"], 1)
        self.assertIsNotNone(idxcol)
        self.assertEqual(idxcol.nelements, self.nrows)

    def test02(self):
        """Checking the indexing of a table in a 4th level hierarchy"""

        # Create an instance of an HDF5 Table
        group = self.h5file.create_group(self.h5file.root, "agroup")
        group = self.h5file.create_group(group, "agroup")
        group = self.h5file.create_group(group, "agroup")

        # Create a table
        title = "This is the IndexArray title"
        table = self.h5file.create_table(
            group, "table", TDescr, title, None, self.nrows
        )
        for i in range(self.nrows):
            # Fill rows with defaults
            table.row.append()
        table.flush()

        # Index some column
        indexrows = table.cols.var1.create_index()
        self.assertIsNotNone(indexrows)
        idxcol = table.cols.var1.index

        # Some sanity checks
        self.assertEqual(table.colindexed["var1"], 1)
        self.assertIsNotNone(idxcol)
        self.assertEqual(idxcol.nelements, self.nrows)

    def test02b(self):
        """Checking the indexing of a table in a 4th level
        (persistent version)"""

        # Create an instance of an HDF5 Table
        group = self.h5file.create_group(self.h5file.root, "agroup")
        group = self.h5file.create_group(group, "agroup")
        group = self.h5file.create_group(group, "agroup")

        # Create a table
        title = "This is the IndexArray title"
        table = self.h5file.create_table(
            group, "table", TDescr, title, None, self.nrows
        )
        for i in range(self.nrows):
            # Fill rows with defaults
            table.row.append()
        table.flush()

        # Index some column
        indexrows = table.cols.var1.create_index()
        self.assertIsNotNone(indexrows)
        idxcol = table.cols.var1.index

        # Close and re-open this file
        self._reopen(mode="a")

        table = self.h5file.root.agroup.agroup.agroup.table
        idxcol = table.cols.var1.index

        # Some sanity checks
        self.assertEqual(table.colindexed["var1"], 1)
        self.assertIsNotNone(idxcol)
        self.assertEqual(idxcol.nelements, self.nrows)

    def test03(self):
        """Checking the indexing of a table in a 100th level hierarchy"""

        # Create an instance of an HDF5 Table
        group = self.h5file.root
        for i in range(100):
            group = self.h5file.create_group(group, "agroup")

        # Create a table
        title = "This is the IndexArray title"
        table = self.h5file.create_table(
            group, "table", TDescr, title, None, self.nrows
        )
        for i in range(self.nrows):
            # Fill rows with defaults
            table.row.append()
        table.flush()

        # Index some column
        indexrows = table.cols.var1.create_index()
        self.assertIsNotNone(indexrows)
        idxcol = table.cols.var1.index

        # Some sanity checks
        self.assertEqual(table.colindexed["var1"], 1)
        self.assertIsNotNone(idxcol)
        self.assertEqual(idxcol.nelements, self.nrows)


class IndexProps:
    def __init__(
        self,
        auto=tb.index.default_auto_index,
        filters=tb.index.default_index_filters,
    ):
        self.auto = auto
        self.filters = filters


DefaultProps = IndexProps()
NoAutoProps = IndexProps(auto=False)
ChangeFiltersProps = IndexProps(
    filters=tb.Filters(
        complevel=6, complib="zlib", shuffle=False, fletcher32=False
    )
)


class AutomaticIndexingTestCase(common.TempFileMixin, common.PyTablesTestCase):
    reopen = 1
    iprops = NoAutoProps
    colsToIndex = ["var1", "var2", "var3"]
    small_blocksizes = (16, 8, 4, 2)

    def setUp(self):
        super().setUp()

        # Create an instance of an HDF5 Table
        title = "This is the IndexArray title"
        root = self.h5file.root

        # Make the chunkshape smaller or equal than small_blocksizes[-1]
        chunkshape = (2,)
        self.table = self.h5file.create_table(
            root,
            "table",
            TDescr,
            title,
            None,
            self.nrows,
            chunkshape=chunkshape,
        )
        self.table.autoindex = self.iprops.auto
        for colname in self.colsToIndex:
            self.table.colinstances[colname].create_index(
                _blocksizes=self.small_blocksizes
            )
        for i in range(self.nrows):
            # Fill rows with defaults
            self.table.row.append()
        self.table.flush()
        if self.reopen:
            self._reopen(mode="a")
            self.table = self.h5file.root.table

    def test01_attrs(self):
        """Checking indexing attributes (part1)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01_attrs..." % self.__class__.__name__)

        table = self.table
        if self.iprops is DefaultProps:
            self.assertEqual(table.indexed, 0)
        else:
            self.assertEqual(table.indexed, 1)
        if self.iprops is DefaultProps:
            self.assertEqual(table.colindexed["var1"], 0)
            self.assertIsNone(table.cols.var1.index)
            self.assertEqual(table.colindexed["var2"], 0)
            self.assertIsNone(table.cols.var2.index)
            self.assertEqual(table.colindexed["var3"], 0)
            self.assertIsNone(table.cols.var3.index)
            self.assertEqual(table.colindexed["var4"], 0)
            self.assertIsNone(table.cols.var4.index)
        else:
            # Check that the var1, var2 and var3 (and only these)
            # has been indexed
            self.assertEqual(table.colindexed["var1"], 1)
            self.assertIsNotNone(table.cols.var1.index)
            self.assertEqual(table.colindexed["var2"], 1)
            self.assertIsNotNone(table.cols.var2.index)
            self.assertEqual(table.colindexed["var3"], 1)
            self.assertIsNotNone(table.cols.var3.index)
            self.assertEqual(table.colindexed["var4"], 0)
            self.assertIsNone(table.cols.var4.index)

    def test02_attrs(self):
        """Checking indexing attributes (part2)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02_attrs..." % self.__class__.__name__)

        table = self.table

        # Check the policy parameters
        if common.verbose:
            if table.indexed:
                print("index props:", table.autoindex)
            else:
                print("Table is not indexed")

        # Check non-default values for index saving policy
        if self.iprops is NoAutoProps:
            self.assertFalse(table.autoindex)
        elif self.iprops is ChangeFiltersProps:
            self.assertTrue(table.autoindex)

        # Check Index() objects exists and are properly placed
        if self.iprops is DefaultProps:
            self.assertEqual(table.cols.var1.index, None)
            self.assertEqual(table.cols.var2.index, None)
            self.assertEqual(table.cols.var3.index, None)
            self.assertEqual(table.cols.var4.index, None)
        else:
            self.assertIsInstance(table.cols.var1.index, tb.index.Index)
            self.assertIsInstance(table.cols.var2.index, tb.index.Index)
            self.assertIsInstance(table.cols.var3.index, tb.index.Index)
            self.assertEqual(table.cols.var4.index, None)

    def test03_counters(self):
        """Checking indexing counters"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03_counters..." % self.__class__.__name__)
        table = self.table

        # Check the counters for indexes
        if common.verbose:
            if table.indexed:
                print("indexedrows:", table._indexedrows)
                print("unsavedindexedrows:", table._unsaved_indexedrows)
                index = table.cols.var1.index
                print("table rows:", table.nrows)
                print("computed indexed rows:", index.nrows * index.slicesize)
            else:
                print("Table is not indexed")
        if self.iprops is not DefaultProps:
            index = table.cols.var1.index
            indexedrows = index.nelements
            self.assertEqual(table._indexedrows, indexedrows)
            indexedrows = index.nelements
            self.assertEqual(
                table._unsaved_indexedrows, self.nrows - indexedrows
            )

    def test04_noauto(self):
        """Checking indexing counters (non-automatic mode)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test04_noauto..." % self.__class__.__name__)
        table = self.table

        # Force a sync in indexes
        table.flush_rows_to_index()

        # Check the counters for indexes
        if common.verbose:
            if table.indexed:
                print("indexedrows:", table._indexedrows)
                print("unsavedindexedrows:", table._unsaved_indexedrows)
                index = table.cols.var1.index
                print("computed indexed rows:", index.nelements)
            else:
                print("Table is not indexed")

        # No unindexated rows should remain
        index = table.cols.var1.index
        if self.iprops is DefaultProps:
            self.assertIsNone(index)
        else:
            indexedrows = index.nelements
            self.assertEqual(table._indexedrows, index.nelements)
            self.assertEqual(
                table._unsaved_indexedrows, self.nrows - indexedrows
            )

        # Check non-default values for index saving policy
        if self.iprops is NoAutoProps:
            self.assertFalse(table.autoindex)
        elif self.iprops is ChangeFiltersProps:
            self.assertTrue(table.autoindex)

    def test05_icounters(self):
        """Checking indexing counters (remove_rows)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test05_icounters..." % self.__class__.__name__)
        table = self.table

        # Force a sync in indexes
        table.flush_rows_to_index()

        # Non indexated rows should remain here
        if self.iprops is not DefaultProps:
            indexedrows = table._indexedrows
            unsavedindexedrows = table._unsaved_indexedrows

        # Now, remove some rows:
        table.remove_rows(2, 4)
        if self.reopen:
            self._reopen(mode="a")
            table = self.h5file.root.table

        # Check the counters for indexes
        if common.verbose:
            if table.indexed:
                print("indexedrows:", table._indexedrows)
                print("original indexedrows:", indexedrows)
                print("unsavedindexedrows:", table._unsaved_indexedrows)
                print("original unsavedindexedrows:", unsavedindexedrows)
                # index = table.cols.var1.index
                print("index dirty:", table.cols.var1.index.dirty)
            else:
                print("Table is not indexed")

        # Check the counters
        self.assertEqual(table.nrows, self.nrows - 2)
        if self.iprops is NoAutoProps:
            self.assertTrue(table.cols.var1.index.dirty)

        # Check non-default values for index saving policy
        if self.iprops is NoAutoProps:
            self.assertFalse(table.autoindex)
        elif self.iprops is ChangeFiltersProps:
            self.assertTrue(table.autoindex)

    def test06_dirty(self):
        """Checking dirty flags (remove_rows action)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test06_dirty..." % self.__class__.__name__)
        table = self.table

        # Force a sync in indexes
        table.flush_rows_to_index()

        # Now, remove some rows:
        table.remove_rows(3, 5)
        if self.reopen:
            self._reopen(mode="a")
            table = self.h5file.root.table

        # Check the dirty flag for indexes
        if common.verbose:
            print("auto flag:", table.autoindex)
            for colname in table.colnames:
                if table.cols._f_col(colname).index:
                    print(
                        "dirty flag col %s: %s"
                        % (colname, table.cols._f_col(colname).index.dirty)
                    )
        # Check the flags
        for colname in table.colnames:
            if table.cols._f_col(colname).index:
                if not table.autoindex:
                    self.assertEqual(
                        table.cols._f_col(colname).index.dirty, True
                    )
                else:
                    self.assertEqual(
                        table.cols._f_col(colname).index.dirty, False
                    )

    def test07_noauto(self):
        """Checking indexing counters (modify_rows, no-auto mode)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test07_noauto..." % self.__class__.__name__)
        table = self.table

        # Force a sync in indexes
        table.flush_rows_to_index()

        # No unindexated rows should remain here
        if self.iprops is not DefaultProps:
            indexedrows = table._indexedrows
            unsavedindexedrows = table._unsaved_indexedrows

        # Now, modify just one row:
        table.modify_rows(3, None, 1, [("asa", 0, 3, 3.1)])
        if self.reopen:
            self._reopen(mode="a")
            table = self.h5file.root.table

        # Check the counters for indexes
        if common.verbose:
            if table.indexed:
                print("indexedrows:", table._indexedrows)
                print("original indexedrows:", indexedrows)
                print("unsavedindexedrows:", table._unsaved_indexedrows)
                print("original unsavedindexedrows:", unsavedindexedrows)
                index = table.cols.var1.index
                print("computed indexed rows:", index.nelements)
            else:
                print("Table is not indexed")

        # Check the counters
        self.assertEqual(table.nrows, self.nrows)
        if self.iprops is NoAutoProps:
            self.assertTrue(table.cols.var1.index.dirty)

        # Check the dirty flag for indexes
        if common.verbose:
            for colname in table.colnames:
                if table.cols._f_col(colname).index:
                    print(
                        "dirty flag col %s: %s"
                        % (colname, table.cols._f_col(colname).index.dirty)
                    )
        for colname in table.colnames:
            if table.cols._f_col(colname).index:
                if not table.autoindex:
                    self.assertEqual(
                        table.cols._f_col(colname).index.dirty, True
                    )
                else:
                    self.assertEqual(
                        table.cols._f_col(colname).index.dirty, False
                    )

    def test07b_noauto(self):
        """Checking indexing queries (modify in iterator, no-auto mode)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test07b_noauto..." % self.__class__.__name__)
        table = self.table

        # Force a sync in indexes
        table.flush_rows_to_index()

        # Do a query that uses indexes
        res = [row.nrow for row in table.where("(var2 == True) & (var3 > 0)")]

        # Now, modify just one row:
        for row in table:
            if row.nrow == 3:
                row["var1"] = "asa"
                row["var2"] = True
                row["var3"] = 3
                row["var4"] = 3.1
                row.update()
        table.flush()
        if self.reopen:
            self._reopen(mode="a")
            table = self.h5file.root.table

        # Do a query that uses indexes
        resq = [row.nrow for row in table.where("(var2 == True) & (var3 > 0)")]
        res_ = res + [3]
        if common.verbose:
            print("AutoIndex?:", table.autoindex)
            print("Query results (original):", res)
            print("Query results (after modifying table):", resq)
            print("Should look like:", res_)
        self.assertEqual(res_, resq)

    def test07c_noauto(self):
        """Checking indexing queries (append, no-auto mode)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test07c_noauto..." % self.__class__.__name__)
        table = self.table

        # Force a sync in indexes
        table.flush_rows_to_index()

        # Do a query that uses indexes
        res = [row.nrow for row in table.where("(var2 == True) & (var3 > 0)")]

        # Now, append three rows
        table.append([("asa", True, 1, 3.1)])
        table.append([("asb", True, 2, 3.1)])
        table.append([("asc", True, 3, 3.1)])
        table.flush()
        if self.reopen:
            self._reopen(mode="a")
            table = self.h5file.root.table

        # Do a query that uses indexes
        resq = [row.nrow for row in table.where("(var2 == True) & (var3 > 0)")]
        res_ = res + [table.nrows - 3, table.nrows - 2, table.nrows - 1]
        if common.verbose:
            print("AutoIndex?:", table.autoindex)
            print("Query results (original):", res)
            print("Query results (after modifying table):", resq)
            print("Should look like:", res_)
        self.assertEqual(res_, resq)

    def test08_dirty(self):
        """Checking dirty flags (modify_columns)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test08_dirty..." % self.__class__.__name__)
        table = self.table

        # Force a sync in indexes
        table.flush_rows_to_index()

        # Non indexed rows should remain here
        if self.iprops is not DefaultProps:
            indexedrows = table._indexedrows
            self.assertIsNotNone(indexedrows)
            unsavedindexedrows = table._unsaved_indexedrows
            self.assertIsNotNone(unsavedindexedrows)

        # Now, modify a couple of rows:
        table.modify_columns(
            1, columns=[["asa", "asb"], [1.0, 2.0]], names=["var1", "var4"]
        )
        if self.reopen:
            self._reopen(mode="a")
            table = self.h5file.root.table

        # Check the counters
        self.assertEqual(table.nrows, self.nrows)
        if self.iprops is NoAutoProps:
            self.assertTrue(table.cols.var1.index.dirty)

        # Check the dirty flag for indexes
        if common.verbose:
            for colname in table.colnames:
                if table.cols._f_col(colname).index:
                    print(
                        "dirty flag col %s: %s"
                        % (colname, table.cols._f_col(colname).index.dirty)
                    )
        for colname in table.colnames:
            if table.cols._f_col(colname).index:
                if not table.autoindex:
                    if colname in ["var1"]:
                        self.assertEqual(
                            table.cols._f_col(colname).index.dirty, True
                        )
                    else:
                        self.assertEqual(
                            table.cols._f_col(colname).index.dirty, False
                        )
                else:
                    self.assertEqual(
                        table.cols._f_col(colname).index.dirty, False
                    )

    def test09a_propIndex(self):
        """Checking propagate Index feature in Table.copy() (attrs)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test09a_propIndex..." % self.__class__.__name__)
        table = self.table

        # Don't force a sync in indexes
        # table.flush_rows_to_index()
        # Non indexed rows should remain here
        if self.iprops is not DefaultProps:
            indexedrows = table._indexedrows
            self.assertIsNotNone(indexedrows)
            unsavedindexedrows = table._unsaved_indexedrows
            self.assertIsNotNone(unsavedindexedrows)

        # Now, remove some rows to make columns dirty
        # table.remove_rows(3,5)
        # Copy a Table to another location
        table2 = table.copy("/", "table2", propindexes=True)
        if self.reopen:
            self._reopen(mode="a")
            table = self.h5file.root.table
            table2 = self.h5file.root.table2

        index1 = table.cols.var1.index
        index2 = table2.cols.var1.index
        if common.verbose:
            print("Copied index:", index2)
            print("Original index:", index1)
            if index1:
                print("Elements in copied index:", index2.nelements)
                print("Elements in original index:", index1.nelements)

        # Check the counters
        self.assertEqual(table.nrows, table2.nrows)
        if table.indexed:
            self.assertTrue(table2.indexed)
        if self.iprops is DefaultProps:
            # No index: the index should not exist
            self.assertIsNone(index1)
            self.assertIsNone(index2)
        elif self.iprops is NoAutoProps:
            self.assertIsNotNone(index2)

        # Check the dirty flag for indexes
        if common.verbose:
            for colname in table2.colnames:
                if table2.cols._f_col(colname).index:
                    print(
                        "dirty flag col %s: %s"
                        % (colname, table2.cols._f_col(colname).index.dirty)
                    )
        for colname in table2.colnames:
            if table2.cols._f_col(colname).index:
                self.assertEqual(
                    table2.cols._f_col(colname).index.dirty, False
                )

    def test09b_propIndex(self):
        """Checking that propindexes=False works"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test09b_propIndex..." % self.__class__.__name__)
        table = self.table

        # Don't force a sync in indexes
        # table.flush_rows_to_index()
        # Non indexed rows should remain here
        if self.iprops is not DefaultProps:
            indexedrows = table._indexedrows
            self.assertIsNotNone(indexedrows)
            unsavedindexedrows = table._unsaved_indexedrows
            self.assertIsNotNone(unsavedindexedrows)

        # Now, remove some rows to make columns dirty
        # table.remove_rows(3,5)
        # Copy a Table to another location
        table2 = table.copy("/", "table2", propindexes=False)
        if self.reopen:
            self._reopen(mode="a")
            table = self.h5file.root.table
            table2 = self.h5file.root.table2

        if common.verbose:
            print("autoindex?:", self.iprops.auto)
            print("Copied index indexed?:", table2.cols.var1.is_indexed)
            print("Original index indexed?:", table.cols.var1.is_indexed)
        if self.iprops is DefaultProps:
            # No index: the index should not exist
            self.assertFalse(table2.cols.var1.is_indexed)
            self.assertFalse(table.cols.var1.is_indexed)
        elif self.iprops is NoAutoProps:
            self.assertFalse(table2.cols.var1.is_indexed)
            self.assertTrue(table.cols.var1.is_indexed)

    def test10_propIndex(self):
        """Checking propagate Index feature in Table.copy() (values)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test10_propIndex..." % self.__class__.__name__)
        table = self.table

        # Don't force a sync in indexes
        # table.flush_rows_to_index()
        # Non indexed rows should remain here
        if self.iprops is not DefaultProps:
            indexedrows = table._indexedrows
            self.assertIsNotNone(indexedrows)
            unsavedindexedrows = table._unsaved_indexedrows
            self.assertIsNotNone(unsavedindexedrows)

        # Now, remove some rows to make columns dirty
        # table.remove_rows(3,5)
        # Copy a Table to another location
        table2 = table.copy("/", "table2", propindexes=True)
        if self.reopen:
            self._reopen(mode="a")
            table = self.h5file.root.table
            table2 = self.h5file.root.table2

        index1 = table.cols.var3.index
        index2 = table2.cols.var3.index
        if common.verbose:
            print("Copied index:", index2)
            print("Original index:", index1)
            if index1:
                print("Elements in copied index:", index2.nelements)
                print("Elements in original index:", index1.nelements)

    def test11_propIndex(self):
        """Checking propagate Index feature in Table.copy() (dirty flags)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test11_propIndex..." % self.__class__.__name__)
        table = self.table

        # Force a sync in indexes
        table.flush_rows_to_index()

        # Non indexed rows should remain here
        if self.iprops is not DefaultProps:
            indexedrows = table._indexedrows
            self.assertIsNotNone(indexedrows)
            unsavedindexedrows = table._unsaved_indexedrows
            self.assertIsNotNone(unsavedindexedrows)

        # Now, modify an indexed column and an unindexed one
        # to make the "var1" dirty
        table.modify_columns(
            1, columns=[["asa", "asb"], [1.0, 2.0]], names=["var1", "var4"]
        )

        # Copy a Table to another location
        table2 = table.copy("/", "table2", propindexes=True)
        if self.reopen:
            self._reopen(mode="a")
            table = self.h5file.root.table
            table2 = self.h5file.root.table2

        index1 = table.cols.var1.index
        index2 = table2.cols.var1.index
        if common.verbose:
            print("Copied index:", index2)
            print("Original index:", index1)
            if index1:
                print("Elements in copied index:", index2.nelements)
                print("Elements in original index:", index1.nelements)

        # Check the dirty flag for indexes
        if common.verbose:
            for colname in table2.colnames:
                if table2.cols._f_col(colname).index:
                    print(
                        "dirty flag col %s: %s"
                        % (colname, table2.cols._f_col(colname).index.dirty)
                    )
        for colname in table2.colnames:
            if table2.cols._f_col(colname).index:
                if table2.autoindex:
                    # All the destination columns should be non-dirty because
                    # the copy removes the dirty state and puts the
                    # index in a sane state
                    self.assertEqual(
                        table2.cols._f_col(colname).index.dirty, False
                    )


# minRowIndex = 10000  # just if one wants more indexed rows to be checked
class AI1TestCase(AutomaticIndexingTestCase):
    # nrows = 10002
    nrows = 102
    reopen = 0
    iprops = NoAutoProps
    colsToIndex = ["var1", "var2", "var3"]


class AI2TestCase(AutomaticIndexingTestCase):
    # nrows = 10002
    nrows = 102
    reopen = 1
    iprops = NoAutoProps
    colsToIndex = ["var1", "var2", "var3"]


class AI4bTestCase(AutomaticIndexingTestCase):
    # nrows = 10012
    nrows = 112
    reopen = 1
    iprops = NoAutoProps
    colsToIndex = ["var1", "var2", "var3"]


class AI5TestCase(AutomaticIndexingTestCase):
    sbs, bs, ss, cs = tb.idxutils.calc_chunksize(minRowIndex, memlevel=1)
    nrows = ss * 11 - 1
    reopen = 0
    iprops = NoAutoProps
    colsToIndex = ["var1", "var2", "var3"]


class AI6TestCase(AutomaticIndexingTestCase):
    sbs, bs, ss, cs = tb.idxutils.calc_chunksize(minRowIndex, memlevel=1)
    nrows = ss * 21 + 1
    reopen = 1
    iprops = NoAutoProps
    colsToIndex = ["var1", "var2", "var3"]


class AI7TestCase(AutomaticIndexingTestCase):
    sbs, bs, ss, cs = tb.idxutils.calc_chunksize(minRowIndex, memlevel=1)
    nrows = ss * 12 - 1
    # nrows = ss * 1-1  # faster test
    reopen = 0
    iprops = NoAutoProps
    colsToIndex = ["var1", "var2", "var3"]


class AI8TestCase(AutomaticIndexingTestCase):
    sbs, bs, ss, cs = tb.idxutils.calc_chunksize(minRowIndex, memlevel=1)
    nrows = ss * 15 + 100
    # nrows = ss * 1 + 100  # faster test
    reopen = 1
    iprops = NoAutoProps
    colsToIndex = ["var1", "var2", "var3"]


class AI9TestCase(AutomaticIndexingTestCase):
    sbs, bs, ss, cs = tb.idxutils.calc_chunksize(minRowIndex, memlevel=1)
    nrows = ss
    reopen = 0
    iprops = DefaultProps
    colsToIndex = []


class AI10TestCase(AutomaticIndexingTestCase):
    # nrows = 10002
    nrows = 102
    reopen = 1
    iprops = DefaultProps
    colsToIndex = []


class AI11TestCase(AutomaticIndexingTestCase):
    # nrows = 10002
    nrows = 102
    reopen = 0
    iprops = ChangeFiltersProps
    colsToIndex = ["var1", "var2", "var3"]


class AI12TestCase(AutomaticIndexingTestCase):
    # nrows = 10002
    nrows = 102
    reopen = 0
    iprops = ChangeFiltersProps
    colsToIndex = ["var1", "var2", "var3"]


class ManyNodesTestCase(common.TempFileMixin, common.PyTablesTestCase):
    opem_kwargs = dict(node_cache_slots=64)

    def test00(self):
        """Indexing many nodes in one single session (based on bug #26)"""

        IdxRecord = {
            "f0": tb.Int8Col(),
            "f1": tb.Int8Col(),
            "f2": tb.Int8Col(),
        }

        for qn in range(5):
            for sn in range(5):
                qchr = "chr" + str(qn)
                name = "chr" + str(sn)
                path = "/at/%s/pt" % (qchr)
                table = self.h5file.create_table(
                    path, name, IdxRecord, createparents=1
                )
                table.cols.f0.create_index()
                table.cols.f1.create_index()
                table.cols.f2.create_index()
                table.row.append()
                table.flush()


class IndexPropsChangeTestCase(common.TempFileMixin, common.PyTablesTestCase):
    """Test case for changing index properties in a table."""

    class MyDescription(tb.IsDescription):
        icol = tb.IntCol()

    oldIndexProps = IndexProps()
    newIndexProps = IndexProps(auto=False, filters=tb.Filters(complevel=9))

    def setUp(self):
        super().setUp()

        table = self.h5file.create_table("/", "test", self.MyDescription)
        table.autoindex = self.oldIndexProps.auto
        row = table.row
        for i in range(100):
            row["icol"] = i % 25
            row.append()
        table.flush()
        self.table = table

    def test_attributes(self):
        """Storing index properties as table attributes."""
        for refprops in [self.oldIndexProps, self.newIndexProps]:
            self.assertEqual(self.table.autoindex, refprops.auto)
            self.table.autoindex = self.newIndexProps.auto

    def test_copyattrs(self):
        """Copying index properties attributes."""
        oldtable = self.table
        newtable = oldtable.copy("/", "test2")
        self.assertEqual(oldtable.autoindex, newtable.autoindex)


class IndexFiltersTestCase(common.TempFileMixin, common.PyTablesTestCase):
    """Test case for setting index filters."""

    def setUp(self):
        super().setUp()
        description = {"icol": tb.IntCol()}
        self.table = self.h5file.create_table("/", "test", description)

    def test_createIndex(self):
        """Checking input parameters in new indexes."""
        # Different from default.
        argfilters = copy.copy(tb.index.default_index_filters)
        argfilters.shuffle = not tb.index.default_index_filters.shuffle

        # Different both from default and the previous one.
        idxfilters = copy.copy(tb.index.default_index_filters)
        idxfilters.shuffle = not tb.index.default_index_filters.shuffle
        idxfilters.fletcher32 = not tb.index.default_index_filters.fletcher32

        icol = self.table.cols.icol

        # First create
        icol.create_index(kind="ultralight", optlevel=4)
        self.assertEqual(icol.index.kind, "ultralight")
        self.assertEqual(icol.index.optlevel, 4)
        self.assertEqual(icol.index.filters, tb.index.default_index_filters)
        icol.remove_index()

        # Second create
        icol.create_index(kind="medium", optlevel=3, filters=argfilters)
        self.assertEqual(icol.index.kind, "medium")
        self.assertEqual(icol.index.optlevel, 3)
        self.assertEqual(icol.index.filters, argfilters)
        icol.remove_index()

    def test_reindex(self):
        """Checking input parameters in recomputed indexes."""
        icol = self.table.cols.icol
        icol.create_index(
            kind="full", optlevel=5, filters=tb.Filters(complevel=3)
        )
        kind = icol.index.kind
        optlevel = icol.index.optlevel
        filters = icol.index.filters
        icol.reindex()
        ni = icol.index
        if common.verbose:
            print(f"Old parameters: {kind}, {optlevel}, {filters}")
            print(
                "New parameters: {}, {}, {}".format(
                    ni.kind, ni.optlevel, ni.filters
                )
            )
        self.assertEqual(ni.kind, kind)
        self.assertEqual(ni.optlevel, optlevel)
        self.assertEqual(ni.filters, filters)


class OldIndexTestCase(common.TestFileMixin, common.PyTablesTestCase):
    h5fname = common.test_filename("idx-std-1.x.h5")

    def test1_x(self):
        """Check that files with 1.x indexes are recognized and warned."""

        self.assertWarns(
            tb.exceptions.OldIndexWarning, self.h5file.get_node, "/table"
        )


# Sensible parameters for indexing with small blocksizes
small_blocksizes = (512, 128, 32, 8)


class CompletelySortedIndexTestCase(
    common.TempFileMixin, common.PyTablesTestCase
):
    """Test case for testing a complete sort in a table."""

    nrows = 100
    nrowsinbuf = 11

    class MyDescription(tb.IsDescription):
        rcol = tb.IntCol(pos=1)
        icol = tb.IntCol(pos=2)

    def setUp(self):
        super().setUp()
        table = self.h5file.create_table("/", "table", self.MyDescription)
        row = table.row
        nrows = self.nrows
        for i in range(nrows):
            row["rcol"] = i
            row["icol"] = nrows - i
            row.append()
        table.flush()
        self.table = table
        self.icol = self.table.cols.icol
        # A full index with maximum optlevel should always be completely sorted
        self.icol.create_csindex(_blocksizes=small_blocksizes)

    def test00_isCompletelySortedIndex(self):
        """Testing the Column.is_csi property."""

        icol = self.icol
        self.assertEqual(icol.index.is_csi, True)
        icol.remove_index()
        # Other kinds than full, should never return a CSI
        icol.create_index(kind="medium", optlevel=9)
        self.assertEqual(icol.index.is_csi, False)
        icol.remove_index()
        # As the table is small, lesser optlevels should be able to
        # create a completely sorted index too.
        icol.create_index(kind="full", optlevel=6)
        self.assertEqual(icol.index.is_csi, True)
        # Checking a CSI in a sorted copy
        self.table.copy("/", "table2", sortby="icol", checkCSI=True)
        self.assertEqual(icol.index.is_csi, True)

    def test01_readSorted1(self):
        """Testing the Index.read_sorted() method with no arguments."""

        icol = self.icol
        sortedcol = np.sort(icol[:])
        sortedcol2 = icol.index.read_sorted()
        if common.verbose:
            print("Original sorted column:", sortedcol)
            print("The values from the index:", sortedcol2)
        self.assertTrue(common.allequal(sortedcol, sortedcol2))

    def test01_readSorted2(self):
        """Testing the Index.read_sorted() method with arguments (I)."""

        icol = self.icol
        sortedcol = np.sort(icol[:])[30:55]
        sortedcol2 = icol.index.read_sorted(30, 55)
        if common.verbose:
            print("Original sorted column:", sortedcol)
            print("The values from the index:", sortedcol2)
        self.assertTrue(common.allequal(sortedcol, sortedcol2))

    def test01_readSorted3(self):
        """Testing the Index.read_sorted() method with arguments (II)."""

        icol = self.icol
        sortedcol = np.sort(icol[:])[33:97]
        sortedcol2 = icol.index.read_sorted(33, 97)
        if common.verbose:
            print("Original sorted column:", sortedcol)
            print("The values from the index:", sortedcol2)
        self.assertTrue(common.allequal(sortedcol, sortedcol2))

    def test02_readIndices1(self):
        """Testing the Index.read_indices() method with no arguments."""

        icol = self.icol
        indicescol = np.argsort(icol[:]).astype("uint64")
        indicescol2 = icol.index.read_indices()
        if common.verbose:
            print("Original indices column:", indicescol)
            print("The values from the index:", indicescol2)
        self.assertTrue(common.allequal(indicescol, indicescol2))

    def test02_readIndices2(self):
        """Testing the Index.read_indices() method with arguments (I)."""

        icol = self.icol
        indicescol = np.argsort(icol[:])[30:55].astype("uint64")
        indicescol2 = icol.index.read_indices(30, 55)
        if common.verbose:
            print("Original indices column:", indicescol)
            print("The values from the index:", indicescol2)
        self.assertTrue(common.allequal(indicescol, indicescol2))

    def test02_readIndices3(self):
        """Testing the Index.read_indices() method with arguments (II)."""

        icol = self.icol
        indicescol = np.argsort(icol[:])[33:97].astype("uint64")
        indicescol2 = icol.index.read_indices(33, 97)
        if common.verbose:
            print("Original indices column:", indicescol)
            print("The values from the index:", indicescol2)
        self.assertTrue(common.allequal(indicescol, indicescol2))

    def test02_readIndices4(self):
        """Testing the Index.read_indices() method with arguments (III)."""

        icol = self.icol
        indicescol = np.argsort(icol[:])[33:97:2].astype("uint64")
        indicescol2 = icol.index.read_indices(33, 97, 2)
        if common.verbose:
            print("Original indices column:", indicescol)
            print("The values from the index:", indicescol2)
        self.assertTrue(common.allequal(indicescol, indicescol2))

    def test02_readIndices5(self):
        """Testing the Index.read_indices() method with arguments (IV)."""

        icol = self.icol
        indicescol = np.argsort(icol[:])[33:55:5].astype("uint64")
        indicescol2 = icol.index.read_indices(33, 55, 5)
        if common.verbose:
            print("Original indices column:", indicescol)
            print("The values from the index:", indicescol2)
        self.assertTrue(common.allequal(indicescol, indicescol2))

    def test02_readIndices6(self):
        """Testing the Index.read_indices() method with step only."""

        icol = self.icol
        indicescol = np.argsort(icol[:])[::3].astype("uint64")
        indicescol2 = icol.index.read_indices(step=3)
        if common.verbose:
            print("Original indices column:", indicescol)
            print("The values from the index:", indicescol2)
        self.assertTrue(common.allequal(indicescol, indicescol2))

    def test03_getitem1(self):
        """Testing the Index.__getitem__() method with no arguments."""

        icol = self.icol
        indicescol = np.argsort(icol[:]).astype("uint64")
        indicescol2 = icol.index[:]
        if common.verbose:
            print("Original indices column:", indicescol)
            print("The values from the index:", indicescol2)
        self.assertTrue(common.allequal(indicescol, indicescol2))

    def test03_getitem2(self):
        """Testing the Index.__getitem__() method with start."""

        icol = self.icol
        indicescol = np.argsort(icol[:])[31].astype("uint64")
        indicescol2 = icol.index[31]
        if common.verbose:
            print("Original indices column:", indicescol)
            print("The values from the index:", indicescol2)
        self.assertTrue(common.allequal(indicescol, indicescol2))

    def test03_getitem3(self):
        """Testing the Index.__getitem__() method with start, stop."""

        icol = self.icol
        indicescol = np.argsort(icol[:])[2:16].astype("uint64")
        indicescol2 = icol.index[2:16]
        if common.verbose:
            print("Original indices column:", indicescol)
            print("The values from the index:", indicescol2)
        self.assertTrue(common.allequal(indicescol, indicescol2))

    def test04_itersorted1(self):
        """Testing the Table.itersorted() method with no arguments."""

        table = self.table
        sortedtable = np.sort(table[:], order="icol")
        sortedtable2 = np.array(
            [row.fetch_all_fields() for row in table.itersorted("icol")],
            dtype=table._v_dtype,
        )
        if common.verbose:
            print("Original sorted table:", sortedtable)
            print("The values from the iterator:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test04_itersorted2(self):
        """Testing the Table.itersorted() method with a start."""

        table = self.table
        sortedtable = np.sort(table[:], order="icol")[15:]
        sortedtable2 = np.array(
            [
                row.fetch_all_fields()
                for row in table.itersorted("icol", start=15)
            ],
            dtype=table._v_dtype,
        )
        if common.verbose:
            print("Original sorted table:", sortedtable)
            print("The values from the iterator:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test04_itersorted3(self):
        """Testing the Table.itersorted() method with a stop."""

        table = self.table
        sortedtable = np.sort(table[:], order="icol")[:20]
        sortedtable2 = np.array(
            [
                row.fetch_all_fields()
                for row in table.itersorted("icol", stop=20)
            ],
            dtype=table._v_dtype,
        )
        if common.verbose:
            print("Original sorted table:", sortedtable)
            print("The values from the iterator:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test04_itersorted4(self):
        """Testing the Table.itersorted() method with a start and stop."""

        table = self.table
        sortedtable = np.sort(table[:], order="icol")[15:20]
        sortedtable2 = np.array(
            [
                row.fetch_all_fields()
                for row in table.itersorted("icol", start=15, stop=20)
            ],
            dtype=table._v_dtype,
        )
        if common.verbose:
            print("Original sorted table:", sortedtable)
            print("The values from the iterator:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test04_itersorted5(self):
        """Testing the Table.itersorted() method with a start, stop and
        step."""

        table = self.table
        sortedtable = np.sort(table[:], order="icol")[15:45:4]
        sortedtable2 = np.array(
            [
                row.fetch_all_fields()
                for row in table.itersorted("icol", start=15, stop=45, step=4)
            ],
            dtype=table._v_dtype,
        )
        if common.verbose:
            print("Original sorted table:", sortedtable)
            print("The values from the iterator:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test04_itersorted6(self):
        """Testing the Table.itersorted() method with a start, stop and
        step."""

        table = self.table
        sortedtable = np.sort(table[:], order="icol")[33:55:5]
        sortedtable2 = np.array(
            [
                row.fetch_all_fields()
                for row in table.itersorted("icol", start=33, stop=55, step=5)
            ],
            dtype=table._v_dtype,
        )
        if common.verbose:
            print("Original sorted table:", sortedtable)
            print("The values from the iterator:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test04_itersorted7(self):
        """Testing the Table.itersorted() method with checkCSI=True."""

        table = self.table
        sortedtable = np.sort(table[:], order="icol")
        sortedtable2 = np.array(
            [
                row.fetch_all_fields()
                for row in table.itersorted("icol", checkCSI=True)
            ],
            dtype=table._v_dtype,
        )
        if common.verbose:
            print("Original sorted table:", sortedtable)
            print("The values from the iterator:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test04_itersorted8(self):
        """Testing the Table.itersorted() method with a start, stop and
        negative step."""

        # see also gh-252
        table = self.table
        sortedtable = np.sort(table[:], order="icol")[55:33:-5]
        sortedtable2 = np.array(
            [
                row.fetch_all_fields()
                for row in table.itersorted("icol", start=55, stop=33, step=-5)
            ],
            dtype=table._v_dtype,
        )
        if common.verbose:
            print("Original sorted table:", sortedtable)
            print("The values from the iterator:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test04_itersorted9(self):
        """Testing the Table.itersorted() method with a negative step -5."""

        # see also gh-252
        table = self.table
        sortedtable = np.sort(table[:], order="icol")[::-5]
        sortedtable2 = np.array(
            [
                row.fetch_all_fields()
                for row in table.itersorted("icol", step=-5)
            ],
            dtype=table._v_dtype,
        )
        if common.verbose:
            print("Original sorted table:", sortedtable)
            print("The values from the iterator:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test04_itersorted10(self):
        """Testing the Table.itersorted() method with a negative step -1."""

        # see also gh-252
        table = self.table
        sortedtable = np.sort(table[:], order="icol")[::-1]
        sortedtable2 = np.array(
            [
                row.fetch_all_fields()
                for row in table.itersorted("icol", step=-1)
            ],
            dtype=table._v_dtype,
        )
        if common.verbose:
            print("Original sorted table:", sortedtable)
            print("The values from the iterator:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test05_readSorted1(self):
        """Testing the Table.read_sorted() method with no arguments."""

        table = self.table
        sortedtable = np.sort(table[:], order="icol")
        sortedtable2 = table.read_sorted("icol")
        if common.verbose:
            print("Original sorted table:", sortedtable)
            print("The values from read_sorted:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test05_readSorted2(self):
        """Testing the Table.read_sorted() method with a start."""

        table = self.table
        sortedtable = np.sort(table[:], order="icol")[16:17]
        sortedtable2 = table.read_sorted("icol", start=16)
        if common.verbose:
            print("Original sorted table:", sortedtable)
            print("The values from read_sorted:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test05_readSorted3(self):
        """Testing the Table.read_sorted() method with a start and stop."""

        table = self.table
        sortedtable = np.sort(table[:], order="icol")[16:33]
        sortedtable2 = table.read_sorted("icol", start=16, stop=33)
        if common.verbose:
            print("Original sorted table:", sortedtable)
            print("The values from read_sorted:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test05_readSorted4(self):
        """Testing the Table.read_sorted() method with a start, stop and
        step."""

        table = self.table
        sortedtable = np.sort(table[:], order="icol")[33:55:5]
        sortedtable2 = table.read_sorted("icol", start=33, stop=55, step=5)
        if common.verbose:
            print("Original sorted table:", sortedtable)
            print("The values from read_sorted:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test05_readSorted5(self):
        """Testing the Table.read_sorted() method with only a step."""

        table = self.table
        sortedtable = np.sort(table[:], order="icol")[::3]
        sortedtable2 = table.read_sorted("icol", step=3)
        if common.verbose:
            print("Original sorted table:", sortedtable)
            print("The values from read_sorted:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test05_readSorted6(self):
        """Testing the Table.read_sorted() method with negative step."""

        table = self.table
        sortedtable = np.sort(table[:], order="icol")[::-1]
        sortedtable2 = table.read_sorted("icol", step=-1)
        if common.verbose:
            print("Original sorted table:", sortedtable)
            print("The values from read_sorted:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test05_readSorted7(self):
        """Testing the Table.read_sorted() method with negative step (II)."""

        table = self.table
        sortedtable = np.sort(table[:], order="icol")[::-2]
        sortedtable2 = table.read_sorted("icol", step=-2)
        if common.verbose:
            print("Original sorted table:", sortedtable)
            print("The values from read_sorted:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test05_readSorted8(self):
        """Testing the Table.read_sorted() method with negative step (III))."""

        table = self.table
        sstart = 100 - 24 - 1
        sstop = 100 - 54 - 1
        sortedtable = np.sort(table[:], order="icol")[sstart:sstop:-1]
        sortedtable2 = table.read_sorted("icol", start=24, stop=54, step=-1)
        if common.verbose:
            print("Original sorted table:", sortedtable)
            print("The values from read_sorted:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test05_readSorted9(self):
        """Testing the Table.read_sorted() method with negative step (IV))."""

        table = self.table
        sstart = 100 - 14 - 1
        sstop = 100 - 54 - 1
        sortedtable = np.sort(table[:], order="icol")[sstart:sstop:-3]
        sortedtable2 = table.read_sorted("icol", start=14, stop=54, step=-3)
        if common.verbose:
            print("Original sorted table:", sortedtable)
            print("The values from read_sorted:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test05_readSorted10(self):
        """Testing the Table.read_sorted() method with negative step (V))."""

        table = self.table
        sstart = 100 - 24 - 1
        sstop = 100 - 25 - 1
        sortedtable = np.sort(table[:], order="icol")[sstart:sstop:-2]
        sortedtable2 = table.read_sorted("icol", start=24, stop=25, step=-2)
        if common.verbose:
            print("Original sorted table:", sortedtable)
            print("The values from read_sorted:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test05_readSorted11(self):
        """Testing the Table.read_sorted() method with start > stop."""

        table = self.table
        sstart = 100 - 137 - 1
        sstop = 100 - 25 - 1
        sortedtable = np.sort(table[:], order="icol")[sstart:sstop:-2]
        sortedtable2 = table.read_sorted("icol", start=137, stop=25, step=-2)
        if common.verbose:
            print("Original sorted table:", sortedtable)
            print("The values from read_sorted:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test05a_readSorted12(self):
        """Testing the Table.read_sorted() method with checkCSI (I)."""

        table = self.table
        sortedtable = np.sort(table[:], order="icol")
        sortedtable2 = table.read_sorted("icol", checkCSI=True)
        if common.verbose:
            print("Original sorted table:", sortedtable)
            print("The values from read_sorted:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test05b_readSorted12(self):
        """Testing the Table.read_sorted() method with checkCSI (II)."""

        table = self.table
        self.assertRaises(
            ValueError, table.read_sorted, "rcol", checkCSI=False
        )

    def test06_copy_sorted1(self):
        """Testing the Table.copy(sortby) method with no arguments."""

        table = self.table
        # Copy to another table
        table.nrowsinbuf = self.nrowsinbuf
        table2 = table.copy("/", "table2", sortby="icol")
        sortedtable = np.sort(table[:], order="icol")
        sortedtable2 = table2[:]
        if common.verbose:
            print("Original sorted table:", sortedtable)
            print("The values from copy:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test06_copy_sorted2(self):
        """Testing the Table.copy(sortby) method with step=-1."""

        table = self.table
        # Copy to another table
        table.nrowsinbuf = self.nrowsinbuf
        table2 = table.copy("/", "table2", sortby="icol", step=-1)
        sortedtable = np.sort(table[:], order="icol")[::-1]
        sortedtable2 = table2[:]
        if common.verbose:
            print("Original sorted table:", sortedtable)
            print("The values from copy:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test06_copy_sorted3(self):
        """Testing the Table.copy(sortby) method with only a start."""

        table = self.table
        # Copy to another table
        table.nrowsinbuf = self.nrowsinbuf
        table2 = table.copy("/", "table2", sortby="icol", start=3)
        sortedtable = np.sort(table[:], order="icol")[3:4]
        sortedtable2 = table2[:]
        if common.verbose:
            print("Original sorted table:", sortedtable)
            print("The values from copy:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test06_copy_sorted4(self):
        """Testing the Table.copy(sortby) method with start, stop."""

        table = self.table
        # Copy to another table
        table.nrowsinbuf = self.nrowsinbuf
        table2 = table.copy("/", "table2", sortby="icol", start=3, stop=40)
        sortedtable = np.sort(table[:], order="icol")[3:40]
        sortedtable2 = table2[:]
        if common.verbose:
            print("Original sorted table:", sortedtable)
            print("The values from copy:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test06_copy_sorted5(self):
        """Testing the Table.copy(sortby) method with start, stop, step."""

        table = self.table
        # Copy to another table
        table.nrowsinbuf = self.nrowsinbuf
        table2 = table.copy(
            "/", "table2", sortby="icol", start=3, stop=33, step=5
        )
        sortedtable = np.sort(table[:], order="icol")[3:33:5]
        sortedtable2 = table2[:]
        if common.verbose:
            print("Original sorted table:", sortedtable)
            print("The values from copy:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test06_copy_sorted6(self):
        """Testing the Table.copy(sortby) method after table re-opening."""

        self._reopen(mode="a")
        table = self.h5file.root.table
        # Copy to another table
        table.nrowsinbuf = self.nrowsinbuf
        table2 = table.copy("/", "table2", sortby="icol")
        sortedtable = np.sort(table[:], order="icol")
        sortedtable2 = table2[:]
        if common.verbose:
            print("Original sorted table:", sortedtable)
            print("The values from copy:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test06_copy_sorted7(self):
        """Testing the `checkCSI` parameter of Table.copy() (I)."""

        table = self.table

        # Copy to another table
        table.nrowsinbuf = self.nrowsinbuf
        table2 = table.copy("/", "table2", sortby="icol")
        self.assertRaises(
            ValueError,
            table2.copy,
            "/",
            "table3",
            sortby="rcol",
            checkCSI=False,
        )

    def test06_copy_sorted8(self):
        """Testing the `checkCSI` parameter of Table.copy() (II)."""

        table = self.table

        # Copy to another table
        table.nrowsinbuf = self.nrowsinbuf
        table2 = table.copy("/", "table2", sortby="icol")
        self.assertRaises(
            ValueError,
            table2.copy,
            "/",
            "table3",
            sortby="rcol",
            checkCSI=True,
        )

    def test07_isCSI_noelements(self):
        """Testing the representation of an index with no elements."""

        t2 = self.h5file.create_table("/", "t2", self.MyDescription)
        irows = t2.cols.rcol.create_csindex()
        if common.verbose:
            print("repr(t2)-->\n", repr(t2))
        self.assertEqual(irows, 0)
        self.assertEqual(t2.colindexes["rcol"].is_csi, False)


class ReadSortedIndexTestCase(common.TempFileMixin, common.PyTablesTestCase):
    """Test case for testing sorted reading in a "full" sorted column."""

    nrows = 100
    nrowsinbuf = 11

    class MyDescription(tb.IsDescription):
        rcol = tb.IntCol(pos=1)
        icol = tb.IntCol(pos=2)

    def setUp(self):
        super().setUp()

        table = self.h5file.create_table("/", "table", self.MyDescription)
        row = table.row
        nrows = self.nrows
        for i in range(nrows):
            row["rcol"] = i
            row["icol"] = nrows - i
            row.append()
        table.flush()
        self.table = table
        self.icol = self.table.cols.icol
        # A full index with maximum optlevel should always be completely sorted
        self.icol.create_index(
            optlevel=self.optlevel, kind="full", _blocksizes=small_blocksizes
        )

    def test01_readSorted1(self):
        """Testing the Table.read_sorted() method with no arguments."""

        table = self.table
        sortedtable = np.sort(table[:], order="icol")
        sortedtable2 = table.read_sorted("icol")
        if common.verbose:
            print("Sorted table:", sortedtable)
            print("The values from read_sorted:", sortedtable2)
        # Compare with the sorted read table because we have no
        # guarantees that read_sorted returns a completely sorted table
        self.assertTrue(
            common.allequal(sortedtable, np.sort(sortedtable2, order="icol"))
        )

    def test01_readSorted2(self):
        """Testing the Table.read_sorted() method with no arguments
        (re-open)."""

        self._reopen()
        table = self.h5file.root.table
        sortedtable = np.sort(table[:], order="icol")
        sortedtable2 = table.read_sorted("icol")
        if common.verbose:
            print("Sorted table:", sortedtable)
            print("The values from read_sorted:", sortedtable2)
        # Compare with the sorted read table because we have no
        # guarantees that read_sorted returns a completely sorted table
        self.assertTrue(
            common.allequal(sortedtable, np.sort(sortedtable2, order="icol"))
        )

    def test02_copy_sorted1(self):
        """Testing the Table.copy(sortby) method."""

        table = self.table
        # Copy to another table
        table.nrowsinbuf = self.nrowsinbuf
        table2 = table.copy("/", "table2", sortby="icol")
        sortedtable = np.sort(table[:], order="icol")
        sortedtable2 = np.sort(table2[:], order="icol")
        if common.verbose:
            print("Original table:", table2[:])
            print("The sorted values from copy:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))

    def test02_copy_sorted2(self):
        """Testing the Table.copy(sortby) method after table re-opening."""

        self._reopen(mode="a")
        table = self.h5file.root.table
        # Copy to another table
        table.nrowsinbuf = self.nrowsinbuf
        table2 = table.copy("/", "table2", sortby="icol")
        sortedtable = np.sort(table[:], order="icol")
        sortedtable2 = np.sort(table2[:], order="icol")
        if common.verbose:
            print("Original table:", table2[:])
            print("The sorted values from copy:", sortedtable2)
        self.assertTrue(common.allequal(sortedtable, sortedtable2))


class ReadSortedIndex0(ReadSortedIndexTestCase):
    optlevel = 0


class ReadSortedIndex3(ReadSortedIndexTestCase):
    optlevel = 3


class ReadSortedIndex6(ReadSortedIndexTestCase):
    optlevel = 6


class ReadSortedIndex9(ReadSortedIndexTestCase):
    optlevel = 9


class Issue156TestBase(common.TempFileMixin, common.PyTablesTestCase):
    # field name in table according to which test_copysort() sorts the table
    sort_field = None

    def setUp(self):
        super().setUp()

        # create nested table
        class Foo(tb.IsDescription):
            frame = tb.UInt16Col()

            class Bar(tb.IsDescription):
                code = tb.UInt16Col()

        table = self.h5file.create_table(
            "/", "foo", Foo, filters=tb.Filters(3, "zlib"), createparents=True
        )

        self.h5file.flush()

        # fill table with 10 random numbers
        for k in range(10):
            row = table.row
            row["frame"] = np.random.randint(0, 2**16 - 1)
            row["Bar/code"] = np.random.randint(0, 2**16 - 1)
            row.append()

        self.h5file.flush()

    def test_copysort(self):
        # copy table
        oldNode = self.h5file.get_node("/foo")

        # create completely sorted index on a main column
        oldNode.colinstances[self.sort_field].create_csindex()

        # this fails on ade2ba123efd267fd31
        # see gh-156
        new_node = oldNode.copy(
            newname="foo2",
            overwrite=True,
            sortby=self.sort_field,
            checkCSI=True,
            propindexes=True,
        )

        # check column is sorted
        self.assertTrue(
            np.all(
                new_node.col(self.sort_field)
                == sorted(oldNode.col(self.sort_field))
            )
        )
        # check index is available
        self.assertIn(self.sort_field, new_node.colindexes)
        # check CSI was propagated
        self.assertTrue(new_node.colindexes[self.sort_field].is_csi)


class Issue156TestCase01(Issue156TestBase):
    # sort by field from non nested entry
    sort_field = "frame"


class Issue156TestCase02(Issue156TestBase):
    # sort by field from nested entry
    sort_field = "Bar/code"


class Issue119Time32ColTestCase(common.TempFileMixin, common.PyTablesTestCase):
    """TimeCol not properly indexing."""

    col_typ = tb.Time32Col
    values = [
        0.93240451618785880,
        0.76322375510776170,
        0.16695030056300875,
        0.91259117097807850,
        0.93977847053454630,
        0.51450406513503090,
        0.24452129962257563,
        0.85475938924825230,
        0.32512326762476930,
        0.75127635627046820,
    ]

    def setUp(self):
        super().setUp()

        class Descr(tb.IsDescription):
            when = self.col_typ(pos=1)
            value = tb.Float32Col(pos=2)

        self.table = self.h5file.create_table("/", "test", Descr)

        self.t = 1321031471.0  # 11/11/11 11:11:11
        data = [(self.t + i, item) for i, item in enumerate(self.values)]
        self.table.append(data)
        self.h5file.flush()

    def test_timecol_issue(self):
        tbl = self.table
        t = self.t

        wherestr = "(when >= %d) & (when < %d)" % (t, t + 5)

        no_index = tbl.read_where(wherestr)

        tbl.cols.when.create_index(_verbose=False)
        with_index = tbl.read_where(wherestr)

        self.assertTrue((no_index == with_index).all())


class Issue119Time64ColTestCase(Issue119Time32ColTestCase):
    col_typ = tb.Time64Col


class TestIndexingNans(common.TempFileMixin, common.PyTablesTestCase):
    def test_issue_282(self):
        trMap = {"index": tb.Int64Col(), "values": tb.FloatCol()}
        table = self.h5file.create_table("/", "table", trMap)

        r = table.row
        for i in range(5):
            r["index"] = i
            r["values"] = np.nan if i == 0 else i
            r.append()
        table.flush()

        table.cols.values.create_index()

        # retrieve
        result = table.read_where("(values >= 0)")
        self.assertEqual(len(result), 4)

    def test_issue_327(self):
        table = self.h5file.create_table(
            "/",
            "table",
            dict(
                index=tb.Int64Col(),
                values=tb.FloatCol(shape=()),
                values2=tb.FloatCol(shape=()),
            ),
        )

        r = table.row
        for i in range(5):
            r["index"] = i
            r["values"] = np.nan if i == 2 or i == 3 else i
            r["values2"] = i
            r.append()
        table.flush()

        table.cols.values.create_index()
        table.cols.values2.create_index()

        results2 = table.read_where("(values2 > 0)")
        self.assertEqual(len(results2), 4)

        results = table.read_where("(values > 0)")
        self.assertEqual(len(results), 2)

    def test_issue_327_b(self):
        table = self.h5file.create_table(
            "/",
            "table",
            dict(
                index=tb.Int64Col(),
                values=tb.FloatCol(shape=()),
                values2=tb.FloatCol(shape=()),
            ),
        )

        r = table.row
        for _ in range(100):
            for i in range(5):
                r["index"] = i
                r["values"] = np.nan if i == 2 or i == 3 else i
                r["values2"] = i
                r.append()
        table.flush()

        table.cols.values.create_index(_blocksizes=small_blocksizes)
        table.cols.values2.create_index(_blocksizes=small_blocksizes)

        results2 = table.read_where("(values2 > 0)")
        self.assertEqual(len(results2), 400)

        results = table.read_where("(values > 0)")
        self.assertEqual(len(results), 200)

    def test_csindex_nans(self):
        table = self.h5file.create_table(
            "/",
            "table",
            dict(
                index=tb.Int64Col(),
                values=tb.FloatCol(shape=()),
                values2=tb.FloatCol(shape=()),
            ),
        )

        r = table.row
        for x in range(100):
            for i in range(5):
                r["index"] = i
                r["values"] = np.nan if i == 2 or i == 3 else i
                r["values2"] = i
                r.append()
        table.flush()

        table.cols.values.create_csindex(_blocksizes=small_blocksizes)
        table.cols.values2.create_csindex(_blocksizes=small_blocksizes)

        results2 = table.read_where("(values2 > 0)")
        self.assertEqual(len(results2), 100 * 4)

        results = table.read_where("(values > 0)")
        self.assertEqual(len(results), 100 * 2)


def suite():
    theSuite = common.unittest.TestSuite()

    niter = 1
    # heavy = 1  # Uncomment this only for testing purposes!

    for n in range(niter):
        theSuite.addTest(common.make_suite(BasicReadTestCase))
        theSuite.addTest(common.make_suite(ZlibReadTestCase))
        theSuite.addTest(common.make_suite(BloscReadTestCase))
        theSuite.addTest(common.make_suite(LZOReadTestCase))
        theSuite.addTest(common.make_suite(Bzip2ReadTestCase))
        theSuite.addTest(common.make_suite(ShuffleReadTestCase))
        theSuite.addTest(common.make_suite(Fletcher32ReadTestCase))
        theSuite.addTest(common.make_suite(ShuffleFletcher32ReadTestCase))
        theSuite.addTest(common.make_suite(OneHalfTestCase))
        theSuite.addTest(common.make_suite(UpperBoundTestCase))
        theSuite.addTest(common.make_suite(LowerBoundTestCase))
        theSuite.addTest(common.make_suite(AI1TestCase))
        theSuite.addTest(common.make_suite(AI2TestCase))
        theSuite.addTest(common.make_suite(AI9TestCase))
        theSuite.addTest(common.make_suite(DeepTableIndexTestCase))
        theSuite.addTest(common.make_suite(IndexPropsChangeTestCase))
        theSuite.addTest(common.make_suite(IndexFiltersTestCase))
        theSuite.addTest(common.make_suite(OldIndexTestCase))
        theSuite.addTest(common.make_suite(CompletelySortedIndexTestCase))
        theSuite.addTest(common.make_suite(ManyNodesTestCase))
        theSuite.addTest(common.make_suite(ReadSortedIndex0))
        theSuite.addTest(common.make_suite(ReadSortedIndex3))
        theSuite.addTest(common.make_suite(ReadSortedIndex6))
        theSuite.addTest(common.make_suite(ReadSortedIndex9))
        theSuite.addTest(common.make_suite(Issue156TestCase01))
        theSuite.addTest(common.make_suite(Issue156TestCase02))
        theSuite.addTest(common.make_suite(Issue119Time32ColTestCase))
        theSuite.addTest(common.make_suite(Issue119Time64ColTestCase))
        theSuite.addTest(common.make_suite(TestIndexingNans))
    if common.heavy:
        # These are too heavy for normal testing
        theSuite.addTest(common.make_suite(AI4bTestCase))
        theSuite.addTest(common.make_suite(AI5TestCase))
        theSuite.addTest(common.make_suite(AI6TestCase))
        theSuite.addTest(common.make_suite(AI7TestCase))
        theSuite.addTest(common.make_suite(AI8TestCase))
        theSuite.addTest(common.make_suite(AI10TestCase))
        theSuite.addTest(common.make_suite(AI11TestCase))
        theSuite.addTest(common.make_suite(AI12TestCase))

    return theSuite


if __name__ == "__main__":
    import sys

    common.parse_argv(sys.argv)
    common.print_versions()
    common.unittest.main(defaultTest="suite")
