import sys
import tempfile
import warnings
from time import perf_counter as clock
from pathlib import Path

import tables as tb
from tables.tests import common


# Test Record class
class Record(tb.IsDescription):
    var1 = tb.StringCol(itemsize=4)  # 4-character String
    var2 = tb.IntCol()  # integer
    var3 = tb.Int16Col()  # short integer
    var4 = tb.FloatCol()  # double (double-precision)
    var5 = tb.Float32Col()  # float  (single-precision)


class TreeTestCase(common.TempFileMixin, common.PyTablesTestCase):
    open_mode = "w"
    title = "This is the table title"
    expectedrows = 10
    appendrows = 5

    def setUp(self):
        super().setUp()

        # Create an instance of HDF5 Table
        self.populateFile()
        self.h5file.close()

    def populateFile(self):
        group = self.h5file.root
        maxshort = 1 << 15
        # maxint   = 2147483647   # (2 ** 31 - 1)
        for j in range(3):
            # Create a table
            table = self.h5file.create_table(
                group,
                "table" + str(j),
                Record,
                title=self.title,
                filters=None,
                expectedrows=self.expectedrows,
            )
            # Get the record object associated with the new table
            d = table.row
            # Fill the table
            for i in range(self.expectedrows):
                d["var1"] = "%04d" % (self.expectedrows - i)
                d["var2"] = i
                d["var3"] = i % maxshort
                d["var4"] = float(i)
                d["var5"] = float(i)
                d.append()  # This injects the Record values
            # Flush the buffer for this table
            table.flush()

            # Create a couple of arrays in each group
            var1List = [x["var1"] for x in table.iterrows()]
            var4List = [x["var4"] for x in table.iterrows()]

            self.h5file.create_array(group, "var1", var1List, "1")
            self.h5file.create_array(group, "var4", var4List, "4")

            # Create a new group (descendant of group)
            group2 = self.h5file.create_group(group, "group" + str(j))
            # Iterate over this new group (group2)
            group = group2

    def test00_getNode(self):
        """Checking the File.get_node() with string node names"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test00_getNode..." % self.__class__.__name__)

        self.h5file = tb.open_file(self.h5fname, "r")
        nodelist = ["/", "/table0", "/group0/var1", "/group0/group1/var4"]
        nodenames = []
        for node in nodelist:
            object = self.h5file.get_node(node)
            nodenames.append(object._v_pathname)

        self.assertEqual(nodenames, nodelist)
        if common.verbose:
            print("get_node(pathname) test passed")
        nodegroups = [
            "/",
            "/group0",
            "/group0/group1",
            "/group0/group1/group2",
        ]
        nodenames = ["var1", "var4"]
        nodepaths = []
        for group in nodegroups:
            for name in nodenames:
                try:
                    object = self.h5file.get_node(group, name)
                except LookupError:
                    pass
                else:
                    nodepaths.append(object._v_pathname)

        self.assertEqual(
            nodepaths,
            [
                "/var1",
                "/var4",
                "/group0/var1",
                "/group0/var4",
                "/group0/group1/var1",
                "/group0/group1/var4",
            ],
        )

        if common.verbose:
            print("get_node(groupname, name) test passed")
        nodelist = [
            "/",
            "/group0",
            "/group0/group1",
            "/group0/group1/group2",
            "/table0",
        ]
        nodenames = []
        groupobjects = []
        # warnings.filterwarnings("error", category=UserWarning)
        for node in nodelist:
            try:
                object = self.h5file.get_node(node, classname="Group")
            except LookupError:
                if common.verbose:
                    (type, value, traceback) = sys.exc_info()
                    print("\nGreat!, the next LookupError was catched!")
                    print(value)
            else:
                nodenames.append(object._v_pathname)
                groupobjects.append(object)

        self.assertEqual(
            nodenames,
            ["/", "/group0", "/group0/group1", "/group0/group1/group2"],
        )
        if common.verbose:
            print("get_node(groupname, classname='Group') test passed")

        # Reset the warning
        # warnings.filterwarnings("default", category=UserWarning)

        nodenames = ["var1", "var4"]
        nodearrays = []
        for group in groupobjects:
            for name in nodenames:
                try:
                    object = self.h5file.get_node(group, name, "Array")
                except Exception:
                    pass
                else:
                    nodearrays.append(object._v_pathname)

        self.assertEqual(
            nodearrays,
            [
                "/var1",
                "/var4",
                "/group0/var1",
                "/group0/var4",
                "/group0/group1/var1",
                "/group0/group1/var4",
            ],
        )
        if common.verbose:
            print("get_node(groupobject, name, classname='Array') test passed")

    def test01_getNodeClass(self):
        """Checking the File.get_node() with instances"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test01_getNodeClass..." % self.__class__.__name__
            )

        self.h5file = tb.open_file(self.h5fname, "r")

        # This tree ways of get_node usage should return a table instance
        table = self.h5file.get_node("/group0/table1")
        self.assertIsInstance(table, tb.Table)
        table = self.h5file.get_node("/group0", "table1")
        self.assertIsInstance(table, tb.Table)
        table = self.h5file.get_node(self.h5file.root.group0, "table1")
        self.assertIsInstance(table, tb.Table)

        # This should return an array instance
        arr = self.h5file.get_node("/group0/var1")
        self.assertIsInstance(arr, tb.Array)
        self.assertIsInstance(arr, tb.Leaf)

        # And this a Group
        group = self.h5file.get_node("/group0", "group1", "Group")
        self.assertIsInstance(group, tb.Group)

    def test02_listNodes(self):
        """Checking the File.list_nodes() method"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02_listNodes..." % self.__class__.__name__)

        # Made the warnings to raise an error
        # warnings.filterwarnings("error", category=UserWarning)
        self.h5file = tb.open_file(self.h5fname, "r")

        self.assertRaises(
            TypeError, self.h5file.list_nodes, "/", "NoSuchClass"
        )

        nodelist = [
            "/",
            "/group0",
            "/group0/table1",
            "/group0/group1/group2",
            "/var1",
        ]
        nodenames = []
        objects = []
        for node in nodelist:
            try:
                objectlist = self.h5file.list_nodes(node)
            except Exception:
                pass
            else:
                objects.extend(objectlist)
                for object in objectlist:
                    nodenames.append(object._v_pathname)

        self.assertEqual(
            nodenames,
            [
                "/group0",
                "/table0",
                "/var1",
                "/var4",
                "/group0/group1",
                "/group0/table1",
                "/group0/var1",
                "/group0/var4",
            ],
        )
        if common.verbose:
            print("list_nodes(pathname) test passed")

        nodenames = []
        for node in objects:
            try:
                objectlist = self.h5file.list_nodes(node)
            except Exception:
                pass
            else:
                for object in objectlist:
                    nodenames.append(object._v_pathname)

        self.assertEqual(
            nodenames,
            [
                "/group0/group1",
                "/group0/table1",
                "/group0/var1",
                "/group0/var4",
                "/group0/group1/group2",
                "/group0/group1/table2",
                "/group0/group1/var1",
                "/group0/group1/var4",
            ],
        )

        if common.verbose:
            print("list_nodes(groupobject) test passed")

        nodenames = []
        for node in objects:
            try:
                objectlist = self.h5file.list_nodes(node, "Leaf")
            except TypeError:
                if common.verbose:
                    (type, value, traceback) = sys.exc_info()
                    print("\nGreat!, the next TypeError was catched!")
                    print(value)
            else:
                for object in objectlist:
                    nodenames.append(object._v_pathname)

        self.assertEqual(
            nodenames,
            [
                "/group0/table1",
                "/group0/var1",
                "/group0/var4",
                "/group0/group1/table2",
                "/group0/group1/var1",
                "/group0/group1/var4",
            ],
        )

        if common.verbose:
            print("list_nodes(groupobject, classname = 'Leaf') test passed")

        nodenames = []
        for node in objects:
            try:
                objectlist = self.h5file.list_nodes(node, "Table")
            except TypeError:
                if common.verbose:
                    (type, value, traceback) = sys.exc_info()
                    print("\nGreat!, the next TypeError was catched!")
                    print(value)
            else:
                for object in objectlist:
                    nodenames.append(object._v_pathname)

        self.assertEqual(
            nodenames, ["/group0/table1", "/group0/group1/table2"]
        )

        if common.verbose:
            print("list_nodes(groupobject, classname = 'Table') test passed")

        # Reset the warning
        # warnings.filterwarnings("default", category=UserWarning)

    def test02b_iterNodes(self):
        """Checking the File.iter_nodes() method"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02b_iterNodes..." % self.__class__.__name__)

        self.h5file = tb.open_file(self.h5fname, "r")

        self.assertRaises(
            TypeError, self.h5file.list_nodes, "/", "NoSuchClass"
        )

        nodelist = [
            "/",
            "/group0",
            "/group0/table1",
            "/group0/group1/group2",
            "/var1",
        ]
        nodenames = []
        objects = []
        for node in nodelist:
            try:
                objectlist = [o for o in self.h5file.iter_nodes(node)]
            except Exception:
                pass
            else:
                objects.extend(objectlist)
                for object in objectlist:
                    nodenames.append(object._v_pathname)

        self.assertEqual(
            nodenames,
            [
                "/group0",
                "/table0",
                "/var1",
                "/var4",
                "/group0/group1",
                "/group0/table1",
                "/group0/var1",
                "/group0/var4",
            ],
        )
        if common.verbose:
            print("iter_nodes(pathname) test passed")

        nodenames = []
        for node in objects:
            try:
                objectlist = [o for o in self.h5file.iter_nodes(node)]
            except Exception:
                pass
            else:
                for object in objectlist:
                    nodenames.append(object._v_pathname)

        self.assertEqual(
            nodenames,
            [
                "/group0/group1",
                "/group0/table1",
                "/group0/var1",
                "/group0/var4",
                "/group0/group1/group2",
                "/group0/group1/table2",
                "/group0/group1/var1",
                "/group0/group1/var4",
            ],
        )

        if common.verbose:
            print("iter_nodes(groupobject) test passed")

        nodenames = []
        for node in objects:
            try:
                objectlist = [o for o in self.h5file.iter_nodes(node, "Leaf")]
            except TypeError:
                if common.verbose:
                    (type, value, traceback) = sys.exc_info()
                    print("\nGreat!, the next TypeError was catched!")
                    print(value)
            else:
                for object in objectlist:
                    nodenames.append(object._v_pathname)

        self.assertEqual(
            nodenames,
            [
                "/group0/table1",
                "/group0/var1",
                "/group0/var4",
                "/group0/group1/table2",
                "/group0/group1/var1",
                "/group0/group1/var4",
            ],
        )

        if common.verbose:
            print("iter_nodes(groupobject, classname = 'Leaf') test passed")

        nodenames = []
        for node in objects:
            try:
                objectlist = [o for o in self.h5file.iter_nodes(node, "Table")]
            except TypeError:
                if common.verbose:
                    (type, value, traceback) = sys.exc_info()
                    print("\nGreat!, the next TypeError was catched!")
                    print(value)
            else:
                for object in objectlist:
                    nodenames.append(object._v_pathname)

        self.assertEqual(
            nodenames, ["/group0/table1", "/group0/group1/table2"]
        )

        if common.verbose:
            print("iter_nodes(groupobject, classname = 'Table') test passed")

        # Reset the warning
        # warnings.filterwarnings("default", category=UserWarning)

    def test03_TraverseTree(self):
        """Checking the File.walk_groups() method"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test03_TraverseTree..." % self.__class__.__name__
            )

        self.h5file = tb.open_file(self.h5fname, "r")
        groups = []
        tables_ = []
        arrays = []
        for group in self.h5file.walk_groups():
            groups.append(group._v_pathname)
            for table in self.h5file.list_nodes(group, "Table"):
                tables_.append(table._v_pathname)
            for arr in self.h5file.list_nodes(group, "Array"):
                arrays.append(arr._v_pathname)

        self.assertEqual(
            groups, ["/", "/group0", "/group0/group1", "/group0/group1/group2"]
        )

        self.assertEqual(
            tables_, ["/table0", "/group0/table1", "/group0/group1/table2"]
        )

        self.assertEqual(
            arrays,
            [
                "/var1",
                "/var4",
                "/group0/var1",
                "/group0/var4",
                "/group0/group1/var1",
                "/group0/group1/var4",
            ],
        )
        if common.verbose:
            print("walk_groups() test passed")

        groups = []
        tables_ = []
        arrays = []
        for group in self.h5file.walk_groups("/group0/group1"):
            groups.append(group._v_pathname)
            for table in self.h5file.list_nodes(group, "Table"):
                tables_.append(table._v_pathname)
            for arr in self.h5file.list_nodes(group, "Array"):
                arrays.append(arr._v_pathname)

        self.assertEqual(groups, ["/group0/group1", "/group0/group1/group2"])

        self.assertEqual(tables_, ["/group0/group1/table2"])

        self.assertEqual(
            arrays, ["/group0/group1/var1", "/group0/group1/var4"]
        )

        if common.verbose:
            print("walk_groups(pathname) test passed")

    def test04_walkNodes(self):
        """Checking File.walk_nodes"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test04_walkNodes..." % self.__class__.__name__)

        self.h5file = tb.open_file(self.h5fname, "r")

        self.assertRaises(
            TypeError, next, self.h5file.walk_nodes("/", "NoSuchClass")
        )

        groups = []
        tables1 = []
        tables2 = []
        arrays = []
        for group in self.h5file.walk_nodes(classname="Group"):
            groups.append(group._v_pathname)
            for table in group._f_iter_nodes(classname="Table"):
                tables1.append(table._v_pathname)

        # Test the recursivity
        for table in self.h5file.root._f_walknodes("Table"):
            tables2.append(table._v_pathname)

        for arr in self.h5file.walk_nodes(classname="Array"):
            arrays.append(arr._v_pathname)

        self.assertEqual(
            groups, ["/", "/group0", "/group0/group1", "/group0/group1/group2"]
        )
        self.assertEqual(
            tables1, ["/table0", "/group0/table1", "/group0/group1/table2"]
        )
        self.assertEqual(
            tables2, ["/table0", "/group0/table1", "/group0/group1/table2"]
        )
        self.assertEqual(
            arrays,
            [
                "/var1",
                "/var4",
                "/group0/var1",
                "/group0/var4",
                "/group0/group1/var1",
                "/group0/group1/var4",
            ],
        )

        if common.verbose:
            print("File.__iter__() and Group.__iter__ test passed")

        groups = []
        tables_ = []
        arrays = []
        for group in self.h5file.walk_nodes(
            "/group0/group1", classname="Group"
        ):
            groups.append(group._v_pathname)
            for table in group._f_walknodes("Table"):
                tables_.append(table._v_pathname)
            for arr in self.h5file.walk_nodes(group, "Array"):
                arrays.append(arr._v_pathname)

        self.assertEqual(groups, ["/group0/group1", "/group0/group1/group2"])

        self.assertEqual(tables_, ["/group0/group1/table2"])

        self.assertEqual(
            arrays, ["/group0/group1/var1", "/group0/group1/var4"]
        )

        if common.verbose:
            print("walk_nodes(pathname, classname) test passed")

    def test05_dir(self):
        """Checking Group.__dir__"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test05_dir..." % self.__class__.__name__)

        self.h5file = tb.open_file(self.h5fname, "r")

        """
        h5file nodes:
        '/table0', '/var1', '/var4'
        '/group0/table1', '/group0/var1', '/group0/var4',
        '/group0/group1/table2', '/group0/group1/var1', '/group0/group1/var4'
        """
        root_dir = dir(self.h5file.root)

        # Check some regular attributes.

        self.assertIn("_v_children", root_dir)
        self.assertIn("_v_attrs", root_dir)
        self.assertIn("_v_groups", root_dir)
        self.assertIn("_g_get_child_group_class", root_dir)
        self.assertIn("_g_get_child_group_class", root_dir)
        self.assertIn("_f_close", root_dir)

        # Check children nodes.

        self.assertIn("group0", root_dir)
        self.assertIn("table0", root_dir)
        self.assertIn("var1", root_dir)
        self.assertNotIn("table1", root_dir)
        self.assertNotIn("table2", root_dir)
        self.assertSequenceEqual(
            sorted(set(root_dir)), sorted(root_dir)
        )  # Check for no duplicates.

        root_group0_dir = dir(self.h5file.root.group0)
        self.assertIn("group1", root_group0_dir)
        self.assertIn("table1", root_group0_dir)
        self.assertNotIn("table0", root_group0_dir)
        self.assertNotIn("table2", root_group0_dir)
        self.assertSequenceEqual(
            sorted(set(root_group0_dir)), sorted(root_group0_dir)
        )

        root_group0_group1_dir = dir(self.h5file.root.group0.group1)
        self.assertIn("group2", root_group0_group1_dir)
        self.assertIn("table2", root_group0_group1_dir)
        self.assertNotIn("table0", root_group0_group1_dir)
        self.assertNotIn("table1", root_group0_group1_dir)
        self.assertNotIn("group0", root_group0_group1_dir)
        self.assertNotIn("group1", root_group0_group1_dir)
        self.assertSequenceEqual(
            sorted(set(root_group0_group1_dir)), sorted(root_group0_group1_dir)
        )

        if common.verbose:
            print("Group.__dir__ test passed")

    def test06_v_groups(self):
        """Checking Group._v_groups"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test06_v_groups..." % self.__class__.__name__)

        self.h5file = tb.open_file(self.h5fname, "r")

        """
        h5file nodes:
        '/table0', '/var1', '/var4'
        '/group0/table1', '/group0/var1', '/group0/var4',
        '/group0/group1/table2', '/group0/group1/var1', '/group0/group1/var4'
        """
        self.assertIsInstance(self.h5file.root._v_groups, dict)
        group_names = {"group0"}
        names = {k for k, v in self.h5file.root._v_groups.iteritems()}
        self.assertEqual(group_names, names)
        groups = list(self.h5file.root._v_groups.itervalues())
        self.assertEqual(len(groups), len(group_names))

        for group in groups:
            with self.subTest(name=group._v_name):
                self.assertIn(group._v_name, group_names)

        if common.verbose:
            print("Group._v_groups test passed")


class DeepTreeTestCase(common.TempFileMixin, common.PyTablesTestCase):
    """Checks for deep hierarchy levels in PyTables trees."""

    def setUp(self):
        super().setUp()

        # Here we put a more conservative limit to deal with more platforms
        # With maxdepth = 64 this test would take less than 40 MB
        # of main memory to run, which is quite reasonable nowadays.
        # With maxdepth = 1024 this test will take around 300 MB.
        if common.heavy:
            self.maxdepth = 256  # Takes around 60 MB of memory!
        else:
            self.maxdepth = 64  # This should be safe for most machines
        if common.verbose:
            print("Maximum depth tested :", self.maxdepth)

        # Open a new empty HDF5 file
        group = self.h5file.root
        if common.verbose:
            print("Depth writing progress: ", end=" ")

        # Iterate until maxdepth
        for depth in range(self.maxdepth):
            # Save it on the HDF5 file
            if common.verbose:
                print("%3d," % (depth), end=" ")
            # Create a couple of arrays here
            self.h5file.create_array(
                group, "array", [1, 1], "depth: %d" % depth
            )
            self.h5file.create_array(
                group, "array2", [1, 1], "depth: %d" % depth
            )
            # And also a group
            self.h5file.create_group(group, "group2_" + str(depth))
            # Finally, iterate over a new group
            group = self.h5file.create_group(group, "group" + str(depth))

        # Close the file
        self.h5file.close()

    def _check_tree(self, filename):
        # Open the previous HDF5 file in read-only mode

        with tb.open_file(filename, mode="r") as h5file:
            group = h5file.root
            if common.verbose:
                print("\nDepth reading progress: ", end=" ")

            # Get the metadata on the previosly saved arrays
            for depth in range(self.maxdepth):
                if common.verbose:
                    print("%3d," % (depth), end=" ")

                # Check the contents
                self.assertEqual(group.array[:], [1, 1])
                self.assertIn("array2", group)
                self.assertIn("group2_" + str(depth), group)

                # Iterate over the next group
                group = h5file.get_node(group, "group" + str(depth))

            if common.verbose:
                print()  # This flush the stdout buffer

    def test00_deepTree(self):
        """Creation of a large depth object tree."""

        self._check_tree(self.h5fname)

    def test01a_copyDeepTree(self):
        """Copy of a large depth object tree."""

        self.h5file = tb.open_file(self.h5fname, mode="r")
        h5fname2 = tempfile.mktemp(".h5")
        try:
            with tb.open_file(h5fname2, mode="w") as h5file2:
                if common.verbose:
                    print("\nCopying deep tree...")

                self.h5file.copy_node(
                    self.h5file.root, h5file2.root, recursive=True
                )
                self.h5file.close()

            self._check_tree(h5fname2)
        finally:
            if Path(h5fname2).is_file():
                Path(h5fname2).unlink()

    def test01b_copyDeepTree(self):
        """Copy of a large depth object tree with small node cache."""

        self.h5file = tb.open_file(self.h5fname, mode="r", node_cache_slots=10)
        h5fname2 = tempfile.mktemp(".h5")
        try:
            with tb.open_file(
                h5fname2, mode="w", node_cache_slots=10
            ) as h5file2:
                if common.verbose:
                    print("\nCopying deep tree...")

                self.h5file.copy_node(
                    self.h5file.root, h5file2.root, recursive=True
                )
                self.h5file.close()

            self._check_tree(h5fname2)
        finally:
            if Path(h5fname2).is_file():
                Path(h5fname2).unlink()

    def test01c_copyDeepTree(self):
        """Copy of a large depth object tree with no node cache."""

        self.h5file = tb.open_file(self.h5fname, mode="r", node_cache_slots=0)
        h5fname2 = tempfile.mktemp(".h5")
        try:
            with tb.open_file(
                h5fname2, mode="w", node_cache_slots=0
            ) as h5file2:
                if common.verbose:
                    print("\nCopying deep tree...")

                self.h5file.copy_node(
                    self.h5file.root, h5file2.root, recursive=True
                )
                self.h5file.close()

            self._check_tree(h5fname2)
        finally:
            if Path(h5fname2).is_file():
                Path(h5fname2).unlink()

    @common.unittest.skipUnless(common.heavy, "only in heavy mode")
    def test01d_copyDeepTree(self):
        """Copy of a large depth object tree with static node cache."""

        self.h5file = tb.open_file(
            self.h5fname, mode="r", node_cache_slots=-256
        )
        h5fname2 = tempfile.mktemp(".h5")
        try:
            with tb.open_file(
                h5fname2, mode="w", node_cache_slots=-256
            ) as h5file2:
                if common.verbose:
                    print("\nCopying deep tree...")

                self.h5file.copy_node(
                    self.h5file.root, h5file2.root, recursive=True
                )
                self.h5file.close()

            self._check_tree(h5fname2)
        finally:
            if Path(h5fname2).is_file():
                Path(h5fname2).unlink()


class WideTreeTestCase(common.TempFileMixin, common.PyTablesTestCase):
    """Checks for maximum number of children for a Group."""

    def test00_Leafs(self):
        """Checking creation of large number of leafs (1024) per group.

        Variable 'maxchildren' controls this check. PyTables support up
        to 4096 children per group, but this would take too much memory
        (up to 64 MB) for testing purposes (maybe we can add a test for
        big platforms). A 1024 children run takes up to 30 MB. A 512
        children test takes around 25 MB.

        """

        if common.heavy:
            maxchildren = 4096
        else:
            maxchildren = 256

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test00_wideTree..." % self.__class__.__name__)
            print("Maximum number of children tested :", maxchildren)

        a = [1, 1]
        if common.verbose:
            print("Children writing progress: ", end=" ")
        for child in range(maxchildren):
            if common.verbose:
                print("%3d," % (child), end=" ")
            self.h5file.create_array(
                self.h5file.root, "array" + str(child), a, "child: %d" % child
            )
        if common.verbose:
            print()

        t1 = clock()
        a = [1, 1]

        # Open the previous HDF5 file in read-only mode
        self._reopen()
        if common.verbose:
            print(
                "\nTime spent opening a file with %d arrays: %s s"
                % (maxchildren, clock() - t1)
            )
            print("\nChildren reading progress: ", end=" ")

        # Get the metadata on the previosly saved arrays
        for child in range(maxchildren):
            if common.verbose:
                print("%3d," % (child), end=" ")

            # Create an array for later comparison
            # Get the actual array
            array_ = getattr(self.h5file.root, "array" + str(child))
            b = array_.read()

            # Arrays a and b must be equal
            self.assertEqual(a, b)
        if common.verbose:
            print()  # This flush the stdout buffer

    def test01_wideTree(self):
        """Checking creation of large number of groups (1024) per group.

        Variable 'maxchildren' controls this check. PyTables support up
        to 4096 children per group, but this would take too much memory
        (up to 64 MB) for testing purposes (maybe we can add a test for
        big platforms). A 1024 children run takes up to 30 MB. A 512
        children test takes around 25 MB.

        """

        if common.heavy:
            # for big platforms!
            maxchildren = 4096
        else:
            # for standard platforms
            maxchildren = 256

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test00_wideTree..." % self.__class__.__name__)
            print("Maximum number of children tested :", maxchildren)

        if common.verbose:
            print("Children writing progress: ", end=" ")
        for child in range(maxchildren):
            if common.verbose:
                print("%3d," % (child), end=" ")
            self.h5file.create_group(
                self.h5file.root, "group" + str(child), "child: %d" % child
            )
        if common.verbose:
            print()

        t1 = clock()

        # Open the previous HDF5 file in read-only mode
        self._reopen()
        if common.verbose:
            print(
                "\nTime spent opening a file with %d groups: %s s"
                % (maxchildren, clock() - t1)
            )
            print("\nChildren reading progress: ", end=" ")

        # Get the metadata on the previosly saved arrays
        for child in range(maxchildren):
            if common.verbose:
                print("%3d," % (child), end=" ")
            # Get the actual group
            group = getattr(self.h5file.root, "group" + str(child))
            # Arrays a and b must be equal
            self.assertEqual(group._v_title, "child: %d" % child)

        if common.verbose:
            print()  # This flush the stdout buffer


class HiddenTreeTestCase(common.TempFileMixin, common.PyTablesTestCase):
    """Check for hidden groups, leaves and hierarchies."""

    def setUp(self):
        super().setUp()

        self.visible = []  # list of visible object paths
        self.hidden = []  # list of hidden object paths

        # Create some visible nodes: a, g, g/a1, g/a2, g/g, g/g/a.
        h5f = self.h5file
        h5f.create_array("/", "a", [0])
        g = h5f.create_group("/", "g")
        h5f.create_array(g, "a1", [0])
        h5f.create_array(g, "a2", [0])
        g_g = h5f.create_group(g, "g")
        h5f.create_array(g_g, "a", [0])

        self.visible.extend(["/a", "/g", "/g/a1", "/g/a2", "/g/g", "/g/g/a"])

        # Create some hidden nodes: _p_a, _p_g, _p_g/a, _p_g/_p_a, g/_p_a.
        h5f.create_array("/", "_p_a", [0])
        hg = h5f.create_group("/", "_p_g")
        h5f.create_array(hg, "a", [0])
        h5f.create_array(hg, "_p_a", [0])
        h5f.create_array(g, "_p_a", [0])

        self.hidden.extend(
            ["/_p_a", "/_p_g", "/_p_g/a", "/_p_g/_p_a", "/g/_p_a"]
        )

    # The test behind commented out because the .objects dictionary
    # has been removed (as well as .leaves and .groups)
    def _test00_objects(self):
        """Absence of hidden nodes in `File.objects`."""

        objects = self.h5file.objects

        warnings.filterwarnings("ignore", category=DeprecationWarning)

        for vpath in self.visible:
            self.assertTrue(
                vpath in objects,
                "Missing visible node ``%s`` from ``File.objects``." % vpath,
            )
        for hpath in self.hidden:
            self.assertTrue(
                hpath not in objects,
                "Found hidden node ``%s`` in ``File.objects``." % hpath,
            )

        warnings.filterwarnings("default", category=DeprecationWarning)

    # The test behind commented out because the .objects dictionary
    # has been removed (as well as .leaves and .groups)
    def _test00b_objects(self):
        """Object dictionaries conformance with ``walk_nodes()``."""

        def dictCheck(dictName, classname):
            file_ = self.h5file

            objects = getattr(file_, dictName)
            walkPaths = [
                node._v_pathname for node in file_.walk_nodes("/", classname)
            ]
            dictPaths = [path for path in objects]
            walkPaths.sort()
            dictPaths.sort()
            self.assertEqual(
                walkPaths,
                dictPaths,
                "nodes in ``%s`` do not match those from ``walk_nodes()``"
                % dictName,
            )
            self.assertEqual(
                len(walkPaths),
                len(objects),
                "length of ``%s`` differs from that of ``walk_nodes()``"
                % dictName,
            )

        warnings.filterwarnings("ignore", category=DeprecationWarning)

        dictCheck("objects", None)
        dictCheck("groups", "Group")
        dictCheck("leaves", "Leaf")

        warnings.filterwarnings("default", category=DeprecationWarning)

    def test01_getNode(self):
        """Node availability via `File.get_node()`."""

        h5f = self.h5file

        for vpath in self.visible:
            h5f.get_node(vpath)
        for hpath in self.hidden:
            h5f.get_node(hpath)

    def test02_walkGroups(self):
        """Hidden group absence in `File.walk_groups()`."""

        hidden = self.hidden

        for group in self.h5file.walk_groups("/"):
            pathname = group._v_pathname
            self.assertNotIn(
                pathname, hidden, f"Walked across hidden group ``{pathname}``."
            )

    def test03_walkNodes(self):
        """Hidden node absence in `File.walk_nodes()`."""

        hidden = self.hidden

        for node in self.h5file.walk_nodes("/"):
            pathname = node._v_pathname
            self.assertNotIn(
                pathname, hidden, f"Walked across hidden node ``{pathname}``."
            )

    def test04_listNodesVisible(self):
        """Listing visible nodes under a visible group (list_nodes)."""

        hidden = self.hidden

        for node in self.h5file.list_nodes("/g"):
            pathname = node._v_pathname
            self.assertNotIn(
                pathname, hidden, f"Listed hidden node ``{pathname}``."
            )

    def test04b_listNodesVisible(self):
        """Listing visible nodes under a visible group (iter_nodes)."""

        hidden = self.hidden

        for node in self.h5file.iter_nodes("/g"):
            pathname = node._v_pathname
            self.assertNotIn(
                pathname, hidden, f"Listed hidden node ``{pathname}``."
            )

    def test05_listNodesHidden(self):
        """Listing visible nodes under a hidden group (list_nodes)."""

        hidden = self.hidden

        node_to_find = "/_p_g/a"
        found_node = False
        for node in self.h5file.list_nodes("/_p_g"):
            pathname = node._v_pathname
            if pathname == node_to_find:
                found_node = True
            self.assertIn(
                pathname, hidden, f"Listed hidden node ``{pathname}``."
            )

        self.assertTrue(
            found_node, "Hidden node ``%s`` was not listed." % node_to_find
        )

    def test05b_iterNodesHidden(self):
        """Listing visible nodes under a hidden group (iter_nodes)."""

        hidden = self.hidden

        node_to_find = "/_p_g/a"
        found_node = False
        for node in self.h5file.iter_nodes("/_p_g"):
            pathname = node._v_pathname
            if pathname == node_to_find:
                found_node = True
            self.assertIn(
                pathname, hidden, f"Listed hidden node ``{pathname}``."
            )

        self.assertTrue(
            found_node, "Hidden node ``%s`` was not listed." % node_to_find
        )

    # The test behind commented out because the .objects dictionary
    # has been removed (as well as .leaves and .groups)
    def _test06_reopen(self):
        """Reopening a file with hidden nodes."""

        self.h5file.close()
        self.h5file = tb.open_file(self.h5fname)
        self.test00_objects()

    def test07_move(self):
        """Moving a node between hidden and visible groups."""

        is_visible_node = self.h5file.is_visible_node

        self.assertFalse(is_visible_node("/_p_g/a"))
        self.h5file.move_node("/_p_g/a", "/g", "a")
        self.assertTrue(is_visible_node("/g/a"))
        self.h5file.move_node("/g/a", "/_p_g", "a")
        self.assertFalse(is_visible_node("/_p_g/a"))

    def test08_remove(self):
        """Removing a visible group with hidden children."""

        self.assertIn("/g/_p_a", self.h5file)
        self.h5file.root.g._f_remove(recursive=True)
        self.assertNotIn("/g/_p_a", self.h5file)


class CreateParentsTestCase(common.TempFileMixin, common.PyTablesTestCase):
    """Test the ``createparents`` flag.

    These are mainly for the user interface.  More thorough tests on the
    workings of the flag can be found in the ``test_do_undo.py`` module.

    """

    filters = tb.Filters(complevel=4)  # simply non-default

    def setUp(self):
        super().setUp()
        self.h5file.create_array("/", "array", [1])
        self.h5file.create_group("/", "group", filters=self.filters)

    def test01_inside(self):
        """Placing a node inside a nonexistent child of itself."""
        self.assertRaises(
            tb.NodeError,
            self.h5file.move_node,
            "/group",
            "/group/foo/bar",
            createparents=True,
        )
        self.assertNotIn("/group/foo", self.h5file)
        self.assertRaises(
            tb.NodeError,
            self.h5file.copy_node,
            "/group",
            "/group/foo/bar",
            recursive=True,
            createparents=True,
        )
        self.assertNotIn("/group/foo", self.h5fname)

    def test02_filters(self):
        """Propagating the filters of created parent groups."""

        self.h5file.create_group("/group/foo/bar", "baz", createparents=True)
        self.assertIn("/group/foo/bar/baz", self.h5file)
        for group in self.h5file.walk_groups("/group"):
            self.assertEqual(self.filters, group._v_filters)


def suite():
    theSuite = common.unittest.TestSuite()
    # This counter is useful when detecting memory leaks
    niter = 1

    for i in range(niter):
        theSuite.addTest(common.make_suite(TreeTestCase))
        theSuite.addTest(common.make_suite(DeepTreeTestCase))
        theSuite.addTest(common.make_suite(WideTreeTestCase))
        theSuite.addTest(common.make_suite(HiddenTreeTestCase))
        theSuite.addTest(common.make_suite(CreateParentsTestCase))

    return theSuite


if __name__ == "__main__":
    common.parse_argv(sys.argv)
    common.print_versions()
    common.unittest.main(defaultTest="suite")
