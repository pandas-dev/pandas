import os
import sys
import queue
import shutil
import platform
import tempfile
import warnings
import threading
import subprocess
from pathlib import Path

import tables

try:
    import multiprocessing as mp

    multiprocessing_imported = True
except ImportError:
    multiprocessing_imported = False

import numpy as np

import tables as tb
from tables.tests import common


class OpenFileFailureTestCase(common.PyTablesTestCase):
    def setUp(self):
        super().setUp()

        self.N = len(tb.file._open_files)
        self.open_files = tb.file._open_files

    def test01_open_file(self):
        """Checking opening of a non existing file."""

        h5fname = tempfile.mktemp(".h5")
        with self.assertRaises(IOError):
            h5file = tb.open_file(h5fname)
            h5file.close()

        self.assertEqual(self.N, len(self.open_files))

    def test02_open_file(self):
        """Checking opening of an existing non HDF5 file."""

        # create a dummy file
        h5fname = tempfile.mktemp(".h5")
        Path(h5fname).write_text("")

        # Try to open the dummy file
        try:
            with self.assertRaises(tb.HDF5ExtError):
                h5file = tb.open_file(h5fname)
                h5file.close()

            self.assertEqual(self.N, len(self.open_files))
        finally:
            Path(h5fname).unlink()

    def test03_open_file(self):
        """Checking opening of an existing file with invalid mode."""

        # See gh-318

        # create a dummy file
        h5fname = tempfile.mktemp(".h5")
        h5file = tb.open_file(h5fname, "w")
        h5file.close()

        try:
            # Try to open the dummy file
            self.assertRaises(ValueError, tb.open_file, h5fname, "ab")
        finally:
            Path(h5fname).unlink()


class OpenFileTestCase(common.TempFileMixin, common.PyTablesTestCase):

    def setUp(self):
        super().setUp()
        self.populateFile()

    def populateFile(self):
        root = self.h5file.root

        # Create an array
        self.h5file.create_array(root, "array", [1, 2], title="Array example")
        self.h5file.create_table(
            root, "table", {"var1": tb.IntCol()}, "Table example"
        )
        root._v_attrs.testattr = 41

        # Create another array object
        self.h5file.create_array(root, "anarray", [1], "Array title")
        self.h5file.create_table(
            root, "atable", {"var1": tb.IntCol()}, "Table title"
        )

        # Create a group object
        group = self.h5file.create_group(root, "agroup", "Group title")
        group._v_attrs.testattr = 42

        # Create a some objects there
        array1 = self.h5file.create_array(
            group, "anarray1", [1, 2, 3, 4, 5, 6, 7], "Array title 1"
        )
        array1.attrs.testattr = 42
        self.h5file.create_array(group, "anarray2", [2], "Array title 2")
        self.h5file.create_table(
            group, "atable1", {"var1": tb.IntCol()}, "Table title 1"
        )
        ra = np.rec.array([(1, 11, "a")], formats="u1,f4,S1")
        self.h5file.create_table(group, "atable2", ra, "Table title 2")

        # Create a lonely group in first level
        self.h5file.create_group(root, "agroup2", "Group title 2")

        # Create a new group in the second level
        group3 = self.h5file.create_group(group, "agroup3", "Group title 3")

        # Create a new group in the third level
        self.h5file.create_group(group3, "agroup4", "Group title 4")

        # Create an array in the root with the same name as one in 'agroup'
        self.h5file.create_array(
            root, "anarray1", [1, 2], title="Array example"
        )

    def test00_newFile(self):
        """Checking creation of a new file."""

        self.h5file.create_array(
            self.h5file.root, "array_new", [1, 2], title="Array example"
        )

        # Get the CLASS attribute of the arr object
        class_ = self.h5file.root.array.attrs.CLASS

        self.assertEqual(class_.capitalize(), "Array")

    def test00_newFile_unicode_filename(self):
        temp_dir = tempfile.mkdtemp()
        try:
            h5fname = str(Path(temp_dir) / "test.h5")
            with tb.open_file(h5fname, "w") as h5file:
                self.assertTrue(h5file, tb.File)
        finally:
            shutil.rmtree(temp_dir)

    def test00_newFile_numpy_str_filename(self):
        temp_dir = tempfile.mkdtemp()
        try:
            h5fname = np.str_(Path(temp_dir) / "test.h5")
            with tb.open_file(h5fname, "w") as h5file:
                self.assertTrue(h5file, tb.File)
        finally:
            shutil.rmtree(temp_dir)

    def test00_newFile_numpy_unicode_filename(self):
        temp_dir = tempfile.mkdtemp()
        try:
            h5fname = np.str_(Path(temp_dir) / "test.h5")
            with tb.open_file(h5fname, "w") as h5file:
                self.assertTrue(h5file, tb.File)
        finally:
            shutil.rmtree(temp_dir)

    def test01_openFile(self):
        """Checking opening of an existing file."""

        # Open the old HDF5 file
        self._reopen(node_cache_slots=self.node_cache_slots)

        # Get the CLASS attribute of the arr object
        title = self.h5file.root.array.get_attr("TITLE")

        self.assertEqual(title, "Array example")

    def test01_open_file_pathlib(self):
        """Checking opening of an existing file."""
        self.h5file.close()
        h5fname = Path(self.h5fname)
        with tables.open_file(h5fname) as h5file:
            title = h5file.root.array.get_attr("TITLE")
            self.assertEqual(title, "Array example")

    def test02_appendFile(self):
        """Checking appending objects to an existing file."""

        # Append a new array to the existing file
        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        self.h5file.create_array(
            self.h5file.root, "array2", [3, 4], title="Title example 2"
        )

        # Open this file in read-only mode
        self._reopen(node_cache_slots=self.node_cache_slots)

        # Get the CLASS attribute of the arr object
        title = self.h5file.root.array2.get_attr("TITLE")

        self.assertEqual(title, "Title example 2")

    def test02b_appendFile2(self):
        """Checking appending objects to an existing file ("a" version)"""

        # Append a new array to the existing file
        self._reopen(mode="a", node_cache_slots=self.node_cache_slots)
        self.h5file.create_array(
            self.h5file.root, "array2", [3, 4], title="Title example 2"
        )

        # Open this file in read-only mode
        self._reopen(node_cache_slots=self.node_cache_slots)

        # Get the CLASS attribute of the arr object
        title = self.h5file.root.array2.get_attr("TITLE")

        self.assertEqual(title, "Title example 2")

    # Begin to raise errors...

    def test03_appendErrorFile(self):
        """Checking appending objects to an existing file in "w" mode."""

        # Append a new array to the existing file but in write mode
        # so, the existing file should be deleted!
        self._reopen(mode="w", node_cache_slots=self.node_cache_slots)
        self.h5file.create_array(
            self.h5file.root, "array2", [3, 4], title="Title example 2"
        )

        # Open this file in read-only mode
        self._reopen(node_cache_slots=self.node_cache_slots)

        with self.assertRaises(LookupError):
            # Try to get the 'array' object in the old existing file
            self.h5file.root.array

    def test04a_openErrorFile(self):
        """Checking opening a non-existing file for reading"""

        with self.assertRaises(IOError):
            tb.open_file(
                "nonexistent.h5",
                mode="r",
                node_cache_slots=self.node_cache_slots,
            )

    def test04b_alternateRootFile(self):
        """Checking alternate root access to the object tree."""

        # Open the existent HDF5 file
        self._reopen(
            root_uep="/agroup", node_cache_slots=self.node_cache_slots
        )

        # Get the CLASS attribute of the arr object
        if common.verbose:
            print("\nFile tree dump:", self.h5file)
        title = self.h5file.root.anarray1.get_attr("TITLE")

        # Get the node again, as this can trigger errors in some situations
        anarray1 = self.h5file.root.anarray1
        self.assertIsNotNone(anarray1)

        self.assertEqual(title, "Array title 1")

    # This test works well, but HDF5 emits a series of messages that
    # may lose the user. It is better to deactivate it.
    def notest04c_alternateRootFile(self):
        """Checking non-existent alternate root access to the object tree"""

        with self.assertRaises(RuntimeError):
            self._reopen(
                root_uep="/nonexistent", node_cache_slots=self.node_cache_slots
            )

    def test05a_removeGroupRecursively(self):
        """Checking removing a group recursively."""

        # Delete a group with leafs
        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        with self.assertRaises(tb.NodeError):
            self.h5file.remove_node(self.h5file.root.agroup)

        # This should work now
        self.h5file.remove_node(self.h5file.root, "agroup", recursive=1)

        # Open this file in read-only mode
        self._reopen(node_cache_slots=self.node_cache_slots)

        # Try to get the removed object
        with self.assertRaises(LookupError):
            self.h5file.root.agroup

        # Try to get a child of the removed object
        with self.assertRaises(LookupError):
            self.h5file.get_node("/agroup/agroup3")

    def test05b_removeGroupRecursively(self):
        """Checking removing a group recursively and access to it
        immediately."""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test05b_removeGroupRecursively..."
                % self.__class__.__name__
            )

        # Delete a group with leafs
        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        with self.assertRaises(tb.NodeError):
            self.h5file.remove_node(self.h5file.root, "agroup")

        # This should work now
        self.h5file.remove_node(self.h5file.root, "agroup", recursive=1)

        # Try to get the removed object
        with self.assertRaises(LookupError):
            self.h5file.root.agroup

        # Try to get a child of the removed object
        with self.assertRaises(LookupError):
            self.h5file.get_node("/agroup/agroup3")

    def test06_removeNodeWithDel(self):
        """Checking removing a node using ``__delattr__()``"""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        with self.assertRaises(AttributeError):
            # This should fail because there is no *Python attribute*
            # called ``agroup``.
            del self.h5file.root.agroup

    def test06a_removeGroup(self):
        """Checking removing a lonely group from an existing file."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        self.h5file.remove_node(self.h5file.root, "agroup2")

        # Open this file in read-only mode
        self._reopen(node_cache_slots=self.node_cache_slots)

        # Try to get the removed object
        with self.assertRaises(LookupError):
            self.h5file.root.agroup2

    def test06b_removeLeaf(self):
        """Checking removing Leaves from an existing file."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        self.h5file.remove_node(self.h5file.root, "anarray")

        # Open this file in read-only mode
        self._reopen(node_cache_slots=self.node_cache_slots)

        # Try to get the removed object
        with self.assertRaises(LookupError):
            self.h5file.root.anarray

    def test06c_removeLeaf(self):
        """Checking removing Leaves and access it immediately."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        self.h5file.remove_node(self.h5file.root, "anarray")

        # Try to get the removed object
        with self.assertRaises(LookupError):
            self.h5file.root.anarray

    def test06d_removeLeaf(self):
        """Checking removing a non-existent node"""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        # Try to get the removed object
        with self.assertRaises(LookupError):
            self.h5file.remove_node(self.h5file.root, "nonexistent")

    def test06e_removeTable(self):
        """Checking removing Tables from an existing file."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        self.h5file.remove_node(self.h5file.root, "atable")

        # Open this file in read-only mode
        self._reopen(node_cache_slots=self.node_cache_slots)

        # Try to get the removed object
        with self.assertRaises(LookupError):
            self.h5file.root.atable

    def test07_renameLeaf(self):
        """Checking renaming a leave and access it after a close/open."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        self.h5file.rename_node(self.h5file.root.anarray, "anarray2")

        # Open this file in read-only mode
        self._reopen(node_cache_slots=self.node_cache_slots)

        # Ensure that the new name exists
        array_ = self.h5file.root.anarray2
        self.assertEqual(array_.name, "anarray2")
        self.assertEqual(array_._v_pathname, "/anarray2")
        self.assertEqual(array_._v_depth, 1)

        # Try to get the previous object with the old name
        with self.assertRaises(LookupError):
            self.h5file.root.anarray

    def test07b_renameLeaf(self):
        """Checking renaming Leaves and accesing them immediately."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        self.h5file.rename_node(self.h5file.root.anarray, "anarray2")

        # Ensure that the new name exists
        array_ = self.h5file.root.anarray2
        self.assertEqual(array_.name, "anarray2")
        self.assertEqual(array_._v_pathname, "/anarray2")
        self.assertEqual(array_._v_depth, 1)

        # Try to get the previous object with the old name
        with self.assertRaises(LookupError):
            self.h5file.root.anarray

    def test07c_renameLeaf(self):
        """Checking renaming Leaves and modify attributes after that."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        self.h5file.rename_node(self.h5file.root.anarray, "anarray2")
        array_ = self.h5file.root.anarray2
        array_.attrs.TITLE = "hello"

        # Ensure that the new attribute has been written correctly
        self.assertEqual(array_.title, "hello")
        self.assertEqual(array_.attrs.TITLE, "hello")

    def test07d_renameLeaf(self):
        """Checking renaming a Group under a nested group."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        self.h5file.rename_node(self.h5file.root.agroup.anarray2, "anarray3")

        # Ensure that we can access n attributes in the new group
        node = self.h5file.root.agroup.anarray3
        self.assertEqual(node._v_title, "Array title 2")

    def test08_renameToExistingLeaf(self):
        """Checking renaming a node to an existing name."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        # Try to get the previous object with the old name
        with self.assertRaises(tb.NodeError):
            self.h5file.rename_node(self.h5file.root.anarray, "array")

        # Now overwrite the destination node.
        anarray = self.h5file.root.anarray
        self.h5file.rename_node(anarray, "array", overwrite=True)
        self.assertNotIn("/anarray", self.h5file)
        self.assertIs(self.h5file.root.array, anarray)

    def test08b_renameToNotValidNaturalName(self):
        """Checking renaming a node to a non-valid natural name"""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        with warnings.catch_warnings():
            warnings.filterwarnings("error", category=tb.NaturalNameWarning)

            # Try to get the previous object with the old name
            with self.assertRaises(tb.NaturalNameWarning):
                self.h5file.rename_node(self.h5file.root.anarray, "array 2")

    def test09_renameGroup(self):
        """Checking renaming a Group and access it after a close/open."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        self.h5file.rename_node(self.h5file.root.agroup, "agroup3")

        # Open this file in read-only mode
        self._reopen(node_cache_slots=self.node_cache_slots)

        # Ensure that the new name exists
        group = self.h5file.root.agroup3
        self.assertEqual(group._v_name, "agroup3")
        self.assertEqual(group._v_pathname, "/agroup3")

        # The children of this group also must be accessible through the
        # new name path
        group2 = self.h5file.get_node("/agroup3/agroup3")
        self.assertEqual(group2._v_name, "agroup3")
        self.assertEqual(group2._v_pathname, "/agroup3/agroup3")

        # Try to get the previous object with the old name
        with self.assertRaises(LookupError):
            self.h5file.root.agroup

        # Try to get a child with the old pathname
        with self.assertRaises(LookupError):
            self.h5file.get_node("/agroup/agroup3")

    def test09b_renameGroup(self):
        """Checking renaming a Group and access it immediately."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        self.h5file.rename_node(self.h5file.root.agroup, "agroup3")

        # Ensure that the new name exists
        group = self.h5file.root.agroup3
        self.assertEqual(group._v_name, "agroup3")
        self.assertEqual(group._v_pathname, "/agroup3")

        # The children of this group also must be accessible through the
        # new name path
        group2 = self.h5file.get_node("/agroup3/agroup3")
        self.assertEqual(group2._v_name, "agroup3")
        self.assertEqual(group2._v_pathname, "/agroup3/agroup3")

        # Try to get the previous object with the old name
        with self.assertRaises(LookupError):
            self.h5file.root.agroup

        # Try to get a child with the old pathname
        with self.assertRaises(LookupError):
            self.h5file.get_node("/agroup/agroup3")

    def test09c_renameGroup(self):
        """Checking renaming a Group and modify attributes afterwards."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        self.h5file.rename_node(self.h5file.root.agroup, "agroup3")

        # Ensure that we can modify attributes in the new group
        group = self.h5file.root.agroup3
        group._v_attrs.TITLE = "Hello"
        self.assertEqual(group._v_title, "Hello")
        self.assertEqual(group._v_attrs.TITLE, "Hello")

    def test09d_renameGroup(self):
        """Checking renaming a Group under a nested group."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        self.h5file.rename_node(self.h5file.root.agroup.agroup3, "agroup4")

        # Ensure that we can access n attributes in the new group
        group = self.h5file.root.agroup.agroup4
        self.assertEqual(group._v_title, "Group title 3")

    def test09e_renameGroup(self):
        """Checking renaming a Group with nested groups in the LRU cache."""
        # This checks for ticket #126.

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        # Load intermediate groups and keep a nested one alive.
        g = self.h5file.root.agroup.agroup3.agroup4
        self.assertIsNotNone(g)
        self.h5file.rename_node("/", name="agroup", newname="agroup_")

        # see ticket #126
        self.assertNotIn("/agroup_/agroup4", self.h5file)

        self.assertNotIn("/agroup", self.h5file)
        for newpath in [
            "/agroup_",
            "/agroup_/agroup3",
            "/agroup_/agroup3/agroup4",
        ]:
            self.assertIn(newpath, self.h5file)
            self.assertEqual(
                newpath, self.h5file.get_node(newpath)._v_pathname
            )

    def test10_moveLeaf(self):
        """Checking moving a leave and access it after a close/open."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        newgroup = self.h5file.create_group("/", "newgroup")
        self.h5file.move_node(self.h5file.root.anarray, newgroup, "anarray2")

        # Open this file in read-only mode
        self._reopen(node_cache_slots=self.node_cache_slots)

        # Ensure that the new name exists
        array_ = self.h5file.root.newgroup.anarray2
        self.assertEqual(array_.name, "anarray2")
        self.assertEqual(array_._v_pathname, "/newgroup/anarray2")
        self.assertEqual(array_._v_depth, 2)

        # Try to get the previous object with the old name
        with self.assertRaises(LookupError):
            self.h5file.root.anarray

    def test10b_moveLeaf(self):
        """Checking moving a leave and access it without a close/open."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        newgroup = self.h5file.create_group("/", "newgroup")
        self.h5file.move_node(self.h5file.root.anarray, newgroup, "anarray2")

        # Ensure that the new name exists
        array_ = self.h5file.root.newgroup.anarray2
        self.assertEqual(array_.name, "anarray2")
        self.assertEqual(array_._v_pathname, "/newgroup/anarray2")
        self.assertEqual(array_._v_depth, 2)

        # Try to get the previous object with the old name
        with self.assertRaises(LookupError):
            self.h5file.root.anarray

    def test10c_moveLeaf(self):
        """Checking moving Leaves and modify attributes after that."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        newgroup = self.h5file.create_group("/", "newgroup")
        self.h5file.move_node(self.h5file.root.anarray, newgroup, "anarray2")
        array_ = self.h5file.root.newgroup.anarray2
        array_.attrs.TITLE = "hello"

        # Ensure that the new attribute has been written correctly
        self.assertEqual(array_.title, "hello")
        self.assertEqual(array_.attrs.TITLE, "hello")

    def test10d_moveToExistingLeaf(self):
        """Checking moving a leaf to an existing name."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        # Try to get the previous object with the old name
        with self.assertRaises(tb.NodeError):
            self.h5file.move_node(
                self.h5file.root.anarray, self.h5file.root, "array"
            )

    def test10_2_moveTable(self):
        """Checking moving a table and access it after a close/open."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        newgroup = self.h5file.create_group("/", "newgroup")
        self.h5file.move_node(self.h5file.root.atable, newgroup, "atable2")

        # Open this file in read-only mode
        self._reopen(node_cache_slots=self.node_cache_slots)

        # Ensure that the new name exists
        table_ = self.h5file.root.newgroup.atable2
        self.assertEqual(table_.name, "atable2")
        self.assertEqual(table_._v_pathname, "/newgroup/atable2")
        self.assertEqual(table_._v_depth, 2)

        # Try to get the previous object with the old name
        with self.assertRaises(LookupError):
            self.h5file.root.atable

    def test10_2b_moveTable(self):
        """Checking moving a table and access it without a close/open."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        newgroup = self.h5file.create_group("/", "newgroup")
        self.h5file.move_node(self.h5file.root.atable, newgroup, "atable2")

        # Ensure that the new name exists
        table_ = self.h5file.root.newgroup.atable2
        self.assertEqual(table_.name, "atable2")
        self.assertEqual(table_._v_pathname, "/newgroup/atable2")
        self.assertEqual(table_._v_depth, 2)

        # Try to get the previous object with the old name
        with self.assertRaises(LookupError):
            self.h5file.root.atable

    def test10_2b_bis_moveTable(self):
        """Checking moving a table and use cached row without a close/open."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        newgroup = self.h5file.create_group("/", "newgroup")

        # Cache the Row attribute prior to the move
        row = self.h5file.root.atable.row
        self.h5file.move_node(self.h5file.root.atable, newgroup, "atable2")

        # Ensure that the new name exists
        table_ = self.h5file.root.newgroup.atable2
        self.assertEqual(table_.name, "atable2")
        self.assertEqual(table_._v_pathname, "/newgroup/atable2")
        self.assertEqual(table_._v_depth, 2)

        # Ensure that cache Row attribute has been updated
        row = table_.row
        self.assertEqual(table_._v_pathname, row.table._v_pathname)
        nrows = table_.nrows

        # Add a new row just to make sure that this works
        row.append()
        table_.flush()
        self.assertEqual(table_.nrows, nrows + 1)

    def test10_2c_moveTable(self):
        """Checking moving tables and modify attributes after that."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        newgroup = self.h5file.create_group("/", "newgroup")
        self.h5file.move_node(self.h5file.root.atable, newgroup, "atable2")
        table_ = self.h5file.root.newgroup.atable2
        table_.attrs.TITLE = "hello"

        # Ensure that the new attribute has been written correctly
        self.assertEqual(table_.title, "hello")
        self.assertEqual(table_.attrs.TITLE, "hello")

    def test10_2d_moveToExistingTable(self):
        """Checking moving a table to an existing name."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        # Try to get the previous object with the old name
        with self.assertRaises(tb.NodeError):
            self.h5file.move_node(
                self.h5file.root.atable, self.h5file.root, "table"
            )

    def test10_2e_moveToExistingTableOverwrite(self):
        """Checking moving a table to an existing name, overwriting it."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        srcNode = self.h5file.root.atable
        self.h5file.move_node(
            srcNode, self.h5file.root, "table", overwrite=True
        )
        dstNode = self.h5file.root.table

        self.assertIs(srcNode, dstNode)

    def test11_moveGroup(self):
        """Checking moving a Group and access it after a close/open."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        newgroup = self.h5file.create_group(self.h5file.root, "newgroup")
        self.h5file.move_node(self.h5file.root.agroup, newgroup, "agroup3")

        # Open this file in read-only mode
        self._reopen(node_cache_slots=self.node_cache_slots)

        # Ensure that the new name exists
        group = self.h5file.root.newgroup.agroup3
        self.assertEqual(group._v_name, "agroup3")
        self.assertEqual(group._v_pathname, "/newgroup/agroup3")
        self.assertEqual(group._v_depth, 2)

        # The children of this group must also be accessible through the
        # new name path
        group2 = self.h5file.get_node("/newgroup/agroup3/agroup3")
        self.assertEqual(group2._v_name, "agroup3")
        self.assertEqual(group2._v_pathname, "/newgroup/agroup3/agroup3")
        self.assertEqual(group2._v_depth, 3)

        # Try to get the previous object with the old name
        with self.assertRaises(LookupError):
            self.h5file.root.agroup

        # Try to get a child with the old pathname
        with self.assertRaises(LookupError):
            self.h5file.get_node("/agroup/agroup3")

    def test11b_moveGroup(self):
        """Checking moving a Group and access it immediately."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        newgroup = self.h5file.create_group(self.h5file.root, "newgroup")
        self.h5file.move_node(self.h5file.root.agroup, newgroup, "agroup3")

        # Ensure that the new name exists
        group = self.h5file.root.newgroup.agroup3
        self.assertEqual(group._v_name, "agroup3")
        self.assertEqual(group._v_pathname, "/newgroup/agroup3")
        self.assertEqual(group._v_depth, 2)

        # The children of this group must also be accessible through the
        # new name path
        group2 = self.h5file.get_node("/newgroup/agroup3/agroup3")
        self.assertEqual(group2._v_name, "agroup3")
        self.assertEqual(group2._v_pathname, "/newgroup/agroup3/agroup3")
        self.assertEqual(group2._v_depth, 3)

        # Try to get the previous object with the old name
        with self.assertRaises(LookupError):
            self.h5file.root.agroup

        # Try to get a child with the old pathname
        with self.assertRaises(LookupError):
            self.h5file.get_node("/agroup/agroup3")

    def test11c_moveGroup(self):
        """Checking moving a Group and modify attributes afterwards."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        newgroup = self.h5file.create_group(self.h5file.root, "newgroup")
        self.h5file.move_node(self.h5file.root.agroup, newgroup, "agroup3")

        # Ensure that we can modify attributes in the new group
        group = self.h5file.root.newgroup.agroup3
        group._v_attrs.TITLE = "Hello"
        group._v_attrs.hola = "Hello"
        self.assertEqual(group._v_title, "Hello")
        self.assertEqual(group._v_attrs.TITLE, "Hello")
        self.assertEqual(group._v_attrs.hola, "Hello")

    def test11d_moveToExistingGroup(self):
        """Checking moving a group to an existing name."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        # Try to get the previous object with the old name
        with self.assertRaises(tb.NodeError):
            self.h5file.move_node(
                self.h5file.root.agroup, self.h5file.root, "agroup2"
            )

    def test11e_moveToExistingGroupOverwrite(self):
        """Checking moving a group to an existing name, overwriting it."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        # agroup2 -> agroup
        srcNode = self.h5file.root.agroup2
        self.h5file.move_node(
            srcNode, self.h5file.root, "agroup", overwrite=True
        )
        dstNode = self.h5file.root.agroup

        self.assertIs(srcNode, dstNode)

    def test12a_moveNodeOverItself(self):
        """Checking moving a node over itself."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        # array -> array
        srcNode = self.h5file.root.array
        self.h5file.move_node(srcNode, self.h5file.root, "array")
        dstNode = self.h5file.root.array

        self.assertIs(srcNode, dstNode)

    def test12b_moveGroupIntoItself(self):
        """Checking moving a group into itself."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        with self.assertRaises(tb.NodeError):
            # agroup2 -> agroup2/
            self.h5file.move_node(
                self.h5file.root.agroup2, self.h5file.root.agroup2
            )

    def test13a_copyLeaf(self):
        """Copying a leaf."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        # array => agroup2/
        new_node = self.h5file.copy_node(
            self.h5file.root.array, self.h5file.root.agroup2
        )
        dstNode = self.h5file.root.agroup2.array

        self.assertIs(new_node, dstNode)

    def test13b_copyGroup(self):
        """Copying a group."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        # agroup2 => agroup/
        new_node = self.h5file.copy_node(
            self.h5file.root.agroup2, self.h5file.root.agroup
        )
        dstNode = self.h5file.root.agroup.agroup2

        self.assertIs(new_node, dstNode)

    def test13c_copyGroupSelf(self):
        """Copying a group into itself."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        # agroup2 => agroup2/
        new_node = self.h5file.copy_node(
            self.h5file.root.agroup2, self.h5file.root.agroup2
        )
        dstNode = self.h5file.root.agroup2.agroup2

        self.assertIs(new_node, dstNode)

    def test13d_copyGroupRecursive(self):
        """Recursively copying a group."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        # agroup => agroup2/
        new_node = self.h5file.copy_node(
            self.h5file.root.agroup, self.h5file.root.agroup2, recursive=True
        )
        dstNode = self.h5file.root.agroup2.agroup

        self.assertIs(new_node, dstNode)
        dstChild1 = dstNode.anarray1
        self.assertIsNotNone(dstChild1)
        dstChild2 = dstNode.anarray2
        self.assertIsNotNone(dstChild2)
        dstChild3 = dstNode.agroup3
        self.assertIsNotNone(dstChild3)

    def test13e_copyRootRecursive(self):
        """Recursively copying the root group into the root of another file."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        h5fname2 = tempfile.mktemp(".h5")
        h5file2 = tb.open_file(
            h5fname2, mode="w", node_cache_slots=self.node_cache_slots
        )
        try:
            # h5file.root => h5file2.root
            new_node = self.h5file.copy_node(
                self.h5file.root, h5file2.root, recursive=True
            )
            dstNode = h5file2.root

            self.assertIs(new_node, dstNode)
            self.assertIn("/agroup", h5file2)
            self.assertIn("/agroup/anarray1", h5file2)
            self.assertIn("/agroup/agroup3", h5file2)

        finally:
            h5file2.close()
            Path(h5fname2).unlink()

    def test13f_copyRootRecursive(self):
        """Recursively copying the root group into a group in another file."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        h5fname2 = tempfile.mktemp(".h5")
        h5file2 = tb.open_file(
            h5fname2, mode="w", node_cache_slots=self.node_cache_slots
        )
        try:
            h5file2.create_group("/", "agroup2")

            # fileh.root => h5file2.root.agroup2
            new_node = self.h5file.copy_node(
                self.h5file.root, h5file2.root.agroup2, recursive=True
            )
            dstNode = h5file2.root.agroup2

            self.assertIs(new_node, dstNode)
            self.assertIn("/agroup2/agroup", h5file2)
            self.assertIn("/agroup2/agroup/anarray1", h5file2)
            self.assertIn("/agroup2/agroup/agroup3", h5file2)

        finally:
            h5file2.close()
            Path(h5fname2).unlink()

    def test13g_copyRootItself(self):
        """Recursively copying the root group into itself."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        agroup2 = self.h5file.root
        self.assertIsNotNone(agroup2)

        # h5file.root => h5file.root
        self.assertRaises(
            IOError,
            self.h5file.copy_node,
            self.h5file.root,
            self.h5file.root,
            recursive=True,
        )

    def test14a_copyNodeExisting(self):
        """Copying over an existing node."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        with self.assertRaises(tb.NodeError):
            # agroup2 => agroup
            self.h5file.copy_node(self.h5file.root.agroup2, newname="agroup")

    def test14b_copyNodeExistingOverwrite(self):
        """Copying over an existing node, overwriting it."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        # agroup2 => agroup
        new_node = self.h5file.copy_node(
            self.h5file.root.agroup2, newname="agroup", overwrite=True
        )
        dstNode = self.h5file.root.agroup

        self.assertIs(new_node, dstNode)

    def test14b2_copyNodeExistingOverwrite(self):
        """Copying over an existing node in other file, overwriting it."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        h5fname2 = tempfile.mktemp(".h5")
        h5file2 = tb.open_file(
            h5fname2, mode="w", node_cache_slots=self.node_cache_slots
        )

        try:
            # file1:/anarray1 => h5fname2:/anarray1
            new_node = self.h5file.copy_node(
                self.h5file.root.agroup.anarray1, newparent=h5file2.root
            )
            # file1:/ => h5fname2:/
            new_node = self.h5file.copy_node(
                self.h5file.root, h5file2.root, overwrite=True, recursive=True
            )
            dstNode = h5file2.root

            self.assertIs(new_node, dstNode)
        finally:
            h5file2.close()
            Path(h5fname2).unlink()

    def test14c_copyNodeExistingSelf(self):
        """Copying over self."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        with self.assertRaises(tb.NodeError):
            # agroup => agroup
            self.h5file.copy_node(self.h5file.root.agroup, newname="agroup")

    def test14d_copyNodeExistingOverwriteSelf(self):
        """Copying over self, trying to overwrite."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        with self.assertRaises(tb.NodeError):
            # agroup => agroup
            self.h5file.copy_node(
                self.h5file.root.agroup, newname="agroup", overwrite=True
            )

    def test14e_copyGroupSelfRecursive(self):
        """Recursively copying a group into itself."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        with self.assertRaises(tb.NodeError):
            # agroup => agroup/
            self.h5file.copy_node(
                self.h5file.root.agroup,
                self.h5file.root.agroup,
                recursive=True,
            )

    def test15a_oneStepMove(self):
        """Moving and renaming a node in a single action."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        # anarray1 -> agroup/array
        srcNode = self.h5file.root.anarray1
        self.h5file.move_node(srcNode, self.h5file.root.agroup, "array")
        dstNode = self.h5file.root.agroup.array

        self.assertIs(srcNode, dstNode)

    def test15b_oneStepCopy(self):
        """Copying and renaming a node in a single action."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        # anarray1 => agroup/array
        new_node = self.h5file.copy_node(
            self.h5file.root.anarray1, self.h5file.root.agroup, "array"
        )
        dstNode = self.h5file.root.agroup.array

        self.assertIs(new_node, dstNode)

    def test16a_fullCopy(self):
        """Copying full data and user attributes."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        # agroup => groupcopy
        srcNode = self.h5file.root.agroup
        new_node = self.h5file.copy_node(
            srcNode, newname="groupcopy", recursive=True
        )
        dstNode = self.h5file.root.groupcopy

        self.assertIs(new_node, dstNode)
        self.assertEqual(srcNode._v_attrs.testattr, dstNode._v_attrs.testattr)
        self.assertEqual(
            srcNode.anarray1.attrs.testattr, dstNode.anarray1.attrs.testattr
        )
        self.assertEqual(srcNode.anarray1.read(), dstNode.anarray1.read())

    def test16b_partialCopy(self):
        """Copying partial data and no user attributes."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        # agroup => groupcopy
        srcNode = self.h5file.root.agroup
        new_node = self.h5file.copy_node(
            srcNode,
            newname="groupcopy",
            recursive=True,
            copyuserattrs=False,
            start=0,
            stop=5,
            step=2,
        )
        dstNode = self.h5file.root.groupcopy

        self.assertIs(new_node, dstNode)
        self.assertFalse(hasattr(dstNode._v_attrs, "testattr"))
        self.assertFalse(hasattr(dstNode.anarray1.attrs, "testattr"))
        self.assertEqual(
            srcNode.anarray1.read()[0:5:2], dstNode.anarray1.read()
        )

    def test16c_fullCopy(self):
        """Copying full data and user attributes (from file to file)."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)

        h5fname2 = tempfile.mktemp(".h5")
        h5file2 = tb.open_file(
            h5fname2, mode="w", node_cache_slots=self.node_cache_slots
        )

        try:
            # file1:/ => h5fname2:groupcopy
            srcNode = self.h5file.root
            new_node = self.h5file.copy_node(
                srcNode, h5file2.root, newname="groupcopy", recursive=True
            )
            dstNode = h5file2.root.groupcopy

            self.assertIs(new_node, dstNode)
            self.assertEqual(
                srcNode._v_attrs.testattr, dstNode._v_attrs.testattr
            )
            self.assertEqual(
                srcNode.agroup.anarray1.attrs.testattr,
                dstNode.agroup.anarray1.attrs.testattr,
            )
            self.assertEqual(
                srcNode.agroup.anarray1.read(), dstNode.agroup.anarray1.read()
            )
        finally:
            h5file2.close()
            Path(h5fname2).unlink()

    def test17a_CopyChunkshape(self):
        """Copying dataset with a chunkshape."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        srcTable = self.h5file.root.table
        newTable = self.h5file.copy_node(
            srcTable, newname="tablecopy", chunkshape=11
        )

        self.assertEqual(newTable.chunkshape, (11,))
        self.assertNotEqual(srcTable.chunkshape, newTable.chunkshape)

    def test17b_CopyChunkshape(self):
        """Copying dataset with a chunkshape with 'keep' value."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        srcTable = self.h5file.root.table
        newTable = self.h5file.copy_node(
            srcTable, newname="tablecopy", chunkshape="keep"
        )

        self.assertEqual(srcTable.chunkshape, newTable.chunkshape)

    def test17c_CopyChunkshape(self):
        """Copying dataset with a chunkshape with 'auto' value."""

        self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
        srcTable = self.h5file.root.table
        newTable = self.h5file.copy_node(
            srcTable, newname="tablecopy", chunkshape=11
        )
        newTable2 = self.h5file.copy_node(
            newTable, newname="tablecopy2", chunkshape="auto"
        )

        self.assertEqual(srcTable.chunkshape, newTable2.chunkshape)

    def test18_closedRepr(self):
        """Representing a closed node as a string."""

        self._reopen(node_cache_slots=self.node_cache_slots)

        for node in [self.h5file.root.agroup, self.h5file.root.anarray]:
            node._f_close()
            self.assertIn("closed", str(node))
            self.assertIn("closed", repr(node))

    def test19_fileno(self):
        """Checking that the 'fileno()' method works."""

        # Open the old HDF5 file
        self._reopen(mode="r", node_cache_slots=self.node_cache_slots)

        # Get the file descriptor for this file
        fd = self.h5file.fileno()
        if common.verbose:
            print("Value of fileno():", fd)
        self.assertGreaterEqual(fd, 0)


class NodeCacheOpenFile(OpenFileTestCase):
    node_cache_slots = tb.parameters.NODE_CACHE_SLOTS
    open_kwargs = dict(node_cache_slots=node_cache_slots)


class NoNodeCacheOpenFile(OpenFileTestCase):
    node_cache_slots = 0
    open_kwargs = dict(node_cache_slots=node_cache_slots)


class DictNodeCacheOpenFile(OpenFileTestCase):
    node_cache_slots = -tb.parameters.NODE_CACHE_SLOTS
    open_kwargs = dict(node_cache_slots=node_cache_slots)


class CheckFileTestCase(common.TempFileMixin, common.PyTablesTestCase):
    def setUp(self):
        super().setUp()

        # Create a regular (text) file
        self.txtfile = tempfile.mktemp(".h5")
        self.fileh = open(self.txtfile, "w")
        self.fileh.write("Hello!")
        self.fileh.close()

    def tearDown(self):
        self.fileh.close()
        Path(self.txtfile).unlink()
        super().tearDown()

    def test00_isHDF5File(self):
        """Checking  tables.is_hdf5_file function (TRUE case)"""

        # Create a PyTables file (and by so, an HDF5 file)
        self.h5file.create_array(
            self.h5file.root, "array", [1, 2], title="Title example"
        )

        # For this method to run, it needs a closed file
        self.h5file.close()

        # When file has an HDF5 format, always returns 1
        if common.verbose:
            print(
                "\nisHDF5File(%s) ==> %d"
                % (self.h5fname, tb.is_hdf5_file(self.h5fname))
            )
        self.assertEqual(tb.is_hdf5_file(self.h5fname), 1)

    def test01_isHDF5File(self):
        """Checking  tables.is_hdf5_file function (FALSE case)"""

        version = tb.is_hdf5_file(self.txtfile)

        # When file is not an HDF5 format, always returns 0 or
        # negative value
        self.assertLessEqual(version, 0)

    def test01x_isHDF5File_nonexistent(self):
        """Identifying a nonexistent HDF5 file."""
        self.assertRaises(IOError, tb.is_hdf5_file, "nonexistent")

    @common.unittest.skipUnless(
        hasattr(os, "getuid") and os.getuid() != 0, "no UID"
    )
    def test01x_isHDF5File_unreadable(self):
        """Identifying an unreadable HDF5 file."""

        self.h5file.close()
        Path(self.h5fname).chmod(0)  # no permissions at all
        self.assertRaises(IOError, tb.is_hdf5_file, self.h5fname)

    def test02_isPyTablesFile(self):
        """Checking is_pytables_file function (TRUE case)"""

        # Create a PyTables h5fname
        self.h5file.create_array(
            self.h5file.root, "array", [1, 2], title="Title example"
        )

        # For this method to run, it needs a closed h5fname
        self.h5file.close()

        version = tb.is_pytables_file(self.h5fname)

        # When h5fname has a PyTables format, always returns "1.0" string or
        # greater
        if common.verbose:
            print()
            print("\nPyTables format version number ==> %s" % version)
        self.assertGreaterEqual(version, "1.0")

    def test03_isPyTablesFile(self):
        """Checking is_pytables_file function (FALSE case)"""

        version = tb.is_pytables_file(self.txtfile)

        # When file is not a PyTables format, always returns 0 or
        # negative value
        if common.verbose:
            print()
            print("\nPyTables format version number ==> %s" % version)
        self.assertIsNone(version)

    def test04_openGenericHDF5File(self):
        """Checking opening of a generic HDF5 file."""

        # Open an existing generic HDF5 file
        h5fname = common.test_filename("ex-noattr.h5")
        with tb.open_file(h5fname, mode="r") as h5file:
            # Check for some objects inside

            # A group
            columns = h5file.get_node("/columns", classname="Group")
            self.assertEqual(columns._v_name, "columns")

            # An Array
            array_ = h5file.get_node(columns, "TDC", classname="Array")
            self.assertEqual(array_._v_name, "TDC")

            # The new LRU code defers the appearance of a warning to this point

            # Here comes an Array of H5T_ARRAY type
            ui = h5file.get_node(columns, "pressure", classname="Array")
            self.assertEqual(ui._v_name, "pressure")
            if common.verbose:
                print("Array object with type H5T_ARRAY -->", repr(ui))
                print("Array contents -->", ui[:])

            # A Table
            table = h5file.get_node("/detector", "table", classname="Table")
            self.assertEqual(table._v_name, "table")

    def test04b_UnImplementedOnLoading(self):
        """Checking failure loading resulting in an ``UnImplemented`` node."""

        # ############## Note for developers ##############################
        # This test fails if you have the line:                           #
        # ##return ChildClass(self, childname)  # uncomment for debugging #
        # uncommented in Group.py!                                        #
        # #################################################################

        h5fname = common.test_filename("smpl_unsupptype.h5")
        with tb.open_file(h5fname) as h5file:
            with self.assertWarns(UserWarning):
                node = h5file.get_node("/CompoundChunked")
            self.assertIsInstance(node, tb.UnImplemented)

    def test04c_UnImplementedScalar(self):
        """Checking opening of HDF5 files containing scalar dataset of
        UnImplemented type."""

        with tb.open_file(common.test_filename("scalar.h5")) as h5file:
            with self.assertWarns(UserWarning):
                node = h5file.get_node("/variable length string")
            self.assertIsInstance(node, tb.UnImplemented)

    def test05_copyUnimplemented(self):
        """Checking that an UnImplemented object cannot be copied."""

        # Open an existing generic HDF5 file
        h5fname = common.test_filename("smpl_unsupptype.h5")
        with tb.open_file(h5fname, mode="r") as h5file:
            self.assertWarns(UserWarning, h5file.get_node, "/CompoundChunked")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ui = h5file.get_node("/CompoundChunked")
            self.assertEqual(ui._v_name, "CompoundChunked")
            if common.verbose:
                print("UnImplement object -->", repr(ui))

            # Check that it cannot be copied to another file:
            self.assertWarns(UserWarning, ui.copy, self.h5file.root, "newui")

    # The next can be used to check the copy of Array objects with H5T_ARRAY
    # in the future
    def _test05_copyUnimplemented(self):
        """Checking that an UnImplemented object cannot be copied."""

        # Open an existing generic HDF5 file
        # We don't need to wrap this in a try clause because
        # it has already been tried and the warning will not happen again
        h5fname2 = common.test_filename("ex-noattr.h5")
        with tb.open_file(h5fname2, mode="r") as h5file2:
            # An unsupported object (the deprecated H5T_ARRAY type in
            # Array, from pytables 0.8 on)
            ui = h5file2.get_node(h5file2.root.columns, "pressure")
            self.assertEqual(ui._v_name, "pressure")
            if common.verbose:
                print("UnImplement object -->", repr(ui))

            # Check that it cannot be copied to another file
            with warnings.catch_warnings():
                # Force the userwarning to issue an error
                warnings.filterwarnings("error", category=UserWarning)
                with self.assertRaises(UserWarning):
                    ui.copy(self.h5file.root, "newui")


@common.unittest.skipIf(
    tb.file._FILE_OPEN_POLICY == "strict", 'FILE_OPEN_POLICY = "strict"'
)
class ThreadingTestCase(common.TempFileMixin, common.PyTablesTestCase):
    def setUp(self):
        super().setUp()
        self.h5file.create_carray(
            "/", "test_array", tb.Int64Atom(), (200, 300)
        )
        self.h5file.close()

    def test(self):
        lock = threading.Lock()

        def syncronized_open_file(*args, **kwargs):
            with lock:
                return tb.open_file(*args, **kwargs)

        def syncronized_close_file(self, *args, **kwargs):
            with lock:
                return self.close(*args, **kwargs)

        filename = self.h5fname

        def run(filename, q):
            try:
                f = syncronized_open_file(filename, mode="r")
                arr = f.root.test_array[8:12, 18:22]
                assert arr.max() == arr.min() == 0
                syncronized_close_file(f)
            except Exception:
                q.put(sys.exc_info())
            else:
                q.put("OK")

        threads = []
        q = queue.Queue()
        for i in range(10):
            t = threading.Thread(target=run, args=(filename, q))
            t.start()
            threads.append(t)

        for i in range(10):
            self.assertEqual(q.get(), "OK")

        for t in threads:
            t.join()


class PythonAttrsTestCase(common.TempFileMixin, common.PyTablesTestCase):
    """Test interactions of Python attributes and child nodes."""

    def test00_attrOverChild(self):
        """Setting a Python attribute over a child node."""

        root = self.h5file.root

        # Create ``/test`` and overshadow it with ``root.test``.
        child = self.h5file.create_array(root, "test", [1])
        attr = "foobar"
        self.assertWarns(tb.NaturalNameWarning, setattr, root, "test", attr)

        self.assertIs(root.test, attr)
        self.assertIs(root._f_get_child("test"), child)

        # Now bring ``/test`` again to light.
        del root.test

        self.assertIs(root.test, child)

        # Now there is no *attribute* named ``test``.
        self.assertRaises(AttributeError, delattr, root, "test")

    def test01_childUnderAttr(self):
        """Creating a child node under a Python attribute."""

        h5file = self.h5file
        root = h5file.root

        # Create ``root.test`` and an overshadowed ``/test``.
        attr = "foobar"
        root.test = attr
        self.assertWarns(
            tb.NaturalNameWarning, h5file.create_array, root, "test", [1]
        )
        child = h5file.get_node("/test")

        self.assertIs(root.test, attr)
        self.assertIs(root._f_get_child("test"), child)

        # Now bring ``/test`` again to light.
        del root.test

        self.assertIs(root.test, child)

        # Now there is no *attribute* named ``test``.
        self.assertRaises(AttributeError, delattr, root, "test")

    def test02_nodeAttrInLeaf(self):
        """Assigning a ``Node`` value as an attribute to a ``Leaf``."""

        h5file = self.h5file

        array1 = h5file.create_array("/", "array1", [1])
        array2 = h5file.create_array("/", "array2", [1])

        # This may make the garbage collector work a little.
        array1.array2 = array2
        array2.array1 = array1

        # Check the assignments.
        self.assertIs(array1.array2, array2)
        self.assertIs(array2.array1, array1)
        self.assertRaises(
            tb.NoSuchNodeError,  # ``/array1`` is not a group
            h5file.get_node,
            "/array1/array2",
        )
        self.assertRaises(
            tb.NoSuchNodeError,  # ``/array2`` is not a group
            h5file.get_node,
            "/array2/array3",
        )

    def test03_nodeAttrInGroup(self):
        """Assigning a ``Node`` value as an attribute to a ``Group``."""

        h5file = self.h5file
        root = h5file.root

        array = h5file.create_array("/", "array", [1])

        # Assign the array to a pair of attributes,
        # one of them overshadowing the original.
        root.arrayAlias = array
        self.assertWarns(tb.NaturalNameWarning, setattr, root, "array", array)

        # Check the assignments.
        self.assertIs(root.arrayAlias, array)
        self.assertIs(root.array, array)
        self.assertRaises(tb.NoSuchNodeError, h5file.get_node, "/arrayAlias")
        self.assertIs(h5file.get_node("/array"), array)

        # Remove the attribute overshadowing the child.
        del root.array

        # Now there is no *attribute* named ``array``.
        self.assertRaises(AttributeError, delattr, root, "array")


class StateTestCase(common.TempFileMixin, common.PyTablesTestCase):
    """Test that ``File`` and ``Node`` operations check their state (open or
    closed, readable or writable) before proceeding."""

    def test00_fileCopyFileClosed(self):
        """Test copying a closed file."""

        self.h5file.close()
        h5cfname = tempfile.mktemp(suffix=".h5")

        try:
            self.assertRaises(
                tb.ClosedFileError, self.h5file.copy_file, h5cfname
            )
        finally:
            if Path(h5cfname).is_file():
                Path(h5cfname).unlink()

    def test01_fileCloseClosed(self):
        """Test closing an already closed file."""

        self.h5file.close()

        try:
            self.h5file.close()
        except tb.ClosedFileError:
            self.fail("could not close an already closed file")

    def test02_fileFlushClosed(self):
        """Test flushing a closed file."""

        self.h5file.close()
        self.assertRaises(tb.ClosedFileError, self.h5file.flush)

    def test03_fileFlushRO(self):
        """Flushing a read-only file."""

        self._reopen("r")

        try:
            self.h5file.flush()
        except tb.FileModeError:
            self.fail("could not flush a read-only file")

    def test04_fileCreateNodeClosed(self):
        """Test creating a node in a closed file."""

        self.h5file.close()
        self.assertRaises(
            tb.ClosedFileError, self.h5file.create_group, "/", "test"
        )

    def test05_fileCreateNodeRO(self):
        """Test creating a node in a read-only file."""

        self._reopen("r")
        self.assertRaises(
            tb.FileModeError, self.h5file.create_group, "/", "test"
        )

    def test06_fileRemoveNodeClosed(self):
        """Test removing a node from a closed file."""

        self.h5file.create_group("/", "test")
        self.h5file.close()
        self.assertRaises(
            tb.ClosedFileError, self.h5file.remove_node, "/", "test"
        )

    def test07_fileRemoveNodeRO(self):
        """Test removing a node from a read-only file."""

        self.h5file.create_group("/", "test")
        self._reopen("r")
        self.assertRaises(
            tb.FileModeError, self.h5file.remove_node, "/", "test"
        )

    def test08_fileMoveNodeClosed(self):
        """Test moving a node in a closed file."""

        self.h5file.create_group("/", "test1")
        self.h5file.create_group("/", "test2")
        self.h5file.close()
        self.assertRaises(
            tb.ClosedFileError, self.h5file.move_node, "/test1", "/", "test2"
        )

    def test09_fileMoveNodeRO(self):
        """Test moving a node in a read-only file."""

        self.h5file.create_group("/", "test1")
        self.h5file.create_group("/", "test2")
        self._reopen("r")
        self.assertRaises(
            tb.FileModeError, self.h5file.move_node, "/test1", "/", "test2"
        )

    def test10_fileCopyNodeClosed(self):
        """Test copying a node in a closed file."""

        self.h5file.create_group("/", "test1")
        self.h5file.create_group("/", "test2")
        self.h5file.close()
        self.assertRaises(
            tb.ClosedFileError, self.h5file.copy_node, "/test1", "/", "test2"
        )

    def test11_fileCopyNodeRO(self):
        """Test copying a node in a read-only file."""

        self.h5file.create_group("/", "test1")
        self._reopen("r")
        self.assertRaises(
            tb.FileModeError, self.h5file.copy_node, "/test1", "/", "test2"
        )

    def test13_fileGetNodeClosed(self):
        """Test getting a node from a closed file."""

        self.h5file.create_group("/", "test")
        self.h5file.close()
        self.assertRaises(tb.ClosedFileError, self.h5file.get_node, "/test")

    def test14_fileWalkNodesClosed(self):
        """Test walking a closed file."""

        self.h5file.create_group("/", "test1")
        self.h5file.create_group("/", "test2")
        self.h5file.close()
        self.assertRaises(tb.ClosedFileError, next, self.h5file.walk_nodes())

    def test15_fileAttrClosed(self):
        """Test setting and deleting a node attribute in a closed file."""

        self.h5file.create_group("/", "test")
        self.h5file.close()
        self.assertRaises(
            tb.ClosedFileError,
            self.h5file.set_node_attr,
            "/test",
            "foo",
            "bar",
        )
        self.assertRaises(
            tb.ClosedFileError, self.h5file.del_node_attr, "/test", "foo"
        )

    def test16_fileAttrRO(self):
        """Test setting and deleting a node attribute in a read-only file."""

        self.h5file.create_group("/", "test")
        self.h5file.set_node_attr("/test", "foo", "foo")
        self._reopen("r")
        self.assertRaises(
            tb.FileModeError, self.h5file.set_node_attr, "/test", "foo", "bar"
        )
        self.assertRaises(
            tb.FileModeError, self.h5file.del_node_attr, "/test", "foo"
        )

    def test17_fileUndoClosed(self):
        """Test undo operations in a closed file."""

        self.h5file.enable_undo()
        self.h5file.create_group("/", "test2")
        self.h5file.close()
        self.assertRaises(tb.ClosedFileError, self.h5file.is_undo_enabled)
        self.assertRaises(tb.ClosedFileError, self.h5file.get_current_mark)
        self.assertRaises(tb.ClosedFileError, self.h5file.undo)
        self.assertRaises(tb.ClosedFileError, self.h5file.disable_undo)

    def test18_fileUndoRO(self):
        """Test undo operations in a read-only file."""

        self.h5file.enable_undo()
        self.h5file.create_group("/", "test")
        self._reopen("r")
        self.assertEqual(self.h5file._undoEnabled, False)
        # self.assertRaises(FileModeError, self.h5file.undo)
        # self.assertRaises(FileModeError, self.h5file.disable_undo)

    def test19a_getNode(self):
        """Test getting a child of a closed node."""

        g1 = self.h5file.create_group("/", "g1")
        g2 = self.h5file.create_group("/g1", "g2")

        # Close this *object* so that it should not be used.
        g1._f_close()
        self.assertRaises(tb.ClosedNodeError, g1._f_get_child, "g2")

        # Getting a node by its closed object is not allowed.
        self.assertRaises(tb.ClosedNodeError, self.h5file.get_node, g1)

        # Going through that *node* should reopen it automatically.
        try:
            g2_ = self.h5file.get_node("/g1/g2")
        except tb.ClosedNodeError:
            self.fail("closed parent group has not been reopened")

        # Already open nodes should be closed now, but not the new ones.
        self.assertIs(
            g2._v_isopen,
            False,
            "open child of closed group has not been closed",
        )
        self.assertIs(
            g2_._v_isopen,
            True,
            "open child of closed group has not been closed",
        )

        # And existing closed ones should remain closed, but not the new ones.
        g1_ = self.h5file.get_node("/g1")
        self.assertIs(
            g1._v_isopen, False, "already closed group is not closed anymore"
        )
        self.assertIs(
            g1_._v_isopen, True, "newly opened group is still closed"
        )

    def test19b_getNode(self):
        """Test getting a node that does not start with a slash ('/')."""

        # Create an array in the root
        self.h5file.create_array("/", "array", [1, 2], title="Title example")

        # Get the array without specifying a leading slash
        self.assertRaises(NameError, self.h5file.get_node, "array")

    def test20_removeNode(self):
        """Test removing a closed node."""

        # This test is a little redundant once we know that ``File.get_node()``
        # will reload a closed node, but anyway...

        group = self.h5file.create_group("/", "group")
        array = self.h5file.create_array("/group", "array", [1])

        # The closed *object* can not be used.
        group._f_close()
        self.assertRaises(tb.ClosedNodeError, group._f_remove)
        self.assertRaises(tb.ClosedNodeError, self.h5file.remove_node, group)

        # Still, the *node* is reloaded when necessary.
        try:
            self.h5file.remove_node("/group", recursive=True)
        except tb.ClosedNodeError:
            self.fail("closed node has not been reloaded")

        # Objects of descendent removed nodes
        # should have been automatically closed when removed.
        self.assertRaises(tb.ClosedNodeError, array._f_remove)

        self.assertNotIn("/group/array", self.h5file)  # just in case
        self.assertNotIn("/group", self.h5file)  # just in case

    def test21_attrsOfNode(self):
        """Test manipulating the attributes of a closed node."""

        node = self.h5file.create_group("/", "test")
        nodeAttrs = node._v_attrs

        nodeAttrs.test = attr = "foo"

        node._f_close()
        self.assertRaises(tb.ClosedNodeError, getattr, node, "_v_attrs")
        # The design of ``AttributeSet`` does not yet allow this test.
        # self.assertRaises(ClosedNodeError, getattr, nodeAttrs, 'test')

        self.assertEqual(self.h5file.get_node_attr("/test", "test"), attr)

    def test21b_attrsOfNode(self):
        """Test manipulating the attributes of a node in a read-only file."""

        self.h5file.create_group("/", "test")
        self.h5file.set_node_attr("/test", "test", "foo")

        self._reopen("r")
        self.assertRaises(
            tb.FileModeError, self.h5file.set_node_attr, "/test", "test", "bar"
        )

    def test22_fileClosesNode(self):
        """Test node closing because of file closing."""

        node = self.h5file.create_group("/", "test")

        self.h5file.close()
        self.assertRaises(tb.ClosedNodeError, getattr, node, "_v_attrs")

    def test23_reopenFile(self):
        """Testing reopening a file and closing it several times."""

        self.h5file.create_array("/", "test", [1, 2, 3])
        self.h5file.close()

        with tb.open_file(self.h5fname, "r") as h5file1:
            if tb.file._FILE_OPEN_POLICY == "strict":
                self.assertRaises(ValueError, tb.open_file, self.h5fname, "r")
            else:
                with tb.open_file(self.h5fname, "r") as h5file2:
                    if common.verbose:
                        print("(h5file1) test[1]:", h5file1.root.test[1])
                    self.assertEqual(h5file1.root.test[1], 2)
                    h5file1.close()

                    if common.verbose:
                        print("(h5file2) test[1]:", h5file2.root.test[1])
                    self.assertEqual(h5file2.root.test[1], 2)


class FlavorTestCase(common.TempFileMixin, common.PyTablesTestCase):
    """Test that setting, getting and changing the ``flavor`` attribute of a
    leaf works as expected."""

    array_data = np.arange(10)
    scalar_data = np.int32(10)

    def _reopen(self, mode="r"):
        super()._reopen(mode)
        self.array = self.h5file.get_node("/array")
        self.scalar = self.h5file.get_node("/scalar")
        return True

    def setUp(self):
        super().setUp()
        self.array = self.h5file.create_array("/", "array", self.array_data)
        self.scalar = self.h5file.create_array("/", "scalar", self.scalar_data)

    def test00_invalid(self):
        """Setting an invalid flavor."""

        self.assertRaises(tb.FlavorError, setattr, self.array, "flavor", "foo")

    def test01_readonly(self):
        """Setting a flavor in a read-only file."""

        self._reopen(mode="r")
        self.assertRaises(
            tb.FileModeError,
            setattr,
            self.array,
            "flavor",
            tb.flavor.internal_flavor,
        )

    def test02_change(self):
        """Changing the flavor and reading data."""

        for flavor in tb.flavor.all_flavors:
            self.array.flavor = flavor
            self.assertEqual(self.array.flavor, flavor)
            idata = tb.flavor.array_of_flavor(self.array_data, flavor)
            odata = self.array[:]
            self.assertTrue(common.allequal(odata, idata, flavor))

    def test03_store(self):
        """Storing a changed flavor."""

        for flavor in tb.flavor.all_flavors:
            self.array.flavor = flavor
            self.assertEqual(self.array.flavor, flavor)
            self._reopen(mode="r+")
            self.assertEqual(self.array.flavor, flavor)

    def test04_missing(self):
        """Reading a dataset of a missing flavor."""

        flavor = self.array.flavor  # default is internal
        self.array._v_attrs.FLAVOR = "foobar"  # breaks flavor
        self._reopen(mode="r")
        idata = tb.flavor.array_of_flavor(self.array_data, flavor)
        with self.assertWarns(tb.FlavorWarning):
            odata = self.array.read()
        self.assertTrue(common.allequal(odata, idata, flavor))

    def test05_delete(self):
        """Deleting the flavor of a dataset."""

        self.array.flavor = "python"  # non-default
        self.assertEqual(self.array.flavor, "python")
        self.assertEqual(self.array.attrs.FLAVOR, "python")
        del self.array.flavor
        self.assertEqual(self.array.flavor, tb.flavor.internal_flavor)
        self.assertRaises(AttributeError, getattr, self.array.attrs, "FLAVOR")

    def test06_copyDeleted(self):
        """Copying a node with a deleted flavor (see #100)."""

        snames = [node._v_name for node in [self.array, self.scalar]]
        dnames = ["%s_copy" % name for name in snames]
        for name in snames:
            node = self.h5file.get_node("/", name)
            del node.flavor
        # Check the copied flavors right after copying and after reopening.
        for fmode in ["r+", "r"]:
            self._reopen(fmode)
            for sname, dname in zip(snames, dnames):
                if fmode == "r+":
                    snode = self.h5file.get_node("/", sname)
                    node = snode.copy("/", dname)
                elif fmode == "r":
                    node = self.h5file.get_node("/", dname)
                self.assertEqual(
                    node.flavor,
                    tb.flavor.internal_flavor,
                    "flavor of node ``%s`` is not internal: %r"
                    % (node._v_pathname, node.flavor),
                )

    def test07_restrict_flavors(self):
        # regression test for gh-163

        all_flavors = list(tb.flavor.all_flavors)
        alias_map = tb.flavor.alias_map.copy()
        converter_map = tb.flavor.converter_map.copy()
        identifier_map = tb.flavor.identifier_map.copy()
        description_map = tb.flavor.description_map.copy()

        try:
            tb.flavor.restrict_flavors(keep=[])
            self.assertLess(len(tb.flavor.alias_map), len(alias_map))
            self.assertLess(len(tb.flavor.converter_map), len(converter_map))
        finally:
            tb.flavor.all_flavors[:] = all_flavors[:]
            tb.flavor.alias_map.update(alias_map)
            tb.flavor.converter_map.update(converter_map)
            tb.flavor.identifier_map.update(identifier_map)
            tb.flavor.description_map.update(description_map)


# @common.unittest.skipIf(
#     sys.getfilesystemencoding() != "utf-8",
#     "need utf-8 file-system encoding",
# )
@common.unittest.skipIf(
    sys.platform == "win32", "no unicode filenames on windows"
)
class UnicodeFilename(common.TempFileMixin, common.PyTablesTestCase):
    unicode_prefix = "para\u0140lel"

    def _getTempFileName(self):
        return tempfile.mktemp(prefix=self.unicode_prefix, suffix=".h5")

    def setUp(self):
        super().setUp()

        self.test = self.h5file.create_array("/", "test", [1, 2])

        # So as to check the reading
        self._reopen()

    def test01(self):
        """Checking creating a filename with Unicode chars."""

        test = self.h5file.root.test
        if common.verbose:
            print("Filename:", self.h5fname)
            print("Array:", test[:])
            print("Should look like:", [1, 2])
        self.assertEqual(test[:], [1, 2], "Values does not match.")

    def test02(self):
        """Checking  tables.is_hdf5_file with a Unicode filename."""

        self.h5file.close()
        if common.verbose:
            print("Filename:", self.h5fname)
            print(" tables.is_hdf5_file?:", tb.is_hdf5_file(self.h5fname))
        self.assertTrue(tb.is_hdf5_file(self.h5fname))

    def test03(self):
        """Checking is_pytables_file with a Unicode filename."""

        self.h5file.close()
        if common.verbose:
            print("Filename:", self.h5fname)
            print("is_pytables_file?:", tb.is_pytables_file(self.h5fname))
        self.assertNotEqual(tb.is_pytables_file(self.h5fname), False)

    @staticmethod
    def _store_carray(name, data, group):
        atom = tb.Atom.from_dtype(data.dtype)
        node = tb.CArray(group, name, shape=data.shape, atom=atom)
        node[:] = data

    def test_store_and_load_with_non_ascii_attributes(self):
        self.h5file.close()
        self.h5file = tb.open_file(self.h5fname, "a")
        root = self.h5file.root
        group = self.h5file.create_group(root, "face_data")
        array_name = "data at 40\N{DEGREE SIGN}C"
        data = np.sinh(np.linspace(-1.4, 1.4, 500))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", tb.NaturalNameWarning)
            self._store_carray(array_name, data, group)
        group = self.h5file.create_group(root, "vertex_data")


@common.unittest.skipIf(
    sys.version_info < (3, 6), "PEP 519 was implemented in Python 3.6"
)
class PathLikeFilename(common.TempFileMixin, common.PyTablesTestCase):

    def _getTempFileName(self):
        from pathlib import Path

        return Path(tempfile.mktemp(suffix=".h5"))

    def setUp(self):
        super().setUp()

        self.test = self.h5file.create_array("/", "test", [1, 2])

        # So as to check the reading
        self._reopen()

    def test01(self):
        """Checking creating a file with a PathLike object as the filename."""

        test = self.h5file.root.test
        if common.verbose:
            print("Filename:", self.h5fname)
            print("Array:", test[:])
            print("Should look like:", [1, 2])
        self.assertEqual(test[:], [1, 2], "Values does not match.")

    def test02(self):
        """Checking tables.is_hdf5_file with a PathLike object as the
        filename."""

        self.h5file.close()
        if common.verbose:
            print("Filename:", self.h5fname)
            print(" tables.is_hdf5_file?:", tb.is_hdf5_file(self.h5fname))
        self.assertTrue(tb.is_hdf5_file(self.h5fname))

    def test03(self):
        """Checking is_pytables_file with a PathLike object as the filename."""

        self.h5file.close()
        if common.verbose:
            print("Filename:", self.h5fname)
            print("is_pytables_file?:", tb.is_pytables_file(self.h5fname))
        self.assertNotEqual(tb.is_pytables_file(self.h5fname), False)

    def test04_str(self):
        str(self.h5file)


class FilePropertyTestCase(common.PyTablesTestCase):
    def setUp(self):
        super().setUp()
        self.h5fname = tempfile.mktemp(".h5")
        self.h5file = None

    def tearDown(self):
        if self.h5file:
            self.h5file.close()

        if Path(self.h5fname).is_file():
            Path(self.h5fname).unlink()
        super().tearDown()

    def test_get_filesize(self):
        data = np.zeros((2000, 2000))
        datasize = np.prod(data.shape) * data.dtype.itemsize

        self.h5file = tb.open_file(self.h5fname, mode="w")
        self.h5file.create_array(self.h5file.root, "array", data)
        h5_filesize = self.h5file.get_filesize()
        self.h5file.close()

        fs_filesize = Path(self.h5fname).stat().st_size

        self.assertGreaterEqual(h5_filesize, datasize)
        self.assertEqual(h5_filesize, fs_filesize)

    def test01_null_userblock_size(self):
        self.h5file = tb.open_file(self.h5fname, mode="w")
        self.h5file.create_array(self.h5file.root, "array", [1, 2])
        self.assertEqual(self.h5file.get_userblock_size(), 0)

    def test02_null_userblock_size(self):
        self.h5file = tb.open_file(self.h5fname, mode="w")
        self.h5file.create_array(self.h5file.root, "array", [1, 2])
        self.h5file.close()
        self.h5file = tb.open_file(self.h5fname, mode="r")
        self.assertEqual(self.h5file.get_userblock_size(), 0)

    def test03_null_userblock_size(self):
        USER_BLOCK_SIZE = 0
        self.h5file = tb.open_file(
            self.h5fname, mode="w", user_block_size=USER_BLOCK_SIZE
        )
        self.h5file.create_array(self.h5file.root, "array", [1, 2])
        self.assertEqual(self.h5file.get_userblock_size(), 0)

    def test01_userblock_size(self):
        USER_BLOCK_SIZE = 512
        self.h5file = tb.open_file(
            self.h5fname, mode="w", user_block_size=USER_BLOCK_SIZE
        )
        self.h5file.create_array(self.h5file.root, "array", [1, 2])
        self.assertEqual(self.h5file.get_userblock_size(), USER_BLOCK_SIZE)

    def test02_userblock_size(self):
        USER_BLOCK_SIZE = 512
        self.h5file = tb.open_file(
            self.h5fname, mode="w", user_block_size=USER_BLOCK_SIZE
        )
        self.h5file.create_array(self.h5file.root, "array", [1, 2])
        self.h5file.close()
        self.h5file = tb.open_file(self.h5fname, mode="r")
        self.assertEqual(self.h5file.get_userblock_size(), USER_BLOCK_SIZE)

    def test_small_userblock_size(self):
        USER_BLOCK_SIZE = 12
        self.assertRaises(
            ValueError,
            tb.open_file,
            self.h5fname,
            mode="w",
            user_block_size=USER_BLOCK_SIZE,
        )

    def test_invalid_userblock_size(self):
        USER_BLOCK_SIZE = 1025
        self.assertRaises(
            ValueError,
            tb.open_file,
            self.h5fname,
            mode="w",
            user_block_size=USER_BLOCK_SIZE,
        )


# Test for reading a file that uses Blosc and created on a big-endian platform
@common.unittest.skipIf(not common.blosc_avail, "Blosc not available")
class BloscBigEndian(common.TestFileMixin, common.PyTablesTestCase):
    h5fname = common.test_filename("blosc_bigendian.h5")

    def test00_bigendian(self):
        """Checking compatibility with Blosc on big-endian machines."""

        # Check that we can read the contents without problems (nor warnings!)
        for dset_name in ("i1", "i2", "i4", "i8"):
            a = np.arange(10, dtype=dset_name)
            dset = self.h5file.get_node("/" + dset_name)
            self.assertTrue(
                common.allequal(a, dset[:]), "Error in big-endian data!"
            )


# Case test for Blosc and subprocesses (via multiprocessing module)


# The worker function for the subprocess (needs to be here because Windows
# has problems pickling nested functions with the multiprocess module :-/)
def _worker(fn, qout=None):
    fp = tb.open_file(fn)
    if common.verbose:
        print("About to load: ", fn)
    rows = fp.root.table.where("(f0 < 10)")
    if common.verbose:
        print("Got the iterator, about to iterate")
    next(rows)
    if common.verbose:
        print("Succeeded in one iteration\n")
    fp.close()

    if qout is not None:
        qout.put("Done")


# From: Yaroslav Halchenko <debian@onerussian.com>
# Subject: Skip the unittest on kFreeBSD and Hurd -- locking seems to
#         be N/A
#
#  on kfreebsd /dev/shm is N/A
#  on Hurd -- inter-process semaphore locking is N/A
@common.unittest.skipIf(
    not multiprocessing_imported, "multiprocessing module not available"
)
@common.unittest.skipIf(
    platform.system().lower() in ("gnu", "gnu/kfreebsd"),
    "multiprocessing module is not " "supported on Hurd/kFreeBSD",
)
@common.unittest.skipIf(not common.blosc_avail, "Blosc not available")
class BloscSubprocess(common.PyTablesTestCase):
    def test_multiprocess(self):
        # Create a relatively large table with Blosc level 9 (large blocks)
        h5fname = tempfile.mktemp(prefix="multiproc-blosc9-", suffix=".h5")
        try:
            size = 300_000
            sa = np.fromiter(
                ((i, i**2, i // 3) for i in range(size)), "i4,i8,f8"
            )
            with tb.open_file(h5fname, "w") as h5file:
                h5file.create_table(
                    h5file.root,
                    "table",
                    sa,
                    filters=tb.Filters(complevel=9, complib="blosc"),
                    chunkshape=(size // 3,),
                )

            if common.verbose:
                print("**** Running from main process:")
            _worker(h5fname)

            if common.verbose:
                print("**** Running from subprocess:")

            try:
                qout = mp.Queue()
            except OSError:
                print("Permission denied due to /dev/shm settings")
            else:
                ps = mp.Process(
                    target=_worker,
                    args=(
                        h5fname,
                        qout,
                    ),
                )
                ps.daemon = True
                ps.start()

                result = qout.get()
                if common.verbose:
                    print(result)

                ps.join()
                # Avoid warnings with later tests forking subprocesses.
                ps.terminate()
        finally:
            Path(h5fname).unlink()


class HDF5ErrorHandling(common.PyTablesTestCase):
    def setUp(self):
        super().setUp()
        self._old_policy = tb.HDF5ExtError.DEFAULT_H5_BACKTRACE_POLICY

    def tearDown(self):
        tb.HDF5ExtError.DEFAULT_H5_BACKTRACE_POLICY = self._old_policy
        super().tearDown()

    def test_silence_messages(self):
        code = """
import tables as tb
tb.silence_hdf5_messages(False)
tb.silence_hdf5_messages()
try:
    tb.open_file(r'%s')
except tb.HDF5ExtError, e:
    pass
"""

        filename = tempfile.mktemp(prefix="hdf5-error-handling-", suffix=".py")
        try:
            with open(filename, "w") as fp:
                fp.write(code % filename)

            p = subprocess.Popen(
                [sys.executable, filename],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            (stdout, stderr) = p.communicate()

            self.assertNotIn("HDF5-DIAG", stderr.decode("ascii"))
        finally:
            Path(filename).unlink()

    # This test is a bit flaky and in some situations it fails
    # E.g. on Mac OSX (arm64), I am getting this:
    # FAIL: None (tables.tests.test_basics.HDF5ErrorHandling)
    # ----------------------------------------------------------------------
    # Traceback (most recent call last):
    #   File "/Users/faltet/software/PyTables-upstream/tables/tests/test_basics.py",
    #   line 2231, in test_enable_messages
    #     self.assertIn("HDF5-DIAG", stderr.decode('ascii'))
    # AssertionError: 'HDF5-DIAG' not found in 'Traceback (most recent call last):\n
    # symbol not found in flat namespace \'_blosc2_cbuffer_sizes\'\n'
    # As the fix is not clear to me, I prefer to disable it until a more robust
    # path is found.
    def _test_enable_messages(self):
        code = """
import tables as tb
tb.silence_hdf5_messages()
tb.silence_hdf5_messages(False)
try:
    tb.open_file(r'%s')
except tb.HDF5ExtError as e:
    pass
"""

        filename = tempfile.mktemp(prefix="hdf5-error-handling-", suffix=".py")
        try:
            Path(filename).write_text(code % filename)

            p = subprocess.Popen(
                [sys.executable, filename],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            (stdout, stderr) = p.communicate()

            self.assertIn("HDF5-DIAG", stderr.decode("ascii"))
        finally:
            Path(filename).unlink()

    def _raise_exterror(self):
        h5fname = tempfile.mktemp(".h5")
        Path(h5fname).write_text("")

        try:
            h5file = tb.open_file(h5fname)
            h5file.close()
        finally:
            Path(h5fname).unlink()

    def test_h5_backtrace_quiet(self):
        tb.HDF5ExtError.DEFAULT_H5_BACKTRACE_POLICY = True

        with self.assertRaises(tb.HDF5ExtError) as cm:
            self._raise_exterror()

        self.assertIsNotNone(cm.exception.h5backtrace)

    def test_h5_backtrace_verbose(self):
        tb.HDF5ExtError.DEFAULT_H5_BACKTRACE_POLICY = "VERBOSE"

        with self.assertRaises(tb.HDF5ExtError) as cm:
            self._raise_exterror()

        self.assertIsNotNone(cm.exception.h5backtrace)
        msg = str(cm.exception)
        self.assertIn(cm.exception.h5backtrace[-1][-1], msg)

    def test_h5_backtrace_ignore(self):
        tb.HDF5ExtError.DEFAULT_H5_BACKTRACE_POLICY = False

        with self.assertRaises(tb.HDF5ExtError) as cm:
            self._raise_exterror()

        self.assertIsNone(cm.exception.h5backtrace)


class TestDescription(common.PyTablesTestCase):
    def test_isdescription_inheritance(self):
        # Regression test for gh-65
        class TestDescParent(tb.IsDescription):
            c = tb.Int32Col()

        class TestDesc(TestDescParent):
            pass

        self.assertIn("c", TestDesc.columns)

    def test_descr_from_dtype(self):
        t = np.dtype([("col1", "int16"), ("col2", float)])
        descr, byteorder = tb.description.descr_from_dtype(t)

        self.assertIn("col1", descr._v_colobjects)
        self.assertIn("col2", descr._v_colobjects)
        self.assertEqual(len(descr._v_colobjects), 2)
        self.assertIsInstance(descr._v_colobjects["col1"], tb.Col)
        self.assertIsInstance(descr._v_colobjects["col2"], tb.Col)
        self.assertEqual(descr._v_colobjects["col1"].dtype, np.int16)
        self.assertEqual(descr._v_colobjects["col2"].dtype, float)

    def test_descr_from_dtype_rich_dtype(self):
        header = [
            (("timestamp", "t"), "u4"),
            (("unit (cluster) id", "unit"), "u2"),
        ]
        t = np.dtype(header)

        descr, byteorder = tb.description.descr_from_dtype(t)
        self.assertEqual(len(descr._v_names), 2)
        self.assertEqual(sorted(descr._v_names), ["t", "unit"])

    def test_descr_from_dtype_comp_01(self):
        d1 = np.dtype([("x", "int16"), ("y", "int16")])
        d_comp = np.dtype([("time", "float64"), ("value", d1)])

        descr, byteorder = tb.description.descr_from_dtype(d_comp)

        self.assertTrue(descr._v_is_nested)
        self.assertIn("time", descr._v_colobjects)
        self.assertIn("value", descr._v_colobjects)
        self.assertEqual(len(descr._v_colobjects), 2)
        self.assertIsInstance(descr._v_colobjects["time"], tb.Col)
        self.assertTrue(
            isinstance(descr._v_colobjects["value"], tb.Description)
        )
        self.assertEqual(descr._v_colobjects["time"].dtype, np.float64)

    def test_descr_from_dtype_comp_02(self):
        d1 = np.dtype([("x", "int16"), ("y", "int16")])

        d_comp = np.dtype([("time", "float64"), ("value", (d1, (1,)))])

        with self.assertWarns(UserWarning):
            descr, byteorder = tb.description.descr_from_dtype(d_comp)

        self.assertTrue(descr._v_is_nested)
        self.assertIn("time", descr._v_colobjects)
        self.assertIn("value", descr._v_colobjects)
        self.assertEqual(len(descr._v_colobjects), 2)
        self.assertIsInstance(descr._v_colobjects["time"], tb.Col)
        self.assertTrue(
            isinstance(descr._v_colobjects["value"], tb.Description)
        )
        self.assertEqual(descr._v_colobjects["time"].dtype, np.float64)

    def test_dtype_from_descr_is_description(self):
        # See gh-152
        class TestDescParent(tb.IsDescription):
            col1 = tb.Int16Col()
            col2 = tb.FloatCol()

        dtype = np.dtype([("col1", "int16"), ("col2", float)])
        t = tb.description.dtype_from_descr(TestDescParent)

        self.assertEqual(t, dtype)

    def test_dtype_from_descr_is_description_instance(self):
        # See gh-152
        class TestDescParent(tb.IsDescription):
            col1 = tb.Int16Col()
            col2 = tb.FloatCol()

        dtype = np.dtype([("col1", "int16"), ("col2", float)])
        t = tb.description.dtype_from_descr(TestDescParent())

        self.assertEqual(t, dtype)

    def test_dtype_from_descr_description_instance(self):
        # See gh-152
        class TestDescParent(tb.IsDescription):
            col1 = tb.Int16Col()
            col2 = tb.FloatCol()

        dtype = np.dtype([("col1", "int16"), ("col2", float)])
        desctiption = tb.Description(TestDescParent().columns)
        t = tb.description.dtype_from_descr(desctiption)

        self.assertEqual(t, dtype)

    def test_dtype_from_descr_dict(self):
        # See gh-152
        dtype = np.dtype([("col1", "int16"), ("col2", float)])
        t = tb.description.dtype_from_descr(
            {"col1": tb.Int16Col(), "col2": tb.FloatCol()}
        )

        self.assertEqual(t, dtype)

    def test_dtype_from_descr_invalid_type(self):
        # See gh-152
        self.assertRaises(ValueError, tb.description.dtype_from_descr, [])

    def test_dtype_from_descr_byteorder(self):
        # See gh-152
        class TestDescParent(tb.IsDescription):
            col1 = tb.Int16Col()
            col2 = tb.FloatCol()

        t = tb.description.dtype_from_descr(TestDescParent, byteorder=">")

        self.assertEqual(t["col1"].byteorder, ">")
        self.assertEqual(t["col2"].byteorder, ">")

    def test_str_names(self):
        # see gh-42
        d = {"name": tb.Int16Col()}
        descr = tb.Description(d)
        self.assertEqual(sorted(descr._v_names), sorted(d.keys()))
        self.assertIsInstance(descr._v_dtype, np.dtype)
        self.assertTrue(sorted(descr._v_dtype.fields), sorted(d.keys()))


class TestAtom(common.PyTablesTestCase):
    def test_atom_attributes01(self):
        shape = (10, 10)
        a = tb.Float64Atom(shape=shape)

        self.assertEqual(a.dflt, 0.0)
        self.assertEqual(a.dtype, np.dtype((np.float64, shape)))
        self.assertEqual(a.itemsize, a.dtype.base.itemsize)
        self.assertEqual(a.kind, "float")
        self.assertEqual(a.ndim, len(shape))
        # self.assertEqual(a.recarrtype, )
        self.assertEqual(a.shape, shape)
        self.assertEqual(a.size, a.itemsize * np.prod(shape))
        self.assertEqual(a.type, "float64")

    def test_atom_copy01(self):
        shape = (10, 10)
        a = tb.Float64Atom(shape=shape)
        aa = a.copy()
        self.assertEqual(aa.shape, shape)

    def test_atom_copy02(self):
        dflt = 2.0
        a = tb.Float64Atom(dflt=dflt)
        aa = a.copy()
        self.assertEqual(aa.dflt, dflt)

    def test_atom_copy_override(self):
        shape = (10, 10)
        dflt = 2.0
        a = tb.Float64Atom(shape=shape, dflt=dflt)
        aa = a.copy(dflt=-dflt)
        self.assertEqual(aa.shape, shape)
        self.assertNotEqual(aa.dflt, dflt)
        self.assertEqual(aa.dflt, -dflt)


class TestCol(common.PyTablesTestCase):
    def test_col_copy01(self):
        shape = (10, 10)
        c = tb.Float64Col(shape=shape)
        cc = c.copy()
        self.assertEqual(cc.shape, shape)

    def test_col_copy02(self):
        dflt = 2.0
        c = tb.Float64Col(dflt=dflt)
        cc = c.copy()
        self.assertEqual(cc.dflt, dflt)

    def test_col_copy_override(self):
        shape = (10, 10)
        dflt = 2.0
        pos = 3
        c = tb.Float64Col(shape=shape, dflt=dflt, pos=pos)
        cc = c.copy(pos=2)
        self.assertEqual(cc.shape, shape)
        self.assertEqual(cc.dflt, dflt)
        self.assertNotEqual(cc._v_pos, pos)
        self.assertEqual(cc._v_pos, 2)


class TestSysattrCompatibility(common.PyTablesTestCase):
    def test_open_python2(self):
        h5fname = common.test_filename("python2.h5")
        with tb.open_file(h5fname, "r") as h5file:
            self.assertTrue(h5file.isopen)

    def test_open_python3(self):
        h5fname = common.test_filename("python3.h5")
        with tb.open_file(h5fname, "r") as h5file:
            self.assertTrue(h5file.isopen)


def suite():
    theSuite = common.unittest.TestSuite()
    niter = 1

    for i in range(niter):
        theSuite.addTest(common.make_suite(OpenFileFailureTestCase))
        theSuite.addTest(common.make_suite(NodeCacheOpenFile))
        theSuite.addTest(common.make_suite(NoNodeCacheOpenFile))
        theSuite.addTest(common.make_suite(DictNodeCacheOpenFile))
        theSuite.addTest(common.make_suite(CheckFileTestCase))
        theSuite.addTest(common.make_suite(ThreadingTestCase))
        theSuite.addTest(common.make_suite(PythonAttrsTestCase))
        theSuite.addTest(common.make_suite(StateTestCase))
        theSuite.addTest(common.make_suite(FlavorTestCase))
        theSuite.addTest(common.make_suite(UnicodeFilename))
        theSuite.addTest(common.make_suite(PathLikeFilename))
        theSuite.addTest(common.make_suite(FilePropertyTestCase))
        theSuite.addTest(common.make_suite(BloscBigEndian))
        theSuite.addTest(common.make_suite(BloscSubprocess))
        theSuite.addTest(common.make_suite(HDF5ErrorHandling))
        theSuite.addTest(common.make_suite(TestDescription))
        theSuite.addTest(common.make_suite(TestAtom))
        theSuite.addTest(common.make_suite(TestCol))
        theSuite.addTest(common.make_suite(TestSysattrCompatibility))

    return theSuite


if __name__ == "__main__":
    common.parse_argv(sys.argv)
    common.print_versions()
    common.unittest.main(defaultTest="suite")
