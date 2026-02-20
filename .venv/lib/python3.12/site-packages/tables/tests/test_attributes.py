"""This test unit checks node attributes that are persistent (AttributeSet)."""

import sys
import datetime
import warnings

import numpy as np
from packaging.version import Version

import tables as tb
from tables.tests import common


class Record(tb.IsDescription):
    var1 = tb.StringCol(itemsize=4)  # 4-character String
    var2 = tb.IntCol()  # integer
    var3 = tb.Int16Col()  # short integer
    var4 = tb.FloatCol()  # double (double-precision)
    var5 = tb.Float32Col()  # float  (single-precision)


class CreateTestCase(common.TempFileMixin, common.PyTablesTestCase):
    def setUp(self):
        super().setUp()
        self.root = self.h5file.root

        # Create a table object
        self.table = self.h5file.create_table(
            self.root, "atable", Record, "Table title"
        )
        # Create an array object
        self.array = self.h5file.create_array(
            self.root, "anarray", [1], "Array title"
        )
        # Create a group object
        self.group = self.h5file.create_group(
            self.root, "agroup", "Group title"
        )

    def test01_setAttributes(self):
        """Checking setting large string attributes (File methods)"""

        attrlength = 2048
        # Try to put a long string attribute on a group object
        self.h5file.set_node_attr(self.root.agroup, "attr1", "p" * attrlength)

        # Now, try with a Table object
        self.h5file.set_node_attr(self.root.atable, "attr1", "a" * attrlength)

        # Finally, try with an Array object
        self.h5file.set_node_attr(self.root.anarray, "attr1", "n" * attrlength)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
            self.root = self.h5file.root

        self.assertEqual(
            self.h5file.get_node_attr(self.root.agroup, "attr1"),
            "p" * attrlength,
        )
        self.assertEqual(
            self.h5file.get_node_attr(self.root.atable, "attr1"),
            "a" * attrlength,
        )
        self.assertEqual(
            self.h5file.get_node_attr(self.root.anarray, "attr1"),
            "n" * attrlength,
        )

    def reopen(self):
        # Reopen
        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
            self.root = self.h5file.root

    def check_missing(self, name):
        self.reopen()
        self.assertNotIn(name, self.root.agroup._v_attrs)
        self.assertNotIn(name, self.root.atable.attrs)
        self.assertNotIn(name, self.root.anarray.attrs)

    def check_name(self, name, val=""):
        """Check validity of attribute name filtering"""
        self.check_missing(name)
        # Using File methods
        self.h5file.set_node_attr(self.root.agroup, name, val)
        self.h5file.set_node_attr(self.root.atable, name, val)
        self.h5file.set_node_attr(self.root.anarray, name, val)
        # Check File methods
        self.reopen()
        self.assertEqual(
            self.h5file.get_node_attr(self.root.agroup, name), val
        )
        self.assertEqual(
            self.h5file.get_node_attr(self.root.atable, name), val
        )
        self.assertEqual(
            self.h5file.get_node_attr(self.root.anarray, name), val
        )
        # Remove, file methods
        self.h5file.del_node_attr(self.root.agroup, name)
        self.h5file.del_node_attr(self.root.atable, name)
        self.h5file.del_node_attr(self.root.anarray, name)
        self.check_missing(name)

        # Using Node methods
        self.root.agroup._f_setattr(name, val)
        self.root.atable.set_attr(name, val)
        self.root.anarray.set_attr(name, val)
        # Check Node methods
        self.reopen()
        self.assertEqual(self.root.agroup._f_getattr(name), val)
        self.assertEqual(self.root.atable.get_attr(name), val)
        self.assertEqual(self.root.anarray.get_attr(name), val)
        self.root.agroup._f_delattr(name)
        self.root.atable.del_attr(name)
        self.root.anarray.del_attr(name)
        self.check_missing(name)

        # Using AttributeSet methods
        setattr(self.root.agroup._v_attrs, name, val)
        setattr(self.root.atable.attrs, name, val)
        setattr(self.root.anarray.attrs, name, val)
        # Check AttributeSet methods
        self.reopen()
        self.assertEqual(getattr(self.root.agroup._v_attrs, name), val)
        self.assertEqual(getattr(self.root.atable.attrs, name), val)
        self.assertEqual(getattr(self.root.anarray.attrs, name), val)
        delattr(self.root.agroup._v_attrs, name)
        delattr(self.root.atable.attrs, name)
        delattr(self.root.anarray.attrs, name)
        self.check_missing(name)

        # Using dict []
        self.root.agroup._v_attrs[name] = val
        self.root.atable.attrs[name] = val
        self.root.anarray.attrs[name] = val
        # Check dict []
        self.reopen()
        self.assertEqual(self.root.agroup._v_attrs[name], val)
        self.assertEqual(self.root.atable.attrs[name], val)
        self.assertEqual(self.root.anarray.attrs[name], val)
        del self.root.agroup._v_attrs[name]
        del self.root.atable.attrs[name]
        del self.root.anarray.attrs[name]
        self.check_missing(name)

    def test01a_setAttributes(self):
        """Checking attribute names validity"""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", tb.NaturalNameWarning)
            self.check_name("a")
            self.check_name("a:b")
            self.check_name("/a/b")
            self.check_name(".")
        self.assertRaises(ValueError, self.check_name, "")
        self.assertRaises(ValueError, self.check_name, "__members__")
        self.assertRaises(TypeError, self.check_name, 0)

    def test02_setAttributes(self):
        """Checking setting large string attributes (Node methods)"""

        attrlength = 2048
        # Try to put a long string attribute on a group object
        self.root.agroup._f_setattr("attr1", "p" * attrlength)
        # Now, try with a Table object
        self.root.atable.set_attr("attr1", "a" * attrlength)

        # Finally, try with an Array object
        self.root.anarray.set_attr("attr1", "n" * attrlength)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
            self.root = self.h5file.root

        self.assertEqual(
            self.root.agroup._f_getattr("attr1"), "p" * attrlength
        )
        self.assertEqual(self.root.atable.get_attr("attr1"), "a" * attrlength)
        self.assertEqual(self.root.anarray.get_attr("attr1"), "n" * attrlength)

    def test03_setAttributes(self):
        """Checking setting large string attributes (AttributeSet methods)"""

        attrlength = 2048
        # Try to put a long string attribute on a group object
        self.group._v_attrs.attr1 = "p" * attrlength
        # Now, try with a Table object
        self.table.attrs.attr1 = "a" * attrlength
        # Finally, try with an Array object
        self.array.attrs.attr1 = "n" * attrlength

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
            self.root = self.h5file.root

        # This should work even when the node cache is disabled
        self.assertEqual(self.root.agroup._v_attrs.attr1, "p" * attrlength)
        self.assertEqual(self.root.atable.attrs.attr1, "a" * attrlength)
        self.assertEqual(self.root.anarray.attrs.attr1, "n" * attrlength)

    def test04_listAttributes(self):
        """Checking listing attributes."""

        # With a Group object
        self.group._v_attrs.pq = "1"
        self.group._v_attrs.qr = "2"
        self.group._v_attrs.rs = "3"
        if common.verbose:
            print("Attribute list:", self.group._v_attrs._f_list())

        # Now, try with a Table object
        self.table.attrs.a = "1"
        self.table.attrs.c = "2"
        self.table.attrs.b = "3"
        if common.verbose:
            print("Attribute list:", self.table.attrs._f_list())

        # Finally, try with an Array object
        self.array.attrs.k = "1"
        self.array.attrs.j = "2"
        self.array.attrs.i = "3"
        if common.verbose:
            print("Attribute list:", self.array.attrs._f_list())

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
            self.root = self.h5file.root

        agroup = self.root.agroup
        self.assertEqual(agroup._v_attrs._f_list("user"), ["pq", "qr", "rs"])
        self.assertEqual(
            agroup._v_attrs._f_list("sys"), ["CLASS", "TITLE", "VERSION"]
        )
        self.assertEqual(
            agroup._v_attrs._f_list("all"),
            ["CLASS", "TITLE", "VERSION", "pq", "qr", "rs"],
        )

        atable = self.root.atable
        self.assertEqual(atable.attrs._f_list(), ["a", "b", "c"])
        self.assertEqual(
            atable.attrs._f_list("sys"),
            [
                "CLASS",
                "FIELD_0_FILL",
                "FIELD_0_NAME",
                "FIELD_1_FILL",
                "FIELD_1_NAME",
                "FIELD_2_FILL",
                "FIELD_2_NAME",
                "FIELD_3_FILL",
                "FIELD_3_NAME",
                "FIELD_4_FILL",
                "FIELD_4_NAME",
                "NROWS",
                "TITLE",
                "VERSION",
            ],
        )
        self.assertEqual(
            atable.attrs._f_list("all"),
            [
                "CLASS",
                "FIELD_0_FILL",
                "FIELD_0_NAME",
                "FIELD_1_FILL",
                "FIELD_1_NAME",
                "FIELD_2_FILL",
                "FIELD_2_NAME",
                "FIELD_3_FILL",
                "FIELD_3_NAME",
                "FIELD_4_FILL",
                "FIELD_4_NAME",
                "NROWS",
                "TITLE",
                "VERSION",
                "a",
                "b",
                "c",
            ],
        )

        anarray = self.root.anarray
        self.assertEqual(anarray.attrs._f_list(), ["i", "j", "k"])
        self.assertEqual(
            anarray.attrs._f_list("sys"),
            ["CLASS", "FLAVOR", "TITLE", "VERSION"],
        )
        self.assertEqual(
            anarray.attrs._f_list("all"),
            ["CLASS", "FLAVOR", "TITLE", "VERSION", "i", "j", "k"],
        )

    def test05_removeAttributes(self):
        """Checking removing attributes."""

        # With a Group object
        self.group._v_attrs.pq = "1"
        self.group._v_attrs.qr = "2"
        self.group._v_attrs.rs = "3"
        # delete an attribute
        del self.group._v_attrs.pq

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
            self.root = self.h5file.root

        agroup = self.root.agroup
        if common.verbose:
            print("Attribute list:", agroup._v_attrs._f_list())
        # Check the local attributes names
        self.assertEqual(agroup._v_attrs._f_list(), ["qr", "rs"])
        if common.verbose:
            print("Attribute list in disk:", agroup._v_attrs._f_list("all"))
        # Check the disk attribute names
        self.assertEqual(
            agroup._v_attrs._f_list("all"),
            ["CLASS", "TITLE", "VERSION", "qr", "rs"],
        )

        # delete an attribute (__delattr__ method)
        del agroup._v_attrs.qr
        if common.verbose:
            print("Attribute list:", agroup._v_attrs._f_list())
        # Check the local attributes names
        self.assertEqual(agroup._v_attrs._f_list(), ["rs"])
        if common.verbose:
            print("Attribute list in disk:", agroup._v_attrs._f_list())
        # Check the disk attribute names
        self.assertEqual(
            agroup._v_attrs._f_list("all"), ["CLASS", "TITLE", "VERSION", "rs"]
        )

    def test05b_removeAttributes(self):
        """Checking removing attributes (using File.del_node_attr())"""

        # With a Group object
        self.group._v_attrs.pq = "1"
        self.group._v_attrs.qr = "2"
        self.group._v_attrs.rs = "3"
        # delete an attribute
        self.h5file.del_node_attr(self.group, "pq")

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
            self.root = self.h5file.root

        agroup = self.root.agroup
        if common.verbose:
            print("Attribute list:", agroup._v_attrs._f_list())
        # Check the local attributes names
        self.assertEqual(agroup._v_attrs._f_list(), ["qr", "rs"])
        if common.verbose:
            print("Attribute list in disk:", agroup._v_attrs._f_list("all"))
        # Check the disk attribute names
        self.assertEqual(
            agroup._v_attrs._f_list("all"),
            ["CLASS", "TITLE", "VERSION", "qr", "rs"],
        )

        # delete an attribute (File.del_node_attr method)
        self.h5file.del_node_attr(self.root, "qr", "agroup")
        if common.verbose:
            print("Attribute list:", agroup._v_attrs._f_list())
        # Check the local attributes names
        self.assertEqual(agroup._v_attrs._f_list(), ["rs"])
        if common.verbose:
            print("Attribute list in disk:", agroup._v_attrs._f_list())
        # Check the disk attribute names
        self.assertEqual(
            agroup._v_attrs._f_list("all"), ["CLASS", "TITLE", "VERSION", "rs"]
        )

    def test06_removeAttributes(self):
        """Checking removing system attributes."""

        # remove a system attribute
        if common.verbose:
            print("Before removing CLASS attribute")
            print("System attrs:", self.group._v_attrs._v_attrnamessys)
        del self.group._v_attrs.CLASS
        self.assertEqual(
            self.group._v_attrs._f_list("sys"), ["TITLE", "VERSION"]
        )
        if common.verbose:
            print("After removing CLASS attribute")
            print("System attrs:", self.group._v_attrs._v_attrnamessys)

    def test07_renameAttributes(self):
        """Checking renaming attributes."""

        # With a Group object
        self.group._v_attrs.pq = "1"
        self.group._v_attrs.qr = "2"
        self.group._v_attrs.rs = "3"
        # rename an attribute
        self.group._v_attrs._f_rename("pq", "op")

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
            self.root = self.h5file.root

        agroup = self.root.agroup
        if common.verbose:
            print("Attribute list:", agroup._v_attrs._f_list())
        # Check the local attributes names (alphabetically sorted)
        self.assertEqual(agroup._v_attrs._f_list(), ["op", "qr", "rs"])
        if common.verbose:
            print("Attribute list in disk:", agroup._v_attrs._f_list("all"))
        # Check the disk attribute names (not sorted)
        self.assertEqual(
            agroup._v_attrs._f_list("all"),
            ["CLASS", "TITLE", "VERSION", "op", "qr", "rs"],
        )

    def test08_renameAttributes(self):
        """Checking renaming system attributes."""

        if common.verbose:
            print("Before renaming CLASS attribute")
            print("All attrs:", self.group._v_attrs._v_attrnames)
        # rename a system attribute
        self.group._v_attrs._f_rename("CLASS", "op")
        if common.verbose:
            print("After renaming CLASS attribute")
            print("All attrs:", self.group._v_attrs._v_attrnames)

        # Check the disk attribute names (not sorted)
        agroup = self.root.agroup
        self.assertEqual(
            agroup._v_attrs._f_list("all"), ["TITLE", "VERSION", "op"]
        )

    def test09_overwriteAttributes(self):
        """Checking overwriting attributes."""

        # With a Group object
        self.group._v_attrs.pq = "1"
        self.group._v_attrs.qr = "2"
        self.group._v_attrs.rs = "3"
        # overwrite attributes
        self.group._v_attrs.pq = "4"
        self.group._v_attrs.qr = 2
        self.group._v_attrs.rs = [1, 2, 3]

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
            self.root = self.h5file.root

        agroup = self.root.agroup
        if common.verbose:
            print("Value of Attribute pq:", agroup._v_attrs.pq)
        # Check the local attributes names (alphabetically sorted)
        self.assertEqual(agroup._v_attrs.pq, "4")
        self.assertEqual(agroup._v_attrs.qr, 2)
        self.assertEqual(agroup._v_attrs.rs, [1, 2, 3])
        if common.verbose:
            print("Attribute list in disk:", agroup._v_attrs._f_list("all"))
        # Check the disk attribute names (not sorted)
        self.assertEqual(
            agroup._v_attrs._f_list("all"),
            ["CLASS", "TITLE", "VERSION", "pq", "qr", "rs"],
        )

    def test10a_copyAttributes(self):
        """Checking copying attributes."""

        # With a Group object
        self.group._v_attrs.pq = "1"
        self.group._v_attrs.qr = "2"
        self.group._v_attrs.rs = "3"
        # copy all attributes from "/agroup" to "/atable"
        self.group._v_attrs._f_copy(self.root.atable)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
            self.root = self.h5file.root

        atable = self.root.atable
        if common.verbose:
            print("Attribute list:", atable._v_attrs._f_list())
        # Check the local attributes names (alphabetically sorted)
        self.assertEqual(atable._v_attrs._f_list(), ["pq", "qr", "rs"])
        if common.verbose:
            print("Complete attribute list:", atable._v_attrs._f_list("all"))
        # Check the disk attribute names (not sorted)
        self.assertEqual(
            atable._v_attrs._f_list("all"),
            [
                "CLASS",
                "FIELD_0_FILL",
                "FIELD_0_NAME",
                "FIELD_1_FILL",
                "FIELD_1_NAME",
                "FIELD_2_FILL",
                "FIELD_2_NAME",
                "FIELD_3_FILL",
                "FIELD_3_NAME",
                "FIELD_4_FILL",
                "FIELD_4_NAME",
                "NROWS",
                "TITLE",
                "VERSION",
                "pq",
                "qr",
                "rs",
            ],
        )

    def test10b_copyAttributes(self):
        """Checking copying attributes (copy_node_attrs)"""

        # With a Group object
        self.group._v_attrs.pq = "1"
        self.group._v_attrs.qr = "2"
        self.group._v_attrs.rs = "3"
        # copy all attributes from "/agroup" to "/atable"
        self.h5file.copy_node_attrs(self.group, self.root.atable)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
            self.root = self.h5file.root

        atable = self.root.atable
        if common.verbose:
            print("Attribute list:", atable._v_attrs._f_list())
        # Check the local attributes names (alphabetically sorted)
        self.assertEqual(atable._v_attrs._f_list(), ["pq", "qr", "rs"])
        if common.verbose:
            print("Complete attribute list:", atable._v_attrs._f_list("all"))
        # Check the disk attribute names (not sorted)
        self.assertEqual(
            atable._v_attrs._f_list("all"),
            [
                "CLASS",
                "FIELD_0_FILL",
                "FIELD_0_NAME",
                "FIELD_1_FILL",
                "FIELD_1_NAME",
                "FIELD_2_FILL",
                "FIELD_2_NAME",
                "FIELD_3_FILL",
                "FIELD_3_NAME",
                "FIELD_4_FILL",
                "FIELD_4_NAME",
                "NROWS",
                "TITLE",
                "VERSION",
                "pq",
                "qr",
                "rs",
            ],
        )

    def test10c_copyAttributes(self):
        """Checking copying attributes during group copies."""

        # With a Group object
        self.group._v_attrs["CLASS"] = "GROUP2"
        self.group._v_attrs["VERSION"] = "1.3"
        # copy "/agroup" to "/agroup2"
        self.h5file.copy_node(self.group, self.root, "agroup2")

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
            self.root = self.h5file.root

        agroup2 = self.root.agroup2
        if common.verbose:
            print("Complete attribute list:", agroup2._v_attrs._f_list("all"))
        self.assertEqual(agroup2._v_attrs["CLASS"], "GROUP2")
        self.assertEqual(agroup2._v_attrs["VERSION"], "1.3")

    def test10d_copyAttributes(self):
        """Checking copying attributes during leaf copies."""

        # With a Group object
        atable = self.root.atable
        atable._v_attrs["CLASS"] = "TABLE2"
        atable._v_attrs["VERSION"] = "1.3"
        # copy "/agroup" to "/agroup2"
        self.h5file.copy_node(atable, self.root, "atable2")

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+", node_cache_slots=self.node_cache_slots)
            self.root = self.h5file.root

        atable2 = self.root.atable2
        if common.verbose:
            print("Complete attribute list:", atable2._v_attrs._f_list("all"))
        self.assertEqual(atable2._v_attrs["CLASS"], "TABLE2")
        self.assertEqual(atable2._v_attrs["VERSION"], "1.3")

    def test11a_getitem(self):
        """Checking the __getitem__ interface."""

        attrs = self.group._v_attrs
        attrs.pq = "1"
        self.assertEqual(attrs["pq"], "1")

    def test11b_setitem(self):
        """Checking the __setitem__ interface."""

        attrs = self.group._v_attrs
        attrs["pq"] = "2"
        self.assertEqual(attrs["pq"], "2")

    def test11c_delitem(self):
        """Checking the __delitem__ interface."""

        attrs = self.group._v_attrs
        attrs.pq = "1"
        del attrs["pq"]
        self.assertNotIn("pq", attrs._f_list())

    def test11d_KeyError(self):
        """Checking that KeyError is raised in __getitem__/__delitem__."""

        attrs = self.group._v_attrs
        self.assertRaises(KeyError, attrs.__getitem__, "pq")
        self.assertRaises(KeyError, attrs.__delitem__, "pq")

    def test_2d_non_contiguous(self):
        """Checking setting 2D and non-contiguous NumPy attributes"""

        # Regression for gh-176 numpy.
        # In the views old implementation PyTAbles performs a copy of the
        # array:
        #
        #     value = np.array(value)
        #
        # in order to get a contiguous array.
        # Unfortunately array with swapped axis are copied as they are so
        # they are stored in to HDF5 attributes without being actually
        # contiguous and ths causes an error whn they are restored.

        data = np.array([[0, 1], [2, 3]])

        self.array.attrs["a"] = data
        self.array.attrs["b"] = data.T.copy()
        self.array.attrs["c"] = data.T

        np.testing.assert_array_equal(self.array.attrs["a"], data)
        np.testing.assert_array_equal(self.array.attrs["b"], data.T)
        # AssertionError:
        np.testing.assert_array_equal(self.array.attrs["c"], data.T)

    def test12_dir(self):
        """Checking AttributeSet.__dir__"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test12_dir..." % self.__class__.__name__)

        attrset = self.group._v_attrs

        user_attr = "good_attr"
        sys_attr = "BETTER_ATTR"
        for a in [user_attr, sys_attr]:
            attrset[a] = 1

        bad_user = "5bad"
        bad_sys = "SYS%"
        for a in [bad_user, bad_sys]:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", tb.NaturalNameWarning)
                attrset[a] = 1

        completions = dir(attrset)

        # Check some regular attributes.
        self.assertIn("__class__", completions)
        self.assertIn("_f_copy", completions)
        self.assertEqual(completions.count("_f_copy"), 1)

        # Check SYS attrs.
        self.assertNotIn(bad_sys, completions)
        self.assertIn(sys_attr, completions)
        self.assertEqual(completions.count(sys_attr), 1)

        # Check USER attrs.
        self.assertIn(user_attr, completions)
        self.assertNotIn(bad_user, completions)
        self.assertEqual(completions.count(user_attr), 1)

        # Now check all for no duplicates.
        self.assertSequenceEqual(sorted(set(completions)), sorted(completions))


class NotCloseCreate(CreateTestCase):
    close = False
    node_cache_slots = tb.parameters.NODE_CACHE_SLOTS
    open_kwargs = dict(node_cache_slots=node_cache_slots)


class CloseCreate(CreateTestCase):
    close = True
    node_cache_slots = tb.parameters.NODE_CACHE_SLOTS
    open_kwargs = dict(node_cache_slots=node_cache_slots)


class NoCacheNotCloseCreate(CreateTestCase):
    close = False
    node_cache_slots = 0
    open_kwargs = dict(node_cache_slots=node_cache_slots)


class NoCacheCloseCreate(CreateTestCase):
    close = True
    node_cache_slots = 0
    open_kwargs = dict(node_cache_slots=node_cache_slots)


class DictCacheNotCloseCreate(CreateTestCase):
    close = False
    node_cache_slots = -tb.parameters.NODE_CACHE_SLOTS
    open_kwargs = dict(node_cache_slots=node_cache_slots)


class DictCacheCloseCreate(CreateTestCase):
    close = True
    node_cache_slots = -tb.parameters.NODE_CACHE_SLOTS
    open_kwargs = dict(node_cache_slots=node_cache_slots)


class TypesTestCase(common.TempFileMixin, common.PyTablesTestCase):

    def setUp(self):
        self.open_kwargs = {"allow_padding": self.allow_padding}
        super().setUp()
        self.root = self.h5file.root

        # Create an array object
        self.array = self.h5file.create_array(
            self.root, "anarray", [1], "Array title"
        )
        # Create a group object
        self.group = self.h5file.create_group(
            self.root, "agroup", "Group title"
        )

    def test00a_setBoolAttributes(self):
        """Checking setting Bool attributes (scalar, Python case)"""

        self.array.attrs.pq = True
        self.array.attrs.qr = False
        self.array.attrs.rs = True

        # Check the results
        if common.verbose:
            print("pq -->", self.array.attrs.pq)
            print("qr -->", self.array.attrs.qr)
            print("rs -->", self.array.attrs.rs)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        self.assertEqual(self.root.anarray.attrs.pq, True)
        self.assertEqual(self.root.anarray.attrs.qr, False)
        self.assertEqual(self.root.anarray.attrs.rs, True)

    def test00b_setBoolAttributes(self):
        """Checking setting Bool attributes (scalar, NumPy case)"""

        self.array.attrs.pq = np.bool_(True)
        self.array.attrs.qr = np.bool_(False)
        self.array.attrs.rs = np.bool_(True)

        # Check the results
        if common.verbose:
            print("pq -->", self.array.attrs.pq)
            print("qr -->", self.array.attrs.qr)
            print("rs -->", self.array.attrs.rs)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        self.assertIsInstance(self.root.anarray.attrs.pq, np.bool_)
        self.assertIsInstance(self.root.anarray.attrs.qr, np.bool_)
        self.assertIsInstance(self.root.anarray.attrs.rs, np.bool_)
        self.assertEqual(self.root.anarray.attrs.pq, True)
        self.assertEqual(self.root.anarray.attrs.qr, False)
        self.assertEqual(self.root.anarray.attrs.rs, True)

    def test00c_setBoolAttributes(self):
        """Checking setting Bool attributes (NumPy, 0-dim case)"""

        self.array.attrs.pq = np.array(True)
        self.array.attrs.qr = np.array(False)
        self.array.attrs.rs = np.array(True)

        # Check the results
        if common.verbose:
            print("pq -->", self.array.attrs.pq)
            print("qr -->", self.array.attrs.qr)
            print("rs -->", self.array.attrs.rs)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        self.assertEqual(self.root.anarray.attrs.pq, True)
        self.assertEqual(self.root.anarray.attrs.qr, False)
        self.assertEqual(self.root.anarray.attrs.rs, True)

    def test00d_setBoolAttributes(self):
        """Checking setting Bool attributes (NumPy, multidim case)"""

        self.array.attrs.pq = np.array([True])
        self.array.attrs.qr = np.array([[False]])
        self.array.attrs.rs = np.array([[True, False], [True, False]])

        # Check the results
        if common.verbose:
            print("pq -->", self.array.attrs.pq)
            print("qr -->", self.array.attrs.qr)
            print("rs -->", self.array.attrs.rs)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        np.testing.assert_array_equal(
            self.root.anarray.attrs.pq, np.array([True])
        )
        np.testing.assert_array_equal(
            self.root.anarray.attrs.qr, np.array([[False]])
        )
        np.testing.assert_array_equal(
            self.root.anarray.attrs.rs,
            np.array([[True, False], [True, False]]),
        )

    def test01a_setIntAttributes(self):
        """Checking setting Int attributes (scalar, Python case)"""

        self.array.attrs.pq = 1
        self.array.attrs.qr = 2
        self.array.attrs.rs = 3

        # Check the results
        if common.verbose:
            print("pq -->", self.array.attrs.pq)
            print("qr -->", self.array.attrs.qr)
            print("rs -->", self.array.attrs.rs)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        self.assertIsInstance(self.root.anarray.attrs.pq, np.int_)
        self.assertIsInstance(self.root.anarray.attrs.qr, np.int_)
        self.assertIsInstance(self.root.anarray.attrs.rs, np.int_)
        self.assertEqual(self.root.anarray.attrs.pq, 1)
        self.assertEqual(self.root.anarray.attrs.qr, 2)
        self.assertEqual(self.root.anarray.attrs.rs, 3)

    def test01b_setIntAttributes(self):
        """Checking setting Int attributes (scalar, NumPy case)"""

        # 'UInt64' not supported on Win
        checktypes = [
            "int8",
            "int16",
            "int32",
            "int64",
            "uint8",
            "uint16",
            "uint32",
        ]

        for dtype in checktypes:
            setattr(self.array.attrs, dtype, np.array(1, dtype=dtype))

        # Check the results
        if common.verbose:
            for dtype in checktypes:
                print(
                    "type, value-->", dtype, getattr(self.array.attrs, dtype)
                )

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        for dtype in checktypes:
            np.testing.assert_array_equal(
                getattr(self.array.attrs, dtype), np.array(1, dtype=dtype)
            )

    def test01c_setIntAttributes(self):
        """Checking setting Int attributes (unidimensional NumPy case)"""

        # 'UInt64' not supported on Win
        checktypes = [
            "int8",
            "int16",
            "int32",
            "int64",
            "uint8",
            "uint16",
            "uint32",
        ]

        for dtype in checktypes:
            setattr(self.array.attrs, dtype, np.array([1, 2], dtype=dtype))

        # Check the results
        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        for dtype in checktypes:
            if common.verbose:
                print(
                    "type, value-->", dtype, getattr(self.array.attrs, dtype)
                )
            np.testing.assert_array_equal(
                getattr(self.array.attrs, dtype), np.array([1, 2], dtype=dtype)
            )

    def test01d_setIntAttributes(self):
        """Checking setting Int attributes (unidimensional, non-contiguous)"""

        # 'UInt64' not supported on Win
        checktypes = [
            "int8",
            "int16",
            "int32",
            "int64",
            "uint8",
            "uint16",
            "uint32",
        ]

        for dtype in checktypes:
            arr = np.array([1, 2, 3, 4], dtype=dtype)[::2]
            setattr(self.array.attrs, dtype, arr)

        # Check the results
        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        for dtype in checktypes:
            arr = np.array([1, 2, 3, 4], dtype=dtype)[::2]
            if common.verbose:
                print(
                    "type, value-->", dtype, getattr(self.array.attrs, dtype)
                )
            np.testing.assert_array_equal(
                getattr(self.array.attrs, dtype), arr
            )

    def test01e_setIntAttributes(self):
        """Checking setting Int attributes (bidimensional NumPy case)"""

        # 'UInt64' not supported on Win
        checktypes = [
            "int8",
            "int16",
            "int32",
            "int64",
            "uint8",
            "uint16",
            "uint32",
        ]

        for dtype in checktypes:
            setattr(
                self.array.attrs,
                dtype,
                np.array([[1, 2], [2, 3]], dtype=dtype),
            )

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        # Check the results
        for dtype in checktypes:
            if common.verbose:
                print(
                    "type, value-->", dtype, getattr(self.array.attrs, dtype)
                )
            np.testing.assert_array_equal(
                getattr(self.array.attrs, dtype),
                np.array([[1, 2], [2, 3]], dtype=dtype),
            )

    def test02a_setFloatAttributes(self):
        """Checking setting Float (double) attributes."""

        # Set some attrs
        self.array.attrs.pq = 1.0
        self.array.attrs.qr = 2.0
        self.array.attrs.rs = 3.0

        # Check the results
        if common.verbose:
            print("pq -->", self.array.attrs.pq)
            print("qr -->", self.array.attrs.qr)
            print("rs -->", self.array.attrs.rs)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        self.assertIsInstance(self.root.anarray.attrs.pq, np.float64)
        self.assertIsInstance(self.root.anarray.attrs.qr, np.float64)
        self.assertIsInstance(self.root.anarray.attrs.rs, np.float64)
        self.assertEqual(self.root.anarray.attrs.pq, 1.0)
        self.assertEqual(self.root.anarray.attrs.qr, 2.0)
        self.assertEqual(self.root.anarray.attrs.rs, 3.0)

    def test02b_setFloatAttributes(self):
        """Checking setting Float attributes (scalar, NumPy case)"""

        checktypes = ["float32", "float64"]

        for dtype in checktypes:
            setattr(self.array.attrs, dtype, np.array(1.1, dtype=dtype))

        # Check the results
        if common.verbose:
            for dtype in checktypes:
                print(
                    "type, value-->", dtype, getattr(self.array.attrs, dtype)
                )

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        for dtype in checktypes:
            # assert getattr(self.array.attrs, dtype) == 1.1
            # In order to make Float32 tests pass. This is legal, not a trick.
            np.testing.assert_almost_equal(
                getattr(self.array.attrs, dtype), 1.1
            )

    def test02c_setFloatAttributes(self):
        """Checking setting Float attributes (unidimensional NumPy case)"""

        checktypes = ["float32", "float64"]

        for dtype in checktypes:
            setattr(self.array.attrs, dtype, np.array([1.1, 2.1], dtype=dtype))

        # Check the results
        if common.verbose:
            for dtype in checktypes:
                print(
                    "type, value-->", dtype, getattr(self.array.attrs, dtype)
                )

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        for dtype in checktypes:
            np.testing.assert_array_equal(
                getattr(self.array.attrs, dtype),
                np.array([1.1, 2.1], dtype=dtype),
            )

    def test02d_setFloatAttributes(self):
        """Checking setting Float attributes (unidimensional,
        non-contiguous)"""

        checktypes = ["float32", "float64"]

        for dtype in checktypes:
            arr = np.array([1.1, 2.1, 3.1, 4.1], dtype=dtype)[1::2]
            setattr(self.array.attrs, dtype, arr)

        # Check the results
        if common.verbose:
            for dtype in checktypes:
                print(
                    "type, value-->", dtype, getattr(self.array.attrs, dtype)
                )

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        for dtype in checktypes:
            arr = np.array([1.1, 2.1, 3.1, 4.1], dtype=dtype)[1::2]
            np.testing.assert_array_equal(
                getattr(self.array.attrs, dtype), arr
            )

    def test02e_setFloatAttributes(self):
        """Checking setting Int attributes (bidimensional NumPy case)"""

        checktypes = ["float32", "float64"]

        for dtype in checktypes:
            setattr(
                self.array.attrs,
                dtype,
                np.array([[1.1, 2.1], [2.1, 3.1]], dtype=dtype),
            )

        # Check the results
        if common.verbose:
            for dtype in checktypes:
                print(
                    "type, value-->", dtype, getattr(self.array.attrs, dtype)
                )

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        for dtype in checktypes:
            np.testing.assert_array_equal(
                getattr(self.array.attrs, dtype),
                np.array([[1.1, 2.1], [2.1, 3.1]], dtype=dtype),
            )

    def test03_setObjectAttributes(self):
        """Checking setting Object attributes."""

        # Set some attrs
        self.array.attrs.pq = [1.0, 2]
        self.array.attrs.qr = (1, 2)
        self.array.attrs.rs = {"ddf": 32.1, "dsd": 1}

        # Check the results
        if common.verbose:
            print("pq -->", self.array.attrs.pq)
            print("qr -->", self.array.attrs.qr)
            print("rs -->", self.array.attrs.rs)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        self.assertEqual(self.root.anarray.attrs.pq, [1.0, 2])
        self.assertEqual(self.root.anarray.attrs.qr, (1, 2))
        self.assertEqual(self.root.anarray.attrs.rs, {"ddf": 32.1, "dsd": 1})

    def test04a_setStringAttributes(self):
        """Checking setting string attributes (scalar case)"""

        self.array.attrs.pq = "foo"
        self.array.attrs.qr = "bar"
        self.array.attrs.rs = "baz"

        # Check the results
        if common.verbose:
            print("pq -->", self.array.attrs.pq)
            print("qr -->", self.array.attrs.qr)
            print("rs -->", self.array.attrs.rs)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        self.assertIsInstance(self.root.anarray.attrs.pq, np.str_)
        self.assertIsInstance(self.root.anarray.attrs.qr, np.str_)
        self.assertIsInstance(self.root.anarray.attrs.rs, np.str_)
        self.assertEqual(self.root.anarray.attrs.pq, "foo")
        self.assertEqual(self.root.anarray.attrs.qr, "bar")
        self.assertEqual(self.root.anarray.attrs.rs, "baz")

    def test04b_setStringAttributes(self):
        """Checking setting string attributes (unidimensional 1-elem case)"""

        self.array.attrs.pq = np.array(["foo"])

        # Check the results
        if common.verbose:
            print("pq -->", self.array.attrs.pq)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        np.testing.assert_array_equal(
            self.root.anarray.attrs.pq, np.array(["foo"])
        )

    def test04c_setStringAttributes(self):
        """Checking setting string attributes (empty unidimensional
        1-elem case)"""

        self.array.attrs.pq = np.array([""])

        # Check the results
        if common.verbose:
            print("pq -->", self.array.attrs.pq)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray
            if common.verbose:
                print("pq -->", self.array.attrs.pq)

        np.testing.assert_array_equal(
            self.root.anarray.attrs.pq, np.array([""])
        )

    def test04d_setStringAttributes(self):
        """Checking setting string attributes (unidimensional 2-elem case)"""

        self.array.attrs.pq = np.array(["foo", "bar3"])

        # Check the results
        if common.verbose:
            print("pq -->", self.array.attrs.pq)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        np.testing.assert_array_equal(
            self.root.anarray.attrs.pq, np.array(["foo", "bar3"])
        )

    def test04e_setStringAttributes(self):
        """Checking setting string attributes (empty unidimensional
        2-elem case)"""

        self.array.attrs.pq = np.array(["", ""])

        # Check the results
        if common.verbose:
            print("pq -->", self.array.attrs.pq)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        np.testing.assert_array_equal(
            self.root.anarray.attrs.pq, np.array(["", ""])
        )

    def test04f_setStringAttributes(self):
        """Checking setting string attributes (bidimensional 4-elem case)"""

        self.array.attrs.pq = np.array([["foo", "foo2"], ["foo3", "foo4"]])

        # Check the results
        if common.verbose:
            print("pq -->", self.array.attrs.pq)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        np.testing.assert_array_equal(
            self.root.anarray.attrs.pq,
            np.array([["foo", "foo2"], ["foo3", "foo4"]]),
        )

    def test05a_setComplexAttributes(self):
        """Checking setting Complex (python) attributes."""

        # Set some attrs
        self.array.attrs.pq = 1.0 + 2j
        self.array.attrs.qr = 2.0 + 3j
        self.array.attrs.rs = 3.0 + 4j

        # Check the results
        if common.verbose:
            print("pq -->", self.array.attrs.pq)
            print("qr -->", self.array.attrs.qr)
            print("rs -->", self.array.attrs.rs)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        self.assertIsInstance(self.root.anarray.attrs.pq, np.complex128)
        self.assertIsInstance(self.root.anarray.attrs.qr, np.complex128)
        self.assertIsInstance(self.root.anarray.attrs.rs, np.complex128)
        self.assertEqual(self.root.anarray.attrs.pq, 1.0 + 2j)
        self.assertEqual(self.root.anarray.attrs.qr, 2.0 + 3j)
        self.assertEqual(self.root.anarray.attrs.rs, 3.0 + 4j)

    def test05b_setComplexAttributes(self):
        """Checking setting Complex attributes (scalar, NumPy case)"""

        checktypes = ["complex64", "complex128"]

        for dtype in checktypes:
            setattr(self.array.attrs, dtype, np.array(1.1 + 2j, dtype=dtype))

        # Check the results
        if common.verbose:
            for dtype in checktypes:
                print(
                    "type, value-->", dtype, getattr(self.array.attrs, dtype)
                )

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        for dtype in checktypes:
            # assert getattr(self.array.attrs, dtype) == 1.1 + 2j
            # In order to make Complex32 tests pass.
            np.testing.assert_almost_equal(
                getattr(self.array.attrs, dtype), 1.1 + 2j
            )

    def test05c_setComplexAttributes(self):
        """Checking setting Complex attributes (unidimensional NumPy case)"""

        checktypes = ["complex64", "complex128"]

        for dtype in checktypes:
            setattr(self.array.attrs, dtype, np.array([1.1, 2.1], dtype=dtype))

        # Check the results
        if common.verbose:
            for dtype in checktypes:
                print(
                    "type, value-->", dtype, getattr(self.array.attrs, dtype)
                )

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        for dtype in checktypes:
            np.testing.assert_array_equal(
                getattr(self.array.attrs, dtype),
                np.array([1.1, 2.1], dtype=dtype),
            )

    def test05d_setComplexAttributes(self):
        """Checking setting Int attributes (bidimensional NumPy case)"""

        checktypes = ["complex64", "complex128"]

        for dtype in checktypes:
            setattr(
                self.array.attrs,
                dtype,
                np.array([[1.1, 2.1], [2.1, 3.1]], dtype=dtype),
            )

        # Check the results
        if common.verbose:
            for dtype in checktypes:
                print(
                    "type, value-->", dtype, getattr(self.array.attrs, dtype)
                )

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        for dtype in checktypes:
            np.testing.assert_array_equal(
                getattr(self.array.attrs, dtype),
                np.array([[1.1, 2.1], [2.1, 3.1]], dtype=dtype),
            )

    def test06a_setUnicodeAttributes(self):
        """Checking setting unicode attributes (scalar case)"""

        self.array.attrs.pq = "para\u0140lel"
        self.array.attrs.qr = ""  # check #213 or gh-64
        self.array.attrs.rs = "baz"

        # Check the results
        if common.verbose:
            if sys.platform != "win32":
                # It seems that Windows cannot print this
                print("pq -->", repr(self.array.attrs.pq))
                # XXX: try to use repr instead
                # print("pq -->", repr(self.array.attrs.pq))
            print("qr -->", self.array.attrs.qr)
            print("rs -->", self.array.attrs.rs)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        self.assertIsInstance(self.array.attrs.pq, np.str_)
        self.assertIsInstance(self.array.attrs.qr, np.str_)
        self.assertIsInstance(self.array.attrs.rs, np.str_)
        self.assertEqual(self.array.attrs.pq, "para\u0140lel")
        self.assertEqual(self.array.attrs.qr, "")
        self.assertEqual(self.array.attrs.rs, "baz")

    def test06b_setUnicodeAttributes(self):
        """Checking setting unicode attributes (unidimensional 1-elem case)"""

        self.array.attrs.pq = np.array(["para\u0140lel"])

        # Check the results
        if common.verbose:
            print("pq -->", self.array.attrs.pq)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        np.testing.assert_array_equal(
            self.array.attrs.pq, np.array(["para\u0140lel"])
        )

    def test06c_setUnicodeAttributes(self):
        """Checking setting unicode attributes (empty unidimensional
        1-elem case)"""

        # The next raises a `TypeError` when unpickled. See:
        # http://projects.scipy.org/numpy/ticket/1037
        # self.array.attrs.pq = np.array([''])
        self.array.attrs.pq = np.array([""], dtype="U1")

        # Check the results
        if common.verbose:
            print("pq -->", self.array.attrs.pq)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray
            if common.verbose:
                print("pq -->", repr(self.array.attrs.pq))

        np.testing.assert_array_equal(
            self.array.attrs.pq, np.array([""], dtype="U1")
        )

    def test06d_setUnicodeAttributes(self):
        """Checking setting unicode attributes (unidimensional 2-elem case)"""

        self.array.attrs.pq = np.array(["para\u0140lel", "bar3"])

        # Check the results
        if common.verbose:
            print("pq -->", self.array.attrs.pq)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        np.testing.assert_array_equal(
            self.array.attrs.pq, np.array(["para\u0140lel", "bar3"])
        )

    def test06e_setUnicodeAttributes(self):
        """Checking setting unicode attributes (empty unidimensional
        2-elem case)"""

        self.array.attrs.pq = np.array(["", ""], dtype="U1")

        # Check the results
        if common.verbose:
            print("pq -->", self.array.attrs.pq)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        np.testing.assert_array_equal(
            self.array.attrs.pq, np.array(["", ""], dtype="U1")
        )

    def test06f_setUnicodeAttributes(self):
        """Checking setting unicode attributes (bidimensional 4-elem case)"""

        self.array.attrs.pq = np.array(
            [["para\u0140lel", "foo2"], ["foo3", "para\u0140lel4"]]
        )

        # Check the results
        if common.verbose:
            print("pq -->", self.array.attrs.pq)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        np.testing.assert_array_equal(
            self.array.attrs.pq,
            np.array([["para\u0140lel", "foo2"], ["foo3", "para\u0140lel4"]]),
        )

    def test07a_setRecArrayAttributes(self):
        """Checking setting RecArray (NumPy) attributes."""

        dt = np.dtype("i4,f8", align=self.aligned)
        # Set some attrs
        self.array.attrs.pq = np.zeros(2, dt)
        self.array.attrs.qr = np.ones((2, 2), dt)
        self.array.attrs.rs = np.array([(1, 2.0)], dt)

        # Check the results
        if common.verbose:
            print("pq -->", self.array.attrs.pq)
            print("qr -->", self.array.attrs.qr)
            print("rs -->", self.array.attrs.rs)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        self.assertIsInstance(self.array.attrs.pq, np.ndarray)
        self.assertIsInstance(self.array.attrs.qr, np.ndarray)
        self.assertIsInstance(self.array.attrs.rs, np.ndarray)
        np.testing.assert_array_equal(self.array.attrs.pq, np.zeros(2, dt))
        np.testing.assert_array_equal(self.array.attrs.qr, np.ones((2, 2), dt))
        np.testing.assert_array_equal(
            self.array.attrs.rs, np.array([(1, 2.0)], dt)
        )

    def test07b_setRecArrayAttributes(self):
        """Checking setting nested RecArray (NumPy) attributes."""

        # Build a nested dtype
        dt = np.dtype([("f1", [("f1", "i2"), ("f2", "f8")])])
        # Set some attrs
        self.array.attrs.pq = np.zeros(2, dt)
        self.array.attrs.qr = np.ones((2, 2), dt)
        self.array.attrs.rs = np.array([((1, 2.0),)], dt)

        # Check the results
        if common.verbose:
            print("pq -->", self.array.attrs.pq)
            print("qr -->", self.array.attrs.qr)
            print("rs -->", self.array.attrs.rs)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        self.assertIsInstance(self.array.attrs.pq, np.ndarray)
        self.assertIsInstance(self.array.attrs.qr, np.ndarray)
        self.assertIsInstance(self.array.attrs.rs, np.ndarray)
        np.testing.assert_array_equal(self.array.attrs.pq, np.zeros(2, dt))
        np.testing.assert_array_equal(self.array.attrs.qr, np.ones((2, 2), dt))
        np.testing.assert_array_equal(
            self.array.attrs.rs, np.array([((1, 2),)], dt)
        )

    def test07c_setRecArrayAttributes(self):
        """Checking setting multidim nested RecArray (NumPy) attributes."""

        # Build a nested dtype
        dt = np.dtype([("f1", [("f1", "i2", (2,)), ("f2", "f8")])], align=True)

        # Set some attrs
        self.array.attrs.pq = np.zeros(2, dt)
        self.array.attrs.qr = np.ones((2, 2), dt)
        self.array.attrs.rs = np.array([(([1, 3], 2.0),)], dt)

        # Check the results
        if common.verbose:
            print("pq -->", self.array.attrs.pq)
            print("qr -->", self.array.attrs.qr)
            print("rs -->", self.array.attrs.rs)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        self.assertIsInstance(self.array.attrs.pq, np.ndarray)
        self.assertIsInstance(self.array.attrs.qr, np.ndarray)
        self.assertIsInstance(self.array.attrs.rs, np.ndarray)
        np.testing.assert_array_equal(self.array.attrs.pq, np.zeros(2, dt))
        np.testing.assert_array_equal(self.array.attrs.qr, np.ones((2, 2), dt))
        np.testing.assert_array_equal(
            self.array.attrs.rs, np.array([(([1, 3], 2),)], dt)
        )

    def test08_setRecArrayNotAllowPadding(self):
        """Checking setting aligned RecArray (NumPy) attributes with
        `allow_aligned` param set to False when reopened."""

        dt = np.dtype("i4,f8", align=self.aligned)
        # Set some attrs
        self.array.attrs.pq = np.zeros(2, dt)
        self.array.attrs.qr = np.ones((2, 2), dt)
        self.array.attrs.rs = np.array([(1, 2.0)], dt)

        # Check the results
        if common.verbose:
            print("pq -->", self.array.attrs.pq)
            print("qr -->", self.array.attrs.qr)
            print("rs -->", self.array.attrs.rs)

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+", allow_align=False)
            self.root = self.h5file.root
            self.array = self.h5file.root.anarray

        self.assertIsInstance(self.array.attrs.pq, np.ndarray)
        self.assertIsInstance(self.array.attrs.qr, np.ndarray)
        self.assertIsInstance(self.array.attrs.rs, np.ndarray)
        np.testing.assert_array_equal(self.array.attrs.pq, np.zeros(2, dt))
        np.testing.assert_array_equal(self.array.attrs.qr, np.ones((2, 2), dt))
        np.testing.assert_array_equal(
            self.array.attrs.rs, np.array([(1, 2.0)], dt)
        )


class NotCloseTypesTestCase(TypesTestCase):
    allow_padding = False
    aligned = False
    close = False


class NoCloseAlignedTypesTestCase(TypesTestCase):
    allow_padding = True
    aligned = True
    close = False


class CloseNotAlignedPaddedTypesTestCase(TypesTestCase):
    allow_padding = False
    aligned = False
    close = True


class CloseTypesTestCase(TypesTestCase):
    allow_padding = True
    aligned = False
    close = True


class CloseAlignedTypesTestCase(TypesTestCase):
    allow_padding = False
    aligned = True
    close = True


class CloseAlignedPaddedTypesTestCase(TypesTestCase):
    allow_padding = True
    aligned = True
    close = True


class NoSysAttrsTestCase(common.TempFileMixin, common.PyTablesTestCase):
    open_kwargs = dict(pytables_sys_attrs=False)

    def setUp(self):
        super().setUp()
        self.root = self.h5file.root

        # Create a table object
        self.table = self.h5file.create_table(
            self.root, "atable", Record, "Table title"
        )
        # Create an array object
        self.array = self.h5file.create_array(
            self.root, "anarray", [1], "Array title"
        )
        # Create a group object
        self.group = self.h5file.create_group(
            self.root, "agroup", "Group title"
        )

    def test00_listAttributes(self):
        """Checking listing attributes (no system attrs version)."""

        # With a Group object
        self.group._v_attrs.pq = "1"
        self.group._v_attrs.qr = "2"
        self.group._v_attrs.rs = "3"
        if common.verbose:
            print("Attribute list:", self.group._v_attrs._f_list())

        # Now, try with a Table object
        self.table.attrs.a = "1"
        self.table.attrs.c = "2"
        self.table.attrs.b = "3"
        if common.verbose:
            print("Attribute list:", self.table.attrs._f_list())

        # Finally, try with an Array object
        self.array.attrs.k = "1"
        self.array.attrs.j = "2"
        self.array.attrs.i = "3"
        if common.verbose:
            print("Attribute list:", self.array.attrs._f_list())

        if self.close:
            if common.verbose:
                print("(closing file version)")
            self._reopen(mode="r+")
            self.root = self.h5file.root

        agroup = self.root.agroup
        self.assertEqual(agroup._v_attrs._f_list("user"), ["pq", "qr", "rs"])
        self.assertEqual(agroup._v_attrs._f_list("sys"), [])
        self.assertEqual(agroup._v_attrs._f_list("all"), ["pq", "qr", "rs"])

        atable = self.root.atable
        self.assertEqual(atable.attrs._f_list(), ["a", "b", "c"])
        self.assertEqual(atable.attrs._f_list("sys"), [])
        self.assertEqual(atable.attrs._f_list("all"), ["a", "b", "c"])

        anarray = self.root.anarray
        self.assertEqual(anarray.attrs._f_list(), ["i", "j", "k"])
        self.assertEqual(anarray.attrs._f_list("sys"), [])
        self.assertEqual(anarray.attrs._f_list("all"), ["i", "j", "k"])


class NoSysAttrsNotClose(NoSysAttrsTestCase):
    close = False


class NoSysAttrsClose(NoSysAttrsTestCase):
    close = True


class CompatibilityTestCase(common.TestFileMixin, common.PyTablesTestCase):
    h5fname = common.test_filename("issue_368.h5")

    @common.unittest.skipIf(
        Version(np.__version__) < Version("1.9.0"), "requires numpy >= 1.9"
    )
    def test_pickled_unicode_attrs(self):
        # See also gh-368 and https://github.com/numpy/numpy/issues/4879.
        #
        # This is a compatibility test. In PyTables < 3.0 unicode
        # attributes were stored as pickled unicode strings.
        # In PyTables >= 3.0 unicode strings are stored as encoded utf-8
        # strings (the utf-8 marker is set at HDF5 level).
        #
        # In any case PyTables (>= 3.0) should be able to handle correctly
        # also data files generated with older versions of PyTables.
        # Unfortunately a bug in numpy < 1.9
        # (https://github.com/numpy/numpy/issues/4879) makes it impossible
        # unpickle numpy arrays with dtype "U" resulting in an incorrect
        # behaviour of PyTables.

        self.assertEqual(
            self.h5file.get_node_attr("/", "py2_pickled_unicode"), "abc"
        )


class PicklePy2UnpicklePy3TestCase(
    common.TestFileMixin, common.PyTablesTestCase
):
    h5fname = common.test_filename("issue_560.h5")

    def test_pickled_datetime_object(self):
        # See also gh-560
        #
        # Objects (classes) that are pickled using python 2 may contain
        # non-ascii characters in the pickled string. This will cause
        # a UnicodeDecodeError when unpickling on python 3.
        # Python 3.4 adds encoding='bytes' to fix this
        # http://bugs.python.org/issue6784
        # Objects pickled in the testfile have non-ascii chars in the
        # picklestring and will throw UnicodeDecodeError when unpickled
        # on python 3.

        # datetime will be unpickled with encoding='bytes'
        self.assertIsInstance(
            self.h5file.get_node_attr("/", "py2_pickled_datetime"),
            datetime.datetime,
        )
        # dict will be unpickled with encoding='latin1'
        d = self.h5file.get_node_attr("/", "py2_pickled_dict")
        self.assertIsInstance(d, dict)
        self.assertEqual(d["s"], "just a string")


class SegFaultPythonTestCase(common.TempFileMixin, common.PyTablesTestCase):

    def test00_segfault(self):
        """Checking workaround for Python unpickle problem (see #253)."""

        self.h5file.root._v_attrs.trouble1 = "0"
        self.assertEqual(self.h5file.root._v_attrs.trouble1, "0")
        self.h5file.root._v_attrs.trouble2 = "0."
        self.assertEqual(self.h5file.root._v_attrs.trouble2, "0.")
        # Problem happens after reopening
        self._reopen()
        self.assertEqual(self.h5file.root._v_attrs.trouble1, "0")
        self.assertEqual(self.h5file.root._v_attrs.trouble2, "0.")
        if common.verbose:
            print("Great! '0' and '0.' values can be safely retrieved.")


class EmbeddedNullsTestCase(common.TempFileMixin, common.PyTablesTestCase):
    # See laso gh-371 (https://github.com/PyTables/PyTables/issues/371)

    def test_unicode(self):
        value = "string with a null byte \x00 in it"

        self.h5file.root._v_attrs.name = value
        self.assertEqual(self.h5file.root._v_attrs.name, value)

        self._reopen()

        self.assertEqual(self.h5file.root._v_attrs.name, value)

    def test_bytes(self):
        value = b"string with a null byte \x00 in it"

        self.h5file.root._v_attrs.name = value
        self.assertEqual(self.h5file.root._v_attrs.name, value)

        self._reopen()

        self.assertEqual(self.h5file.root._v_attrs.name, value)


class VlenStrAttrTestCase(common.PyTablesTestCase):
    def setUp(self):
        super().setUp()
        self.h5fname = common.test_filename("vlstr_attr.h5")
        self.h5file = tb.open_file(self.h5fname)

    def tearDown(self):
        self.h5file.close()
        super().tearDown()

    def test01_vlen_str_scalar(self):
        """Checking file with variable length string attributes."""

        attr = "vlen_str_scalar"
        self.assertEqual(
            self.h5file.get_node_attr("/", attr), attr.encode("ascii")
        )

    def test02_vlen_str_array(self):
        """Checking file with variable length string attributes (1d)."""

        attr = "vlen_str_array"
        v = self.h5file.get_node_attr("/", attr)
        self.assertEqual(v.ndim, 1)
        for idx, item in enumerate(v):
            value = "%s_%d" % (attr, idx)
            self.assertEqual(item, value.encode("ascii"))

    def test03_vlen_str_matrix(self):
        """Checking file with variable length string attributes (2d)."""

        attr = "vlen_str_matrix"
        m = self.h5file.get_node_attr("/", attr)
        self.assertEqual(m.ndim, 2)
        for row, rowdata in enumerate(m):
            for col, item in enumerate(rowdata):
                value = "%s_%d%d" % (attr, row, col)
                self.assertEqual(item, value.encode("ascii"))


class UnsupportedAttrTypeTestCase(
    common.TestFileMixin, common.PyTablesTestCase
):
    h5fname = common.test_filename("attr-u16.h5")

    def test00_unsupportedType(self):
        """Checking file with unsupported type."""

        self.assertWarns(tb.exceptions.DataTypeWarning, repr, self.h5file)


# Test for specific system attributes
class SpecificAttrsTestCase(common.TempFileMixin, common.PyTablesTestCase):

    def test00_earray(self):
        """Testing EArray specific attrs (create)."""

        ea = self.h5file.create_earray("/", "ea", tb.Int32Atom(), (2, 0, 4))
        if common.verbose:
            print("EXTDIM-->", ea.attrs.EXTDIM)
        self.assertEqual(ea.attrs.EXTDIM, 1)

    def test01_earray(self):
        """Testing EArray specific attrs (open)."""

        ea = self.h5file.create_earray("/", "ea", tb.Int32Atom(), (0, 1, 4))
        self._reopen("r")
        ea = self.h5file.root.ea
        if common.verbose:
            print("EXTDIM-->", ea.attrs.EXTDIM)
        self.assertEqual(ea.attrs.EXTDIM, 0)


def suite():
    theSuite = common.unittest.TestSuite()
    niter = 1

    for i in range(niter):
        theSuite.addTest(common.make_suite(NotCloseCreate))
        theSuite.addTest(common.make_suite(CloseCreate))
        theSuite.addTest(common.make_suite(NoCacheNotCloseCreate))
        theSuite.addTest(common.make_suite(NoCacheCloseCreate))
        theSuite.addTest(common.make_suite(DictCacheNotCloseCreate))
        theSuite.addTest(common.make_suite(DictCacheCloseCreate))
        theSuite.addTest(common.make_suite(NotCloseTypesTestCase))
        theSuite.addTest(common.make_suite(CloseTypesTestCase))
        theSuite.addTest(common.make_suite(CloseNotAlignedPaddedTypesTestCase))
        theSuite.addTest(common.make_suite(NoCloseAlignedTypesTestCase))
        theSuite.addTest(common.make_suite(CloseAlignedTypesTestCase))
        theSuite.addTest(common.make_suite(CloseAlignedPaddedTypesTestCase))
        theSuite.addTest(common.make_suite(NoSysAttrsNotClose))
        theSuite.addTest(common.make_suite(NoSysAttrsClose))
        theSuite.addTest(common.make_suite(CompatibilityTestCase))
        theSuite.addTest(common.make_suite(PicklePy2UnpicklePy3TestCase))
        theSuite.addTest(common.make_suite(SegFaultPythonTestCase))
        theSuite.addTest(common.make_suite(EmbeddedNullsTestCase))
        theSuite.addTest(common.make_suite(VlenStrAttrTestCase))
        theSuite.addTest(common.make_suite(UnsupportedAttrTypeTestCase))
        theSuite.addTest(common.make_suite(SpecificAttrsTestCase))

    return theSuite


if __name__ == "__main__":
    common.parse_argv(sys.argv)
    common.print_versions()
    common.unittest.main(defaultTest="suite")
