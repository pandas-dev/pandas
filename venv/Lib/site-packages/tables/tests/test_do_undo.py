import warnings

import tables as tb
from tables.tests import common


class BasicTestCase(common.TempFileMixin, common.PyTablesTestCase):
    """Test for basic Undo/Redo operations."""

    _reopen_flag = False
    """Whether to reopen the file at certain points."""

    def _do_reopen(self):
        if self._reopen_flag:
            self._reopen("r+")

    def setUp(self):
        super().setUp()

        h5file = self.h5file
        root = h5file.root

        # Create an array
        h5file.create_array(root, "array", [1, 2], title="Title example")

        # Create another array object
        h5file.create_array(root, "anarray", [1], "Array title")

        # Create a group object
        group = h5file.create_group(root, "agroup", "Group title")

        # Create a couple of objects there
        h5file.create_array(group, "anarray1", [2], "Array title 1")
        h5file.create_array(group, "anarray2", [2], "Array title 2")

        # Create a lonely group in first level
        h5file.create_group(root, "agroup2", "Group title 2")

        # Create a new group in the second level
        h5file.create_group(group, "agroup3", "Group title 3")

    def test00_simple(self):
        """Checking simple do/undo."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test00_simple..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.create_array("/", "otherarray", [3, 4], "Another array")

        # Now undo the past operation
        self.h5file.undo()

        # Check that otherarray does not exist in the object tree
        self.assertNotIn("/otherarray", self.h5file)
        self.assertEqual(self.h5file._curaction, 0)
        self.assertEqual(self.h5file._curmark, 0)

        # Redo the operation
        self._do_reopen()
        self.h5file.redo()
        if common.verbose:
            print("Object tree after redo:", self.h5file)

        # Check that otherarray has come back to life in a sane state
        self.assertIn("/otherarray", self.h5file)
        self.assertEqual(self.h5file.root.otherarray.read(), [3, 4])
        self.assertEqual(self.h5file.root.otherarray.title, "Another array")
        self.assertEqual(self.h5file._curaction, 1)
        self.assertEqual(self.h5file._curmark, 0)

    def test01_twice(self):
        """Checking do/undo (twice operations intertwined)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01_twice..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.create_array("/", "otherarray", [3, 4], "Another array")
        self.h5file.create_array("/", "otherarray2", [4, 5], "Another array 2")

        # Now undo the past operations
        self._do_reopen()
        self.h5file.undo()
        self.assertNotIn("/otherarray", self.h5file)
        self.assertNotIn("/otherarray2", self.h5file)
        self.assertEqual(self.h5file._curaction, 0)
        self.assertEqual(self.h5file._curmark, 0)

        # Redo the operation
        self.h5file.redo()

        # Check that otherarray has come back to life in a sane state
        self.assertIn("/otherarray", self.h5file)
        self.assertIn("/otherarray2", self.h5file)
        self.assertEqual(self.h5file.root.otherarray.read(), [3, 4])
        self.assertEqual(self.h5file.root.otherarray2.read(), [4, 5])
        self.assertEqual(self.h5file.root.otherarray.title, "Another array")
        self.assertEqual(self.h5file.root.otherarray2.title, "Another array 2")
        self.assertEqual(self.h5file._curaction, 2)
        self.assertEqual(self.h5file._curmark, 0)

    def test02_twice2(self):
        """Checking twice ops and two marks."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02_twice2..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.create_array("/", "otherarray", [3, 4], "Another array")

        # Put a mark
        self._do_reopen()
        self.h5file.mark()
        self.h5file.create_array("/", "otherarray2", [4, 5], "Another array 2")
        self.assertEqual(self.h5file._curaction, 3)
        self.assertEqual(self.h5file._curmark, 1)

        # Unwind just one mark
        self.h5file.undo()
        self.assertIn("/otherarray", self.h5file)
        self.assertNotIn("/otherarray2", self.h5file)
        self.assertEqual(self.h5file._curaction, 2)
        self.assertEqual(self.h5file._curmark, 1)

        # Unwind another mark
        self.h5file.undo()
        self.assertEqual(self.h5file._curaction, 0)
        self.assertEqual(self.h5file._curmark, 0)
        self.assertNotIn("/otherarray", self.h5file)
        self.assertNotIn("/otherarray2", self.h5file)

        # Redo until the next mark
        self.h5file.redo()
        self.assertIn("/otherarray", self.h5file)
        self.assertNotIn("/otherarray2", self.h5file)
        self._do_reopen()
        self.assertEqual(self.h5file._curaction, 2)
        self.assertEqual(self.h5file._curmark, 1)

        # Redo until the end
        self.h5file.redo()
        self.assertIn("/otherarray", self.h5file)
        self.assertIn("/otherarray2", self.h5file)
        self.assertEqual(self.h5file.root.otherarray.read(), [3, 4])
        self.assertEqual(self.h5file.root.otherarray2.read(), [4, 5])
        self.assertEqual(self.h5file.root.otherarray.title, "Another array")
        self.assertEqual(self.h5file.root.otherarray2.title, "Another array 2")
        self.assertEqual(self.h5file._curaction, 3)
        self.assertEqual(self.h5file._curmark, 1)

    def test03_6times3marks(self):
        """Checking with six ops and three marks."""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test03_6times3marks..." % self.__class__.__name__
            )

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.create_array("/", "otherarray1", [3, 4], "Another array 1")
        self.h5file.create_array("/", "otherarray2", [4, 5], "Another array 2")

        # Put a mark
        self.h5file.mark()
        self.h5file.create_array("/", "otherarray3", [5, 6], "Another array 3")
        self.h5file.create_array("/", "otherarray4", [6, 7], "Another array 4")

        # Put a mark
        self._do_reopen()
        self.h5file.mark()
        self.h5file.create_array("/", "otherarray5", [7, 8], "Another array 5")
        self.h5file.create_array("/", "otherarray6", [8, 9], "Another array 6")

        # Unwind just one mark
        self.h5file.undo()
        self.assertIn("/otherarray1", self.h5file)
        self.assertIn("/otherarray2", self.h5file)
        self.assertIn("/otherarray3", self.h5file)
        self.assertIn("/otherarray4", self.h5file)
        self.assertNotIn("/otherarray5", self.h5file)
        self.assertNotIn("/otherarray6", self.h5file)

        # Unwind another mark
        self.h5file.undo()
        self.assertIn("/otherarray1", self.h5file)
        self.assertIn("/otherarray2", self.h5file)
        self.assertNotIn("/otherarray3", self.h5file)
        self.assertNotIn("/otherarray4", self.h5file)
        self.assertNotIn("/otherarray5", self.h5file)
        self.assertNotIn("/otherarray6", self.h5file)

        # Unwind all marks
        self.h5file.undo()
        self.assertNotIn("/otherarray1", self.h5file)
        self.assertNotIn("/otherarray2", self.h5file)
        self.assertNotIn("/otherarray3", self.h5file)
        self.assertNotIn("/otherarray4", self.h5file)
        self.assertNotIn("/otherarray5", self.h5file)
        self.assertNotIn("/otherarray6", self.h5file)

        # Redo until the next mark
        self._do_reopen()
        self.h5file.redo()
        self.assertIn("/otherarray1", self.h5file)
        self.assertIn("/otherarray2", self.h5file)
        self.assertNotIn("/otherarray3", self.h5file)
        self.assertNotIn("/otherarray4", self.h5file)
        self.assertNotIn("/otherarray5", self.h5file)
        self.assertNotIn("/otherarray6", self.h5file)

        # Redo until the next mark
        self.h5file.redo()
        self.assertIn("/otherarray1", self.h5file)
        self.assertIn("/otherarray2", self.h5file)
        self.assertIn("/otherarray3", self.h5file)
        self.assertIn("/otherarray4", self.h5file)
        self.assertNotIn("/otherarray5", self.h5file)
        self.assertNotIn("/otherarray6", self.h5file)

        # Redo until the end
        self.h5file.redo()
        self.assertIn("/otherarray1", self.h5file)
        self.assertIn("/otherarray2", self.h5file)
        self.assertIn("/otherarray3", self.h5file)
        self.assertIn("/otherarray4", self.h5file)
        self.assertIn("/otherarray5", self.h5file)
        self.assertIn("/otherarray6", self.h5file)
        self.assertEqual(self.h5file.root.otherarray1.read(), [3, 4])
        self.assertEqual(self.h5file.root.otherarray2.read(), [4, 5])
        self.assertEqual(self.h5file.root.otherarray3.read(), [5, 6])
        self.assertEqual(self.h5file.root.otherarray4.read(), [6, 7])
        self.assertEqual(self.h5file.root.otherarray5.read(), [7, 8])
        self.assertEqual(self.h5file.root.otherarray6.read(), [8, 9])
        self.assertEqual(self.h5file.root.otherarray1.title, "Another array 1")
        self.assertEqual(self.h5file.root.otherarray2.title, "Another array 2")
        self.assertEqual(self.h5file.root.otherarray3.title, "Another array 3")
        self.assertEqual(self.h5file.root.otherarray4.title, "Another array 4")
        self.assertEqual(self.h5file.root.otherarray5.title, "Another array 5")
        self.assertEqual(self.h5file.root.otherarray6.title, "Another array 6")

    def test04_6times3marksro(self):
        """Checking with six operations, three marks and do/undo in random
        order."""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test04_6times3marksro..." % self.__class__.__name__
            )

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.create_array("/", "otherarray1", [3, 4], "Another array 1")
        self.h5file.create_array("/", "otherarray2", [4, 5], "Another array 2")

        # Put a mark
        self.h5file.mark()
        self._do_reopen()
        self.h5file.create_array("/", "otherarray3", [5, 6], "Another array 3")
        self.h5file.create_array("/", "otherarray4", [6, 7], "Another array 4")

        # Unwind the previous mark
        self.h5file.undo()
        self.assertIn("/otherarray1", self.h5file)
        self.assertIn("/otherarray2", self.h5file)
        self.assertNotIn("/otherarray3", self.h5file)
        self.assertNotIn("/otherarray4", self.h5file)

        # Put a mark in the middle of stack
        if common.verbose:
            print("All nodes:", self.h5file.walk_nodes())
        self.h5file.mark()
        self._do_reopen()
        self.h5file.create_array("/", "otherarray5", [7, 8], "Another array 5")
        self.h5file.create_array("/", "otherarray6", [8, 9], "Another array 6")
        self.assertIn("/otherarray1", self.h5file)
        self.assertIn("/otherarray2", self.h5file)
        self.assertNotIn("/otherarray3", self.h5file)
        self.assertNotIn("/otherarray4", self.h5file)
        self.assertIn("/otherarray5", self.h5file)
        self.assertIn("/otherarray6", self.h5file)

        # Unwind previous mark
        self.h5file.undo()
        self.assertIn("/otherarray1", self.h5file)
        self.assertIn("/otherarray2", self.h5file)
        self.assertNotIn("/otherarray3", self.h5file)
        self.assertNotIn("/otherarray4", self.h5file)
        self.assertNotIn("/otherarray5", self.h5file)
        self.assertNotIn("/otherarray6", self.h5file)

        # Redo until the last mark
        self.h5file.redo()
        self.assertIn("/otherarray1", self.h5file)
        self.assertIn("/otherarray2", self.h5file)
        self.assertNotIn("/otherarray3", self.h5file)
        self.assertNotIn("/otherarray4", self.h5file)
        self.assertIn("/otherarray5", self.h5file)
        self.assertIn("/otherarray6", self.h5file)

        # Redo until the next mark (non-existent, so no action)
        self._do_reopen()
        self.h5file.redo()
        self.assertIn("/otherarray1", self.h5file)
        self.assertIn("/otherarray2", self.h5file)
        self.assertNotIn("/otherarray3", self.h5file)
        self.assertNotIn("/otherarray4", self.h5file)
        self.assertIn("/otherarray5", self.h5file)
        self.assertIn("/otherarray6", self.h5file)
        self.assertEqual(self.h5file.root.otherarray1.read(), [3, 4])
        self.assertEqual(self.h5file.root.otherarray2.read(), [4, 5])
        self.assertEqual(self.h5file.root.otherarray5.read(), [7, 8])
        self.assertEqual(self.h5file.root.otherarray6.read(), [8, 9])
        self.assertEqual(self.h5file.root.otherarray1.title, "Another array 1")
        self.assertEqual(self.h5file.root.otherarray2.title, "Another array 2")
        self.assertEqual(self.h5file.root.otherarray5.title, "Another array 5")
        self.assertEqual(self.h5file.root.otherarray6.title, "Another array 6")

    def test05_destructive(self):
        """Checking with a destructive action during undo."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test05_destructive..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.create_array("/", "otherarray1", [3, 4], "Another array 1")

        # Put a mark
        self.h5file.mark()
        self._do_reopen()
        self.h5file.create_array("/", "otherarray2", [4, 5], "Another array 2")

        # Now undo the past operation
        self.h5file.undo()

        # Do the destructive operation
        self._do_reopen()
        self.h5file.create_array("/", "otherarray3", [5, 6], "Another array 3")

        # Check objects
        self.assertIn("/otherarray1", self.h5file)
        self.assertEqual(self.h5file.root.otherarray1.read(), [3, 4])
        self.assertEqual(self.h5file.root.otherarray1.title, "Another array 1")
        self.assertNotIn("/otherarray2", self.h5file)
        self.assertIn("/otherarray3", self.h5file)
        self.assertEqual(self.h5file.root.otherarray3.read(), [5, 6])
        self.assertEqual(self.h5file.root.otherarray3.title, "Another array 3")

    def test05b_destructive(self):
        """Checking with a destructive action during undo (II)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test05b_destructive..." % self.__class__.__name__
            )

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.create_array("/", "otherarray1", [3, 4], "Another array 1")

        # Put a mark
        self._do_reopen()
        self.h5file.mark()
        self.h5file.create_array("/", "otherarray2", [4, 5], "Another array 2")

        # Now undo the past operation
        self.h5file.undo()

        # Do the destructive operation
        self.h5file.create_array("/", "otherarray3", [5, 6], "Another array 3")

        # Put a mark
        self._do_reopen()
        self.h5file.mark()
        self.h5file.create_array("/", "otherarray4", [6, 7], "Another array 4")
        self.assertIn("/otherarray4", self.h5file)

        # Now undo the past operation
        self.h5file.undo()

        # Check objects
        self.assertIn("/otherarray1", self.h5file)
        self.assertEqual(self.h5file.root.otherarray1.read(), [3, 4])
        self.assertEqual(self.h5file.root.otherarray1.title, "Another array 1")
        self.assertNotIn("/otherarray2", self.h5file)
        self.assertIn("/otherarray3", self.h5file)
        self.assertEqual(self.h5file.root.otherarray3.read(), [5, 6])
        self.assertEqual(self.h5file.root.otherarray3.title, "Another array 3")
        self.assertNotIn("/otherarray4", self.h5file)

    def test05c_destructive(self):
        """Checking with a destructive action during undo (III)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test05c_destructive..." % self.__class__.__name__
            )

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.create_array("/", "otherarray1", [3, 4], "Another array 1")

        # Put a mark
        self.h5file.mark()
        self._do_reopen()
        self.h5file.create_array("/", "otherarray2", [4, 5], "Another array 2")

        # Now undo the past operation
        self.h5file.undo()

        # Do the destructive operation
        self.h5file.create_array("/", "otherarray3", [5, 6], "Another array 3")

        # Put a mark
        self.h5file.mark()
        self._do_reopen()
        self.h5file.create_array("/", "otherarray4", [6, 7], "Another array 4")
        self.assertIn("/otherarray4", self.h5file)

        # Now unwind twice
        self.h5file.undo()
        self._do_reopen()
        self.h5file.undo()

        # Check objects
        self.assertIn("/otherarray1", self.h5file)
        self.assertNotIn("/otherarray2", self.h5file)
        self.assertNotIn("/otherarray3", self.h5file)
        self.assertNotIn("/otherarray4", self.h5file)

    def test05d_destructive(self):
        """Checking with a destructive action during undo (IV)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test05d_destructive..." % self.__class__.__name__
            )

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.create_array("/", "otherarray1", [3, 4], "Another array 1")

        # Put a mark
        self._do_reopen()
        self.h5file.mark()
        self.h5file.create_array("/", "otherarray2", [4, 5], "Another array 2")

        # Now undo the past operation
        self.h5file.undo()

        # Do the destructive operation
        self.h5file.create_array("/", "otherarray3", [5, 6], "Another array 3")

        # Put a mark
        self.h5file.mark()
        self.h5file.create_array("/", "otherarray4", [6, 7], "Another array 4")
        self.assertIn("/otherarray4", self.h5file)

        # Now, go to the first mark
        self._do_reopen()
        self.h5file.undo(0)

        # Check objects
        self.assertNotIn("/otherarray1", self.h5file)
        self.assertNotIn("/otherarray2", self.h5file)
        self.assertNotIn("/otherarray3", self.h5file)
        self.assertNotIn("/otherarray4", self.h5file)

    def test05e_destructive(self):
        """Checking with a destructive action during undo (V)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test05e_destructive..." % self.__class__.__name__
            )

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.create_array("/", "otherarray1", [3, 4], "Another array 1")

        # Put a mark
        self.h5file.mark()
        self.h5file.create_array("/", "otherarray2", [4, 5], "Another array 2")

        # Now undo the past operation
        self.h5file.undo()
        self._do_reopen()

        # Do the destructive operation
        self.h5file.create_array("/", "otherarray3", [5, 6], "Another array 3")

        # Now, unwind the actions
        self.h5file.undo(0)
        self._do_reopen()

        # Check objects
        self.assertNotIn("/otherarray1", self.h5file)
        self.assertNotIn("/otherarray2", self.h5file)
        self.assertNotIn("/otherarray3", self.h5file)

    def test05f_destructive(self):
        """Checking with a destructive creation of existing node during undo"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test05f_destructive..." % self.__class__.__name__
            )

        self.h5file.enable_undo()
        self.h5file.create_array("/", "newarray", [1])
        self.h5file.undo()
        self._do_reopen()
        self.assertNotIn("/newarray", self.h5file)
        newarr = self.h5file.create_array("/", "newarray", [1])
        self.h5file.undo()
        self.assertNotIn("/newarray", self.h5file)
        self._do_reopen()
        self.h5file.redo()
        self.assertIn("/newarray", self.h5file)
        if not self._reopen_flag:
            self.assertIs(self.h5file.root.newarray, newarr)

    def test06_totalunwind(self):
        """Checking do/undo (total unwind)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test06_totalunwind..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.create_array("/", "otherarray", [3, 4], "Another array")
        self.h5file.mark()
        self.h5file.create_array("/", "otherarray2", [4, 5], "Another array 2")

        # Now undo the past operations
        self._do_reopen()
        self.h5file.undo(0)
        self.assertNotIn("/otherarray", self.h5file)
        self.assertNotIn("/otherarray2", self.h5file)

    def test07_totalrewind(self):
        """Checking do/undo (total rewind)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test07_totalunwind..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.create_array("/", "otherarray", [3, 4], "Another array")
        self.h5file.mark()
        self.h5file.create_array("/", "otherarray2", [4, 5], "Another array 2")

        # Now undo the past operations
        self.h5file.undo(0)

        # Redo all the operations
        self._do_reopen()
        self.h5file.redo(-1)

        # Check that objects has come back to life in a sane state
        self.assertIn("/otherarray", self.h5file)
        self.assertIn("/otherarray2", self.h5file)
        self.assertEqual(self.h5file.root.otherarray.read(), [3, 4])
        self.assertEqual(self.h5file.root.otherarray2.read(), [4, 5])
        self.assertEqual(self.h5file.root.otherarray.title, "Another array")
        self.assertEqual(self.h5file.root.otherarray2.title, "Another array 2")

    def test08_marknames(self):
        """Checking mark names."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test08_marknames..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.create_array("/", "otherarray1", [3, 4], "Another array 1")
        self.h5file.mark("first")
        self.h5file.create_array("/", "otherarray2", [4, 5], "Another array 2")
        self.h5file.mark("second")
        self.h5file.create_array("/", "otherarray3", [5, 6], "Another array 3")
        self.h5file.mark("third")
        self.h5file.create_array("/", "otherarray4", [6, 7], "Another array 4")

        # Now go to mark "first"
        self.h5file.undo("first")
        self._do_reopen()
        self.assertIn("/otherarray1", self.h5file)
        self.assertNotIn("/otherarray2", self.h5file)
        self.assertNotIn("/otherarray3", self.h5file)
        self.assertNotIn("/otherarray4", self.h5file)

        # Go to mark "third"
        self.h5file.redo("third")
        self.assertIn("/otherarray1", self.h5file)
        self.assertIn("/otherarray2", self.h5file)
        self.assertIn("/otherarray3", self.h5file)
        self.assertNotIn("/otherarray4", self.h5file)

        # Now go to mark "second"
        self.h5file.undo("second")
        self.assertIn("/otherarray1", self.h5file)
        self.assertIn("/otherarray2", self.h5file)
        self.assertNotIn("/otherarray3", self.h5file)
        self.assertNotIn("/otherarray4", self.h5file)

        # Go to the end
        self._do_reopen()
        self.h5file.redo(-1)
        self.assertIn("/otherarray1", self.h5file)
        self.assertIn("/otherarray2", self.h5file)
        self.assertIn("/otherarray3", self.h5file)
        self.assertIn("/otherarray4", self.h5file)

        # Check that objects has come back to life in a sane state
        self.assertEqual(self.h5file.root.otherarray1.read(), [3, 4])
        self.assertEqual(self.h5file.root.otherarray2.read(), [4, 5])
        self.assertEqual(self.h5file.root.otherarray3.read(), [5, 6])
        self.assertEqual(self.h5file.root.otherarray4.read(), [6, 7])

    def test08_initialmark(self):
        """Checking initial mark."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test08_initialmark..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()
        initmid = self.h5file.get_current_mark()

        # Create a new array
        self.h5file.create_array("/", "otherarray", [3, 4], "Another array")
        self.h5file.mark()
        self._do_reopen()
        self.h5file.create_array("/", "otherarray2", [4, 5], "Another array 2")

        # Now undo the past operations
        self.h5file.undo(initmid)
        self.assertNotIn("/otherarray", self.h5file)
        self.assertNotIn("/otherarray2", self.h5file)

        # Redo all the operations
        self.h5file.redo(-1)
        self._do_reopen()

        # Check that objects has come back to life in a sane state
        self.assertIn("/otherarray", self.h5file)
        self.assertIn("/otherarray2", self.h5file)
        self.assertEqual(self.h5file.root.otherarray.read(), [3, 4])
        self.assertEqual(self.h5file.root.otherarray2.read(), [4, 5])
        self.assertEqual(self.h5file.root.otherarray.title, "Another array")
        self.assertEqual(self.h5file.root.otherarray2.title, "Another array 2")

    def test09_marknames(self):
        """Checking mark names (wrong direction)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test09_marknames..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.create_array("/", "otherarray1", [3, 4], "Another array 1")
        self.h5file.mark("first")
        self.h5file.create_array("/", "otherarray2", [4, 5], "Another array 2")
        self.h5file.mark("second")
        self._do_reopen()
        self.h5file.create_array("/", "otherarray3", [5, 6], "Another array 3")
        self.h5file.mark("third")
        self.h5file.create_array("/", "otherarray4", [6, 7], "Another array 4")

        # Now go to mark "first"
        self.h5file.undo("first")

        # Try to undo up to mark "third"
        with self.assertRaises(tb.UndoRedoError):
            self.h5file.undo("third")

        # Now go to mark "third"
        self.h5file.redo("third")
        self._do_reopen()

        # Try to redo up to mark "second"
        with self.assertRaises(tb.UndoRedoError):
            self.h5file.redo("second")

        # Final checks
        self.assertIn("/otherarray1", self.h5file)
        self.assertIn("/otherarray2", self.h5file)
        self.assertIn("/otherarray3", self.h5file)
        self.assertNotIn("/otherarray4", self.h5file)

    def test10_goto(self):
        """Checking mark names (goto)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test10_goto..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.create_array("/", "otherarray1", [3, 4], "Another array 1")
        self._do_reopen()
        self.h5file.mark("first")
        self.h5file.create_array("/", "otherarray2", [4, 5], "Another array 2")
        self.h5file.mark("second")
        self.h5file.create_array("/", "otherarray3", [5, 6], "Another array 3")
        self._do_reopen()
        self.h5file.mark("third")
        self.h5file.create_array("/", "otherarray4", [6, 7], "Another array 4")

        # Now go to mark "first"
        self.h5file.goto("first")
        self.assertIn("/otherarray1", self.h5file)
        self.assertNotIn("/otherarray2", self.h5file)
        self.assertNotIn("/otherarray3", self.h5file)
        self.assertNotIn("/otherarray4", self.h5file)

        # Go to mark "third"
        self.h5file.goto("third")
        self.assertIn("/otherarray1", self.h5file)
        self.assertIn("/otherarray2", self.h5file)
        self.assertIn("/otherarray3", self.h5file)
        self.assertNotIn("/otherarray4", self.h5file)

        # Now go to mark "second"
        self._do_reopen()
        self.h5file.goto("second")
        self.assertIn("/otherarray1", self.h5file)
        self.assertIn("/otherarray2", self.h5file)
        self.assertNotIn("/otherarray3", self.h5file)
        self.assertNotIn("/otherarray4", self.h5file)

        # Go to the end
        self.h5file.goto(-1)
        self.assertIn("/otherarray1", self.h5file)
        self.assertIn("/otherarray2", self.h5file)
        self.assertIn("/otherarray3", self.h5file)
        self.assertIn("/otherarray4", self.h5file)

        # Check that objects has come back to life in a sane state
        self.assertIn("/otherarray2", self.h5file)
        self.assertEqual(self.h5file.root.otherarray1.read(), [3, 4])
        self.assertEqual(self.h5file.root.otherarray2.read(), [4, 5])
        self.assertEqual(self.h5file.root.otherarray3.read(), [5, 6])
        self.assertEqual(self.h5file.root.otherarray4.read(), [6, 7])

    def test10_gotoint(self):
        """Checking mark sequential ids (goto)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test10_gotoint..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.create_array("/", "otherarray1", [3, 4], "Another array 1")
        self.h5file.mark("first")
        self.h5file.create_array("/", "otherarray2", [4, 5], "Another array 2")
        self.h5file.mark("second")
        self._do_reopen()
        self.h5file.create_array("/", "otherarray3", [5, 6], "Another array 3")
        self.h5file.mark("third")
        self.h5file.create_array("/", "otherarray4", [6, 7], "Another array 4")

        # Now go to mark "first"
        self.h5file.goto(1)
        self._do_reopen()
        self.assertIn("/otherarray1", self.h5file)
        self.assertNotIn("/otherarray2", self.h5file)
        self.assertNotIn("/otherarray3", self.h5file)
        self.assertNotIn("/otherarray4", self.h5file)

        # Go to beginning
        self.h5file.goto(0)
        self.assertNotIn("/otherarray1", self.h5file)
        self.assertNotIn("/otherarray2", self.h5file)
        self.assertNotIn("/otherarray3", self.h5file)
        self.assertNotIn("/otherarray4", self.h5file)

        # Go to mark "third"
        self._do_reopen()
        self.h5file.goto(3)
        self.assertIn("/otherarray1", self.h5file)
        self.assertIn("/otherarray2", self.h5file)
        self.assertIn("/otherarray3", self.h5file)
        self.assertNotIn("/otherarray4", self.h5file)

        # Now go to mark "second"
        self.h5file.goto(2)
        self.assertIn("/otherarray1", self.h5file)
        self.assertIn("/otherarray2", self.h5file)
        self.assertNotIn("/otherarray3", self.h5file)
        self.assertNotIn("/otherarray4", self.h5file)

        # Go to the end
        self._do_reopen()
        self.h5file.goto(-1)
        self.assertIn("/otherarray1", self.h5file)
        self.assertIn("/otherarray2", self.h5file)
        self.assertIn("/otherarray3", self.h5file)
        self.assertIn("/otherarray4", self.h5file)

        # Check that objects has come back to life in a sane state
        self.assertIn("/otherarray2", self.h5file)
        self.assertEqual(self.h5file.root.otherarray1.read(), [3, 4])
        self.assertEqual(self.h5file.root.otherarray2.read(), [4, 5])
        self.assertEqual(self.h5file.root.otherarray3.read(), [5, 6])
        self.assertEqual(self.h5file.root.otherarray4.read(), [6, 7])

    def test11_contiguous(self):
        """Creating contiguous marks"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test11_contiguous..." % self.__class__.__name__)

        self.h5file.enable_undo()
        m1 = self.h5file.mark()
        m2 = self.h5file.mark()
        self.assertNotEqual(m1, m2)
        self._do_reopen()
        self.h5file.undo(m1)
        self.assertEqual(self.h5file.get_current_mark(), m1)
        self.h5file.redo(m2)
        self.assertEqual(self.h5file.get_current_mark(), m2)
        self.h5file.goto(m1)
        self.assertEqual(self.h5file.get_current_mark(), m1)
        self.h5file.goto(m2)
        self.assertEqual(self.h5file.get_current_mark(), m2)
        self.h5file.goto(-1)
        self._do_reopen()
        self.assertEqual(self.h5file.get_current_mark(), m2)
        self.h5file.goto(0)
        self.assertEqual(self.h5file.get_current_mark(), 0)

    def test12_keepMark(self):
        """Ensuring the mark is kept after an UNDO operation"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test12_keepMark..." % self.__class__.__name__)

        self.h5file.enable_undo()
        self.h5file.create_array("/", "newarray1", [1])

        mid = self.h5file.mark()
        self.assertIsNotNone(mid)
        self._do_reopen()
        self.h5file.undo()

        # We should have moved to the initial mark.
        self.assertEqual(self.h5file.get_current_mark(), 0)

        # So /newarray1 should not be there.
        self.assertNotIn("/newarray1", self.h5file)

    def test13_severalEnableDisable(self):
        """Checking that successive enable/disable Undo works"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test13_severalEnableDisable..."
                % self.__class__.__name__
            )

        self.h5file.enable_undo()
        self.h5file.create_array("/", "newarray1", [1])
        self.h5file.undo()
        self._do_reopen()

        # We should have moved to 'mid' mark, not the initial mark.
        self.assertEqual(self.h5file.get_current_mark(), 0)

        # So /newarray1 should still be there.
        self.assertNotIn("/newarray1", self.h5file)

        # Close this do/undo session
        self.h5file.disable_undo()

        # Do something
        self.h5file.create_array("/", "newarray2", [1])

        # Enable again do/undo
        self.h5file.enable_undo()
        self.h5file.create_array("/", "newarray3", [1])
        mid = self.h5file.mark()
        self.h5file.create_array("/", "newarray4", [1])
        self.h5file.undo()

        # We should have moved to 'mid' mark, not the initial mark.
        self.assertEqual(self.h5file.get_current_mark(), mid)

        # So /newarray2 and /newarray3 should still be there.
        self.assertNotIn("/newarray1", self.h5file)
        self.assertIn("/newarray2", self.h5file)
        self.assertIn("/newarray3", self.h5file)
        self.assertNotIn("/newarray4", self.h5file)

        # Close this do/undo session
        self._do_reopen()
        self.h5file.disable_undo()

        # Enable again do/undo
        self.h5file.enable_undo()
        self.h5file.create_array("/", "newarray1", [1])
        self.h5file.create_array("/", "newarray4", [1])

        # So /newarray2 and /newarray3 should still be there.
        self.assertIn("/newarray1", self.h5file)
        self.assertIn("/newarray2", self.h5file)
        self.assertIn("/newarray3", self.h5file)
        self.assertIn("/newarray4", self.h5file)
        self.h5file.undo()
        self._do_reopen()
        self.assertNotIn("/newarray1", self.h5file)
        self.assertIn("/newarray2", self.h5file)
        self.assertIn("/newarray3", self.h5file)
        self.assertNotIn("/newarray4", self.h5file)

        # Close this do/undo session
        self.h5file.disable_undo()


class PersistenceTestCase(BasicTestCase):
    """Test for basic Undo/Redo operations with persistence."""

    _reopen_flag = True


class CreateArrayTestCase(common.TempFileMixin, common.PyTablesTestCase):
    """Test for create_array operations"""

    def setUp(self):
        super().setUp()

        h5file = self.h5file
        root = h5file.root

        # Create an array
        h5file.create_array(root, "array", [1, 2], title="Title example")

        # Create another array object
        h5file.create_array(root, "anarray", [1], "Array title")

        # Create a group object
        group = h5file.create_group(root, "agroup", "Group title")

        # Create a couple of objects there
        h5file.create_array(group, "anarray1", [2], "Array title 1")
        h5file.create_array(group, "anarray2", [2], "Array title 2")

        # Create a lonely group in first level
        h5file.create_group(root, "agroup2", "Group title 2")

        # Create a new group in the second level
        h5file.create_group(group, "agroup3", "Group title 3")

    def test00(self):
        """Checking one action."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test00..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.create_array("/", "otherarray1", [1, 2], "Another array 1")

        # Now undo the past operation
        self.h5file.undo()

        # Check that otherarray does not exist in the object tree
        self.assertNotIn("/otherarray1", self.h5file)

        # Redo the operation
        self.h5file.redo()

        # Check that otherarray has come back to life in a sane state
        self.assertIn("/otherarray1", self.h5file)
        self.assertEqual(self.h5file.root.otherarray1.title, "Another array 1")
        self.assertEqual(self.h5file.root.otherarray1.read(), [1, 2])

    def test01(self):
        """Checking two actions."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.create_array("/", "otherarray1", [1, 2], "Another array 1")
        self.h5file.create_array("/", "otherarray2", [2, 3], "Another array 2")

        # Now undo the past operation
        self.h5file.undo()

        # Check that otherarray does not exist in the object tree
        self.assertNotIn("/otherarray1", self.h5file)
        self.assertNotIn("/otherarray2", self.h5file)

        # Redo the operation
        self.h5file.redo()

        # Check that otherarray has come back to life in a sane state
        self.assertIn("/otherarray1", self.h5file)
        self.assertIn("/otherarray2", self.h5file)
        self.assertEqual(self.h5file.root.otherarray1.title, "Another array 1")
        self.assertEqual(self.h5file.root.otherarray2.title, "Another array 2")
        self.assertEqual(self.h5file.root.otherarray1.read(), [1, 2])
        self.assertEqual(self.h5file.root.otherarray2.read(), [2, 3])

    def test02(self):
        """Checking three actions."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.create_array("/", "otherarray1", [1, 2], "Another array 1")
        self.h5file.create_array("/", "otherarray2", [2, 3], "Another array 2")
        self.h5file.create_array("/", "otherarray3", [3, 4], "Another array 3")

        # Now undo the past operation
        self.h5file.undo()

        # Check that otherarray does not exist in the object tree
        self.assertNotIn("/otherarray1", self.h5file)
        self.assertNotIn("/otherarray2", self.h5file)
        self.assertNotIn("/otherarray3", self.h5file)

        # Redo the operation
        self.h5file.redo()

        # Check that otherarray has come back to life in a sane state
        self.assertIn("/otherarray1", self.h5file)
        self.assertIn("/otherarray2", self.h5file)
        self.assertIn("/otherarray3", self.h5file)
        self.assertEqual(self.h5file.root.otherarray1.title, "Another array 1")
        self.assertEqual(self.h5file.root.otherarray2.title, "Another array 2")
        self.assertEqual(self.h5file.root.otherarray3.title, "Another array 3")
        self.assertEqual(self.h5file.root.otherarray1.read(), [1, 2])
        self.assertEqual(self.h5file.root.otherarray2.read(), [2, 3])
        self.assertEqual(self.h5file.root.otherarray3.read(), [3, 4])

    def test03(self):
        """Checking three actions in different depth levels."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.create_array("/", "otherarray1", [1, 2], "Another array 1")
        self.h5file.create_array(
            "/agroup", "otherarray2", [2, 3], "Another array 2"
        )
        self.h5file.create_array(
            "/agroup/agroup3", "otherarray3", [3, 4], "Another array 3"
        )

        # Now undo the past operation
        self.h5file.undo()

        # Check that otherarray does not exist in the object tree
        self.assertNotIn("/otherarray1", self.h5file)
        self.assertNotIn("/agroup/otherarray2", self.h5file)
        self.assertNotIn("/agroup/agroup3/otherarray3", self.h5file)

        # Redo the operation
        self.h5file.redo()

        # Check that otherarray has come back to life in a sane state
        self.assertIn("/otherarray1", self.h5file)
        self.assertIn("/agroup/otherarray2", self.h5file)
        self.assertIn("/agroup/agroup3/otherarray3", self.h5file)
        self.assertEqual(self.h5file.root.otherarray1.title, "Another array 1")
        self.assertEqual(
            self.h5file.root.agroup.otherarray2.title, "Another array 2"
        )
        self.assertEqual(
            self.h5file.root.agroup.agroup3.otherarray3.title,
            "Another array 3",
        )
        self.assertEqual(self.h5file.root.otherarray1.read(), [1, 2])
        self.assertEqual(self.h5file.root.agroup.otherarray2.read(), [2, 3])
        self.assertEqual(
            self.h5file.root.agroup.agroup3.otherarray3.read(), [3, 4]
        )


class CreateGroupTestCase(common.TempFileMixin, common.PyTablesTestCase):
    """Test for create_group operations"""

    def setUp(self):
        super().setUp()

        h5file = self.h5file
        root = h5file.root

        # Create an array
        h5file.create_array(root, "array", [1, 2], title="Title example")

        # Create another array object
        h5file.create_array(root, "anarray", [1], "Array title")

        # Create a group object
        group = h5file.create_group(root, "agroup", "Group title")

        # Create a couple of objects there
        h5file.create_array(group, "anarray1", [2], "Array title 1")
        h5file.create_array(group, "anarray2", [2], "Array title 2")

        # Create a lonely group in first level
        h5file.create_group(root, "agroup2", "Group title 2")

        # Create a new group in the second level
        h5file.create_group(group, "agroup3", "Group title 3")

    def test00(self):
        """Checking one action."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test00..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new group
        self.h5file.create_group("/", "othergroup1", "Another group 1")

        # Now undo the past operation
        self.h5file.undo()

        # Check that othergroup1 does not exist in the object tree
        self.assertNotIn("/othergroup1", self.h5file)

        # Redo the operation
        self.h5file.redo()

        # Check that othergroup1 has come back to life in a sane state
        self.assertIn("/othergroup1", self.h5file)
        self.assertEqual(
            self.h5file.root.othergroup1._v_title, "Another group 1"
        )

    def test01(self):
        """Checking two actions."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new group
        self.h5file.create_group("/", "othergroup1", "Another group 1")
        self.h5file.create_group("/", "othergroup2", "Another group 2")

        # Now undo the past operation
        self.h5file.undo()

        # Check that othergroup does not exist in the object tree
        self.assertNotIn("/othergroup1", self.h5file)
        self.assertNotIn("/othergroup2", self.h5file)

        # Redo the operation
        self.h5file.redo()

        # Check that othergroup* has come back to life in a sane state
        self.assertIn("/othergroup1", self.h5file)
        self.assertIn("/othergroup2", self.h5file)
        self.assertEqual(
            self.h5file.root.othergroup1._v_title, "Another group 1"
        )
        self.assertEqual(
            self.h5file.root.othergroup2._v_title, "Another group 2"
        )

    def test02(self):
        """Checking three actions."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new group
        self.h5file.create_group("/", "othergroup1", "Another group 1")
        self.h5file.create_group("/", "othergroup2", "Another group 2")
        self.h5file.create_group("/", "othergroup3", "Another group 3")

        # Now undo the past operation
        self.h5file.undo()

        # Check that othergroup* does not exist in the object tree
        self.assertNotIn("/othergroup1", self.h5file)
        self.assertNotIn("/othergroup2", self.h5file)
        self.assertNotIn("/othergroup3", self.h5file)

        # Redo the operation
        self.h5file.redo()

        # Check that othergroup* has come back to life in a sane state
        self.assertIn("/othergroup1", self.h5file)
        self.assertIn("/othergroup2", self.h5file)
        self.assertIn("/othergroup3", self.h5file)
        self.assertEqual(
            self.h5file.root.othergroup1._v_title, "Another group 1"
        )
        self.assertEqual(
            self.h5file.root.othergroup2._v_title, "Another group 2"
        )
        self.assertEqual(
            self.h5file.root.othergroup3._v_title, "Another group 3"
        )

    def test03(self):
        """Checking three actions in different depth levels."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new group
        self.h5file.create_group("/", "othergroup1", "Another group 1")
        self.h5file.create_group(
            "/othergroup1", "othergroup2", "Another group 2"
        )
        self.h5file.create_group(
            "/othergroup1/othergroup2", "othergroup3", "Another group 3"
        )

        # Now undo the past operation
        self.h5file.undo()

        # Check that othergroup* does not exist in the object tree
        self.assertNotIn("/othergroup1", self.h5file)
        self.assertNotIn("/othergroup1/othergroup2", self.h5file)
        self.assertTrue(
            "/othergroup1/othergroup2/othergroup3" not in self.h5file
        )

        # Redo the operation
        self.h5file.redo()

        # Check that othergroup* has come back to life in a sane state
        self.assertIn("/othergroup1", self.h5file)
        self.assertIn("/othergroup1/othergroup2", self.h5file)
        self.assertIn("/othergroup1/othergroup2/othergroup3", self.h5file)
        self.assertEqual(
            self.h5file.root.othergroup1._v_title, "Another group 1"
        )
        self.assertEqual(
            self.h5file.root.othergroup1.othergroup2._v_title,
            "Another group 2",
        )
        self.assertEqual(
            self.h5file.root.othergroup1.othergroup2.othergroup3._v_title,
            "Another group 3",
        )


minRowIndex = 10


def populateTable(where, name):
    """Create a table under where with name name"""

    class Indexed(tb.IsDescription):
        var1 = tb.StringCol(itemsize=4, dflt=b"", pos=1)
        var2 = tb.BoolCol(dflt=0, pos=2)
        var3 = tb.IntCol(dflt=0, pos=3)
        var4 = tb.FloatCol(dflt=0, pos=4)

    nrows = minRowIndex
    table = where._v_file.create_table(
        where, name, Indexed, "Indexed", None, nrows
    )
    for i in range(nrows):
        table.row["var1"] = str(i)

        # table.row['var2'] = i > 2
        table.row["var2"] = i % 2
        table.row["var3"] = i
        table.row["var4"] = float(nrows - i - 1)
        table.row.append()
    table.flush()

    # Index all entries:
    indexrows = table.cols.var1.create_index()
    indexrows = table.cols.var2.create_index()
    indexrows = table.cols.var3.create_index()

    # Do not index the var4 column
    # indexrows = table.cols.var4.create_index()
    if common.verbose:
        print("Number of written rows:", nrows)
        print("Number of indexed rows:", table.cols.var1.index.nelements)
        print("Number of indexed rows(2):", indexrows)


class RenameNodeTestCase(common.TempFileMixin, common.PyTablesTestCase):
    """Test for rename_node operations"""

    def setUp(self):
        super().setUp()

        h5file = self.h5file
        root = h5file.root

        # Create an array
        h5file.create_array(root, "array", [1, 2], title="Title example")

        # Create another array object
        h5file.create_array(root, "anarray", [1], "Array title")

        # Create a group object
        group = h5file.create_group(root, "agroup", "Group title")

        # Create a couple of objects there
        h5file.create_array(group, "anarray1", [2], "Array title 1")
        h5file.create_array(group, "anarray2", [2], "Array title 2")

        # Create a lonely group in first level
        h5file.create_group(root, "agroup2", "Group title 2")

        # Create a new group in the second level
        h5file.create_group(group, "agroup3", "Group title 3")

        # Create a table in root
        populateTable(self.h5file.root, "table")

    def test00(self):
        """Checking rename_node (over Groups without children)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test00..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.rename_node("/agroup2", "agroup3")

        # Now undo the past operation
        self.h5file.undo()

        # Check that it does not exist in the object tree
        self.assertIn("/agroup2", self.h5file)
        self.assertNotIn("/agroup3", self.h5file)
        self.assertEqual(self.h5file.root.agroup2._v_title, "Group title 2")

        # Redo the operation
        self.h5file.redo()

        # Check that otherarray has come back to life in a sane state
        self.assertNotIn("/agroup2", self.h5file)
        self.assertIn("/agroup3", self.h5file)
        self.assertEqual(self.h5file.root.agroup3._v_title, "Group title 2")

    def test01(self):
        """Checking rename_node (over Groups with children)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.rename_node("/agroup", "agroup3")

        # Now undo the past operation
        self.h5file.undo()

        # Check that it does not exist in the object tree
        self.assertIn("/agroup", self.h5file)
        self.assertNotIn("/agroup3", self.h5file)

        # Check that children are reachable
        self.assertIn("/agroup/anarray1", self.h5file)
        self.assertIn("/agroup/anarray2", self.h5file)
        self.assertIn("/agroup/agroup3", self.h5file)
        self.assertEqual(self.h5file.root.agroup._v_title, "Group title")

        # Redo the operation
        self.h5file.redo()

        # Check that otherarray has come back to life in a sane state
        self.assertNotIn("/agroup", self.h5file)
        self.assertIn("/agroup3", self.h5file)
        self.assertEqual(self.h5file.root.agroup3._v_title, "Group title")

        # Check that children are reachable
        self.assertIn("/agroup3/anarray1", self.h5file)
        self.assertIn("/agroup3/anarray2", self.h5file)
        self.assertIn("/agroup3/agroup3", self.h5file)

    def test01b(self):
        """Checking rename_node (over Groups with children 2)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01b..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.rename_node("/agroup", "agroup3")
        self.h5file.rename_node("/agroup3", "agroup4")

        # Now undo the past operation
        self.h5file.undo()

        # Check that it does not exist in the object tree
        self.assertIn("/agroup", self.h5file)
        self.assertNotIn("/agroup4", self.h5file)

        # Check that children are reachable
        self.assertIn("/agroup/anarray1", self.h5file)
        self.assertIn("/agroup/anarray2", self.h5file)
        self.assertIn("/agroup/agroup3", self.h5file)
        self.assertEqual(self.h5file.root.agroup._v_title, "Group title")

        # Redo the operation
        self.h5file.redo()

        # Check that otherarray has come back to life in a sane state
        self.assertNotIn("/agroup", self.h5file)
        self.assertIn("/agroup4", self.h5file)
        self.assertEqual(self.h5file.root.agroup4._v_title, "Group title")

        # Check that children are reachable
        self.assertIn("/agroup4/anarray1", self.h5file)
        self.assertIn("/agroup4/anarray2", self.h5file)
        self.assertIn("/agroup4/agroup3", self.h5file)

    def test02(self):
        """Checking rename_node (over Leaves)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.rename_node("/anarray", "anarray2")

        # Now undo the past operation
        self.h5file.undo()

        # Check that otherarray does not exist in the object tree
        self.assertIn("/anarray", self.h5file)
        self.assertNotIn("/anarray2", self.h5file)
        self.assertEqual(self.h5file.root.anarray.title, "Array title")

        # Redo the operation
        self.h5file.redo()

        # Check that otherarray has come back to life in a sane state
        self.assertNotIn("/anarray", self.h5file)
        self.assertIn("/anarray2", self.h5file)
        self.assertEqual(self.h5file.root.anarray2.title, "Array title")

    def test03(self):
        """Checking rename_node (over Tables)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.rename_node("/table", "table2")

        # Now undo the past operation
        self.h5file.undo()

        # Check that table2 does not exist in the object tree
        self.assertIn("/table", self.h5file)
        table = self.h5file.root.table
        self.assertIsNotNone(table.cols.var1.index)
        self.assertIsNotNone(table.cols.var2.index)
        self.assertIsNotNone(table.cols.var3.index)
        self.assertIsNone(table.cols.var4.index)
        self.assertEqual(table.cols.var1.index.nelements, minRowIndex)
        self.assertEqual(table.cols.var2.index.nelements, minRowIndex)
        self.assertEqual(table.cols.var3.index.nelements, minRowIndex)
        self.assertNotIn("/table2", self.h5file)
        self.assertEqual(self.h5file.root.table.title, "Indexed")

        # Redo the operation
        self.h5file.redo()

        # Check that table2 has come back to life in a sane state
        self.assertNotIn("/table", self.h5file)
        self.assertIn("/table2", self.h5file)
        self.assertEqual(self.h5file.root.table2.title, "Indexed")
        table = self.h5file.root.table2
        self.assertIsNotNone(table.cols.var1.index)
        self.assertIsNotNone(table.cols.var2.index)
        self.assertIsNotNone(table.cols.var3.index)
        self.assertEqual(table.cols.var1.index.nelements, minRowIndex)
        self.assertEqual(table.cols.var2.index.nelements, minRowIndex)
        self.assertEqual(table.cols.var3.index.nelements, minRowIndex)
        self.assertIsNone(table.cols.var4.index)


class MoveNodeTestCase(common.TempFileMixin, common.PyTablesTestCase):
    """Tests for move_node operations"""

    def setUp(self):
        super().setUp()

        h5file = self.h5file
        root = h5file.root

        # Create an array
        h5file.create_array(root, "array", [1, 2], title="Title example")

        # Create another array object
        h5file.create_array(root, "anarray", [1], "Array title")

        # Create a group object
        group = h5file.create_group(root, "agroup", "Group title")

        # Create a couple of objects there
        h5file.create_array(group, "anarray1", [2], "Array title 1")
        h5file.create_array(group, "anarray2", [2], "Array title 2")

        # Create a lonely group in first level
        h5file.create_group(root, "agroup2", "Group title 2")

        # Create a new group in the second level
        h5file.create_group(group, "agroup3", "Group title 3")

        # Create a table in root
        populateTable(self.h5file.root, "table")

    def test00(self):
        """Checking move_node (over Leaf)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test00..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.move_node("/anarray", "/agroup/agroup3")

        # Now undo the past operation
        self.h5file.undo()

        # Check that it does not exist in the object tree
        self.assertIn("/anarray", self.h5file)
        self.assertNotIn("/agroup/agroup3/anarray", self.h5file)
        self.assertEqual(self.h5file.root.anarray.title, "Array title")

        # Redo the operation
        self.h5file.redo()

        # Check that otherarray has come back to life in a sane state
        self.assertNotIn("/anarray", self.h5file)
        self.assertIn("/agroup/agroup3/anarray", self.h5file)
        self.assertEqual(
            self.h5file.root.agroup.agroup3.anarray.title, "Array title"
        )

    def test01(self):
        """Checking move_node (over Groups with children)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.move_node("/agroup", "/agroup2", "agroup3")

        # Now undo the past operation
        self.h5file.undo()

        # Check that it does not exist in the object tree
        self.assertIn("/agroup", self.h5file)
        self.assertNotIn("/agroup2/agroup3", self.h5file)

        # Check that children are reachable
        self.assertIn("/agroup/anarray1", self.h5file)
        self.assertIn("/agroup/anarray2", self.h5file)
        self.assertIn("/agroup/agroup3", self.h5file)
        self.assertEqual(self.h5file.root.agroup._v_title, "Group title")

        # Redo the operation
        self.h5file.redo()

        # Check that otherarray has come back to life in a sane state
        self.assertNotIn("/agroup", self.h5file)
        self.assertIn("/agroup2/agroup3", self.h5file)
        self.assertEqual(
            self.h5file.root.agroup2.agroup3._v_title, "Group title"
        )

        # Check that children are reachable
        self.assertIn("/agroup2/agroup3/anarray1", self.h5file)
        self.assertIn("/agroup2/agroup3/anarray2", self.h5file)
        self.assertIn("/agroup2/agroup3/agroup3", self.h5file)

    def test01b(self):
        """Checking move_node (over Groups with children 2)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01b..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.move_node("/agroup", "/", "agroup3")
        self.h5file.move_node("/agroup3", "/agroup2", "agroup4")

        # Now undo the past operation
        self.h5file.undo()

        # Check that it does not exist in the object tree
        self.assertIn("/agroup", self.h5file)
        self.assertNotIn("/agroup2/agroup4", self.h5file)

        # Check that children are reachable
        self.assertIn("/agroup/anarray1", self.h5file)
        self.assertIn("/agroup/anarray2", self.h5file)
        self.assertIn("/agroup/agroup3", self.h5file)
        self.assertEqual(self.h5file.root.agroup._v_title, "Group title")

        # Redo the operation
        self.h5file.redo()

        # Check that otherarray has come back to life in a sane state
        self.assertNotIn("/agroup", self.h5file)
        self.assertIn("/agroup2/agroup4", self.h5file)
        self.assertEqual(
            self.h5file.root.agroup2.agroup4._v_title, "Group title"
        )

        # Check that children are reachable
        self.assertIn("/agroup2/agroup4/anarray1", self.h5file)
        self.assertIn("/agroup2/agroup4/anarray2", self.h5file)
        self.assertIn("/agroup2/agroup4/agroup3", self.h5file)

    def test02(self):
        """Checking move_node (over Leaves)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.move_node("/anarray", "/agroup2", "anarray2")

        # Now undo the past operation
        self.h5file.undo()

        # Check that otherarray does not exist in the object tree
        self.assertIn("/anarray", self.h5file)
        self.assertNotIn("/agroup2/anarray2", self.h5file)
        self.assertEqual(self.h5file.root.anarray.title, "Array title")

        # Redo the operation
        self.h5file.redo()

        # Check that otherarray has come back to life in a sane state
        self.assertNotIn("/anarray", self.h5file)
        self.assertIn("/agroup2/anarray2", self.h5file)
        self.assertEqual(
            self.h5file.root.agroup2.anarray2.title, "Array title"
        )

    def test03(self):
        """Checking move_node (over Tables)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.move_node("/table", "/agroup2", "table2")

        # Now undo the past operation
        self.h5file.undo()

        # Check that table2 does not exist in the object tree
        self.assertIn("/table", self.h5file)
        self.assertNotIn("/agroup2/table2", self.h5file)
        table = self.h5file.root.table
        self.assertIsNotNone(table.cols.var1.index)
        self.assertIsNotNone(table.cols.var2.index)
        self.assertIsNotNone(table.cols.var3.index)
        self.assertIsNone(table.cols.var4.index)
        self.assertEqual(table.cols.var1.index.nelements, minRowIndex)
        self.assertEqual(table.cols.var2.index.nelements, minRowIndex)
        self.assertEqual(table.cols.var3.index.nelements, minRowIndex)
        self.assertEqual(self.h5file.root.table.title, "Indexed")

        # Redo the operation
        self.h5file.redo()

        # Check that table2 has come back to life in a sane state
        self.assertNotIn("/table", self.h5file)
        self.assertIn("/agroup2/table2", self.h5file)
        self.assertEqual(self.h5file.root.agroup2.table2.title, "Indexed")
        table = self.h5file.root.agroup2.table2
        self.assertIsNotNone(table.cols.var1.index)
        self.assertIsNotNone(table.cols.var2.index)
        self.assertIsNotNone(table.cols.var3.index)
        self.assertEqual(table.cols.var1.index.nelements, minRowIndex)
        self.assertEqual(table.cols.var2.index.nelements, minRowIndex)
        self.assertEqual(table.cols.var3.index.nelements, minRowIndex)
        self.assertIsNone(table.cols.var4.index)


class RemoveNodeTestCase(common.TempFileMixin, common.PyTablesTestCase):
    """Test for remove_node operations"""

    def setUp(self):
        super().setUp()

        h5file = self.h5file
        root = h5file.root

        # Create an array
        h5file.create_array(root, "array", [1, 2], title="Title example")

        # Create another array object
        h5file.create_array(root, "anarray", [1], "Array title")

        # Create a group object
        group = h5file.create_group(root, "agroup", "Group title")

        # Create a couple of objects there
        h5file.create_array(group, "anarray1", [2], "Array title 1")
        h5file.create_array(group, "anarray2", [2], "Array title 2")

        # Create a lonely group in first level
        h5file.create_group(root, "agroup2", "Group title 2")

        # Create a new group in the second level
        h5file.create_group(group, "agroup3", "Group title 3")

        # Create a table in root
        populateTable(self.h5file.root, "table")

    def test00(self):
        """Checking remove_node (over Leaf)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test00..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Delete an existing array
        self.h5file.remove_node("/anarray")

        # Now undo the past operation
        self.h5file.undo()

        # Check that it does exist in the object tree
        self.assertIn("/anarray", self.h5file)
        self.assertEqual(self.h5file.root.anarray.title, "Array title")

        # Redo the operation
        self.h5file.redo()

        # Check that array has gone again
        self.assertNotIn("/anarray", self.h5file)

    def test00b(self):
        """Checking remove_node (over several Leaves)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test00b..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Delete a couple of arrays
        self.h5file.remove_node("/anarray")
        self.h5file.remove_node("/agroup/anarray2")

        # Now undo the past operation
        self.h5file.undo()

        # Check that arrays has come into life
        self.assertIn("/anarray", self.h5file)
        self.assertIn("/agroup/anarray2", self.h5file)
        self.assertEqual(self.h5file.root.anarray.title, "Array title")
        self.assertEqual(
            self.h5file.root.agroup.anarray2.title, "Array title 2"
        )

        # Redo the operation
        self.h5file.redo()

        # Check that arrays has disappeared again
        self.assertNotIn("/anarray", self.h5file)
        self.assertNotIn("/agroup/anarray2", self.h5file)

    def test00c(self):
        """Checking remove_node (over Tables)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test00c..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Create a new array
        self.h5file.remove_node("/table")

        # Now undo the past operation
        self.h5file.undo()

        # Check that table2 does not exist in the object tree
        self.assertIn("/table", self.h5file)
        table = self.h5file.root.table
        self.assertIsNotNone(table.cols.var1.index)
        self.assertIsNotNone(table.cols.var2.index)
        self.assertIsNotNone(table.cols.var3.index)
        self.assertIsNone(table.cols.var4.index)
        self.assertEqual(table.cols.var1.index.nelements, minRowIndex)
        self.assertEqual(table.cols.var2.index.nelements, minRowIndex)
        self.assertEqual(table.cols.var3.index.nelements, minRowIndex)
        self.assertEqual(self.h5file.root.table.title, "Indexed")

        # Redo the operation
        self.h5file.redo()

        # Check that table2 has come back to life in a sane state
        self.assertNotIn("/table", self.h5file)

    def test01(self):
        """Checking remove_node (over Groups with children)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Delete a group recursively
        self.h5file.remove_node("/agroup", recursive=1)

        # Now undo the past operation
        self.h5file.undo()

        # Check that parent and children has come into life in a sane state
        self.assertIn("/agroup", self.h5file)
        self.assertIn("/agroup/anarray1", self.h5file)
        self.assertIn("/agroup/anarray2", self.h5file)
        self.assertIn("/agroup/agroup3", self.h5file)
        self.assertEqual(self.h5file.root.agroup._v_title, "Group title")

        # Redo the operation
        self.h5file.redo()

        # Check that parent and children are not reachable
        self.assertNotIn("/agroup", self.h5file)
        self.assertNotIn("/agroup/anarray1", self.h5file)
        self.assertNotIn("/agroup/anarray2", self.h5file)
        self.assertNotIn("/agroup/agroup3", self.h5file)

    def test01b(self):
        """Checking remove_node (over Groups with children 2)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01b..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # Remove a couple of groups
        self.h5file.remove_node("/agroup", recursive=1)
        self.h5file.remove_node("/agroup2")

        # Now undo the past operation
        self.h5file.undo()

        # Check that they does exist in the object tree
        self.assertIn("/agroup", self.h5file)
        self.assertIn("/agroup2", self.h5file)

        # Check that children are reachable
        self.assertIn("/agroup/anarray1", self.h5file)
        self.assertIn("/agroup/anarray2", self.h5file)
        self.assertIn("/agroup/agroup3", self.h5file)
        self.assertEqual(self.h5file.root.agroup._v_title, "Group title")

        # Redo the operation
        self.h5file.redo()

        # Check that groups does not exist again
        self.assertNotIn("/agroup", self.h5file)
        self.assertNotIn("/agroup2", self.h5file)

        # Check that children are not reachable
        self.assertNotIn("/agroup/anarray1", self.h5file)
        self.assertNotIn("/agroup/anarray2", self.h5file)
        self.assertNotIn("/agroup/agroup3", self.h5file)


class CopyNodeTestCase(common.TempFileMixin, common.PyTablesTestCase):
    """Tests for copy_node and copy_children operations"""

    def setUp(self):
        super().setUp()

        h5file = self.h5file
        root = h5file.root

        # Create an array
        h5file.create_array(root, "array", [1, 2], title="Title example")

        # Create another array object
        h5file.create_array(root, "anarray", [1], "Array title")

        # Create a group object
        group = h5file.create_group(root, "agroup", "Group title")

        # Create a couple of objects there
        h5file.create_array(group, "anarray1", [2], "Array title 1")
        h5file.create_array(group, "anarray2", [2], "Array title 2")

        # Create a lonely group in first level
        h5file.create_group(root, "agroup2", "Group title 2")

        # Create a new group in the second level
        h5file.create_group(group, "agroup3", "Group title 3")

        # Create a table in root
        populateTable(self.h5file.root, "table")

    def test00_copyLeaf(self):
        """Checking copy_node (over Leaves)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test00_copyLeaf..." % self.__class__.__name__)

        # Enable undo/redo.
        self.h5file.enable_undo()

        # /anarray => /agroup/agroup3/
        new_node = self.h5file.copy_node("/anarray", "/agroup/agroup3")

        # Undo the copy.
        self.h5file.undo()

        # Check that the copied node does not exist in the object tree.
        self.assertNotIn("/agroup/agroup3/anarray", self.h5file)

        # Redo the copy.
        self.h5file.redo()

        # Check that the copied node exists again in the object tree.
        self.assertIn("/agroup/agroup3/anarray", self.h5file)
        self.assertIs(self.h5file.root.agroup.agroup3.anarray, new_node)

    def test00b_copyTable(self):
        """Checking copy_node (over Tables)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test00b_copyTable..." % self.__class__.__name__)

        # open the do/undo
        self.h5file.enable_undo()

        # /table => /agroup/agroup3/
        warnings.filterwarnings("ignore", category=UserWarning)
        table = self.h5file.copy_node(
            "/table", "/agroup/agroup3", propindexes=True
        )
        warnings.filterwarnings("default", category=UserWarning)
        self.assertIn("/agroup/agroup3/table", self.h5file)

        table = self.h5file.root.agroup.agroup3.table
        self.assertEqual(table.title, "Indexed")
        self.assertIsNotNone(table.cols.var1.index)
        self.assertIsNotNone(table.cols.var2.index)
        self.assertIsNotNone(table.cols.var3.index)
        self.assertEqual(table.cols.var1.index.nelements, minRowIndex)
        self.assertEqual(table.cols.var2.index.nelements, minRowIndex)
        self.assertEqual(table.cols.var3.index.nelements, minRowIndex)
        self.assertIsNone(table.cols.var4.index)

        # Now undo the past operation
        self.h5file.undo()
        table = self.h5file.root.table
        self.assertIsNotNone(table.cols.var1.index)
        self.assertIsNotNone(table.cols.var2.index)
        self.assertIsNotNone(table.cols.var3.index)
        self.assertIsNone(table.cols.var4.index)
        self.assertEqual(table.cols.var1.index.nelements, minRowIndex)
        self.assertEqual(table.cols.var2.index.nelements, minRowIndex)
        self.assertEqual(table.cols.var3.index.nelements, minRowIndex)

        # Check that the copied node does not exist in the object tree.
        self.assertNotIn("/agroup/agroup3/table", self.h5file)

        # Redo the operation
        self.h5file.redo()

        # Check that table has come back to life in a sane state
        self.assertIn("/table", self.h5file)
        self.assertIn("/agroup/agroup3/table", self.h5file)
        table = self.h5file.root.agroup.agroup3.table
        self.assertEqual(table.title, "Indexed")
        self.assertIsNotNone(table.cols.var1.index)
        self.assertIsNotNone(table.cols.var2.index)
        self.assertIsNotNone(table.cols.var3.index)
        self.assertEqual(table.cols.var1.index.nelements, minRowIndex)
        self.assertEqual(table.cols.var2.index.nelements, minRowIndex)
        self.assertEqual(table.cols.var3.index.nelements, minRowIndex)
        self.assertIsNone(table.cols.var4.index)

    def test01_copyGroup(self):
        """Copying a group (recursively)."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01_copyGroup..." % self.__class__.__name__)

        # Enable undo/redo.
        self.h5file.enable_undo()

        # /agroup => /acopy
        new_node = self.h5file.copy_node(
            "/agroup", newname="acopy", recursive=True
        )

        # Undo the copy.
        self.h5file.undo()

        # Check that the copied node does not exist in the object tree.
        self.assertNotIn("/acopy", self.h5file)
        self.assertNotIn("/acopy/anarray1", self.h5file)
        self.assertNotIn("/acopy/anarray2", self.h5file)
        self.assertNotIn("/acopy/agroup3", self.h5file)

        # Redo the copy.
        self.h5file.redo()

        # Check that the copied node exists again in the object tree.
        self.assertIn("/acopy", self.h5file)
        self.assertIn("/acopy/anarray1", self.h5file)
        self.assertIn("/acopy/anarray2", self.h5file)
        self.assertIn("/acopy/agroup3", self.h5file)
        self.assertIs(self.h5file.root.acopy, new_node)

    def test02_copyLeafOverwrite(self):
        """Copying a leaf, overwriting destination."""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test02_copyLeafOverwrite..."
                % self.__class__.__name__
            )

        # Enable undo/redo.
        self.h5file.enable_undo()

        # /anarray => /agroup/agroup
        oldNode = self.h5file.root.agroup
        new_node = self.h5file.copy_node(
            "/anarray", newname="agroup", overwrite=True
        )

        # Undo the copy.
        self.h5file.undo()

        # Check that the copied node does not exist in the object tree.
        # Check that the overwritten node exists again in the object tree.
        self.assertIs(self.h5file.root.agroup, oldNode)

        # Redo the copy.
        self.h5file.redo()

        # Check that the copied node exists again in the object tree.
        # Check that the overwritten node does not exist in the object tree.
        self.assertIs(self.h5file.root.agroup, new_node)

    def test03_copyChildren(self):
        """Copying the children of a group"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test03_copyChildren..." % self.__class__.__name__
            )

        # Enable undo/redo.
        self.h5file.enable_undo()

        # /agroup/* => /agroup/
        self.h5file.copy_children("/agroup", "/agroup2", recursive=True)

        # Undo the copy.
        self.h5file.undo()

        # Check that the copied nodes do not exist in the object tree.
        self.assertNotIn("/agroup2/anarray1", self.h5file)
        self.assertNotIn("/agroup2/anarray2", self.h5file)
        self.assertNotIn("/agroup2/agroup3", self.h5file)

        # Redo the copy.
        self.h5file.redo()

        # Check that the copied nodes exist again in the object tree.
        self.assertIn("/agroup2/anarray1", self.h5file)
        self.assertIn("/agroup2/anarray2", self.h5file)
        self.assertIn("/agroup2/agroup3", self.h5file)


class ComplexTestCase(common.TempFileMixin, common.PyTablesTestCase):
    """Tests for a mix of all operations"""

    def setUp(self):
        super().setUp()

        h5file = self.h5file
        root = h5file.root

        # Create an array
        h5file.create_array(root, "array", [1, 2], title="Title example")

        # Create another array object
        h5file.create_array(root, "anarray", [1], "Array title")

        # Create a group object
        group = h5file.create_group(root, "agroup", "Group title")

        # Create a couple of objects there
        h5file.create_array(group, "anarray1", [2], "Array title 1")
        h5file.create_array(group, "anarray2", [2], "Array title 2")

        # Create a lonely group in first level
        h5file.create_group(root, "agroup2", "Group title 2")

        # Create a new group in the second level
        h5file.create_group(group, "agroup3", "Group title 3")

    def test00(self):
        """Mix of create_array, create_group, renameNone, move_node,
        remove_node, copy_node and copy_children."""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test00..." % self.__class__.__name__)

        # Enable undo/redo.
        self.h5file.enable_undo()

        # Create an array
        self.h5file.create_array(
            self.h5file.root, "anarray3", [1], "Array title 3"
        )
        # Create a group
        self.h5file.create_group(self.h5file.root, "agroup3", "Group title 3")

        # /anarray => /agroup/agroup3/
        new_node = self.h5file.copy_node("/anarray3", "/agroup/agroup3")
        new_node = self.h5file.copy_children(
            "/agroup", "/agroup3", recursive=1
        )

        # rename anarray
        self.h5file.rename_node("/anarray", "anarray4")

        # Move anarray
        new_node = self.h5file.copy_node("/anarray3", "/agroup")

        # Remove anarray4
        self.h5file.remove_node("/anarray4")

        # Undo the actions
        self.h5file.undo()
        self.assertNotIn("/anarray4", self.h5file)
        self.assertNotIn("/anarray3", self.h5file)
        self.assertNotIn("/agroup/agroup3/anarray3", self.h5file)
        self.assertNotIn("/agroup3", self.h5file)
        self.assertNotIn("/anarray4", self.h5file)
        self.assertIn("/anarray", self.h5file)

        # Redo the actions
        self.h5file.redo()

        # Check that the copied node exists again in the object tree.
        self.assertIn("/agroup/agroup3/anarray3", self.h5file)
        self.assertIn("/agroup/anarray3", self.h5file)
        self.assertIn("/agroup3/agroup3/anarray3", self.h5file)
        self.assertNotIn("/agroup3/anarray3", self.h5file)
        self.assertIs(self.h5file.root.agroup.anarray3, new_node)
        self.assertNotIn("/anarray", self.h5file)
        self.assertNotIn("/anarray4", self.h5file)

    def test01(self):
        """Test with multiple generations (Leaf case)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test01..." % self.__class__.__name__)

        # Enable undo/redo.
        self.h5file.enable_undo()

        # remove /anarray
        self.h5file.remove_node("/anarray")

        # Create an array in the same place
        self.h5file.create_array(
            self.h5file.root, "anarray", [2], "Array title 2"
        )
        # remove the array again
        self.h5file.remove_node("/anarray")

        # Create an array
        self.h5file.create_array(
            self.h5file.root, "anarray", [3], "Array title 3"
        )
        # remove the array again
        self.h5file.remove_node("/anarray")

        # Create an array
        self.h5file.create_array(
            self.h5file.root, "anarray", [4], "Array title 4"
        )
        # Undo the actions
        self.h5file.undo()

        # Check that /anarray is in the correct state before redoing
        self.assertEqual(self.h5file.root.anarray.title, "Array title")
        self.assertEqual(self.h5file.root.anarray[:], [1])

        # Redo the actions
        self.h5file.redo()
        self.assertEqual(self.h5file.root.anarray.title, "Array title 4")
        self.assertEqual(self.h5file.root.anarray[:], [4])

    def test02(self):
        """Test with multiple generations (Group case)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02..." % self.__class__.__name__)

        # Enable undo/redo.
        self.h5file.enable_undo()

        # remove /agroup
        self.h5file.remove_node("/agroup2")

        # Create a group in the same place
        self.h5file.create_group(self.h5file.root, "agroup2", "Group title 22")

        # remove the group
        self.h5file.remove_node("/agroup2")

        # Create a group
        self.h5file.create_group(self.h5file.root, "agroup2", "Group title 3")

        # remove the group
        self.h5file.remove_node("/agroup2")

        # Create a group
        self.h5file.create_group(self.h5file.root, "agroup2", "Group title 4")

        # Create a child group
        self.h5file.create_group(
            self.h5file.root.agroup2, "agroup5", "Group title 5"
        )

        # Undo the actions
        self.h5file.undo()

        # Check that /agroup is in the state before enabling do/undo
        self.assertEqual(self.h5file.root.agroup2._v_title, "Group title 2")
        self.assertIn("/agroup2", self.h5file)

        # Redo the actions
        self.h5file.redo()
        self.assertEqual(self.h5file.root.agroup2._v_title, "Group title 4")
        self.assertEqual(
            self.h5file.root.agroup2.agroup5._v_title, "Group title 5"
        )

    def test03(self):
        """Test with multiple generations (Group case, recursive remove)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03..." % self.__class__.__name__)

        # Enable undo/redo.
        self.h5file.enable_undo()

        # remove /agroup
        self.h5file.remove_node("/agroup", recursive=1)

        # Create a group in the same place
        self.h5file.create_group(self.h5file.root, "agroup", "Group title 2")

        # remove the group
        self.h5file.remove_node("/agroup")

        # Create a group
        self.h5file.create_group(self.h5file.root, "agroup", "Group title 3")

        # remove the group
        self.h5file.remove_node("/agroup")

        # Create a group
        self.h5file.create_group(self.h5file.root, "agroup", "Group title 4")

        # Create a child group
        self.h5file.create_group(
            self.h5file.root.agroup, "agroup5", "Group title 5"
        )
        # Undo the actions
        self.h5file.undo()

        # Check that /agroup is in the state before enabling do/undo
        self.assertIn("/agroup", self.h5file)
        self.assertEqual(self.h5file.root.agroup._v_title, "Group title")
        self.assertIn("/agroup/anarray1", self.h5file)
        self.assertIn("/agroup/anarray2", self.h5file)
        self.assertIn("/agroup/agroup3", self.h5file)
        self.assertNotIn("/agroup/agroup5", self.h5file)

        # Redo the actions
        self.h5file.redo()
        self.assertIn("/agroup", self.h5file)
        self.assertEqual(self.h5file.root.agroup._v_title, "Group title 4")
        self.assertIn("/agroup/agroup5", self.h5file)
        self.assertEqual(
            self.h5file.root.agroup.agroup5._v_title, "Group title 5"
        )

    def test03b(self):
        """Test with multiple generations (Group case, recursive remove,
        case 2)"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test03b..." % self.__class__.__name__)

        # Enable undo/redo.
        self.h5file.enable_undo()

        # Create a new group with a child
        self.h5file.create_group(self.h5file.root, "agroup3", "Group title 3")
        self.h5file.create_group(
            self.h5file.root.agroup3, "agroup4", "Group title 4"
        )

        # remove /agroup3
        self.h5file.remove_node("/agroup3", recursive=1)

        # Create a group in the same place
        self.h5file.create_group(self.h5file.root, "agroup3", "Group title 4")

        # Undo the actions
        self.h5file.undo()

        # Check that /agroup is in the state before enabling do/undo
        self.assertNotIn("/agroup3", self.h5file)

        # Redo the actions
        self.h5file.redo()
        self.assertEqual(self.h5file.root.agroup3._v_title, "Group title 4")
        self.assertIn("/agroup3", self.h5file)
        self.assertNotIn("/agroup/agroup4", self.h5file)


class AttributesTestCase(common.TempFileMixin, common.PyTablesTestCase):
    """Tests for operation on attributes"""

    def setUp(self):
        super().setUp()

        # Create an array.
        array = self.h5file.create_array("/", "array", [1, 2])

        # Set some attributes on it.
        attrs = array.attrs
        attrs.attr_1 = 10
        attrs.attr_2 = 20
        attrs.attr_3 = 30

    def test00_setAttr(self):
        """Setting a nonexistent attribute"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test00_setAttr..." % self.__class__.__name__)

        array = self.h5file.root.array
        attrs = array.attrs

        self.h5file.enable_undo()
        setattr(attrs, "attr_0", 0)
        self.assertIn("attr_0", attrs)
        self.assertEqual(attrs.attr_0, 0)
        self.h5file.undo()
        self.assertNotIn("attr_0", attrs)
        self.h5file.redo()
        self.assertIn("attr_0", attrs)
        self.assertEqual(attrs.attr_0, 0)

    def test01_setAttrExisting(self):
        """Setting an existing attribute"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test01_setAttrExisting..."
                % self.__class__.__name__
            )

        array = self.h5file.root.array
        attrs = array.attrs

        self.h5file.enable_undo()
        setattr(attrs, "attr_1", 11)
        self.assertIn("attr_1", attrs)
        self.assertEqual(attrs.attr_1, 11)
        self.h5file.undo()
        self.assertIn("attr_1", attrs)
        self.assertEqual(attrs.attr_1, 10)
        self.h5file.redo()
        self.assertIn("attr_1", attrs)
        self.assertEqual(attrs.attr_1, 11)

    def test02_delAttr(self):
        """Removing an attribute"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test02_delAttr..." % self.__class__.__name__)

        array = self.h5file.root.array
        attrs = array.attrs

        self.h5file.enable_undo()
        delattr(attrs, "attr_1")
        self.assertNotIn("attr_1", attrs)
        self.h5file.undo()
        self.assertIn("attr_1", attrs)
        self.assertEqual(attrs.attr_1, 10)
        self.h5file.redo()
        self.assertNotIn("attr_1", attrs)

    def test03_copyNodeAttrs(self):
        """Copying an attribute set"""

        if common.verbose:
            print("\n", "-=" * 30)
            print(
                "Running %s.test03_copyNodeAttrs..." % self.__class__.__name__
            )

        rattrs = self.h5file.root._v_attrs
        rattrs.attr_0 = 0
        rattrs.attr_1 = 100

        array = self.h5file.root.array
        attrs = array.attrs

        self.h5file.enable_undo()
        attrs._f_copy(self.h5file.root)
        self.assertEqual(rattrs.attr_0, 0)
        self.assertEqual(rattrs.attr_1, 10)
        self.assertEqual(rattrs.attr_2, 20)
        self.assertEqual(rattrs.attr_3, 30)
        self.h5file.undo()
        self.assertEqual(rattrs.attr_0, 0)
        self.assertEqual(rattrs.attr_1, 100)
        self.assertNotIn("attr_2", rattrs)
        self.assertNotIn("attr_3", rattrs)
        self.h5file.redo()
        self.assertEqual(rattrs.attr_0, 0)
        self.assertEqual(rattrs.attr_1, 10)
        self.assertEqual(rattrs.attr_2, 20)
        self.assertEqual(rattrs.attr_3, 30)

    def test04_replaceNode(self):
        """Replacing a node with a rewritten attribute"""

        if common.verbose:
            print("\n", "-=" * 30)
            print("Running %s.test04_replaceNode..." % self.__class__.__name__)

        array = self.h5file.root.array
        attrs = array.attrs

        self.h5file.enable_undo()
        attrs.attr_1 = 11
        self.h5file.remove_node("/array")
        arr = self.h5file.create_array("/", "array", [1])
        arr.attrs.attr_1 = 12
        self.h5file.undo()
        self.assertIn("attr_1", self.h5file.root.array.attrs)
        self.assertEqual(self.h5file.root.array.attrs.attr_1, 10)
        self.h5file.redo()
        self.assertIn("attr_1", self.h5file.root.array.attrs)
        self.assertEqual(self.h5file.root.array.attrs.attr_1, 12)


class NotLoggedTestCase(common.TempFileMixin, common.PyTablesTestCase):
    """Test not logged nodes."""

    class NotLoggedArray(tb.node.NotLoggedMixin, tb.Array):
        pass

    def test00_hierarchy(self):
        """Performing hierarchy operations on a not logged node."""

        self.h5file.create_group("/", "tgroup")
        self.h5file.enable_undo()

        # Node creation is not undone.
        arr = self.NotLoggedArray(
            self.h5file.root, "test", [1], self._getMethodName()
        )
        self.h5file.undo()
        self.assertIn("/test", self.h5file)

        # Node movement is not undone.
        arr.move("/tgroup")
        self.h5file.undo()
        self.assertIn("/tgroup/test", self.h5file)

        # Node removal is not undone.
        arr.remove()
        self.h5file.undo()
        self.assertNotIn("/tgroup/test", self.h5file)

    def test01_attributes(self):
        """Performing attribute operations on a not logged node."""

        arr = self.NotLoggedArray(
            self.h5file.root, "test", [1], self._getMethodName()
        )
        self.h5file.enable_undo()

        # Attribute creation is not undone.
        arr._v_attrs.foo = "bar"
        self.h5file.undo()
        self.assertEqual(arr._v_attrs.foo, "bar")

        # Attribute change is not undone.
        arr._v_attrs.foo = "baz"
        self.h5file.undo()
        self.assertEqual(arr._v_attrs.foo, "baz")

        # Attribute removal is not undone.
        del arr._v_attrs.foo
        self.h5file.undo()
        self.assertRaises(AttributeError, getattr, arr._v_attrs, "foo")


class CreateParentsTestCase(common.TempFileMixin, common.PyTablesTestCase):
    """Test the ``createparents`` flag."""

    def setUp(self):
        super().setUp()
        g1 = self.h5file.create_group("/", "g1")
        self.h5file.create_group(g1, "g2")

    def existing(self, paths):
        """Return a set of the existing paths in `paths`."""
        return frozenset(path for path in paths if path in self.h5file)

    def basetest(self, doit, pre, post):
        pre()
        self.h5file.enable_undo()

        paths = ["/g1", "/g1/g2", "/g1/g2/g3", "/g1/g2/g3/g4"]
        for newpath in paths:
            before = self.existing(paths)
            doit(newpath)
            after = self.existing(paths)
            self.assertTrue(after.issuperset(before))

            self.h5file.undo()
            post(newpath)
            after = self.existing(paths)
            self.assertEqual(after, before)

    def test00_create(self):
        """Test creating a node."""

        def pre():
            pass

        def doit(newpath):
            self.h5file.create_array(newpath, "array", [1], createparents=True)
            self.assertIn(tb.path.join_path(newpath, "array"), self.h5file)

        def post(newpath):
            self.assertNotIn(tb.path.join_path(newpath, "array"), self.h5file)

        self.basetest(doit, pre, post)

    def test01_move(self):
        """Test moving a node."""

        def pre():
            self.h5file.create_array("/", "array", [1])

        def doit(newpath):
            self.h5file.move_node("/array", newpath, createparents=True)
            self.assertNotIn("/array", self.h5file)
            self.assertIn(tb.path.join_path(newpath, "array"), self.h5file)

        def post(newpath):
            self.assertIn("/array", self.h5file)
            self.assertNotIn(tb.path.join_path(newpath, "array"), self.h5file)

        self.basetest(doit, pre, post)

    def test02_copy(self):
        """Test copying a node."""

        def pre():
            self.h5file.create_array("/", "array", [1])

        def doit(newpath):
            self.h5file.copy_node("/array", newpath, createparents=True)
            self.assertIn(tb.path.join_path(newpath, "array"), self.h5file)

        def post(newpath):
            self.assertNotIn(tb.path.join_path(newpath, "array"), self.h5file)

        self.basetest(doit, pre, post)

    def test03_copyChildren(self):
        """Test copying the children of a group."""

        def pre():
            g = self.h5file.create_group("/", "group")
            self.h5file.create_array(g, "array1", [1])
            self.h5file.create_array(g, "array2", [1])

        def doit(newpath):
            self.h5file.copy_children("/group", newpath, createparents=True)
            self.assertIn(tb.path.join_path(newpath, "array1"), self.h5file)
            self.assertIn(tb.path.join_path(newpath, "array2"), self.h5file)

        def post(newpath):
            self.assertNotIn(tb.path.join_path(newpath, "array1"), self.h5file)
            self.assertNotIn(tb.path.join_path(newpath, "array2"), self.h5file)

        self.basetest(doit, pre, post)


def suite():
    theSuite = common.unittest.TestSuite()
    niter = 1
    # common.heavy = 1  # uncomment this only for testing purposes

    for n in range(niter):
        theSuite.addTest(common.make_suite(BasicTestCase))
        theSuite.addTest(common.make_suite(PersistenceTestCase))
        theSuite.addTest(common.make_suite(CreateArrayTestCase))
        theSuite.addTest(common.make_suite(CreateGroupTestCase))
        theSuite.addTest(common.make_suite(RenameNodeTestCase))
        theSuite.addTest(common.make_suite(MoveNodeTestCase))
        theSuite.addTest(common.make_suite(RemoveNodeTestCase))
        theSuite.addTest(common.make_suite(CopyNodeTestCase))
        theSuite.addTest(common.make_suite(AttributesTestCase))
        theSuite.addTest(common.make_suite(ComplexTestCase))
        theSuite.addTest(common.make_suite(NotLoggedTestCase))
        theSuite.addTest(common.make_suite(CreateParentsTestCase))
    if common.heavy:
        pass

    return theSuite


if __name__ == "__main__":
    import sys

    common.parse_argv(sys.argv)
    common.print_versions()
    common.unittest.main(defaultTest="suite")
