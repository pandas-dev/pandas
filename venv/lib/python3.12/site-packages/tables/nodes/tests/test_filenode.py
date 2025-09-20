"""Unit test for the filenode module."""

import os
import shutil
import tempfile
import warnings
from pathlib import Path

from ... import open_file, file, NoSuchNodeError
from ...nodes import filenode
from ...tests.common import (
    unittest,
    TempFileMixin,
    parse_argv,
    print_versions,
    make_suite,
)
from ...tests.common import PyTablesTestCase as TestCase


def test_file(name):
    from importlib import resources

    return resources.files("tables.nodes.tests") / name


class NewFileTestCase(TempFileMixin, TestCase):
    """Tests creating a new file node with the new_node() function."""

    def test00_NewFile(self):
        """Creation of a brand new file node."""

        try:
            fnode = filenode.new_node(self.h5file, where="/", name="test")
            node = self.h5file.get_node("/test")
        except LookupError:
            self.fail("filenode.new_node() failed to create a new node.")
        else:
            self.assertEqual(
                fnode.node,
                node,
                "filenode.new_node() created a node in the wrong place.",
            )

    def test01_NewFileTooFewArgs(self):
        """Creation of a new file node without arguments for node creation."""

        self.assertRaises(TypeError, filenode.new_node, self.h5file)

    def test02_NewFileWithExpectedSize(self):
        """Creation of a new file node with 'expectedsize' argument."""

        try:
            filenode.new_node(
                self.h5file, where="/", name="test", expectedsize=100_000
            )
        except TypeError:
            self.fail(
                "filenode.new_node() failed to accept 'expectedsize' argument."
            )

    def test03_NewFileWithExpectedRows(self):
        """Creation of a new file node with illegal 'expectedrows' argument."""

        self.assertRaises(
            TypeError,
            filenode.new_node,
            self.h5file,
            where="/",
            name="test",
            expectedrows=100_000,
        )


class ClosedFileTestCase(TempFileMixin, TestCase):
    """Tests calling several methods on a closed file."""

    def setUp(self):
        """setUp() -> None

        This method sets the following instance attributes:
          * 'h5fname', the name of the temporary HDF5 file
          * 'h5file', the writable, temporary HDF5 file with a '/test' node
          * 'fnode', the closed file node in '/test'

        """

        super().setUp()
        self.fnode = filenode.new_node(self.h5file, where="/", name="test")
        self.fnode.close()

    def tearDown(self):
        """tearDown() -> None

        Closes 'h5file'; removes 'h5fname'.

        """

        self.fnode = None
        super().tearDown()

    # All these tests mey seem odd, but Python (2.3) files
    # do test whether the file is not closed regardless of their mode.
    def test00_Close(self):
        """Closing a closed file."""

        try:
            self.fnode.close()
        except ValueError:
            self.fail("Could not close an already closed file.")

    def test01_Flush(self):
        """Flushing a closed file."""

        self.assertRaises(ValueError, self.fnode.flush)

    def test02_Next(self):
        """Getting the next line of a closed file."""

        self.assertRaises(ValueError, next, self.fnode)

    def test03_Read(self):
        """Reading a closed file."""

        self.assertRaises(ValueError, self.fnode.read)

    def test04_Readline(self):
        """Reading a line from a closed file."""

        self.assertRaises(ValueError, self.fnode.readline)

    def test05_Readlines(self):
        """Reading lines from a closed file."""

        self.assertRaises(ValueError, self.fnode.readlines)

    def test06_Seek(self):
        """Seeking a closed file."""

        self.assertRaises(ValueError, self.fnode.seek, 0)

    def test07_Tell(self):
        """Getting the pointer position in a closed file."""

        self.assertRaises(ValueError, self.fnode.tell)

    def test08_Truncate(self):
        """Truncating a closed file."""

        self.assertRaises(ValueError, self.fnode.truncate)

    def test09_Write(self):
        """Writing a closed file."""

        self.assertRaises(ValueError, self.fnode.write, b"foo")

    def test10_Writelines(self):
        """Writing lines to a closed file."""

        self.assertRaises(ValueError, self.fnode.writelines, [b"foo\n"])


def copyFileToFile(srcfile, dstfile, blocksize=4096):
    """copyFileToFile(srcfile, dstfile[, blocksize]) -> None

    Copies a readable opened file 'srcfile' to a writable opened file
    'destfile' in blocks of 'blocksize' bytes (4 KiB by default).

    """

    data = srcfile.read(blocksize)
    while len(data) > 0:
        dstfile.write(data)
        data = srcfile.read(blocksize)


class WriteFileTestCase(TempFileMixin, TestCase):
    """Tests writing, seeking and truncating a new file node."""

    datafname = "test_filenode.dat"

    def setUp(self):
        """setUp() -> None

        This method sets the following instance attributes:
          * 'h5fname', the name of the temporary HDF5 file
          * 'h5file', the writable, temporary HDF5 file with a '/test' node
          * 'fnode', the writable file node in '/test'

        """

        super().setUp()
        self.fnode = filenode.new_node(self.h5file, where="/", name="test")
        self.datafname = test_file(self.datafname)

    def tearDown(self):
        """tearDown() -> None

        Closes 'fnode' and 'h5file'; removes 'h5fname'.

        """

        self.fnode.close()
        self.fnode = None
        super().tearDown()

    def test00_WriteFile(self):
        """Writing a whole file node."""

        datafile = open(self.datafname, "rb")
        try:
            copyFileToFile(datafile, self.fnode)
        finally:
            datafile.close()

    def test01_SeekFile(self):
        """Seeking and writing file node."""

        self.fnode.write(b"0123")
        self.fnode.seek(8)
        self.fnode.write(b"4567")
        self.fnode.seek(3)
        data = self.fnode.read(6)
        self.assertEqual(
            data,
            b"3\0\0\0\0" b"4",
            "Gap caused by forward seek was not properly filled.",
        )

        self.fnode.seek(0)
        self.fnode.write(b"test")

        self.fnode.seek(0)
        data = self.fnode.read(4)
        self.assertNotEqual(
            data, b"test", "Data was overwritten instead of appended."
        )

        self.fnode.seek(-4, 2)
        data = self.fnode.read(4)
        self.assertEqual(data, b"test", "Written data was not appended.")

        self.fnode.seek(0, 2)
        oldendoff = self.fnode.tell()
        self.fnode.seek(-2, 2)
        self.fnode.write(b"test")
        newendoff = self.fnode.tell()
        self.assertEqual(
            newendoff,
            oldendoff - 2 + 4,
            "Pointer was not correctly moved on append.",
        )

    def test02_TruncateFile(self):
        """Truncating a file node."""

        self.fnode.write(b"test")

        self.fnode.seek(2)
        self.assertRaises(IOError, self.fnode.truncate)

        self.fnode.seek(6)
        self.fnode.truncate()
        self.fnode.seek(0)
        data = self.fnode.read()
        self.assertEqual(
            data, b"test\0\0", "File was not grown to the current offset."
        )

        self.fnode.truncate(8)
        self.fnode.seek(0)
        data = self.fnode.read()
        self.assertEqual(
            data, b"test\0\0\0\0", "File was not grown to an absolute size."
        )


class OpenFileTestCase(TempFileMixin, TestCase):
    """Tests opening an existing file node for reading and writing."""

    def setUp(self):
        """setUp() -> None

        This method sets the following instance attributes:
          * 'h5fname', the name of the temporary HDF5 file
          * 'h5file', the writable, temporary HDF5 file with a '/test' node

        """

        super().setUp()
        fnode = filenode.new_node(self.h5file, where="/", name="test")
        fnode.close()

    def test00_OpenFileRead(self):
        """Opening an existing file node for reading."""

        node = self.h5file.get_node("/test")
        fnode = filenode.open_node(node)
        self.assertEqual(
            fnode.node, node, "filenode.open_node() opened the wrong node."
        )
        self.assertEqual(
            fnode.mode,
            "r",
            f"File was opened with an invalid mode {fnode.mode!r}.",
        )
        self.assertEqual(
            fnode.tell(),
            0,
            "Pointer is not positioned at the beginning of the file.",
        )
        fnode.close()

    def test01_OpenFileReadAppend(self):
        """Opening an existing file node for reading and appending."""

        node = self.h5file.get_node("/test")
        fnode = filenode.open_node(node, "a+")
        self.assertEqual(
            fnode.node, node, "filenode.open_node() opened the wrong node."
        )
        self.assertEqual(
            fnode.mode,
            "a+",
            f"File was opened with an invalid mode {fnode.mode!r}.",
        )

        self.assertEqual(
            fnode.tell(),
            0,
            "Pointer is not positioned at the beginning of the file.",
        )
        fnode.close()

    def test02_OpenFileInvalidMode(self):
        """Opening an existing file node with an invalid mode."""

        self.assertRaises(
            IOError, filenode.open_node, self.h5file.get_node("/test"), "w"
        )

    # This no longer works since type and type version attributes
    # are now system attributes.  ivb(2004-12-29)
    # def test03_OpenFileNoAttrs(self):
    #      "Opening a node with no type attributes."
    #
    #      node = self.h5file.get_node('/test')
    #      self.h5file.del_node_attr('/test', '_type')
    #      # Another way to get the same result is changing the value.
    #      ##self.h5file.set_node_attr('/test', '_type', 'foobar')
    #      self.assertRaises(ValueError, filenode.open_node, node)


class ReadFileTestCase(TempFileMixin, TestCase):
    """Tests reading from an existing file node."""

    datafname = "test_filenode.xbm"

    def setUp(self):
        """setUp() -> None

        This method sets the following instance attributes:
          * 'datafile', the opened data file
          * 'h5fname', the name of the temporary HDF5 file
          * 'h5file', the writable, temporary HDF5 file with a '/test' node
          * 'fnode', the readable file node in '/test', with data in it

        """

        self.datafname = test_file(self.datafname)
        self.datafile = open(self.datafname, "rb")

        super().setUp()

        fnode = filenode.new_node(self.h5file, where="/", name="test")
        copyFileToFile(self.datafile, fnode)
        fnode.close()

        self.datafile.seek(0)
        self.fnode = filenode.open_node(self.h5file.get_node("/test"))

    def tearDown(self):
        """tearDown() -> None

        Closes 'fnode', 'h5file' and 'datafile'; removes 'h5fname'.

        """

        self.fnode.close()
        self.fnode = None

        self.datafile.close()
        self.datafile = None

        super().tearDown()

    def test00_CompareFile(self):
        """Reading and comparing a whole file node."""

        import hashlib

        dfiledigest = hashlib.md5(self.datafile.read()).digest()
        fnodedigest = hashlib.md5(self.fnode.read()).digest()

        self.assertEqual(
            dfiledigest,
            fnodedigest,
            "Data read from file node differs from that in the file on disk.",
        )

    def test01_Write(self):
        """Writing on a read-only file."""

        self.assertRaises(IOError, self.fnode.write, "no way")

    def test02_UseAsImageFile(self):
        """Using a file node with Python Imaging Library."""

        try:
            from PIL import Image

            Image.open(self.fnode)
        except ImportError:
            # PIL not available, nothing to do.
            pass
        except OSError:
            self.fail(
                "PIL was not able to create an image from the file node."
            )

    def test_fileno(self):
        self.assertIsNot(self.fnode.fileno(), None)


class ReadlineTestCase(TempFileMixin, TestCase):
    """Base class for text line-reading test cases.

    It provides a set of tests independent of the line separator string.
    Sub-classes must provide the 'line_separator' attribute.

    """

    def setUp(self):
        """This method sets the following instance attributes:

        * ``h5fname``: the name of the temporary HDF5 file.
        * ``h5file``: the writable, temporary HDF5 file with a ``/test`` node.
        * ``fnode``: the readable file node in ``/test``, with text in it.

        """

        super().setUp()

        linesep = self.line_separator

        # Fill the node file with some text.
        fnode = filenode.new_node(self.h5file, where="/", name="test")
        # fnode.line_separator = linesep
        fnode.write(linesep)
        data = "short line%sshort line%s%s" % ((linesep.decode("ascii"),) * 3)
        data = data.encode("ascii")
        fnode.write(data)
        fnode.write(b"long line " * 20 + linesep)
        fnode.write(b"unterminated")
        fnode.close()

        # Re-open it for reading.
        self.fnode = filenode.open_node(self.h5file.get_node("/test"))
        # self.fnode.line_separator = linesep

    def tearDown(self):
        """tearDown() -> None

        Closes 'fnode' and 'h5file'; removes 'h5fname'.

        """

        self.fnode.close()
        self.fnode = None
        super().tearDown()

    def test00_Readline(self):
        """Reading individual lines."""

        linesep = self.line_separator

        line = self.fnode.readline()
        self.assertEqual(line, linesep)

        line = self.fnode.readline()  # 'short line' + linesep
        line = self.fnode.readline()
        self.assertEqual(line, b"short line" + linesep)
        line = self.fnode.readline()
        self.assertEqual(line, linesep)

        line = self.fnode.readline()
        self.assertEqual(line, b"long line " * 20 + linesep)

        line = self.fnode.readline()
        self.assertEqual(line, b"unterminated")

        line = self.fnode.readline()
        self.assertEqual(line, b"")

        line = self.fnode.readline()
        self.assertEqual(line, b"")

    def test01_ReadlineSeek(self):
        """Reading individual lines and seeking back and forth."""

        linesep = self.line_separator
        lseplen = len(linesep)

        self.fnode.readline()  # linesep
        self.fnode.readline()  # 'short line' + linesep

        self.fnode.seek(-(lseplen + 4), 1)
        line = self.fnode.readline()
        self.assertEqual(
            line, b"line" + linesep, "Seeking back yielded different data."
        )

        self.fnode.seek(lseplen + 20, 1)  # Into the long line.
        line = self.fnode.readline()
        self.assertEqual(
            line[-(lseplen + 10) :],
            b"long line " + linesep,
            "Seeking forth yielded unexpected data.",
        )

    def test02_Iterate(self):
        """Iterating over the lines."""

        linesep = self.line_separator

        # Iterate to the end.
        for line in self.fnode:
            pass

        self.assertRaises(StopIteration, next, self.fnode)

        self.fnode.seek(0)

        line = next(self.fnode)
        self.assertEqual(line, linesep)

        line = next(self.fnode)
        self.assertEqual(line, b"short line" + linesep)

    def test03_Readlines(self):
        """Reading a list of lines."""

        linesep = self.line_separator

        lines = self.fnode.readlines()
        self.assertEqual(
            lines,
            [
                linesep,
                b"short line" + linesep,
                b"short line" + linesep,
                linesep,
                b"long line " * 20 + linesep,
                b"unterminated",
            ],
        )

    def test04_ReadlineSize(self):
        """Reading individual lines of limited size."""

        linesep = self.line_separator
        lseplen = len(linesep)

        line = self.fnode.readline()  # linesep

        line = self.fnode.readline(lseplen + 20)
        self.assertEqual(line, b"short line" + linesep)

        line = self.fnode.readline(5)
        self.assertEqual(line, b"short")

        line = self.fnode.readline(lseplen + 20)
        self.assertEqual(line, b" line" + linesep)

        line = self.fnode.readline(lseplen)
        self.assertEqual(line, linesep)

        self.fnode.seek(-4, 2)
        line = self.fnode.readline(4)
        self.assertEqual(line, b"ated")

        self.fnode.seek(-4, 2)
        line = self.fnode.readline(20)
        self.assertEqual(line, b"ated")

    def test05_ReadlinesSize(self):
        """Reading a list of lines with a limited size."""

        linesep = self.line_separator

        data = "%sshort line%sshort" % ((linesep.decode("ascii"),) * 2)
        data = data.encode("ascii")
        lines = self.fnode.readlines(len(data))
        # self.assertEqual(lines, [linesep, b'short line' + linesep, b'short'])
        #
        # line = self.fnode.readline()
        # self.assertEqual(line, b' line' + linesep)

        # NOTE: the test is relaxed because the *hint* parameter of
        # io.BaseIO.readlines controls the amout of read data in a coarse way
        self.assertEqual(len(lines), len(data.split(b"\n")))
        self.assertEqual(lines[:-1], [linesep, b"short line" + linesep])
        self.assertTrue(lines[-1].startswith(b"short"))


class MonoReadlineTestCase(ReadlineTestCase):
    """Tests reading one-byte-separated text lines from an existing
    file node."""

    line_separator = b"\n"


# class MultiReadlineTestCase(ReadlineTestCase):
#    "Tests reading multibyte-separated text lines from an existing file node."
#
#    line_separator = b'<br/>'


# class LineSeparatorTestCase(TempFileMixin, TestCase):
#    "Tests text line separator manipulation in a file node."
#
#    def setUp(self):
#        """setUp() -> None
#
#        This method sets the following instance attributes:
#          * 'h5fname', the name of the temporary HDF5 file
#          * 'h5file', the writable, temporary HDF5 file with a '/test' node
#          * 'fnode', the writable file node in '/test'
#        """
#        super().setUp()
#        self.fnode = filenode.new_node(self.h5file, where='/', name='test')
#
#    def tearDown(self):
#        """tearDown() -> None
#
#        Closes 'fnode' and 'h5file'; removes 'h5fname'.
#        """
#        self.fnode.close()
#        self.fnode = None
#        super().tearDown()
#
#    def test00_DefaultLineSeparator(self):
#        "Default line separator."
#
#        self.assertEqual(
#            self.fnode.line_separator, os.linesep.encode('ascii'),
#            "Default line separator does not match that in os.linesep.")
#
#    def test01_SetLineSeparator(self):
#        "Setting a valid line separator."
#
#        try:
#            self.fnode.line_separator = b'SEPARATOR'
#        except ValueError:
#            self.fail("Valid line separator was not accepted.")
#        else:
#            self.assertEqual(
#                self.fnode.line_separator, b'SEPARATOR',
#                "Line separator was not correctly set.")
#
#    def test02_SetInvalidLineSeparator(self):
#        "Setting an invalid line separator."
#
#        self.assertRaises(
#            ValueError, setattr, self.fnode, 'line_separator', b'')
#        self.assertRaises(
#            ValueError, setattr, self.fnode, 'line_separator', b'x' * 1024)
#        self.assertRaises(
#            TypeError, setattr, self.fnode, 'line_separator', 'x')


class AttrsTestCase(TempFileMixin, TestCase):
    """Tests setting and getting file node attributes."""

    def setUp(self):
        """setUp() -> None

        This method sets the following instance attributes:
          * 'h5fname', the name of the temporary HDF5 file
          * 'h5file', the writable, temporary HDF5 file with a '/test' node
          * 'fnode', the writable file node in '/test'

        """

        super().setUp()
        self.fnode = filenode.new_node(self.h5file, where="/", name="test")

    def tearDown(self):
        """tearDown() -> None

        Closes 'fnode' and 'h5file'; removes 'h5fname'.

        """

        self.fnode.close()
        self.fnode = None
        super().tearDown()

    # This no longer works since type and type version attributes
    # are now system attributes.  ivb(2004-12-29)
    # def test00_GetTypeAttr(self):
    #      "Getting the type attribute of a file node."
    #
    #      self.assertEqual(
    #          getattr(self.fnode.attrs, '_type', None), filenode.NodeType,
    #          "File node has no '_type' attribute.")
    def test00_MangleTypeAttrs(self):
        """Mangling the type attributes on a file node."""

        nodeType = getattr(self.fnode.attrs, "NODE_TYPE", None)
        self.assertEqual(
            nodeType,
            filenode.NodeType,
            "File node does not have a valid 'NODE_TYPE' attribute.",
        )

        nodeTypeVersion = getattr(self.fnode.attrs, "NODE_TYPE_VERSION", None)
        self.assertTrue(
            nodeTypeVersion in filenode.NodeTypeVersions,
            "File node does not have a valid 'NODE_TYPE_VERSION' attribute.",
        )

        # System attributes are now writable.  ivb(2004-12-30)
        # self.assertRaises(
        #      AttributeError,
        #      setattr, self.fnode.attrs, 'NODE_TYPE', 'foobar')
        # self.assertRaises(
        #      AttributeError,
        #      setattr, self.fnode.attrs, 'NODE_TYPE_VERSION', 'foobar')

        # System attributes are now removables.  F. Alted (2007-03-06)

    #         self.assertRaises(
    #                 AttributeError,
    #                 delattr, self.fnode.attrs, 'NODE_TYPE')
    #         self.assertRaises(
    #                 AttributeError,
    #                 delattr, self.fnode.attrs, 'NODE_TYPE_VERSION')

    # System attributes are now writable.  ivb(2004-12-30)
    # def test01_SetSystemAttr(self):
    #      "Setting a system attribute on a file node."
    #
    #      self.assertRaises(
    # AttributeError, setattr, self.fnode.attrs, 'CLASS', 'foobar')
    def test02_SetGetDelUserAttr(self):
        """Setting a user attribute on a file node."""

        self.assertEqual(
            getattr(self.fnode.attrs, "userAttr", None),
            None,
            "Inexistent attribute has a value that is not 'None'.",
        )

        self.fnode.attrs.userAttr = "foobar"
        self.assertEqual(
            getattr(self.fnode.attrs, "userAttr", None),
            "foobar",
            "User attribute was not correctly set.",
        )

        self.fnode.attrs.userAttr = "bazquux"
        self.assertEqual(
            getattr(self.fnode.attrs, "userAttr", None),
            "bazquux",
            "User attribute was not correctly changed.",
        )

        del self.fnode.attrs.userAttr
        self.assertEqual(
            getattr(self.fnode.attrs, "userAttr", None),
            None,
            "User attribute was not deleted.",
        )
        # Another way is looking up the attribute in the attribute list.
        # if 'userAttr' in self.fnode.attrs._f_list():
        #      self.fail("User attribute was not deleted.")

    def test03_AttrsOnClosedFile(self):
        """Accessing attributes on a closed file node."""

        self.fnode.close()
        self.assertRaises(AttributeError, getattr, self.fnode, "attrs")


class ClosedH5FileTestCase(TempFileMixin, TestCase):
    """Tests accessing a file node in a closed PyTables file."""

    def setUp(self):
        """setUp() -> None

        This method sets the following instance attributes:
          * 'h5fname', the name of the temporary HDF5 file
          * 'h5file', the closed HDF5 file with a '/test' node
          * 'fnode', the writable file node in '/test'

        """

        super().setUp()
        self.fnode = filenode.new_node(self.h5file, where="/", name="test")
        self.h5file.close()

    def tearDown(self):
        """tearDown() -> None

        Closes 'fnode'; removes 'h5fname'.

        """

        # ivilata:  We know that a UserWarning will be raised
        #   because the PyTables file has already been closed.
        #   However, we don't want it to pollute the test output.
        warnings.filterwarnings("ignore", category=UserWarning)
        try:
            self.fnode.close()
        except ValueError:
            pass
        finally:
            warnings.filterwarnings("default", category=UserWarning)

        self.fnode = None
        super().tearDown()

    def test00_Write(self):
        """Writing to a file node in a closed PyTables file."""

        self.assertRaises(ValueError, self.fnode.write, "data")

    def test01_Attrs(self):
        """Accessing the attributes of a file node in a closed
        PyTables file."""

        self.assertRaises(ValueError, getattr, self.fnode, "attrs")


class OldVersionTestCase(TestCase):
    """Base class for old version compatibility test cases.

    It provides some basic tests for file operations and attribute handling.
    Sub-classes must provide the 'oldversion' attribute
    and the 'oldh5fname' attribute.

    """

    def setUp(self):
        """This method sets the following instance attributes:

        * ``h5fname``: the name of the temporary HDF5 file.
        * ``h5file``: the writable, temporary HDF5 file with a ``/test`` node.
        * ``fnode``: the readable file node in ``/test``.

        """

        super().setUp()
        self.h5fname = tempfile.mktemp(suffix=".h5")

        self.oldh5fname = test_file(self.oldh5fname)
        oldh5f = open_file(self.oldh5fname)
        oldh5f.copy_file(self.h5fname)
        oldh5f.close()

        self.h5file = open_file(
            self.h5fname,
            "r+",
            title="Test for file node old version compatibility",
        )
        self.fnode = filenode.open_node(self.h5file.root.test, "a+")

    def tearDown(self):
        """Closes ``fnode`` and ``h5file``; removes ``h5fname``."""

        self.fnode.close()
        self.fnode = None
        self.h5file.close()
        self.h5file = None
        Path(self.h5fname).unlink()
        super().tearDown()

    def test00_Read(self):
        """Reading an old version file node."""

        # self.fnode.line_separator = '\n'

        line = self.fnode.readline()
        self.assertEqual(line, "This is only\n")

        line = self.fnode.readline()
        self.assertEqual(line, "a test file\n")

        line = self.fnode.readline()
        self.assertEqual(line, "for FileNode version %d\n" % self.oldversion)

        line = self.fnode.readline()
        self.assertEqual(line, "")

        self.fnode.seek(0)
        line = self.fnode.readline()
        self.assertEqual(line, "This is only\n")

    def test01_Write(self):
        """Writing an old version file node."""

        # self.fnode.line_separator = '\n'

        self.fnode.write("foobar\n")
        self.fnode.seek(-7, 2)
        line = self.fnode.readline()
        self.assertEqual(line, "foobar\n")

    def test02_Attributes(self):
        """Accessing attributes in an old version file node."""

        self.fnode.attrs.userAttr = "foobar"
        self.assertEqual(
            getattr(self.fnode.attrs, "userAttr", None),
            "foobar",
            "User attribute was not correctly set.",
        )

        self.fnode.attrs.userAttr = "bazquux"
        self.assertEqual(
            getattr(self.fnode.attrs, "userAttr", None),
            "bazquux",
            "User attribute was not correctly changed.",
        )

        del self.fnode.attrs.userAttr
        self.assertEqual(
            getattr(self.fnode.attrs, "userAttr", None),
            None,
            "User attribute was not deleted.",
        )


class Version1TestCase(OldVersionTestCase):
    """Basic test for version 1 format compatibility."""

    oldversion = 1
    oldh5fname = "test_filenode_v1.h5"


class DirectReadWriteTestCase(TempFileMixin, TestCase):

    datafname = "test_filenode.dat"

    def setUp(self):
        """This method sets the following instance attributes:

        * ``h5fname``: the name of the temporary HDF5 file.
        * ``h5file``, the writable, temporary HDF5 file with a '/test' node
        * ``datafname``: the name of the data file to be stored in the
          temporary HDF5 file.
        * ``data``: the contents of the file ``datafname``
        * ``testfname``: the name of a temporary file to be written to.

        """

        super().setUp()
        self.datafname = test_file(self.datafname)
        self.testfname = tempfile.mktemp()
        self.testh5fname = tempfile.mktemp(suffix=".h5")
        self.data = Path(self.datafname).read_bytes()
        self.testdir = tempfile.mkdtemp()

    def tearDown(self):
        """tearDown() -> None

        Closes 'fnode' and 'h5file'; removes 'h5fname'.

        """
        if os.access(self.testfname, os.R_OK):
            Path(self.testfname).unlink()
        if os.access(self.testh5fname, os.R_OK):
            Path(self.testh5fname).unlink()
        shutil.rmtree(self.testdir)
        super().tearDown()

    def test01_WriteToPathlibPath(self):
        testh5fname = Path(self.testh5fname)
        datafname = Path(self.datafname)
        filenode.save_to_filenode(testh5fname, datafname, "/test1")

    def test01_WriteToFilename(self):
        # write contents of datafname to h5 testfile
        filenode.save_to_filenode(self.testh5fname, self.datafname, "/test1")
        # make sure writing to an existing node doesn't work ...
        self.assertRaises(
            IOError,
            filenode.save_to_filenode,
            self.testh5fname,
            self.datafname,
            "/test1",
        )
        # ... except if overwrite is True
        filenode.save_to_filenode(
            self.testh5fname, self.datafname, "/test1", overwrite=True
        )
        # write again, this time specifying a name
        filenode.save_to_filenode(
            self.testh5fname, self.datafname, "/", name="test2"
        )
        # read from test h5file
        filenode.read_from_filenode(self.testh5fname, self.testfname, "/test1")
        # and compare result to what it should be
        self.assertEqual(Path(self.testfname).read_bytes(), self.data)
        # make sure extracting to an existing file doesn't work ...
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.assertRaises(
                IOError,
                filenode.read_from_filenode,
                self.testh5fname,
                self.testfname,
                "/test1",
            )
        # except overwrite is True.  And try reading with a name
        filenode.read_from_filenode(
            self.testh5fname, self.testfname, "/", name="test2", overwrite=True
        )
        # and compare to what it should be
        self.assertEqual(Path(self.testfname).read_bytes(), self.data)
        # cleanup
        Path(self.testfname).unlink()
        Path(self.testh5fname).unlink()

    def test02_WriteToHDF5File(self):
        # write contents of datafname to h5 testfile
        filenode.save_to_filenode(self.h5file, self.datafname, "/test1")
        # make sure writing to an existing node doesn't work ...
        self.assertRaises(
            IOError,
            filenode.save_to_filenode,
            self.h5file,
            self.datafname,
            "/test1",
        )
        # ... except if overwrite is True
        filenode.save_to_filenode(
            self.h5file, self.datafname, "/test1", overwrite=True
        )
        # read from test h5file
        filenode.read_from_filenode(self.h5file, self.testfname, "/test1")
        # and compare result to what it should be
        self.assertEqual(Path(self.testfname).read_bytes(), self.data)
        # make sure extracting to an existing file doesn't work ...
        self.assertRaises(
            IOError,
            filenode.read_from_filenode,
            self.h5file,
            self.testfname,
            "/test1",
        )
        # make sure the original h5file is still alive and kicking
        self.assertEqual(isinstance(self.h5file, file.File), True)
        self.assertEqual(self.h5file.mode, "w")

    def test03_AutomaticNameGuessing(self):
        # write using the filename as node name
        filenode.save_to_filenode(self.testh5fname, self.datafname, "/")
        # and read again
        datafname = Path(self.datafname).name
        filenode.read_from_filenode(
            self.testh5fname,
            self.testdir,
            "/",
            name=datafname.replace(".", "_"),
        )
        # test if the output file really has the expected name
        self.assertEqual(
            os.access(Path(self.testdir) / datafname, os.R_OK), True
        )
        # and compare result to what it should be
        self.assertEqual(
            (Path(self.testdir) / datafname).read_bytes(), self.data
        )

    def test04_AutomaticNameGuessingWithFilenameAttribute(self):
        # write using the filename as node name
        filenode.save_to_filenode(self.testh5fname, self.datafname, "/")
        # and read again
        datafname = Path(self.datafname).name
        filenode.read_from_filenode(
            self.testh5fname, self.testdir, "/", name=datafname
        )
        # test if the output file really has the expected name
        self.assertEqual(
            os.access(Path(self.testdir) / datafname, os.R_OK), True
        )
        # and compare result to what it should be
        self.assertEqual(
            (Path(self.testdir) / datafname).read_bytes(), self.data
        )

    def test05_ReadFromNonexistingNodeRaises(self):
        # write using the filename as node name
        filenode.save_to_filenode(self.testh5fname, self.datafname, "/")
        # and read again
        self.assertRaises(
            NoSuchNodeError,
            filenode.read_from_filenode,
            self.testh5fname,
            self.testdir,
            "/",
            name="THISNODEDOESNOTEXIST",
        )


def suite():
    """suite() -> test suite

    Returns a test suite consisting of all the test cases in the module.

    """

    theSuite = unittest.TestSuite()

    theSuite.addTest(make_suite(NewFileTestCase))
    theSuite.addTest(make_suite(ClosedFileTestCase))
    theSuite.addTest(make_suite(WriteFileTestCase))
    theSuite.addTest(make_suite(OpenFileTestCase))
    theSuite.addTest(make_suite(ReadFileTestCase))
    theSuite.addTest(make_suite(MonoReadlineTestCase))
    # theSuite.addTest(make_suite(MultiReadlineTestCase))
    # theSuite.addTest(make_suite(LineSeparatorTestCase))
    theSuite.addTest(make_suite(AttrsTestCase))
    theSuite.addTest(make_suite(ClosedH5FileTestCase))
    theSuite.addTest(make_suite(DirectReadWriteTestCase))

    return theSuite


if __name__ == "__main__":
    import sys

    parse_argv(sys.argv)
    print_versions()
    unittest.main(defaultTest="suite")
