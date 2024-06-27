import sys
from io import StringIO

from unittest.mock import patch

import tables.scripts.ptrepack as ptrepack
import tables.scripts.ptdump as ptdump
import tables.scripts.pttree as pttree
from tables.tests import common


class ptrepackTestCase(common.PyTablesTestCase):
    """Test ptrepack"""

    @patch.object(ptrepack, 'copy_leaf')
    @patch.object(ptrepack.tb, 'open_file')
    def test_paths_windows(self, mock_open_file, mock_copy_leaf):
        """Checking handling of windows filenames: test gh-616"""

        # this filename has a semi-colon to check for
        # regression of gh-616
        src_fn = 'D:\\window~1\\path\\000\\infile'
        src_path = '/'
        dst_fn = 'another\\path\\'
        dst_path = '/path/in/outfile'

        argv = ['ptrepack', src_fn + ':' + src_path, dst_fn + ':' + dst_path]
        with patch.object(sys, 'argv', argv):
            ptrepack.main()

        args, kwargs = mock_open_file.call_args_list[0]
        self.assertEqual(args, (src_fn, 'r'))

        args, kwargs = mock_copy_leaf.call_args_list[0]
        self.assertEqual(args, (src_fn, dst_fn, src_path, dst_path))


class ptdumpTestCase(common.PyTablesTestCase):
    """Test ptdump"""

    @patch.object(ptdump.tb, 'open_file')
    @patch('sys.stdout', new_callable=StringIO)
    def test_paths_windows(self, _, mock_open_file):
        """Checking handling of windows filenames: test gh-616"""

        # this filename has a semi-colon to check for
        # regression of gh-616 (in ptdump)
        src_fn = 'D:\\window~1\\path\\000\\ptdump'
        src_path = '/'

        argv = ['ptdump', src_fn + ':' + src_path]
        with patch.object(sys, 'argv', argv):
            ptdump.main()

        args, kwargs = mock_open_file.call_args_list[0]
        self.assertEqual(args, (src_fn, 'r'))


class pttreeTestCase(common.PyTablesTestCase):
    """Test ptdump"""

    @patch.object(pttree.tb, 'open_file')
    @patch.object(pttree, 'get_tree_str')
    @patch('sys.stdout', new_callable=StringIO)
    def test_paths_windows(self, _, mock_get_tree_str, mock_open_file):
        """Checking handling of windows filenames: test gh-616"""

        # this filename has a semi-colon to check for
        # regression of gh-616 (in pttree)
        src_fn = 'D:\\window~1\\path\\000\\pttree'
        src_path = '/'

        argv = ['pttree', src_fn + ':' + src_path]
        with patch.object(sys, 'argv', argv):
            pttree.main()

        args, kwargs = mock_open_file.call_args_list[0]
        self.assertEqual(args, (src_fn, 'r'))


def suite():
    theSuite = common.unittest.TestSuite()

    theSuite.addTest(common.unittest.makeSuite(ptrepackTestCase))
    theSuite.addTest(common.unittest.makeSuite(ptdumpTestCase))
    theSuite.addTest(common.unittest.makeSuite(pttreeTestCase))

    return theSuite


if __name__ == '__main__':
    common.parse_argv(sys.argv)
    common.print_versions()
    common.unittest.main(defaultTest='suite')
