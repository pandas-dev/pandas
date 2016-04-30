# -*- coding: utf-8 -*-

"""
Tests that features that are currently unsupported in
either the Python or C parser are actually enforced
and are clearly communicated to the user.

Ultimately, the goal is to remove test cases from this
test suite as new feature support is added to the parsers.
"""

import nose

import pandas.io.parsers as parsers
import pandas.util.testing as tm

from pandas.compat import StringIO
from pandas.io.common import CParserError
from pandas.io.parsers import read_csv, read_table


class TestUnsupportedFeatures(tm.TestCase):
    def test_c_engine(self):
        # see gh-6607
        data = 'a b c\n1 2 3'
        msg = 'does not support'

        # specify C-unsupported options with python-unsupported option
        # (options will be ignored on fallback, raise)
        with tm.assertRaisesRegexp(ValueError, msg):
            read_table(StringIO(data), sep=None,
                       delim_whitespace=False, dtype={'a': float})
        with tm.assertRaisesRegexp(ValueError, msg):
            read_table(StringIO(data), sep='\s', dtype={'a': float})
        with tm.assertRaisesRegexp(ValueError, msg):
            read_table(StringIO(data), skip_footer=1, dtype={'a': float})

        # specify C engine with unsupported options (raise)
        with tm.assertRaisesRegexp(ValueError, msg):
            read_table(StringIO(data), engine='c',
                       sep=None, delim_whitespace=False)
        with tm.assertRaisesRegexp(ValueError, msg):
            read_table(StringIO(data), engine='c', sep='\s')
        with tm.assertRaisesRegexp(ValueError, msg):
            read_table(StringIO(data), engine='c', skip_footer=1)

        # specify C-unsupported options without python-unsupported options
        with tm.assert_produces_warning(parsers.ParserWarning):
            read_table(StringIO(data), sep=None, delim_whitespace=False)
        with tm.assert_produces_warning(parsers.ParserWarning):
            read_table(StringIO(data), sep='\s')
        with tm.assert_produces_warning(parsers.ParserWarning):
            read_table(StringIO(data), skip_footer=1)

        text = """                      A       B       C       D        E
one two three   four
a   b   10.0032 5    -0.5109 -2.3358 -0.4645  0.05076  0.3640
a   q   20      4     0.4473  1.4152  0.2834  1.00661  0.1744
x   q   30      3    -0.6662 -0.5243 -0.3580  0.89145  2.5838"""
        msg = 'Error tokenizing data'

        with tm.assertRaisesRegexp(CParserError, msg):
            read_table(StringIO(text), sep='\s+')
        with tm.assertRaisesRegexp(CParserError, msg):
            read_table(StringIO(text), engine='c', sep='\s+')

        msg = "Only length-1 thousands markers supported"
        data = """A|B|C
1|2,334|5
10|13|10.
"""
        with tm.assertRaisesRegexp(ValueError, msg):
            read_csv(StringIO(data), thousands=',,')
        with tm.assertRaisesRegexp(ValueError, msg):
            read_csv(StringIO(data), thousands='')

        msg = "Only length-1 line terminators supported"
        data = 'a,b,c~~1,2,3~~4,5,6'
        with tm.assertRaisesRegexp(ValueError, msg):
            read_csv(StringIO(data), lineterminator='~~')

    def test_python_engine(self):
        from pandas.io.parsers import _python_unsupported as py_unsupported

        data = """1,2,3,,
1,2,3,4,
1,2,3,4,5
1,2,,,
1,2,3,4,"""
        engines = 'python', 'python-fwf'

        for engine in engines:
            for default in py_unsupported:
                msg = ('The %r option is not supported '
                       'with the %r engine' % (default, engine))

                kwargs = {default: object()}
                with tm.assertRaisesRegexp(ValueError, msg):
                    read_csv(StringIO(data), engine=engine, **kwargs)

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
