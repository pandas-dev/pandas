import pytest

from pandas.compat import StringIO
from pandas import read_sas


class TestSas(object):

    def test_sas_buffer_format(self):
        # see gh-14947
        b = StringIO("")

        msg = ("If this is a buffer object rather than a string "
               "name, you must specify a format string")
        with pytest.raises(ValueError, match=msg):
            read_sas(b)
