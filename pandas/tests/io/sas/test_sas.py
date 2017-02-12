import pandas.util.testing as tm
from pandas.compat import StringIO
from pandas import read_sas


class TestSas(tm.TestCase):

    def test_sas_buffer_format(self):

        # GH14947
        b = StringIO("")
        with self.assertRaises(ValueError):
            read_sas(b)
