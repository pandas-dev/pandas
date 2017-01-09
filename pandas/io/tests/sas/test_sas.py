import pandas.util.testing as tm
from pandas.compat import StringIO
from pandas import read_sas


class TestSas(tm.TestCase):

    def test_sas_buffer_format(self):
        b = StringIO("")
        with self.assertRaises(TypeError):
            read_sas(b)
