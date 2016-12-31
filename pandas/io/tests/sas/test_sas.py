from pandas.io.sas.sasreader import read_sas
import pandas.util.testing as tm
import pandas as pd

class TestSasBuff(tm.TestCase):

    def test_sas_buffer_format(self):
        import StringIO
        from pandas.io.sas.sasreader import read_sas
        b = StringIO.StringIO("")
        with self.assertRaises(TypeError):
            result=read_sas(b)

