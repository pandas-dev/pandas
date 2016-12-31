import pandas.util.testing as tm

class TestSasBuff(tm.TestCase):
    def test_sas_buffer_format(self):
        import StringIO
        from pandas.io.sas.sasreader import read_sas
        b = StringIO.StringIO("")
        with self.assertRaises(TypeError):
            read_sas(b)
