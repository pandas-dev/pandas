def test_sas_buffer_format(self):
    import StringIO
    b = StringIO.StringIO("")
    with self.assertRaises(TypeError):
        result=pd.read_sas(b)
