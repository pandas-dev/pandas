import unittest


class TestExample(unittest.TestCase):

    def test_documentation_example1(self):
        # magictoken.ex_rec_arr_const_index.begin
        import numpy as np
        from numba import njit

        arr = np.array([(1, 2)], dtype=[('a1', 'f8'), ('a2', 'f8')])
        fields_gl = ('a1', 'a2')

        @njit
        def get_field_sum(rec):
            fields_lc = ('a1', 'a2')
            field_name1 = fields_lc[0]
            field_name2 = fields_gl[1]
            return rec[field_name1] + rec[field_name2]

        get_field_sum(arr[0])  # returns 3
        # magictoken.ex_rec_arr_const_index.end
        self.assertEqual(get_field_sum(arr[0]), 3)

    def test_documentation_example2(self):
        # magictoken.ex_rec_arr_lit_unroll_index.begin
        import numpy as np
        from numba import njit, literal_unroll

        arr = np.array([(1, 2)], dtype=[('a1', 'f8'), ('a2', 'f8')])
        fields_gl = ('a1', 'a2')

        @njit
        def get_field_sum(rec):
            out = 0
            for f in literal_unroll(fields_gl):
                out += rec[f]
            return out

        get_field_sum(arr[0])   # returns 3
        # magictoken.ex_rec_arr_lit_unroll_index.end
        self.assertEqual(get_field_sum(arr[0]), 3)


if __name__ == '__main__':
    unittest.main()
