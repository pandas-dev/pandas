import os
import subprocess

import pandas.util.testing as tm

class TestParser(tm.TestCase):
    _multiprocess_can_split_ = True

    def test_parse_trim_buffers(self):
        code_ = """\n
import pandas as pd
from cStringIO import StringIO
record_ = "9999-9,99:99,,,,ZZ,ZZ,,,ZZZ-ZZZZ,.Z-ZZZZ,-9.99,,,9.99,ZZZZZ,,-99,9,ZZZ-ZZZZ,ZZ-ZZZZ,,9.99,ZZZ-ZZZZZ,ZZZ-ZZZZZ,ZZZ-ZZZZ,ZZZ-ZZZZ,ZZZ-ZZZZ,ZZZ-ZZZZ,ZZZ-ZZZZ,ZZZ-ZZZZ,999,ZZZ-ZZZZ,,ZZ-ZZZZ,,,,,ZZZZ,ZZZ-ZZZZZ,ZZZ-ZZZZ,,,9,9,9,9,99,99,999,999,ZZZZZ,ZZZ-ZZZZZ,ZZZ-ZZZZ,9,ZZ-ZZZZ,9.99,ZZ-ZZZZ,ZZ-ZZZZ,,,,ZZZZ,,,ZZ,ZZ,,,,,,,,,,,,,9,,,999.99,999.99,,,ZZZZZ,,,Z9,,,,,,,ZZZ,ZZZ,,,,,,,,,,,ZZZZZ,ZZZZZ,ZZZ-ZZZZZZ,ZZZ-ZZZZZZ,ZZ-ZZZZ,ZZ-ZZZZ,ZZ-ZZZZ,ZZ-ZZZZ,,,999999,999999,ZZZ,ZZZ,,,ZZZ,ZZZ,999.99,999.99,,,,ZZZ-ZZZ,ZZZ-ZZZ,-9.99,-9.99,9,9,,99,,9.99,9.99,9,9,9.99,9.99,,,,9.99,9.99,,99,,99,9.99,9.99,,,ZZZ,ZZZ,,999.99,,999.99,ZZZ,ZZZ-ZZZZ,ZZZ-ZZZZ,,,ZZZZZ,ZZZZZ,ZZZ,ZZZ,9,9,,,,,,ZZZ-ZZZZ,ZZZ999Z,,,999.99,,999.99,ZZZ-ZZZZ,,,9.999,9.999,9.999,9.999,-9.999,-9.999,-9.999,-9.999,9.999,9.999,9.999,9.999,9.999,9.999,9.999,9.999,99999,ZZZ-ZZZZ,,9.99,ZZZ,,,,,,,,ZZZ,,,,,9,,,,9,,,,,,,,,,ZZZ-ZZZZ,ZZZ-ZZZZ,,ZZZZZ,ZZZZZ,ZZZZZ,ZZZZZ,,,9.99,,ZZ-ZZZZ,ZZ-ZZZZ,ZZ,999,,,,ZZ-ZZZZ,ZZZ,ZZZ,ZZZ-ZZZZ,ZZZ-ZZZZ,,,99.99,99.99,,,9.99,9.99,9.99,9.99,ZZZ-ZZZZ,,,ZZZ-ZZZZZ,,,,,-9.99,-9.99,-9.99,-9.99,,,,,,,,,ZZZ-ZZZZ,,9,9.99,9.99,99ZZ,,-9.99,-9.99,ZZZ-ZZZZ,,,,,,,ZZZ-ZZZZ,9.99,9.99,9999,,,,,,,,,,-9.9,Z/Z-ZZZZ,999.99,9.99,,999.99,ZZ-ZZZZ,ZZ-ZZZZ,9.99,9.99,9.99,9.99,9.99,9.99,,ZZZ-ZZZZZ,ZZZ-ZZZZZ,ZZZ-ZZZZZ,ZZZ-ZZZZZ,ZZZ-ZZZZZ,ZZZ,ZZZ,ZZZ,ZZZ,9.99,,,-9.99,ZZ-ZZZZ,-999.99,,-9999,,999.99,,,,999.99,99.99,,,ZZ-ZZZZZZZZ,ZZ-ZZZZ-ZZZZZZZ,,,,ZZ-ZZ-ZZZZZZZZ,ZZZZZZZZ,ZZZ-ZZZZ,9999,999.99,ZZZ-ZZZZ,-9.99,-9.99,ZZZ-ZZZZ,99:99:99,,99,99,,9.99,,-99.99,,,,,,9.99,ZZZ-ZZZZ,-9.99,-9.99,9.99,9.99,,ZZZ,,,,,,,ZZZ,ZZZ,,,,,"
csv_data = "\\n".join([record_]*173) + "\\n"
for n_lines in range(82, 90):
    iterator_ = pd.read_csv(StringIO(csv_data), header=None, engine="c",
                            dtype=object, chunksize=n_lines, na_filter=True)
    for chunk_ in iterator_:
        print n_lines, chunk_.iloc[0, 0], chunk_.iloc[-1, 0]
exit(0)
"""
        expected_ = "".join("%d 9999-9 9999-9\n"%(n_lines,)
                            for n_lines in range(82, 90)
                            for _ in range((173 + n_lines - 1) // n_lines))

        # Run the faulty code via ang explicit argumnet to python
        proc_ = subprocess.Popen(("python", "-c", code_), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Wait until the subprocess finishes and then collect the output
        stdout_, stderr_ = proc_.communicate()
        exit_code = proc_.poll()

        # Check whether a segfault or memory corruption occurred
        # self.assertTrue(exit_code == -11 or (exit_code == 0 and stdout_ != expected_))

        # Check for correct exit code and output
        self.assertTrue(exit_code == 0 and stdout_ == expected_, msg="success")

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
