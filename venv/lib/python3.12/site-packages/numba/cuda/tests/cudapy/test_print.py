from numba.cuda.testing import CUDATestCase, skip_on_cudasim
import subprocess
import sys
import unittest


cuhello_usecase = """\
from numba import cuda

@cuda.jit
def cuhello():
    i = cuda.grid(1)
    print(i, 999)
    print(-42)

cuhello[2, 3]()
cuda.synchronize()
"""


printfloat_usecase = """\
from numba import cuda

@cuda.jit
def printfloat():
    i = cuda.grid(1)
    print(i, 23, 34.75, 321)

printfloat[1, 1]()
cuda.synchronize()
"""


printstring_usecase = """\
from numba import cuda

@cuda.jit
def printstring():
    i = cuda.grid(1)
    print(i, "hop!", 999)

printstring[1, 3]()
cuda.synchronize()
"""

printempty_usecase = """\
from numba import cuda

@cuda.jit
def printempty():
    print()

printempty[1, 1]()
cuda.synchronize()
"""


print_too_many_usecase = """\
from numba import cuda
import numpy as np

@cuda.jit
def print_too_many(r):
    print(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9], r[10],
          r[11], r[12], r[13], r[14], r[15], r[16], r[17], r[18], r[19], r[20],
          r[21], r[22], r[23], r[24], r[25], r[26], r[27], r[28], r[29], r[30],
          r[31], r[32])

print_too_many[1, 1](np.arange(33))
cuda.synchronize()
"""


class TestPrint(CUDATestCase):
    # Note that in these tests we generally strip the output to avoid dealing
    # with platform-specific line ending issues, e.g. '\r\n' vs '\n' etc.

    def run_code(self, code):
        """Runs code in a subprocess and returns the captured output"""
        cmd = [sys.executable, "-c", code]
        cp = subprocess.run(cmd, timeout=60, capture_output=True, check=True)
        return cp.stdout.decode(), cp.stderr.decode()

    def test_cuhello(self):
        output, _ = self.run_code(cuhello_usecase)
        actual = [line.strip() for line in output.splitlines()]
        expected = ['-42'] * 6 + ['%d 999' % i for i in range(6)]
        # The output of GPU threads is intermingled, but each print()
        # call is still atomic
        self.assertEqual(sorted(actual), expected)

    def test_printfloat(self):
        output, _ = self.run_code(printfloat_usecase)
        # CUDA and the simulator use different formats for float formatting
        expected_cases = ["0 23 34.750000 321", "0 23 34.75 321"]
        self.assertIn(output.strip(), expected_cases)

    def test_printempty(self):
        output, _ = self.run_code(printempty_usecase)
        self.assertEqual(output.strip(), "")

    def test_string(self):
        output, _ = self.run_code(printstring_usecase)
        lines = [line.strip() for line in output.splitlines(True)]
        expected = ['%d hop! 999' % i for i in range(3)]
        self.assertEqual(sorted(lines), expected)

    @skip_on_cudasim('cudasim can print unlimited output')
    def test_too_many_args(self):
        # Tests that we emit the format string and warn when there are more
        # than 32 arguments, in common with CUDA C/C++ printf - this is due to
        # a limitation in CUDA vprintf, see:
        # https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#limitations

        output, errors = self.run_code(print_too_many_usecase)

        # Check that the format string was printed instead of formatted garbage
        expected_fmt_string = ' '.join(['%lld' for _ in range(33)])
        self.assertIn(expected_fmt_string, output)

        # Check for the expected warning about formatting more than 32 items
        warn_msg = ('CUDA print() cannot print more than 32 items. The raw '
                    'format string will be emitted by the kernel instead.')
        self.assertIn(warn_msg, errors)


if __name__ == '__main__':
    unittest.main()
