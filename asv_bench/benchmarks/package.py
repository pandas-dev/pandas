"""
Benchmarks for pandas at the package-level.
"""
import subprocess
import sys

from pandas.compat import PY37

class TimeImport:
    def time_import(self):
        if PY37:
            cmd = [sys.executable, "-X", "importtime", "-c", "import pandas as pd"]
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            _, stderr = p.communicate()

            line = stderr.splitlines()[-1]
            field = line.split(b"|")[-2].strip()
            total = int(field)  # microseconds
            return total

        cmd = [sys.executable, "-c", "import pandas as pd"]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p.communicate()
