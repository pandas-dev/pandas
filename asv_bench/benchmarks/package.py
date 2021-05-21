"""
Benchmarks for pandas at the package-level.
"""
import subprocess
import sys


class TimeImport:
    def time_import(self):
        # on py37+ we the "-X importtime" usage gives us a more precise
        #  measurement of the import time we actually care about,
        #  without the subprocess or interpreter overhead
        cmd = [sys.executable, "-X", "importtime", "-c", "import pandas as pd"]
        p = subprocess.run(cmd, stderr=subprocess.PIPE)

        line = p.stderr.splitlines()[-1]
        field = line.split(b"|")[-2].strip()
        total = int(field)  # microseconds
        return total
