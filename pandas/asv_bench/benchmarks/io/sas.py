import os

from pandas import read_sas


class SAS:

    params = ["sas7bdat", "xport"]
    param_names = ["format"]

    def setup(self, format):
        # Read files that are located in 'pandas/tests/io/sas/data'
        files = {"sas7bdat": "test1.sas7bdat", "xport": "paxraw_d_short.xpt"}
        file = files[format]
        paths = [
            os.path.dirname(__file__),
            "..",
            "..",
            "..",
            "pandas",
            "tests",
            "io",
            "sas",
            "data",
            file,
        ]
        self.f = os.path.join(*paths)

    def time_read_sas(self, format):
        read_sas(self.f, format=format)
