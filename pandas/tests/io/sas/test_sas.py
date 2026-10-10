from io import StringIO

import pytest

import pandas as pd
import pandas._testing as tm


class TestSas:
    def test_sas_buffer_format(self):
        # see gh-14947
        b = StringIO("")

        msg = (
            "If this is a buffer object rather than a string "
            "name, you must specify a format string"
        )
        with pytest.raises(ValueError, match=msg):
            pd.read_sas(b)

    def test_sas_read_no_format_or_extension(self, temp_file):
        # see gh-24548
        msg = "unable to infer format of SAS file.+"
        with pytest.raises(ValueError, match=msg):
            pd.read_sas(temp_file)


def test_sas_archive(datapath):
    fname_uncompressed = datapath("io", "sas", "data", "airline.sas7bdat")
    df_uncompressed = pd.read_sas(fname_uncompressed)
    fname_compressed = datapath("io", "sas", "data", "airline.sas7bdat.gz")
    df_compressed = pd.read_sas(fname_compressed, format="sas7bdat")
    tm.assert_frame_equal(df_uncompressed, df_compressed)
