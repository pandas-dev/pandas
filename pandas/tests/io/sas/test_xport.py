import numpy as np
import pytest

from pandas.errors import Pandas4Warning

import pandas as pd
import pandas._testing as tm

from pandas.io.sas.sas_xport import (
    XportReader,
    _parse_float_vec,
)
from pandas.io.sas.sasreader import read_sas

# CSV versions of test xpt files were obtained using the R foreign library

# Numbers in a SAS xport file are always float64, so need to convert
# before making comparisons.


def numeric_as_float(data):
    for v in data.columns:
        if data[v].dtype is np.dtype("int64"):
            data[v] = data[v].astype(np.float64)


class TestXport:
    @pytest.mark.slow
    def test1_basic(self, datapath):
        # Tests with DEMO_G.xpt (all numeric file)

        # Compare to this
        file01 = datapath("io", "sas", "data", "DEMO_G.xpt")
        data_csv = pd.read_csv(file01.replace(".xpt", ".csv"))
        numeric_as_float(data_csv)

        # Read full file
        data = read_sas(file01, format="xport")
        tm.assert_frame_equal(data, data_csv)
        num_rows = data.shape[0]

        # Test reading beyond end of file
        with read_sas(file01, format="xport", iterator=True) as reader:
            data = reader.read(num_rows + 100)
        assert data.shape[0] == num_rows

        # Test incremental read with `read` method.
        with read_sas(file01, format="xport", iterator=True) as reader:
            data = reader.read(10)
        tm.assert_frame_equal(data, data_csv.iloc[0:10, :])

        # Test incremental read with `get_chunk` method.
        with read_sas(file01, format="xport", chunksize=10) as reader:
            data = reader.get_chunk()
        tm.assert_frame_equal(data, data_csv.iloc[0:10, :])

        # Test read in loop
        m = 0
        with read_sas(file01, format="xport", chunksize=1000) as reader:
            for x in reader:
                m += x.shape[0]
        assert m == num_rows

        # Read full file with `read_sas` method
        data = read_sas(file01)
        tm.assert_frame_equal(data, data_csv)

    def test1_index(self, datapath):
        # Tests with DEMO_G.xpt using index (all numeric file)

        # Compare to this
        file01 = datapath("io", "sas", "data", "DEMO_G.xpt")
        data_csv = pd.read_csv(file01.replace(".xpt", ".csv"))
        data_csv = data_csv.set_index("SEQN")
        numeric_as_float(data_csv)

        # Read full file
        data = read_sas(file01, index="SEQN", format="xport")
        tm.assert_frame_equal(data, data_csv, check_index_type=False)

        # Test incremental read with `read` method.
        with read_sas(file01, index="SEQN", format="xport", iterator=True) as reader:
            data = reader.read(10)
        tm.assert_frame_equal(data, data_csv.iloc[0:10, :], check_index_type=False)

        # Test incremental read with `get_chunk` method.
        with read_sas(file01, index="SEQN", format="xport", chunksize=10) as reader:
            data = reader.get_chunk()
        tm.assert_frame_equal(data, data_csv.iloc[0:10, :], check_index_type=False)

    def test1_incremental(self, datapath):
        # Test with DEMO_G.xpt, reading full file incrementally

        file01 = datapath("io", "sas", "data", "DEMO_G.xpt")
        data_csv = pd.read_csv(file01.replace(".xpt", ".csv"))
        data_csv = data_csv.set_index("SEQN")
        numeric_as_float(data_csv)

        with read_sas(file01, index="SEQN", chunksize=1000) as reader:
            all_data = list(reader)
        data = pd.concat(all_data, axis=0)

        tm.assert_frame_equal(data, data_csv, check_index_type=False)

    def test2(self, datapath):
        # Test with SSHSV1_A.xpt

        file02 = datapath("io", "sas", "data", "SSHSV1_A.xpt")
        # Compare to this
        data_csv = pd.read_csv(file02.replace(".xpt", ".csv"))
        numeric_as_float(data_csv)

        data = read_sas(file02)
        tm.assert_frame_equal(data, data_csv)

    def test2_binary(self, datapath):
        # Test with SSHSV1_A.xpt, read as a binary file

        # Compare to this
        file02 = datapath("io", "sas", "data", "SSHSV1_A.xpt")
        data_csv = pd.read_csv(file02.replace(".xpt", ".csv"))
        numeric_as_float(data_csv)

        with open(file02, "rb") as fd:
            # GH#35693 ensure that if we pass an open file, we
            #  dont incorrectly close it in read_sas
            data = read_sas(fd, format="xport")

        tm.assert_frame_equal(data, data_csv)

    def test_multiple_types(self, datapath):
        # Test with DRXFCD_G.xpt (contains text and numeric variables)

        # Compare to this
        file03 = datapath("io", "sas", "data", "DRXFCD_G.xpt")
        data_csv = pd.read_csv(file03.replace(".xpt", ".csv"))

        data = read_sas(file03, encoding="utf-8")
        tm.assert_frame_equal(data, data_csv)

    def test_encoding_default_deprecated(self, datapath):
        # GH#66470
        file03 = datapath("io", "sas", "data", "DRXFCD_G.xpt")
        msg = "The default value of 'encoding' in read_sas is deprecated"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            result = read_sas(file03)
        assert result["DRXFCSD"].iloc[0] == b"MILK, HUMAN"

        for encoding in [None, "infer", "utf-8"]:
            with tm.assert_produces_warning(None):
                read_sas(file03, encoding=encoding)

    def test_encoding_infer(self, datapath):
        # GH#66470 XPORT files record no encoding, so "infer" means the
        #  XportReader default and must not raise
        file03 = datapath("io", "sas", "data", "DRXFCD_G.xpt")
        result = read_sas(file03, encoding="infer")
        expected = read_sas(file03, encoding=XportReader._default_encoding)
        tm.assert_frame_equal(result, expected)

        with XportReader(file03) as reader:
            tm.assert_frame_equal(reader.read(), expected)

    def test_truncated_float_support(self, datapath):
        # Test with paxraw_d_short.xpt, a shortened version of:
        # http://wwwn.cdc.gov/Nchs/Nhanes/2005-2006/PAXRAW_D.ZIP
        # This file has truncated floats (5 bytes in this case).

        # GH 11713
        file04 = datapath("io", "sas", "data", "paxraw_d_short.xpt")
        data_csv = pd.read_csv(file04.replace(".xpt", ".csv"))

        data = read_sas(file04, format="xport")
        tm.assert_frame_equal(data.astype("int64"), data_csv)

    @pytest.mark.parametrize("fname", ["DEMO_G.xpt", "paxraw_d_short.xpt"])
    def test_zero_read_exactly(self, datapath, fname):
        # GH#50670 IBM zero was read as 5.397605e-79, which assert_frame_equal's
        # default tolerance hides; paxraw_d_short.xpt has truncated floats
        path = datapath("io", "sas", "data", fname)
        expected = pd.read_csv(path.replace(".xpt", ".csv"))

        result = read_sas(path, format="xport", encoding="infer")
        assert (expected == 0).any(axis=None)
        tm.assert_frame_equal(result == 0, expected == 0)

    def test_parse_float_vec_zero(self):
        # GH#50670 negative zero and a zero fraction with nonzero exponent
        vec = np.array(
            [b"\x00" * 8, b"\x80" + b"\x00" * 7, b"\x40" + b"\x00" * 7], dtype="S8"
        )
        result = _parse_float_vec(vec)
        tm.assert_numpy_array_equal(result, np.zeros(3))
        assert list(np.signbit(result)) == [False, True, False]

    def test_convert_dates(self, datapath):
        # GH#70876 numeric columns with a SAS date or datetime format are
        # converted like the SAS7BDAT reader does. dates.xpt was written with
        # pyreadstat (file_format_version=5): DT has format DATE9., DTTM has
        # DATETIME20., NUM has the non-date format BEST12., ID has none.
        path = datapath("io", "sas", "data", "dates.xpt")
        expected = pd.DataFrame(
            {
                "ID": [1.0, 2.0, 3.0, 4.0],
                "DT": pd.to_datetime(
                    ["1960-01-01", "2021-02-16", None, "1899-12-31"]
                ).as_unit("s"),
                "DTTM": pd.to_datetime(
                    [
                        "1960-01-01 00:00:00",
                        "2021-02-16 10:07:55",
                        None,
                        "1959-12-31 23:59:59",
                    ]
                ).as_unit("ms"),
                "NUM": [1.5, -2.25, np.nan, 1e10],
                "TXT": ["a", "bb", "", "dddd"],
            }
        )

        result = read_sas(path, encoding="infer")
        tm.assert_frame_equal(result, expected)

        with XportReader(path, encoding="infer") as reader:
            tm.assert_frame_equal(reader.read(), expected)

    def test_convert_dates_false(self, datapath):
        # GH#70876 the raw SAS day and second counts
        path = datapath("io", "sas", "data", "dates.xpt")
        with XportReader(path, encoding="infer", convert_dates=False) as reader:
            result = reader.read()

        expected = pd.DataFrame(
            {
                "ID": [1.0, 2.0, 3.0, 4.0],
                "DT": [0.0, 22327.0, np.nan, -21915.0],
                "DTTM": [0.0, 22327.0 * 86400 + 36475, np.nan, -1.0],
                "NUM": [1.5, -2.25, np.nan, 1e10],
                "TXT": ["a", "bb", "", "dddd"],
            }
        )
        tm.assert_frame_equal(result, expected)

    def test_convert_dates_chunked(self, datapath):
        # GH#70876 every chunk is converted, not only the first
        path = datapath("io", "sas", "data", "dates.xpt")
        expected = read_sas(path, encoding="infer")
        assert expected["DT"].dtype == "M8[s]"
        assert expected["DTTM"].dtype == "M8[ms]"

        with read_sas(path, encoding="infer", chunksize=3) as reader:
            chunks = list(reader)
        assert [len(chunk) for chunk in chunks] == [3, 1]
        for chunk in chunks:
            assert chunk["DT"].dtype == "M8[s]"
            assert chunk["DTTM"].dtype == "M8[ms]"
        tm.assert_frame_equal(pd.concat(chunks), expected)

    def test_cport_header_found_raises(self, datapath):
        # Test with DEMO_PUF.cpt, the beginning of puf2019_1_fall.xpt
        # from https://www.cms.gov/files/zip/puf2019.zip
        # (despite the extension, it's a cpt file)
        msg = "Header record indicates a CPORT file, which is not readable."
        with pytest.raises(ValueError, match=msg):
            read_sas(datapath("io", "sas", "data", "DEMO_PUF.cpt"), format="xport")
