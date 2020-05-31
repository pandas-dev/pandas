"""Test logfmt format support"""

import pytest

import pandas as pd
from pandas import DataFrame, read_logfmt
import pandas._testing as tm


@pytest.fixture
def logfmt_file(datapath):
    """Path to logfmt file"""
    return datapath("io", "parser", "data", "logfmt.log")


@pytest.fixture
def logfmt_gzip_file(datapath):
    """Path to gzipped logfmt file"""
    return datapath("io", "parser", "data", "logfmt.log.gz")


class TestLogfmt:
    def test_read_logfmt(self, logfmt_file, logfmt_gzip_file):
        expected_df = DataFrame(
            [["first", 1, 1.0], ["second line", 2, 2.0], ['"third line"', 3, 3.0]],
            columns=["tag", "foo", "bar"],
        )

        parsed_df = read_logfmt(logfmt_file)
        tm.assert_frame_equal(parsed_df, expected_df)

        parsed_gzip_df = read_logfmt(logfmt_gzip_file)
        tm.assert_frame_equal(parsed_gzip_df, expected_df)

        logfmt_reader = read_logfmt(logfmt_file, chunksize=2)
        tm.assert_frame_equal(logfmt_reader.read(), expected_df)
