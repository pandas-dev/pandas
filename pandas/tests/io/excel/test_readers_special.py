# Tests that don't need or don't work with the autouse fixtures in test_readers.py
import pytest

import pandas as pd


def test_unreadable_bytes():
    with pytest.raises(
        ValueError, match=r"Could not find engine for .+, content was b'rubbish'"
    ):
        pd.read_excel(b"rubbish")


def test_unreadable_file(tmp_path):
    bad = tmp_path / "bad"
    bad.write_bytes(b"rubbish")
    with pytest.raises(
        ValueError, match=r"Could not find engine for .+, content was b'rubbish'"
    ):
        pd.read_excel(bad)
