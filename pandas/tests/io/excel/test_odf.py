import functools

import numpy as np
import pytest

import pandas as pd
import pandas.util.testing as tm

pytest.importorskip("odf")


@pytest.fixture(autouse=True)
def cd_and_set_engine(monkeypatch, datapath):
    func = functools.partial(pd.read_excel, engine="odf")
    monkeypatch.setattr(pd, "read_excel", func)
    monkeypatch.chdir(datapath("io", "data", "excel"))


def test_read_invalid_types_raises():
    # the invalid_value_type.ods required manually editing
    # of the included content.xml file
    with pytest.raises(ValueError, match="Unrecognized type awesome_new_type"):
        pd.read_excel("invalid_value_type.ods")


def test_read_writer_table():
    # Also test reading tables from an text OpenDocument file
    # (.odt)
    index = pd.Index(["Row 1", "Row 2", "Row 3"], name="Header")
    expected = pd.DataFrame(
        [[1, np.nan, 7], [2, np.nan, 8], [3, np.nan, 9]],
        index=index,
        columns=["Column 1", "Unnamed: 2", "Column 3"],
    )

    result = pd.read_excel("writertable.odt", "Table1", index_col=0)

    tm.assert_frame_equal(result, expected)


def test_nonexistent_sheetname_raises(read_ext):
    # GH-27676
    # Specifying a non-existent sheet_name parameter should throw an error
    # with the sheet name.
    with pytest.raises(ValueError, match="sheet xyz not found"):
        pd.read_excel("blank.ods", sheet_name="xyz")
