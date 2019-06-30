import functools
from collections import OrderedDict

import numpy as np
import pytest

import pandas as pd
import pandas.util.testing as tm

pytest.importorskip("odf")


@pytest.fixture(autouse=True)
def cd_and_set_engine(monkeypatch, datapath):
    func = functools.partial(pd.read_excel, engine="odf")
    monkeypatch.setattr(pd, 'read_excel', func)
    monkeypatch.chdir(datapath("io", "data"))


def test_read_invalid_types_raises():
    # the invalid_value_type.ods required manually editing
    # of the included content.xml file
    with pytest.raises(ValueError,
                       match="Unrecognized type awesome_new_type"):
        pd.read_excel("invalid_value_type.ods", header=None)


def test_read_lower_diagonal():
    # Make sure we can parse:
    # 1
    # 2 3
    # 4 5 6
    # 7 8 9 10

    sheet = pd.read_excel("lowerdiagonal.ods", 'Sheet1',
                          index_col=None, header=None)

    assert sheet.shape == (4, 4)


def test_read_writer_table():
    # Also test reading tables from an text OpenDocument file
    # (.odt)
    index = pd.Index(["Row 1", "Row 2", "Row 3"], name="Header")
    expected = pd.DataFrame([
        [1, np.nan, 7],
        [2, np.nan, 8],
        [3, np.nan, 9],
    ], index=index, columns=["Column 1", "Unnamed: 2", "Column 3"])

    result = pd.read_excel("writertable.odt", 'Table1', index_col=0)

    tm.assert_frame_equal(result, expected)


def test_blank_row_repeat():
    table = pd.read_excel("blank-row-repeat.ods", 'Value')

    assert table.shape == (14, 2)
    assert table['value'][7] == 9.0
    assert pd.isnull(table['value'][8])
    assert not pd.isnull(table['value'][11])


def test_runlengthencoding():
    sheet = pd.read_excel("runlengthencoding.ods", 'Sheet1', header=None)
    assert sheet.shape == (5, 3)
    # check by column, not by row.
    assert list(sheet[0]) == [1.0, 1.0, 2.0, 2.0, 2.0]
    assert list(sheet[1]) == [1.0, 2.0, 2.0, 2.0, 2.0]
    assert list(sheet[2]) == [1.0, 2.0, 2.0, 2.0, 2.0]
