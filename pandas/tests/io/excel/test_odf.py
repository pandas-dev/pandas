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


def test_read_types():
    sheet = pd.read_excel("datatypes.ods", header=None)

    expected = pd.DataFrame(
        [[1.0],
         [1.25],
         ['a'],
         [pd.Timestamp(2003, 1, 2)],
         [False],
         [0.35],
         [pd.Timedelta(hours=3, minutes=45),
          pd.Timedelta(hours=17, minutes=53),
          pd.Timedelta(hours=14, minutes=8)],
         # though what should the value of a hyperlink be?
         ['UBERON:0002101']])
    tm.assert_equal(sheet, expected)


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


def test_read_headers():
    sheet = pd.read_excel("headers.ods", 'Sheet1', index_col=0)

    expected = pd.DataFrame.from_dict(OrderedDict([
        ("Header", ["Row 1", "Row 2"]),
        ("Column 1", [1.0, 2.0]),
        ("Column 2", [3.0, 4.0]),
        # Empty Column
        ("Column 4", [7.0, 8.0]),
        # Empty Column 2
        ("Column 6", [11.0, 12.0])]))
    expected.set_index("Header", inplace=True)
    columns = ["Column 1", "Column 2", "Column 4", "Column 6"]
    tm.assert_equal(sheet[columns], expected)
    empties = [None, 'None.1']
    for name in empties:
        for value in sheet[name]:
            assert pd.isnull(value)


def test_read_writer_table():
    # Also test reading tables from an text OpenDocument file
    # (.odt)

    table = pd.read_excel("writertable.odt", 'Table1', index_col=0)

    assert table.shape == (3, 3)
    expected = pd.DataFrame.from_dict(OrderedDict([
        ("Header", ["Row 1", "Row 2", "Row 3"]),
        ("Column 1", [1.0, 2.0, 3.0]),
        ("Unnamed: 2", [np.nan, np.nan, np.nan]),
        ("Column 3", [7.0, 8.0, 9.0])]))
    expected.set_index("Header", inplace=True)
    columns = ["Column 1", "Column 3"]
    tm.assert_equal(table[columns], expected[columns])

    # make sure pandas gives a name to the unnamed column
    for i in range(3):
        assert pd.isnull(table["Unnamed: 2"][i])


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
