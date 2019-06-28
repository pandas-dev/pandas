from collections import OrderedDict

import numpy as np
import pytest

import pandas as pd
import pandas.util.testing as tm

pytest.importorskip("odf")


def test_get_sheet(datapath):
    from pandas.io.excel._odfreader import ODFReader

    pth = datapath("io", "data", "datatypes.ods")
    book = ODFReader(pth)

    assert len(book.sheet_names) == 1
    assert book.sheet_names == ['Sheet1']


def test_get_sheet_raises(datapath):
    from pandas.io.excel._odfreader import ODFReader

    pth = datapath("io", "data", 'datatypes.ods')
    book = ODFReader(pth)

    with pytest.raises(ValueError):
        book._get_sheet(3.14)

    with pytest.raises(ValueError):
        book.get_sheet_by_name("Invalid Sheet 77")

    with pytest.raises(IndexError):
        book.get_sheet_by_index(-33)


def test_read_types(datapath):
    path = datapath("io", "data", "datatypes.ods")
    sheet = pd.read_excel(path, header=None, engine='odf')

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


def test_read_invalid_types_raises(datapath):
    # the invalid_value_type.ods required manually editing
    # of the included content.xml file
    path = datapath("io", "data", "invalid_value_type.ods")
    with pytest.raises(ValueError,
                       match="Unrecognized type awesome_new_type"):
        pd.read_excel(path, header=None, engine='odf')


def test_read_lower_diagonal(datapath):
    # Make sure we can parse:
    # 1
    # 2 3
    # 4 5 6
    # 7 8 9 10
    path = datapath("io", "data", "lowerdiagonal.ods")
    sheet = pd.read_excel(path, 'Sheet1',
                          index_col=None, header=None, engine='odf')

    assert sheet.shape == (4, 4)


def test_read_headers(datapath):
    path = datapath("io", "data", "headers.ods")
    sheet = pd.read_excel(path, 'Sheet1', index_col=0, engine='odf')

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


def test_read_writer_table(datapath):
    # Also test reading tables from an text OpenDocument file
    # (.odt)

    path = datapath("io", "data", "writertable.odt")
    table = pd.read_excel(path, 'Table1', index_col=0, engine='odf')

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


def test_blank_row_repeat(datapath):
    path = datapath("io", "data", "blank-row-repeat.ods")
    table = pd.read_excel(path, 'Value', engine='odf')

    assert table.shape == (14, 2)
    assert table['value'][7] == 9.0
    assert pd.isnull(table['value'][8])
    assert not pd.isnull(table['value'][11])


def test_runlengthencoding(datapath):
    path = datapath("io", "data", "runlengthencoding.ods")
    sheet = pd.read_excel(path, 'Sheet1', header=None, engine='odf')
    assert sheet.shape == (5, 3)
    # check by column, not by row.
    assert list(sheet[0]) == [1.0, 1.0, 2.0, 2.0, 2.0]
    assert list(sheet[1]) == [1.0, 2.0, 2.0, 2.0, 2.0]
    assert list(sheet[2]) == [1.0, 2.0, 2.0, 2.0, 2.0]
