import functools
import uuid

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm

pytest.importorskip("odf")


@pytest.fixture
def ext():
    return ".ods"


@pytest.fixture
def tmp_excel(ext, tmp_path):
    tmp = tmp_path / f"{uuid.uuid4()}{ext}"
    tmp.touch()
    return str(tmp)


from odf.opendocument import OpenDocumentSpreadsheet
from odf.table import (
    CoveredTableCell,
    Table,
    TableCell,
    TableRow,
)


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
    # Also test reading tables from a text OpenDocument file
    # (.odt)
    index = pd.Index(["Row 1", "Row 2", "Row 3"], name="Header")
    expected = pd.DataFrame(
        [[1, np.nan, 7], [2, np.nan, 8], [3, np.nan, 9]],
        index=index,
        columns=["Column 1", "Unnamed: 2", "Column 3"],
    )

    result = pd.read_excel("writertable.odt", sheet_name="Table1", index_col=0)

    tm.assert_frame_equal(result, expected)


def test_read_newlines_between_xml_elements_table():
    # GH#45598
    expected = pd.DataFrame(
        [[1.0, 4.0, 7], [np.nan, np.nan, 8], [3.0, 6.0, 9]],
        columns=["Column 1", "Column 2", "Column 3"],
    )

    result = pd.read_excel("test_newlines.ods")

    tm.assert_frame_equal(result, expected)


def test_read_unempty_cells():
    expected = pd.DataFrame(
        [1, np.nan, 3, np.nan, 5],
        columns=["Column 1"],
    )

    result = pd.read_excel("test_unempty_cells.ods")

    tm.assert_frame_equal(result, expected)


def test_read_cell_annotation():
    expected = pd.DataFrame(
        ["test", np.nan, "test 3"],
        columns=["Column 1"],
    )

    result = pd.read_excel("test_cell_annotation.ods")

    tm.assert_frame_equal(result, expected)


def test_read_covered_table_cell_value(tmp_excel):
    # GH#66579
    doc = OpenDocumentSpreadsheet()
    sheet = Table(name="Sheet1")
    doc.spreadsheet.addElement(sheet)

    row0 = TableRow()
    sheet.addElement(row0)
    row0.addElement(TableCell(valuetype="float", value="1"))
    row0.addElement(TableCell(valuetype="float", value="100"))

    row1 = TableRow()
    sheet.addElement(row1)
    row1.addElement(CoveredTableCell(valuetype="float", value="42"))
    row1.addElement(TableCell(valuetype="float", value="200"))

    row2 = TableRow()
    sheet.addElement(row2)
    row2.addElement(CoveredTableCell(valuetype="float"))
    row2.addElement(TableCell(valuetype="float", value="300"))

    row3 = TableRow()
    sheet.addElement(row3)
    row3.addElement(CoveredTableCell(valuetype="void"))
    row3.addElement(TableCell(valuetype="float", value="400"))

    row4 = TableRow()
    sheet.addElement(row4)
    row4.addElement(CoveredTableCell(valuetype="boolean", booleanvalue="true"))
    row4.addElement(TableCell(valuetype="float", value="500"))

    doc.save(tmp_excel)

    result = pd.read_excel(tmp_excel, header=None)

    expected = pd.DataFrame(
        [
            [1, 100],
            [42, 200],
            [np.nan, 300],
            [np.nan, 400],
            [np.nan, 500],
        ]
    )

    tm.assert_frame_equal(result, expected)
