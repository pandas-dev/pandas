from io import BytesIO

import numpy as np
from odf.opendocument import OpenDocumentSpreadsheet
from odf.table import (
    Table,
    TableCell,
    TableRow,
)
from odf.text import P

from pandas import (
    DataFrame,
    ExcelWriter,
    Index,
    date_range,
    read_excel,
)


def _generate_dataframe():
    N = 2000
    C = 5
    df = DataFrame(
        np.random.randn(N, C),
        columns=[f"float{i}" for i in range(C)],
        index=date_range("20000101", periods=N, freq="h"),
    )
    df["object"] = Index([f"i-{i}" for i in range(N)], dtype=object)
    return df


class WriteExcel:
    params = ["openpyxl", "xlsxwriter"]
    param_names = ["engine"]

    def setup(self, engine):
        self.df = _generate_dataframe()

    def time_write_excel(self, engine):
        bio = BytesIO()
        bio.seek(0)
        with ExcelWriter(bio, engine=engine) as writer:
            self.df.to_excel(writer, sheet_name="Sheet1")


class WriteExcelStyled:
    params = ["openpyxl", "xlsxwriter"]
    param_names = ["engine"]

    def setup(self, engine):
        self.df = _generate_dataframe()

    def time_write_excel_style(self, engine):
        bio = BytesIO()
        bio.seek(0)
        with ExcelWriter(bio, engine=engine) as writer:
            df_style = self.df.style
            df_style.map(lambda x: "border: red 1px solid;")
            df_style.map(lambda x: "color: blue")
            df_style.map(lambda x: "border-color: green black", subset=["float1"])
            df_style.to_excel(writer, sheet_name="Sheet1")


class ReadExcel:
    params = ["openpyxl", "odf"]
    param_names = ["engine"]
    fname_excel = "spreadsheet.xlsx"
    fname_odf = "spreadsheet.ods"

    def _create_odf(self):
        doc = OpenDocumentSpreadsheet()
        table = Table(name="Table1")
        for row in self.df.values:
            tr = TableRow()
            for val in row:
                tc = TableCell(valuetype="string")
                tc.addElement(P(text=val))
                tr.addElement(tc)
            table.addElement(tr)

        doc.spreadsheet.addElement(table)
        doc.save(self.fname_odf)

    def setup_cache(self):
        self.df = _generate_dataframe()

        self.df.to_excel(self.fname_excel, sheet_name="Sheet1")
        self._create_odf()

    def time_read_excel(self, engine):
        if engine == "odf":
            fname = self.fname_odf
        else:
            fname = self.fname_excel
        read_excel(fname, engine=engine)


class ReadExcelNRows(ReadExcel):
    def time_read_excel(self, engine):
        if engine == "odf":
            fname = self.fname_odf
        else:
            fname = self.fname_excel
        read_excel(fname, engine=engine, nrows=10)


from ..pandas_vb_common import setup  # noqa: F401 isort:skip
