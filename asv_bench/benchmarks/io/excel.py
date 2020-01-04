from io import BytesIO

import numpy as np
from odf.opendocument import OpenDocumentSpreadsheet
from odf.table import Table, TableCell, TableRow
from odf.text import P

from pandas import DataFrame, ExcelWriter, date_range, read_excel
import pandas._testing as tm


def _generate_dataframe():
    N = 2000
    C = 5
    df = DataFrame(
        np.random.randn(N, C),
        columns=[f"float{i}" for i in range(C)],
        index=date_range("20000101", periods=N, freq="H"),
    )
    df["object"] = tm.makeStringIndex(N)
    return df


class WriteExcel:

    params = ["openpyxl", "xlsxwriter", "xlwt"]
    param_names = ["engine"]

    def setup(self, engine):
        self.df = _generate_dataframe()

    def time_write_excel(self, engine):
        bio = BytesIO()
        bio.seek(0)
        writer = ExcelWriter(bio, engine=engine)
        self.df.to_excel(writer, sheet_name="Sheet1")
        writer.save()


class ReadExcel:

    params = ["xlrd", "openpyxl", "odf"]
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
        fname = self.fname_odf if engine == "odf" else self.fname_excel
        read_excel(fname, engine=engine)


from ..pandas_vb_common import setup  # noqa: F401 isort:skip
