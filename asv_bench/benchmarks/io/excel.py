from io import BytesIO
import numpy as np
from pandas import DataFrame, date_range, ExcelWriter, read_excel
import pandas.util.testing as tm
from odf.opendocument import OpenDocumentSpreadsheet
from odf.text import P
from odf.table import Table, TableRow, TableCell


def _generate_dataframe():
    N = 2000
    C = 5
    df = DataFrame(
        np.random.randn(N, C),
        columns=["float{}".format(i) for i in range(C)],
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
        bio_write = BytesIO()
        bio_write.seek(0)
        writer_write = ExcelWriter(bio_write, engine=engine)
        self.df.to_excel(writer_write, sheet_name="Sheet1")
        writer_write.save()


class ReadExcel:

    params = ["xlrd", "openpyxl", "odf"]
    param_names = ["engine"]

    def _generate_odf(self):
        doc = OpenDocumentSpreadsheet()
        table = Table(name="Table1")
        for row in self.df.values:
            tr = TableRow()
            for val in row:
                tc = TableCell(valuetype='string')
                tc.addElement(P(text=val))
                tr.addElement(tc)
            table.addElement(tr)

        doc.spreadsheet.addElement(table)

        return doc

    def setup(self, engine):
        self.df = _generate_dataframe()

        self.bio_read = BytesIO()
        self.writer_read = ExcelWriter(self.bio_read)
        self.df.to_excel(self.writer_read, sheet_name="Sheet1")
        self.writer_read.save()
        self.bio_read.seek(0)

        self.bio_read_odf = BytesIO()
        odf_doc = self._generate_odf()
        odf_doc.write(self.bio_read_odf)
        self.bio_read_odf.seek(0)

    def time_read_excel(self, engine):
        bio = self.bio_read_odf if engine == "odf" else self.bio_read
        read_excel(bio, engine=engine)


from ..pandas_vb_common import setup  # noqa: F401
