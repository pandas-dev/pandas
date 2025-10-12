import io
import pytest
import pandas as pd

pytest.importorskip("xlsxwriter")
openpyxl = pytest.importorskip("openpyxl")


def test_to_excel_xlsxwriter_autofilter_and_bold():
    df = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="xlsxwriter") as writer:
        df.to_excel(
            writer,
            index=False,
            engine_kwargs={"autofilter_header": True, "header_bold": True},
        )
    buf.seek(0)
    wb = openpyxl.load_workbook(buf)
    ws = wb.active
    # Autofilter should be set spanning header+data
    assert ws.auto_filter is not None
    assert ws.auto_filter.ref is not None and ws.auto_filter.ref != ""
    # Header row (row 1) should be bold
    assert all(ws.cell(row=1, column=c).font.bold for c in range(1, df.shape[1] + 1))
