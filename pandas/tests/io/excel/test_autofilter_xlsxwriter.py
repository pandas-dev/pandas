import io
import zipfile

import pytest

import pandas as pd

pytest.importorskip("xlsxwriter")
openpyxl = pytest.importorskip("openpyxl")


def test_to_excel_xlsxwriter_autofilter():
    df = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="xlsxwriter") as writer:
        # Test autofilter
        df.to_excel(writer, index=False, autofilter=True)
    buf.seek(0)
    with zipfile.ZipFile(buf) as zf:
        with zf.open("xl/worksheets/sheet1.xml") as f:
            sheet = f.read().decode("utf-8")
    # Check for autofilter
    assert '<autoFilter ref="A1:B3"' in sheet


def test_to_excel_xlsxwriter_styler_bold_header():
    df = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="xlsxwriter") as writer:
        # Test header bold using Styler
        df.style.set_properties(
            **{"font-weight": "bold"}, subset=pd.IndexSlice[0, :]
        ).to_excel(writer, index=False)
    buf.seek(0)
    with zipfile.ZipFile(buf) as zf:
        # Check styles.xml for the bold font definition
        with zf.open("xl/styles.xml") as f:
            styles = f.read().decode("utf-8")
            print("\n===== STYLES XML =====")
            print(styles)
            print("===========================\n")
        # Check sheet1.xml for style references
        with zf.open("xl/worksheets/sheet1.xml") as f:
            sheet = f.read().decode("utf-8")
            print("\n===== SHEET1 XML =====")
            print(sheet)
            print("===========================\n")
    # Check for bold style in the styles.xml and its reference in sheet1.xml
    assert "<b/>" in styles, "Bold style not found in styles.xml"
    # Check that the header row (first row) uses a style with bold
    assert 'r="1"' in sheet, "Header row not found in sheet1.xml"
