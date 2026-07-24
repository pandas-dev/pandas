from openpyxl.worksheet.worksheet import Worksheet

class WorksheetCopy:
    source: Worksheet
    target: Worksheet
    def __init__(self, source_worksheet: Worksheet, target_worksheet: Worksheet) -> None: ...
    def copy_worksheet(self) -> None: ...
