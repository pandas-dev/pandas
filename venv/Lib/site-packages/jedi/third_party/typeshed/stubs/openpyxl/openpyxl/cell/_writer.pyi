from _typeshed import Unused

from openpyxl.cell import _CellOrMergedCell

def etree_write_cell(xf, worksheet: Unused, cell: _CellOrMergedCell, styled=None) -> None: ...
def lxml_write_cell(xf, worksheet: Unused, cell: _CellOrMergedCell, styled: bool = False) -> None: ...

write_cell = lxml_write_cell
write_cell = etree_write_cell
