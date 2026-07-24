import os
import sys
import xml.etree.ElementTree as ET
from . import biff12
from .reader import BIFF12Reader
from collections import namedtuple

if sys.version_info > (3,):
  xrange = range

Cell = namedtuple('Cell', ['r', 'c', 'v'])

class Worksheet(object):
  def __init__(self, name, fp, rels_fp=None, stringtable=None, debug=False):
    super(Worksheet, self).__init__()
    self.name = name
    self._reader = BIFF12Reader(fp=fp, debug=debug)
    self._rels_fp = rels_fp
    self._rels = ET.parse(rels_fp).getroot() if rels_fp is not None else None
    self._stringtable = stringtable
    self._data_offset = 0
    self.dimension = None
    self.cols = []
    self.rels = {}
    self.hyperlinks = {}
    self._parse()

  def __enter__(self):
    return self

  def __exit__(self, type, value, traceback):
    self.close()

  def __iter__(self):
    return self.rows()

  def _parse(self):
    if self._rels is not None:
      for el in self._rels:
        self.rels[el.attrib['Id']] = el.attrib['Target']

    for item in self._reader:
      if item[0] == biff12.DIMENSION:
        self.dimension = item[1]
      elif item[0] == biff12.COL:
        self.cols.append(item[1])
      elif item[0] == biff12.SHEETDATA:
        self._data_offset = self._reader.tell()
        if self._rels is None:
          break
      elif item[0] == biff12.HYPERLINK and self._rels is not None:
        for r in xrange(item[1].h):
          for c in xrange(item[1].w):
            self.hyperlinks[item[1].r + r, item[1].c + c] = item[1].rId

  def rows(self, sparse=False):
    self._reader.seek(self._data_offset, os.SEEK_SET)
    row_num = -1
    row = None
    for item in self._reader:
      if item[0] == biff12.ROW and item[1].r != row_num:
        if row is not None:
          yield row
        if not sparse:
          while row_num < item[1].r - 1:
            row_num += 1
            yield [Cell(row_num, i, None) for i in xrange(self.dimension.c + self.dimension.w)]
        row_num = item[1].r
        row = [Cell(row_num, i, None) for i in xrange(self.dimension.c + self.dimension.w)]
      elif item[0] >= biff12.BLANK and item[0] <= biff12.FORMULA_BOOLERR:
        if item[0] == biff12.STRING and self._stringtable is not None:
          row[item[1].c] = Cell(row_num, item[1].c, self._stringtable[item[1].v])
        else:
          row[item[1].c] = Cell(row_num, item[1].c, item[1].v)
      elif item[0] == biff12.SHEETDATA_END:
        if row is not None:
          yield row
        break

  def close(self):
    self._reader.close()
    if self._rels_fp is not None:
      self._rels_fp.close()
