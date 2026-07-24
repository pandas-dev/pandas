from . import biff12
from collections import namedtuple

class Handler(object):
  def __init__(self):
    super(Handler, self).__init__()

  def read(self, reader, recid, reclen):
    if reclen > 0:
      reader.skip(reclen)


class BasicHandler(Handler):
  def __init__(self, name=None):
    super(BasicHandler, self).__init__()
    self.name = name

  def read(self, reader, recid, reclen):
    super(BasicHandler, self).read(reader, recid, reclen)
    return self.name


class StringTableHandler(Handler):
  cls = namedtuple('sst', ['count', 'uniqueCount'])

  def __init__(self):
    super(StringTableHandler, self).__init__()

  def read(self, reader, recid, reclen):
    count = reader.read_int()
    unique = reader.read_int()
    return self.cls._make([count, unique])


class StringInstanceHandler(Handler):
  cls = namedtuple('si', ['t'])

  def __init__(self):
    super(StringInstanceHandler, self).__init__()

  def read(self, reader, recid, reclen):
    reader.skip(1)
    val = reader.read_string()
    return self.cls._make([val])


class SheetHandler(Handler):
  cls = namedtuple('sheet', ['sheetId', 'rId', 'name'])

  def __init__(self):
    super(SheetHandler, self).__init__()

  def read(self, reader, recid, reclen):
    reader.skip(4)
    sheetid = reader.read_int()
    relid = reader.read_string()
    name = reader.read_string()
    return self.cls._make([sheetid, relid, name])


class DimensionHandler(Handler):
  cls = namedtuple('dimension', ['r', 'c', 'h', 'w'])

  def __init__(self):
    super(DimensionHandler, self).__init__()

  def read(self, reader, recid, reclen):
    r1 = reader.read_int()
    r2 = reader.read_int()
    c1 = reader.read_int()
    c2 = reader.read_int()
    return self.cls._make([r1, c1, r2 - r1 + 1, c2 - c1 + 1])


class ColumnHandler(Handler):
  cls = namedtuple('col', ['c1', 'c2', 'width', 'style'])

  def __init__(self):
    super(ColumnHandler, self).__init__()

  def read(self, reader, recid, reclen):
    c1 = reader.read_int()
    c2 = reader.read_int()
    width = reader.read_int() / 256
    style = reader.read_int()
    return self.cls._make([c1, c2, width, style])


class RowHandler(Handler):
  cls = namedtuple('row', ['r'])

  def __init__(self):
    super(RowHandler, self).__init__()

  def read(self, reader, recid, reclen):
    r = reader.read_int()
    return self.cls._make([r])


class CellHandler(Handler):
  cls = namedtuple('c', ['c', 'v', 'f', 'style'])

  def __init__(self):
    super(CellHandler, self).__init__()

  def read(self, reader, recid, reclen):
    col = reader.read_int()
    style = reader.read_int()
    val = None
    if recid == biff12.NUM:
      val = reader.read_float()
    elif recid == biff12.BOOLERR:
      val = hex(reader.read_byte())
    elif recid == biff12.BOOL:
      val = reader.read_byte() != 0
    elif recid == biff12.FLOAT:
      val = reader.read_double()
    elif recid == biff12.STRING:
      val = reader.read_int()
    elif recid == biff12.FORMULA_STRING:
      val = reader.read_string()
    elif recid == biff12.FORMULA_FLOAT:
      val = reader.read_double()
    elif recid == biff12.FORMULA_BOOL:
      val = reader.read_byte() != 0
    elif recid == biff12.FORMULA_BOOLERR:
      val = hex(reader.read_byte())
    return self.cls._make([col, val, None, style])


class HyperlinkHandler(Handler):
  cls = namedtuple('hyperlink', ['r', 'c', 'h', 'w', 'rId'])

  def __init__(self):
    super(HyperlinkHandler, self).__init__()

  def read(self, reader, recid, reclen):
    r1 = reader.read_int()
    r2 = reader.read_int()
    c1 = reader.read_int()
    c2 = reader.read_int()
    rId = reader.read_string()
    return self.cls._make([r1, c1, r2 - r1 + 1, c2 - c1 + 1, rId])
