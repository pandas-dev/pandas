import io
import os
import struct
from . import biff12
from .handlers import *

uint8_t = struct.Struct('<B')
uint16_t = struct.Struct('<H')
int32_t = struct.Struct('<i')
uint32_t = struct.Struct('<I')
double_t = struct.Struct('<d')

class RecordReader(object):
  def __init__(self, buf, enc='utf-16'):
    self._fp = io.BytesIO(buf)
    self._enc = enc

  def __enter__(self):
    return self

  def __exit__(self, type, value, traceback):
    self._fp.close()

  def tell(self):
    return self._fp.tell()

  def seek(self, offset, whence=os.SEEK_SET):
    self._fp.seek(offset, whence)

  def skip(self, size):
    self._fp.seek(size, os.SEEK_CUR)

  def read(self, size):
    return self._fp.read(size)

  def read_int(self):
    buff = self._fp.read(4)
    if len(buff) < 4:
      return None
    return uint32_t.unpack(buff)[0]

  def read_short(self):
    buff = self._fp.read(2)
    if len(buff) < 2:
      return None
    return uint16_t.unpack(buff)[0]

  def read_byte(self):
    byte = self._fp.read(1)
    if not byte:
      return None
    return uint8_t.unpack(byte)[0]

  def read_float(self):
    buff = self._fp.read(4)
    if len(buff) < 4:
      return None
    v = 0.0
    intval = int32_t.unpack(buff)[0]
    if intval & 0x02 != 0:
      v = float(intval >> 2)
    else:
      v = double_t.unpack(b'\x00\x00\x00\x00' + uint32_t.pack(intval & 0xFFFFFFFC))[0]
    if intval & 0x01 != 0:
      v /= 100
    return v

  def read_double(self):
    buff = self._fp.read(8)
    if len(buff) < 8:
      return None
    return double_t.unpack(buff)[0]

  def read_string(self):
    l = self.read_int()
    if l is None:
      return None
    buff = self.read(l * 2)
    if len(buff) < l * 2:
      return None
    return buff.decode(self._enc, errors='replace')


class BIFF12Reader(object):
  handlers = {
    # Workbook part handlers
    biff12.WORKBOOK:   BasicHandler('workbook'),
    biff12.SHEETS:     BasicHandler('sheets'),
    biff12.SHEETS_END: BasicHandler('/sheets'),
    biff12.SHEET:      SheetHandler(),

    # SharedStrings part handlers
    biff12.SST:     StringTableHandler(),
    biff12.SST_END: BasicHandler('/sst'),
    biff12.SI:      StringInstanceHandler(),

    # Worksheet part handlers
    biff12.WORKSHEET:       BasicHandler('worksheet'),
    biff12.WORKSHEET_END:   BasicHandler('/worksheet'),
    biff12.DIMENSION:       DimensionHandler(),
    biff12.SHEETDATA:       BasicHandler('sheetData'),
    biff12.SHEETDATA_END:   BasicHandler('/sheetData'),
    biff12.COLS:            BasicHandler('cols'),
    biff12.COLS_END:        BasicHandler('/cols'),
    biff12.COL:             ColumnHandler(),
    biff12.ROW:             RowHandler(),
    biff12.BLANK:           CellHandler(),
    biff12.NUM:             CellHandler(),
    biff12.BOOLERR:         CellHandler(),
    biff12.BOOL:            CellHandler(),
    biff12.FLOAT:           CellHandler(),
    biff12.STRING:          CellHandler(),
    biff12.FORMULA_STRING:  CellHandler(),
    biff12.FORMULA_FLOAT:   CellHandler(),
    biff12.FORMULA_BOOL:    CellHandler(),
    biff12.FORMULA_BOOLERR: CellHandler(),
    biff12.HYPERLINK:       HyperlinkHandler()
  }

  def __init__(self, fp, debug=False):
    super(BIFF12Reader, self).__init__()
    self._debug = debug
    self._fp = fp

  def __iter__(self):
    return self

  def __next__(self):
    return self.next()

  def __enter__(self):
    return self

  def __exit__(self, type, value, traceback):
    self.close()

  def tell(self):
    return self._fp.tell()

  def seek(self, offset, whence=os.SEEK_SET):
    self._fp.seek(offset, whence)

  def read_id(self):
    v = 0
    for i in range(4):
      byte = self._fp.read(1)
      if not byte:
        return None
      byte = uint8_t.unpack(byte)[0]
      v += byte << 8 * i
      if byte & 0x80 == 0:
        break
    return v

  def read_len(self):
    v = 0
    for i in range(4):
      byte = self._fp.read(1)
      if not byte:
        return None
      byte = uint8_t.unpack(byte)[0]
      v += (byte & 0x7F) << (7 * i)
      if byte & 0x80 == 0:
        break
    return v

  def register_handler(self, recid, handler):
    self.handlers[recid] = handler

  def next(self):
    ret = None
    while ret is None:
      if self._debug:
        pos = self._fp.tell()
      recid = self.read_id()
      reclen = self.read_len()
      if recid is None or reclen is None:
        raise StopIteration
      recdata = self._fp.read(reclen)
      with RecordReader(recdata) as reader:
        ret = (self.handlers.get(recid) or Handler()).read(reader, recid, reclen)
      if self._debug:
        print('{:08X}  {:04X}  {:<6} {} {}'.format(pos, recid, reclen, ' '.join('{:02X}'.format(b) for b in recdata), ret))
    return (recid, ret)

  def close(self):
    self._fp.close()
