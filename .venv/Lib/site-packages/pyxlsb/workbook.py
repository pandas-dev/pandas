import os
import sys
import xml.etree.ElementTree as ET
from . import biff12
from .reader import BIFF12Reader
from .stringtable import StringTable
from .worksheet import Worksheet
from tempfile import TemporaryFile

if sys.version_info > (3,):
  basestring = (str, bytes)

class Workbook(object):
  def __init__(self, fp, debug=False):
    super(Workbook, self).__init__()
    self._zf = fp
    self._debug = debug
    self._sheets = []
    self.stringtable = None
    self._parse()

  def __enter__(self):
    return self

  def __exit__(self, type, value, traceback):
    self.close()

  @property
  def sheets(self):
    return [v[0] for v in self._sheets]

  def _parse(self):
    rels = {}
    with self._zf.open('xl/_rels/workbook.bin.rels', 'r') as zf:
      for el in ET.parse(zf).getroot():
        rels[el.attrib['Id']] = el.attrib['Target']

    with TemporaryFile() as temp:
      with self._zf.open('xl/workbook.bin', 'r') as zf:
        temp.write(zf.read())
        temp.seek(0, os.SEEK_SET)
      reader = BIFF12Reader(fp=temp, debug=self._debug)
      for item in reader:
        if item[0] == biff12.SHEET:
          self._sheets.append((item[1].name, rels[item[1].rId]))
        elif item[0] == biff12.SHEETS_END:
          break

    try:
      temp = TemporaryFile()
      with self._zf.open('xl/sharedStrings.bin', 'r') as zf:
        temp.write(zf.read())
        temp.seek(0, os.SEEK_SET)
      self.stringtable = StringTable(fp=temp)
    except KeyError:
      temp.close()
    except Exception:
      temp.close()
      raise

  def get_sheet(self, idx, rels=False):
    if isinstance(idx, basestring):
      idx = [s.lower() for s, _ in self._sheets].index(idx.lower()) + 1
    if idx < 1 or idx > len(self._sheets):
      raise IndexError('sheet index out of range')

    name = self._sheets[idx - 1][0]
    target = self._sheets[idx - 1][1].split('/')

    temp = TemporaryFile()
    with self._zf.open('xl/{}/{}'.format(target[0], target[-1]), 'r') as zf:
      temp.write(zf.read())
      temp.seek(0, os.SEEK_SET)

    if rels:
      rels_temp = TemporaryFile()
      with self._zf.open('xl/{}/_rels/{}.rels'.format(target[0], target[-1]), 'r') as zf:
        rels_temp.write(zf.read())
        rels_temp.seek(0, os.SEEK_SET)
    else:
      rels_temp = None

    return Worksheet(name=name, fp=temp, rels_fp=rels_temp, stringtable=self.stringtable, debug=self._debug)

  def close(self):
    self._zf.close()
    if self.stringtable is not None:
      self.stringtable.close()
