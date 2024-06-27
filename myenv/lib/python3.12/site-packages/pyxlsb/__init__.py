from .handlers import Handler
from .reader import BIFF12Reader
from .workbook import Workbook
from .worksheet import Worksheet

__version__ = '1.0.10'

def open_workbook(name, debug=False):
  from zipfile import ZipFile
  zf = ZipFile(name, 'r')
  return Workbook(fp=zf, debug=debug)

def convert_date(date):
  if not isinstance(date, int) and not isinstance(date, float):
    return None

  from datetime import datetime, timedelta
  if int(date) == 0:
    return datetime(1900, 1, 1, 0, 0, 0) + timedelta(seconds=round(date * 24 * 60 * 60))
  elif int(date) >= 61:
    # According to Lotus 1-2-3, Feb 29th 1900 is a real thing, therefore we have to remove one day after that date
    return datetime(1899, 12, 31, 0, 0, 0) + timedelta(days=int(date) - 1, seconds=round((date % 1) * 24 * 60 * 60))
  else:
    # Feb 29th 1900 will show up as Mar 1st 1900 because Python won't handle that date
    return datetime(1899, 12, 31, 0, 0, 0) + timedelta(days=int(date), seconds=round((date % 1) * 24 * 60 * 60))
