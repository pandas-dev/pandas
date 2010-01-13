# pylint: disable-msg=W0614,W0401,W0611,W0622

__docformat__ = 'restructuredtext'

from datetime import datetime

import numpy as np

from pandas.version import __version__
from pandas.info import __doc__

from pandas.core.api import *
from pandas.io.parsers import parseCSV, parseText, parseExcel
from pandas.stats.api import *
