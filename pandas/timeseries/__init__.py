"""TimeSeries

:author: Pierre GF Gerard-Marchant & Matt Knox
:contact: pierregm_at_uga_dot_edu - mattknox_ca_at_hotmail_dot_com
:version: $Id$
"""


__author__ = "Pierre GF Gerard-Marchant  & Matt Knox ($Author$)"
__revision__ = "$Revision$"
__date__     = '$Date$'

import const
import tdates
from tdates import *
import tseries
from tseries import *
import trecords
from trecords import *
_c = const
from extras import tsfromtxt, guess_freq

from scikits.timeseries.version import __version__

__all__ = [
    '_c', 'const', 'tdates','tseries','trecords', 'tsfromtxt', 'guess_freq']
__all__.extend(tdates.__all__)
__all__.extend(tseries.__all__)
__all__.extend(trecords.__all__)

from numpy.testing import Tester
test = Tester().test
