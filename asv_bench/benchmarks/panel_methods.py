import warnings

import numpy as np
from pandas import Panel

from .pandas_vb_common import Panel


class PanelMethods(object):

    goal_time = 0.2
    params = ['items', 'major', 'minor']
    param_names = ['axis']

    def setup(self, axis):
        with warnings.catch_warnings(record=True):
            self.panel = Panel(np.random.randn(100, 1000, 100))

    def time_pct_change(self, axis):
        with warnings.catch_warnings(record=True):
            self.panel.pct_change(1, axis=axis)

    def time_shift(self, axis):
        with warnings.catch_warnings(record=True):
            self.panel.shift(1, axis=axis)


from .pandas_vb_common import setup  # noqa: F401
