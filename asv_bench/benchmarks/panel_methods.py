from pandas_vb_common import *


class panel_pct_change_items(object):
    goal_time = 0.2

    def setup(self):
        self.index = date_range(start='2000', freq='D', periods=1000)
        self.panel = Panel(np.random.randn(100, len(self.index), 1000))

    def time_panel_pct_change_items(self):
        self.panel.pct_change(1, axis='items')


class panel_pct_change_major(object):
    goal_time = 0.2

    def setup(self):
        self.index = date_range(start='2000', freq='D', periods=1000)
        self.panel = Panel(np.random.randn(100, len(self.index), 1000))

    def time_panel_pct_change_major(self):
        self.panel.pct_change(1, axis='major')


class panel_pct_change_minor(object):
    goal_time = 0.2

    def setup(self):
        self.index = date_range(start='2000', freq='D', periods=1000)
        self.panel = Panel(np.random.randn(100, len(self.index), 1000))

    def time_panel_pct_change_minor(self):
        self.panel.pct_change(1, axis='minor')


class panel_shift(object):
    goal_time = 0.2

    def setup(self):
        self.index = date_range(start='2000', freq='D', periods=1000)
        self.panel = Panel(np.random.randn(100, len(self.index), 1000))

    def time_panel_shift(self):
        self.panel.shift(1)


class panel_shift_minor(object):
    goal_time = 0.2

    def setup(self):
        self.index = date_range(start='2000', freq='D', periods=1000)
        self.panel = Panel(np.random.randn(100, len(self.index), 1000))

    def time_panel_shift_minor(self):
        self.panel.shift(1, axis='minor')