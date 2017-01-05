from .pandas_vb_common import *


class Constructors1(object):
    goal_time = 0.2

    def setup(self):
        self.data_frames = {}
        self.start = datetime(1990, 1, 1)
        self.end = datetime(2012, 1, 1)
        for x in range(100):
            self.end += timedelta(days=1)
            self.dr = np.asarray(date_range(self.start, self.end))
            self.df = DataFrame({'a': ([0] * len(self.dr)), 'b': ([1] * len(self.dr)), 'c': ([2] * len(self.dr)), }, index=self.dr)
            self.data_frames[x] = self.df

    def time_panel_from_dict_all_different_indexes(self):
        Panel.from_dict(self.data_frames)


class Constructors2(object):
    goal_time = 0.2

    def setup(self):
        self.data_frames = {}
        for x in range(100):
            self.dr = np.asarray(DatetimeIndex(start=datetime(1990, 1, 1), end=datetime(2012, 1, 1), freq=datetools.Day(1)))
            self.df = DataFrame({'a': ([0] * len(self.dr)), 'b': ([1] * len(self.dr)), 'c': ([2] * len(self.dr)), }, index=self.dr)
            self.data_frames[x] = self.df

    def time_panel_from_dict_equiv_indexes(self):
        Panel.from_dict(self.data_frames)


class Constructors3(object):
    goal_time = 0.2

    def setup(self):
        self.dr = np.asarray(DatetimeIndex(start=datetime(1990, 1, 1), end=datetime(2012, 1, 1), freq=datetools.Day(1)))
        self.data_frames = {}
        for x in range(100):
            self.df = DataFrame({'a': ([0] * len(self.dr)), 'b': ([1] * len(self.dr)), 'c': ([2] * len(self.dr)), }, index=self.dr)
            self.data_frames[x] = self.df

    def time_panel_from_dict_same_index(self):
        Panel.from_dict(self.data_frames)


class Constructors4(object):
    goal_time = 0.2

    def setup(self):
        self.data_frames = {}
        self.start = datetime(1990, 1, 1)
        self.end = datetime(2012, 1, 1)
        for x in range(100):
            if (x == 50):
                self.end += timedelta(days=1)
            self.dr = np.asarray(date_range(self.start, self.end))
            self.df = DataFrame({'a': ([0] * len(self.dr)), 'b': ([1] * len(self.dr)), 'c': ([2] * len(self.dr)), }, index=self.dr)
            self.data_frames[x] = self.df

    def time_panel_from_dict_two_different_indexes(self):
        Panel.from_dict(self.data_frames)
