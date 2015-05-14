
from pandas.core.index import Index


class Interval(object):
    """
    Represents an interval of time defined by two timestamps
    """

    def __init__(self, start, end):
        self.start = start
        self.end = end


class PeriodInterval(object):
    """
    Represents an interval of time defined by two Period objects (time ordinals)
    """

    def __init__(self, start, end):
        self.start = start
        self.end = end


class IntervalIndex(Index):
    """

    """
    def __new__(self, starts, ends):
        pass

    def dtype(self):
        return self.values.dtype

if __name__ == '__main__':
    pass
