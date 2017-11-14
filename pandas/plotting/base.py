"""
Base module for plotting engines to override and register
with pandas.
"""
from pandas.core.base import PandasObject

engines = {}


def register_engine(name, kind, engine):
    # XXX: get rid of the kind parameter
    global engines
    engines[(name, kind)] = engine


def deregister_engine(name, kind):
    # XXX: get rid of the kind parameter
    global engines
    engines.pop((name, kind))


def get_engine(kind):
    # XXX: get rid of the kind parameter
    from pandas import get_option

    active = get_option('plotting.engine')
    if active == 'auto':
        active = 'matplotlib'

    return engines[(active, kind)]


class Dispatcher(object):

    def __init__(self, data):
        self._data = data

    def __call__(self, *args, **kwargs):
        kind = 'frame' if self._data.ndim == 2 else 'series'
        engine = get_engine(kind)
        return engine(self._data)(*args, **kwargs)

    def __getattribute__(self, name):
        if name == '_data':
            return object.__getattribute__(self, name)
        kind = 'frame' if self._data.ndim == 2 else 'series'

        engine = get_engine(kind)(self._data)
        return getattr(engine, name)


class BasePlotMethods(PandasObject):

    def __init__(self, data):
        self._data = data

    def __call__(self, *args, **kwargs):
        """Make a plot"""
        raise NotImplementedError

    def area(self, **kwargs):
        raise NotImplementedError("This backend doesn't support this method")

    def bar(self, **kwargs):
        raise NotImplementedError("This backend doesn't support this method")

    def barh(self, **kwargs):
        raise NotImplementedError("This backend doesn't support this method")

    def box(self, **kwargs):
        raise NotImplementedError("This backend doesn't support this method")

    def density(self, **kwargs):
        raise NotImplementedError("This backend doesn't support this method")

    def hist(self, **kwargs):
        raise NotImplementedError("This backend doesn't support this method")

    def line(self, **kwargs):
        raise NotImplementedError("This backend doesn't support this method")

    def pie(self, **kwargs):
        raise NotImplementedError("This backend doesn't support this method")


class SeriesPlotMethods(BasePlotMethods):
    pass


class FramePlotMethods(BasePlotMethods):

    def hexbin(self, x, y, **kwargs):
        raise NotImplementedError("This backend doesn't support this method")

    def scatter(self, x, y, **kwargs):
        raise NotImplementedError("This backend doesn't support this method")
