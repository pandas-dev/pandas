class PandasError(Exception):
    pass


class DataError(PandasError):
    pass


class AmbiguousIndexError(PandasError, KeyError):
    pass


class RemoteDataError(DataError):
    pass
