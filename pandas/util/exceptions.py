class PandasError(Exception):
    pass


class DataError(PandasError):
    pass


class DataConflictError(DataError):
    pass


class RemoteDataError(DataError):
    pass
