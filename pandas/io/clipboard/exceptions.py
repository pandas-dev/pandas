import ctypes


class PyperclipException(RuntimeError):
    pass


class PyperclipWindowsException(PyperclipException):
    def __init__(self, message):
        # error: Module has no attribute "WinError"
        message += " ({err})".format(
            err=ctypes.WinError()  # type: ignore[attr-defined]
        )
        super().__init__(message)
