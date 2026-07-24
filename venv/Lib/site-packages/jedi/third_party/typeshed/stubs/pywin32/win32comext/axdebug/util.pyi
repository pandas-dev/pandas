from _typeshed import Incomplete

import win32com.server.dispatcher

debugging: int

def trace(*args) -> None: ...
def RaiseNotImpl(who: Incomplete | None = ...) -> None: ...

class Dispatcher(win32com.server.dispatcher.DispatcherWin32trace):
    def __init__(self, policyClass, object) -> None: ...
