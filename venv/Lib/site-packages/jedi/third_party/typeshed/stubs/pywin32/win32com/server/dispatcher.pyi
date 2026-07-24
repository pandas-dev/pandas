from logging import Logger
from typing_extensions import TypeAlias

from win32com.server.policy import BasicWrapPolicy

class DispatcherBase:
    policy: BasicWrapPolicy
    logger: Logger
    def __init__(self, policyClass, object) -> None: ...

class DispatcherTrace(DispatcherBase): ...

class DispatcherWin32trace(DispatcherTrace):
    def __init__(self, policyClass, object) -> None: ...

class DispatcherOutputDebugString(DispatcherTrace): ...

DefaultDebugDispatcher: TypeAlias = DispatcherTrace
