from types import TracebackType

from win32com.server.exception import COMException
from win32comext.axscript.client.debug import DebugManager
from win32comext.axscript.client.framework import AXScriptCodeBlock, COMScript
from win32comext.axscript.server.axsite import AXSite

debugging: int

def FormatForAX(text: str) -> str: ...
def ExpandTabs(text: str) -> str: ...
def AddCR(text: str) -> str: ...

class IActiveScriptError:
    def GetSourceLineText(self) -> str | None: ...
    def GetSourcePosition(self) -> tuple[int, int, int]: ...
    def GetExceptionInfo(self) -> AXScriptException: ...

class AXScriptException(COMException):
    sourceContext: int
    startLineNo: int
    linetext: str
    def __init__(
        self,
        site: COMScript,
        codeBlock: AXScriptCodeBlock | None,
        exc_type: None = None,
        exc_value: BaseException | None = None,
        exc_traceback: None = None,
    ) -> None: ...
    def ExtractTracebackInfo(self, tb: TracebackType, site: COMScript) -> tuple[str, int, str, str | None]: ...

def ProcessAXScriptException(
    scriptingSite: AXSite, debugManager: DebugManager, exceptionInstance: AXScriptException
) -> None | COMException | AXScriptException: ...
