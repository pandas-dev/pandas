import pythoncom

class COMException(pythoncom.com_error):
    scode: int
    description: str
    source: str | None
    helpfile: str | None
    helpcontext: int | None
    def __init__(
        self,
        description: str | None = None,
        scode: int | None = None,
        source: str | None = None,
        helpfile: str | None = None,
        helpContext: int | None = None,
        desc: str | None = None,
        hresult: int | None = None,
    ) -> None: ...

def IsCOMException(t: type[BaseException] | None = None) -> bool: ...
def IsCOMServerException(t: type[BaseException] | None = None) -> bool: ...
