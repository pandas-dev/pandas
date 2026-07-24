from antlr4.InputStream import InputStream as InputStream

class FileStream(InputStream):
    __slots__ = "fileName"
    fileName: str
    def __init__(self, fileName: str, encoding: str = "ascii", errors: str = "strict") -> None: ...
    def readDataFrom(self, fileName: str, encoding: str, errors: str = "strict") -> str: ...
