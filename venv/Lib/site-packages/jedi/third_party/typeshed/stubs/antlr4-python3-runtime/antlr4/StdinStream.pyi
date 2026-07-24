from antlr4.InputStream import InputStream as InputStream

class StdinStream(InputStream):
    def __init__(self, encoding: str = "ascii", errors: str = "strict") -> None: ...
