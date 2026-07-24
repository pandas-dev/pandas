from _typeshed import FileDescriptorOrPath, Unused
from collections.abc import Iterable
from typing import Final, type_check_only

RequirePyRXP: int
simpleparse: int

@type_check_only
class _smartDecode:
    @staticmethod
    def __call__(s: str | bytes | bytearray) -> str: ...

smartDecode: _smartDecode
NONAME: Final = ""
NAMEKEY: Final = 0
CONTENTSKEY: Final = 1
CDATAMARKER: Final = "<![CDATA["
LENCDATAMARKER: Final = 9
CDATAENDMARKER: Final = "]]>"
replacelist: list[tuple[str, str]]

def unEscapeContentList(contentList: Iterable[str]) -> list[str]: ...
def parsexmlSimple(xmltext, oneOutermostTag: int = 0, eoCB: Unused = None, entityReplacer=...): ...

parsexml = parsexmlSimple

def parseFile(filename: FileDescriptorOrPath): ...

verbose: int

def skip_prologue(text, cursor): ...
def parsexml0(xmltext, startingat: int = 0, toplevel: int = 1, entityReplacer=...): ...
def pprettyprint(parsedxml): ...
def testparse(s, dump: int = 0) -> None: ...
def test(dump: int = 0) -> None: ...
