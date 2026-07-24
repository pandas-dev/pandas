from _typeshed import Incomplete
from typing import Final
from typing_extensions import LiteralString

__version__: Final[str]

def asNative(s) -> str: ...

EXTENSIONS: Final = [".ttf", ".ttc", ".otf", ".pfb", ".pfa"]
FF_FIXED: Final = 1
FF_SERIF: Final = 2
FF_SYMBOLIC: Final = 4
FF_SCRIPT: Final = 8
FF_NONSYMBOLIC: Final = 32
FF_ITALIC: Final = 64
FF_ALLCAP: Final = 65536
FF_SMALLCAP: Final = 131072
FF_FORCEBOLD: Final = 262144

class FontDescriptor:
    name: Incomplete
    fullName: Incomplete
    familyName: Incomplete
    styleName: Incomplete
    isBold: bool
    isItalic: bool
    isFixedPitch: bool
    isSymbolic: bool
    typeCode: Incomplete
    fileName: Incomplete
    metricsFileName: Incomplete
    timeModified: int
    def __init__(self) -> None: ...
    def getTag(self) -> LiteralString: ...

class FontFinder:
    useCache: Incomplete
    validate: Incomplete
    verbose: Incomplete
    def __init__(
        self, dirs=[], useCache: bool = True, validate: bool = False, recur: bool = False, fsEncoding=None, verbose: int = 0
    ) -> None: ...
    def addDirectory(self, dirName, recur=None) -> None: ...
    def addDirectories(self, dirNames, recur=None) -> None: ...
    def getFamilyNames(self) -> list[bytes]: ...
    def getFontsInFamily(self, familyName): ...
    def getFamilyXmlReport(self) -> LiteralString: ...
    def getFontsWithAttributes(self, **kwds) -> list[FontDescriptor]: ...
    def getFont(self, familyName, bold: bool = False, italic: bool = False) -> FontDescriptor: ...
    def save(self, fileName) -> None: ...
    def load(self, fileName) -> None: ...
    def search(self) -> None: ...

def test() -> None: ...
