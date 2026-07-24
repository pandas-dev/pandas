import os
from _typeshed import Incomplete
from typing import Final

from reportlab.pdfbase.pdfdoc import PDFObject
from reportlab.platypus.flowables import Flowable

__version__: Final[str]

def xorKey(num: int, key: bytes) -> bytes: ...

CLOBBERID: int
CLOBBERPERMISSIONS: int
DEBUG: int
reserved1: int
reserved2: int
printable: int
modifiable: int
copypastable: int
annotatable: int
higherbits: int

os_urandom = os.urandom

class StandardEncryption:
    prepared: int
    userPassword: Incomplete
    ownerPassword: Incomplete
    revision: int
    canPrint: Incomplete
    canModify: Incomplete
    canCopy: Incomplete
    canAnnotate: Incomplete
    O: Incomplete
    def __init__(
        self,
        userPassword,
        ownerPassword=None,
        canPrint: int = 1,
        canModify: int = 1,
        canCopy: int = 1,
        canAnnotate: int = 1,
        strength=None,
    ) -> None: ...
    def setAllPermissions(self, value) -> None: ...
    def permissionBits(self) -> int: ...
    def encode(self, t) -> bytes: ...
    P: Incomplete
    key: Incomplete
    U: Incomplete
    UE: Incomplete
    OE: Incomplete
    Perms: Incomplete
    objnum: Incomplete
    def prepare(self, document, overrideID=None) -> None: ...
    version: Incomplete
    def register(self, objnum, version) -> None: ...
    def info(self) -> StandardEncryptionDictionary: ...

class StandardEncryptionDictionary(PDFObject):
    __RefOnly__: int
    revision: Incomplete
    def __init__(self, O, OE, U, UE, P, Perms, revision) -> None: ...
    def format(self, document) -> bytes: ...

padding: str

def hexText(text: str | bytes | bytearray) -> str: ...
def unHexText(hexText) -> bytes: ...

PadString: bytes

def checkRevision(revision): ...
def encryptionkey(password, OwnerKey, Permissions, FileId1, revision=None) -> bytes: ...
def computeO(userPassword, ownerPassword, revision) -> bytes: ...
def computeU(
    encryptionkey,
    encodestring=b"(\xbfN^Nu\x8aAd\x00NV\xff\xfa\x01\x08..\x00\xb6\xd0h>\x80/\x0c\xa9\xfedSiz",
    revision=None,
    documentId=None,
) -> bytes: ...
def checkU(encryptionkey, U) -> None: ...
def encodePDF(key, objectNumber, generationNumber, string, revision=None) -> bytes: ...
def equalityCheck(observed, expected, label) -> None: ...
def test() -> None: ...
def encryptCanvas(
    canvas,
    userPassword,
    ownerPassword=None,
    canPrint: int = 1,
    canModify: int = 1,
    canCopy: int = 1,
    canAnnotate: int = 1,
    strength: int = 40,
) -> None: ...

class EncryptionFlowable(StandardEncryption, Flowable):
    def wrap(self, availWidth, availHeight): ...
    def draw(self) -> None: ...

def encryptDocTemplate(
    dt,
    userPassword,
    ownerPassword=None,
    canPrint: int = 1,
    canModify: int = 1,
    canCopy: int = 1,
    canAnnotate: int = 1,
    strength: int = 40,
) -> None: ...
def encryptPdfInMemory(
    inputPDF,
    userPassword,
    ownerPassword=None,
    canPrint: int = 1,
    canModify: int = 1,
    canCopy: int = 1,
    canAnnotate: int = 1,
    strength: int = 40,
) -> bytes: ...
def encryptPdfOnDisk(
    inputFileName,
    outputFileName,
    userPassword,
    ownerPassword=None,
    canPrint: int = 1,
    canModify: int = 1,
    canCopy: int = 1,
    canAnnotate: int = 1,
    strength: int = 40,
) -> int: ...
def scriptInterp() -> None: ...
def main() -> None: ...
