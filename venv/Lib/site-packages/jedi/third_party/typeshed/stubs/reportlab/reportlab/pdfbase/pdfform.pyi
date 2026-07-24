from typing import Literal

from reportlab.pdfbase.pdfdoc import PDFDictionary, PDFObject, PDFStream, PDFString
from reportlab.pdfbase.pdfpattern import PDFPattern, PDFPatternIf

def textFieldAbsolute(
    canvas, title, x, y, width, height, value: str = "", maxlen: int = 1000000, multiline: bool | Literal[0, 1] = 0
) -> None: ...
def textFieldRelative(
    canvas, title, xR, yR, width, height, value: str = "", maxlen: int = 1000000, multiline: bool | Literal[0, 1] = 0
) -> None: ...
def buttonFieldAbsolute(canvas, title, value, x, y, width: float = 16.7704, height: float = 14.907) -> None: ...
def buttonFieldRelative(canvas, title, value, xR, yR, width: float = 16.7704, height: float = 14.907) -> None: ...
def selectFieldAbsolute(canvas, title, value, options, x, y, width, height) -> None: ...
def selectFieldRelative(canvas, title, value, options, xR, yR, width, height) -> None: ...
def getForm(canvas) -> AcroForm: ...

class AcroForm(PDFObject):
    fields: list[PDFPattern]
    def __init__(self) -> None: ...
    def textField(
        self, canvas, title, xmin, ymin, xmax, ymax, value: str = "", maxlen: int = 1000000, multiline: bool | Literal[0, 1] = 0
    ) -> None: ...
    def selectField(self, canvas, title, value, options, xmin, ymin, xmax, ymax) -> None: ...
    def buttonField(self, canvas, title, value, xmin, ymin, width: float = 16.7704, height: float = 14.907) -> None: ...
    def format(self, document) -> bytes: ...

FormPattern: list[str | list[str] | PDFString | PDFPatternIf]

def FormFontsDictionary() -> PDFDictionary: ...
def FormResources() -> PDFPattern: ...

ZaDbPattern: list[str]
FormResourcesDictionaryPattern: list[str | list[str]]
FORMFONTNAMES: dict[str, str]
EncodingPattern: list[str | list[str]]
PDFDocEncodingPattern: list[str]

def FormFont(BaseFont, Name) -> PDFPattern: ...

FormFontPattern: list[str | list[str]]

def resetPdfForm() -> None: ...
def TextField(
    title,
    value,
    xmin,
    ymin,
    xmax,
    ymax,
    page,
    maxlen: int = 1000000,
    font: str = "Helvetica-Bold",
    fontsize: int = 9,
    R: int = 0,
    G: int = 0,
    B: float = 0.627,
    multiline: bool | Literal[0, 1] = 0,
) -> PDFPattern: ...

TextFieldPattern: list[str | list[str]]

def SelectField(
    title,
    value,
    options,
    xmin,
    ymin,
    xmax,
    ymax,
    page,
    font: str = "Helvetica-Bold",
    fontsize: int = 9,
    R: int = 0,
    G: int = 0,
    B: float = 0.627,
) -> PDFPattern: ...

SelectFieldPattern: list[str | list[str]]

def ButtonField(title, value, xmin, ymin, page, width: float = 16.7704, height: float = 14.907) -> PDFPattern: ...

ButtonFieldPattern: list[str | list[str] | PDFString]

def buttonStreamDictionary(width: float = 16.7704, height: float = 14.907) -> PDFDictionary: ...
def ButtonStream(content, width: float = 16.7704, height: float = 14.907) -> PDFStream: ...
