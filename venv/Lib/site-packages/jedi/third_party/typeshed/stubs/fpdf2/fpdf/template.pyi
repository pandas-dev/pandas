from os import PathLike
from typing import Any

__author__: str
__copyright__: str
__license__: str

class FlexTemplate:
    pdf: Any
    splitting_pdf: Any
    handlers: Any
    texts: Any
    def __init__(self, pdf, elements=None) -> None: ...
    elements: Any
    keys: Any
    def load_elements(self, elements) -> None: ...
    def parse_json(self, infile: PathLike[Any], encoding: str = "utf-8") -> None: ...
    def parse_csv(
        self, infile: PathLike[Any], delimiter: str = ",", decimal_sep: str = ".", encoding: str | None = None
    ) -> None: ...
    def __setitem__(self, name, value) -> None: ...
    set: Any
    def __contains__(self, name): ...
    def __getitem__(self, name): ...
    def split_multicell(self, text: str, element_name: str) -> list[str]: ...
    def render(self, offsetx: float = 0.0, offsety: float = 0.0, rotate: float = 0.0, scale: float = 1.0): ...

class Template(FlexTemplate):
    def __init__(
        self,
        infile=None,
        elements=None,
        format: str = "A4",
        orientation: str = "portrait",
        unit: str = "mm",
        title: str = "",
        author: str = "",
        subject: str = "",
        creator: str = "",
        keywords: str = "",
    ) -> None: ...
    def add_page(self) -> None: ...
    def render(self, outfile=None, dest=None) -> None: ...  # type: ignore[override]
