from _typeshed import StrPath

from PyInstaller.building.datastruct import Target, _TOCTuple

# Referenced in https://pyinstaller.org/en/stable/spec-files.html#example-merge-spec-file
# Not to be imported during runtime, but is the type reference for spec files which are executed as python code
class Splash(Target):
    image_file: str
    full_tk: bool
    tcl_lib: str
    tk_lib: str
    name: str
    script_name: StrPath
    minify_script: bool
    max_img_size: tuple[int, int]
    text_pos: tuple[int, int] | None
    text_size: int
    text_font: str
    text_color: str
    text_default: str
    always_on_top: bool
    uses_tkinter: bool
    script: str
    splash_requirements: set[str]
    binaries: list[_TOCTuple]
    def __init__(
        self,
        image_file: StrPath,
        binaries: list[_TOCTuple],
        datas: list[_TOCTuple],
        *,
        text_pos: tuple[int, int] | None = ...,
        text_size: int = 12,
        text_font: str = ...,
        text_color: str = "black",
        text_default: str = "Initializing",
        full_tk: bool = False,
        minify_script: bool = True,
        name: str = ...,
        script_name: StrPath = ...,
        max_img_size: tuple[int, int] | None = (760, 480),
        always_on_top: bool = True,
    ) -> None: ...
    def assemble(self) -> None: ...
    def generate_script(self) -> str: ...
