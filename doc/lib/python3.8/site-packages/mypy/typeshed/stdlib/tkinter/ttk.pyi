import _tkinter
import sys
import tkinter
from tkinter.font import _FontDescription
from typing import Any, Callable, Dict, List, Optional, Tuple, Union, overload
from typing_extensions import Literal

def tclobjs_to_py(adict): ...
def setup_master(master: Optional[Any] = ...): ...

# from ttk_widget (aka ttk::widget) manual page, differs from tkinter._Compound
_TtkCompound = Literal["text", "image", tkinter._Compound]

class Style:
    master: Any
    tk: _tkinter.TkappType
    def __init__(self, master: Optional[Any] = ...): ...
    def configure(self, style, query_opt: Optional[Any] = ..., **kw): ...
    def map(self, style, query_opt: Optional[Any] = ..., **kw): ...
    def lookup(self, style, option, state: Optional[Any] = ..., default: Optional[Any] = ...): ...
    def layout(self, style, layoutspec: Optional[Any] = ...): ...
    def element_create(self, elementname, etype, *args, **kw): ...
    def element_names(self): ...
    def element_options(self, elementname): ...
    def theme_create(self, themename, parent: Optional[Any] = ..., settings: Optional[Any] = ...): ...
    def theme_settings(self, themename, settings): ...
    def theme_names(self): ...
    def theme_use(self, themename: Optional[Any] = ...): ...

class Widget(tkinter.Widget):
    def __init__(self, master: Optional[tkinter.Misc], widgetname, kw: Optional[Any] = ...): ...
    def identify(self, x, y): ...
    def instate(self, statespec, callback: Optional[Any] = ..., *args, **kw): ...
    def state(self, statespec: Optional[Any] = ...): ...

class Button(Widget):
    def __init__(
        self,
        master: Optional[tkinter.Misc] = ...,
        *,
        class_: str = ...,
        command: tkinter._ButtonCommand = ...,
        compound: _TtkCompound = ...,
        cursor: tkinter._Cursor = ...,
        default: Literal["normal", "active", "disabled"] = ...,
        image: tkinter._ImageSpec = ...,
        name: str = ...,
        padding: Any = ...,  # undocumented
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        text: Union[float, str] = ...,
        textvariable: tkinter.Variable = ...,
        underline: int = ...,
        width: Union[int, Literal[""]] = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: Optional[Dict[str, Any]] = ...,
        *,
        command: tkinter._ButtonCommand = ...,
        compound: _TtkCompound = ...,
        cursor: tkinter._Cursor = ...,
        default: Literal["normal", "active", "disabled"] = ...,
        image: tkinter._ImageSpec = ...,
        padding: Any = ...,
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        text: Union[float, str] = ...,
        textvariable: tkinter.Variable = ...,
        underline: int = ...,
        width: Union[int, Literal[""]] = ...,
    ) -> Optional[Dict[str, Tuple[str, str, str, Any, Any]]]: ...
    @overload
    def configure(self, cnf: str) -> Tuple[str, str, str, Any, Any]: ...
    config = configure
    def invoke(self): ...

class Checkbutton(Widget):
    def __init__(
        self,
        master: Optional[tkinter.Misc] = ...,
        *,
        class_: str = ...,
        command: tkinter._ButtonCommand = ...,
        compound: _TtkCompound = ...,
        cursor: tkinter._Cursor = ...,
        image: tkinter._ImageSpec = ...,
        name: str = ...,
        offvalue: Any = ...,
        onvalue: Any = ...,
        padding: Any = ...,  # undocumented
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        text: Union[float, str] = ...,
        textvariable: tkinter.Variable = ...,
        underline: int = ...,
        # Seems like variable can be empty string, but actually setting it to
        # empty string segfaults before Tcl 8.6.9. Search for ttk::checkbutton
        # here: https://sourceforge.net/projects/tcl/files/Tcl/8.6.9/tcltk-release-notes-8.6.9.txt/view
        variable: tkinter.Variable = ...,
        width: Union[int, Literal[""]] = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: Optional[Dict[str, Any]] = ...,
        *,
        command: tkinter._ButtonCommand = ...,
        compound: _TtkCompound = ...,
        cursor: tkinter._Cursor = ...,
        image: tkinter._ImageSpec = ...,
        offvalue: Any = ...,
        onvalue: Any = ...,
        padding: Any = ...,
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        text: Union[float, str] = ...,
        textvariable: tkinter.Variable = ...,
        underline: int = ...,
        variable: tkinter.Variable = ...,
        width: Union[int, Literal[""]] = ...,
    ) -> Optional[Dict[str, Tuple[str, str, str, Any, Any]]]: ...
    @overload
    def configure(self, cnf: str) -> Tuple[str, str, str, Any, Any]: ...
    config = configure
    def invoke(self): ...

class Entry(Widget, tkinter.Entry):
    def __init__(
        self,
        master: Optional[tkinter.Misc] = ...,
        widget: Optional[str] = ...,
        *,
        background: tkinter._Color = ...,  # undocumented
        class_: str = ...,
        cursor: tkinter._Cursor = ...,
        exportselection: bool = ...,
        font: _FontDescription = ...,
        foreground: tkinter._Color = ...,
        invalidcommand: tkinter._EntryValidateCommand = ...,
        justify: Literal["left", "center", "right"] = ...,
        name: str = ...,
        show: str = ...,
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        textvariable: tkinter.Variable = ...,
        validate: Literal["none", "focus", "focusin", "focusout", "key", "all"] = ...,
        validatecommand: tkinter._EntryValidateCommand = ...,
        width: int = ...,
        xscrollcommand: tkinter._XYScrollCommand = ...,
    ) -> None: ...
    @overload  # type: ignore
    def configure(
        self,
        cnf: Optional[Dict[str, Any]] = ...,
        *,
        background: tkinter._Color = ...,
        cursor: tkinter._Cursor = ...,
        exportselection: bool = ...,
        font: _FontDescription = ...,
        foreground: tkinter._Color = ...,
        invalidcommand: tkinter._EntryValidateCommand = ...,
        justify: Literal["left", "center", "right"] = ...,
        show: str = ...,
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        textvariable: tkinter.Variable = ...,
        validate: Literal["none", "focus", "focusin", "focusout", "key", "all"] = ...,
        validatecommand: tkinter._EntryValidateCommand = ...,
        width: int = ...,
        xscrollcommand: tkinter._XYScrollCommand = ...,
    ) -> Optional[Dict[str, Tuple[str, str, str, Any, Any]]]: ...
    @overload
    def configure(self, cnf: str) -> Tuple[str, str, str, Any, Any]: ...
    # config must be copy/pasted, otherwise ttk.Entry().config is mypy error (don't know why)
    @overload  # type: ignore
    def config(
        self,
        cnf: Optional[Dict[str, Any]] = ...,
        *,
        background: tkinter._Color = ...,
        cursor: tkinter._Cursor = ...,
        exportselection: bool = ...,
        font: _FontDescription = ...,
        foreground: tkinter._Color = ...,
        invalidcommand: tkinter._EntryValidateCommand = ...,
        justify: Literal["left", "center", "right"] = ...,
        show: str = ...,
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        textvariable: tkinter.Variable = ...,
        validate: Literal["none", "focus", "focusin", "focusout", "key", "all"] = ...,
        validatecommand: tkinter._EntryValidateCommand = ...,
        width: int = ...,
        xscrollcommand: tkinter._XYScrollCommand = ...,
    ) -> Optional[Dict[str, Tuple[str, str, str, Any, Any]]]: ...
    @overload
    def config(self, cnf: str) -> Tuple[str, str, str, Any, Any]: ...
    def bbox(self, index): ...
    def identify(self, x, y): ...
    def validate(self): ...

class Combobox(Entry):
    def __init__(
        self,
        master: Optional[tkinter.Misc] = ...,
        *,
        background: tkinter._Color = ...,  # undocumented
        class_: str = ...,
        cursor: tkinter._Cursor = ...,
        exportselection: bool = ...,
        font: _FontDescription = ...,  # undocumented
        foreground: tkinter._Color = ...,  # undocumented
        height: int = ...,
        invalidcommand: tkinter._EntryValidateCommand = ...,  # undocumented
        justify: Literal["left", "center", "right"] = ...,
        name: str = ...,
        postcommand: Union[Callable[[], Any], str] = ...,
        show: Any = ...,  # undocumented
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        textvariable: tkinter.Variable = ...,
        validate: Literal["none", "focus", "focusin", "focusout", "key", "all"] = ...,  # undocumented
        validatecommand: tkinter._EntryValidateCommand = ...,  # undocumented
        values: tkinter._TkinterSequence[str] = ...,
        width: int = ...,
        xscrollcommand: tkinter._XYScrollCommand = ...,  # undocumented
    ) -> None: ...
    @overload  # type: ignore
    def configure(
        self,
        cnf: Optional[Dict[str, Any]] = ...,
        *,
        background: tkinter._Color = ...,
        cursor: tkinter._Cursor = ...,
        exportselection: bool = ...,
        font: _FontDescription = ...,
        foreground: tkinter._Color = ...,
        height: int = ...,
        invalidcommand: tkinter._EntryValidateCommand = ...,
        justify: Literal["left", "center", "right"] = ...,
        postcommand: Union[Callable[[], Any], str] = ...,
        show: Any = ...,
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        textvariable: tkinter.Variable = ...,
        validate: Literal["none", "focus", "focusin", "focusout", "key", "all"] = ...,
        validatecommand: tkinter._EntryValidateCommand = ...,
        values: tkinter._TkinterSequence[str] = ...,
        width: int = ...,
        xscrollcommand: tkinter._XYScrollCommand = ...,
    ) -> Optional[Dict[str, Tuple[str, str, str, Any, Any]]]: ...
    @overload
    def configure(self, cnf: str) -> Tuple[str, str, str, Any, Any]: ...
    # config must be copy/pasted, otherwise ttk.Combobox().config is mypy error (don't know why)
    @overload  # type: ignore
    def config(
        self,
        cnf: Optional[Dict[str, Any]] = ...,
        *,
        background: tkinter._Color = ...,
        cursor: tkinter._Cursor = ...,
        exportselection: bool = ...,
        font: _FontDescription = ...,
        foreground: tkinter._Color = ...,
        height: int = ...,
        invalidcommand: tkinter._EntryValidateCommand = ...,
        justify: Literal["left", "center", "right"] = ...,
        postcommand: Union[Callable[[], Any], str] = ...,
        show: Any = ...,
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        textvariable: tkinter.Variable = ...,
        validate: Literal["none", "focus", "focusin", "focusout", "key", "all"] = ...,
        validatecommand: tkinter._EntryValidateCommand = ...,
        values: tkinter._TkinterSequence[str] = ...,
        width: int = ...,
        xscrollcommand: tkinter._XYScrollCommand = ...,
    ) -> Optional[Dict[str, Tuple[str, str, str, Any, Any]]]: ...
    @overload
    def config(self, cnf: str) -> Tuple[str, str, str, Any, Any]: ...
    def current(self, newindex: Optional[Any] = ...): ...
    def set(self, value): ...

class Frame(Widget):
    def __init__(
        self,
        master: Optional[tkinter.Misc] = ...,
        *,
        border: tkinter._ScreenUnits = ...,
        borderwidth: tkinter._ScreenUnits = ...,
        class_: str = ...,
        cursor: tkinter._Cursor = ...,
        height: tkinter._ScreenUnits = ...,
        name: str = ...,
        padding: tkinter._Padding = ...,
        relief: tkinter._Relief = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        width: tkinter._ScreenUnits = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: Optional[Dict[str, Any]] = ...,
        *,
        border: tkinter._ScreenUnits = ...,
        borderwidth: tkinter._ScreenUnits = ...,
        cursor: tkinter._Cursor = ...,
        height: tkinter._ScreenUnits = ...,
        padding: tkinter._Padding = ...,
        relief: tkinter._Relief = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        width: tkinter._ScreenUnits = ...,
    ) -> Optional[Dict[str, Tuple[str, str, str, Any, Any]]]: ...
    @overload
    def configure(self, cnf: str) -> Tuple[str, str, str, Any, Any]: ...
    config = configure

class Label(Widget):
    def __init__(
        self,
        master: Optional[tkinter.Misc] = ...,
        *,
        anchor: tkinter._Anchor = ...,
        background: tkinter._Color = ...,
        border: tkinter._ScreenUnits = ...,  # alias for borderwidth
        borderwidth: tkinter._ScreenUnits = ...,  # undocumented
        class_: str = ...,
        compound: _TtkCompound = ...,
        cursor: tkinter._Cursor = ...,
        font: _FontDescription = ...,
        foreground: tkinter._Color = ...,
        image: tkinter._ImageSpec = ...,
        justify: Literal["left", "center", "right"] = ...,
        name: str = ...,
        padding: tkinter._Padding = ...,
        relief: tkinter._Relief = ...,
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        text: Union[float, str] = ...,
        textvariable: tkinter.Variable = ...,
        underline: int = ...,
        width: Union[int, Literal[""]] = ...,
        wraplength: tkinter._ScreenUnits = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: Optional[Dict[str, Any]] = ...,
        *,
        anchor: tkinter._Anchor = ...,
        background: tkinter._Color = ...,
        border: tkinter._ScreenUnits = ...,
        borderwidth: tkinter._ScreenUnits = ...,
        compound: _TtkCompound = ...,
        cursor: tkinter._Cursor = ...,
        font: _FontDescription = ...,
        foreground: tkinter._Color = ...,
        image: tkinter._ImageSpec = ...,
        justify: Literal["left", "center", "right"] = ...,
        padding: tkinter._Padding = ...,
        relief: tkinter._Relief = ...,
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        text: Union[float, str] = ...,
        textvariable: tkinter.Variable = ...,
        underline: int = ...,
        width: Union[int, Literal[""]] = ...,
        wraplength: tkinter._ScreenUnits = ...,
    ) -> Optional[Dict[str, Tuple[str, str, str, Any, Any]]]: ...
    @overload
    def configure(self, cnf: str) -> Tuple[str, str, str, Any, Any]: ...
    config = configure

class Labelframe(Widget):
    def __init__(
        self,
        master: Optional[tkinter.Misc] = ...,
        *,
        border: tkinter._ScreenUnits = ...,
        borderwidth: tkinter._ScreenUnits = ...,  # undocumented
        class_: str = ...,
        cursor: tkinter._Cursor = ...,
        height: tkinter._ScreenUnits = ...,
        labelanchor: Literal["nw", "n", "ne", "en", "e", "es", "se", "s", "sw", "ws", "w", "wn"] = ...,
        labelwidget: tkinter.Misc = ...,
        name: str = ...,
        padding: tkinter._Padding = ...,
        relief: tkinter._Relief = ...,  # undocumented
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        text: Union[float, str] = ...,
        underline: int = ...,
        width: tkinter._ScreenUnits = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: Optional[Dict[str, Any]] = ...,
        *,
        border: tkinter._ScreenUnits = ...,
        borderwidth: tkinter._ScreenUnits = ...,
        cursor: tkinter._Cursor = ...,
        height: tkinter._ScreenUnits = ...,
        labelanchor: Literal["nw", "n", "ne", "en", "e", "es", "se", "s", "sw", "ws", "w", "wn"] = ...,
        labelwidget: tkinter.Misc = ...,
        padding: tkinter._Padding = ...,
        relief: tkinter._Relief = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        text: Union[float, str] = ...,
        underline: int = ...,
        width: tkinter._ScreenUnits = ...,
    ) -> Optional[Dict[str, Tuple[str, str, str, Any, Any]]]: ...
    @overload
    def configure(self, cnf: str) -> Tuple[str, str, str, Any, Any]: ...
    config = configure

LabelFrame = Labelframe

class Menubutton(Widget):
    def __init__(
        self,
        master: Optional[tkinter.Misc] = ...,
        *,
        class_: str = ...,
        compound: _TtkCompound = ...,
        cursor: tkinter._Cursor = ...,
        direction: Literal["above", "below", "left", "right", "flush"] = ...,
        image: tkinter._ImageSpec = ...,
        menu: tkinter.Menu = ...,
        name: str = ...,
        padding: Any = ...,  # undocumented
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        text: Union[float, str] = ...,
        textvariable: tkinter.Variable = ...,
        underline: int = ...,
        width: Union[int, Literal[""]] = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: Optional[Dict[str, Any]] = ...,
        *,
        compound: _TtkCompound = ...,
        cursor: tkinter._Cursor = ...,
        direction: Literal["above", "below", "left", "right", "flush"] = ...,
        image: tkinter._ImageSpec = ...,
        menu: tkinter.Menu = ...,
        padding: Any = ...,
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        text: Union[float, str] = ...,
        textvariable: tkinter.Variable = ...,
        underline: int = ...,
        width: Union[int, Literal[""]] = ...,
    ) -> Optional[Dict[str, Tuple[str, str, str, Any, Any]]]: ...
    @overload
    def configure(self, cnf: str) -> Tuple[str, str, str, Any, Any]: ...
    config = configure

class Notebook(Widget):
    def __init__(
        self,
        master: Optional[tkinter.Misc] = ...,
        *,
        class_: str = ...,
        cursor: tkinter._Cursor = ...,
        height: int = ...,
        name: str = ...,
        padding: tkinter._Padding = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        width: int = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: Optional[Dict[str, Any]] = ...,
        *,
        cursor: tkinter._Cursor = ...,
        height: int = ...,
        padding: tkinter._Padding = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        width: int = ...,
    ) -> Optional[Dict[str, Tuple[str, str, str, Any, Any]]]: ...
    @overload
    def configure(self, cnf: str) -> Tuple[str, str, str, Any, Any]: ...
    config = configure
    def add(self, child, **kw): ...
    def forget(self, tab_id): ...
    def hide(self, tab_id): ...
    def identify(self, x, y): ...
    def index(self, tab_id): ...
    def insert(self, pos, child, **kw): ...
    def select(self, tab_id: Optional[Any] = ...): ...
    def tab(self, tab_id, option: Optional[Any] = ..., **kw): ...
    def tabs(self): ...
    def enable_traversal(self): ...

class Panedwindow(Widget, tkinter.PanedWindow):
    def __init__(
        self,
        master: Optional[tkinter.Misc] = ...,
        *,
        class_: str = ...,
        cursor: tkinter._Cursor = ...,
        # width and height for tkinter.ttk.Panedwindow are int but for tkinter.PanedWindow they are screen units
        height: int = ...,
        name: str = ...,
        orient: Literal["vertical", "horizontal"] = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        width: int = ...,
    ) -> None: ...
    @overload  # type: ignore
    def configure(
        self,
        cnf: Optional[Dict[str, Any]] = ...,
        *,
        cursor: tkinter._Cursor = ...,
        height: int = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        width: int = ...,
    ) -> Optional[Dict[str, Tuple[str, str, str, Any, Any]]]: ...
    @overload
    def configure(self, cnf: str) -> Tuple[str, str, str, Any, Any]: ...
    # config must be copy/pasted, otherwise ttk.Panedwindow().config is mypy error (don't know why)
    @overload  # type: ignore
    def config(
        self,
        cnf: Optional[Dict[str, Any]] = ...,
        *,
        cursor: tkinter._Cursor = ...,
        height: int = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        width: int = ...,
    ) -> Optional[Dict[str, Tuple[str, str, str, Any, Any]]]: ...
    @overload
    def config(self, cnf: str) -> Tuple[str, str, str, Any, Any]: ...
    forget: Any
    def insert(self, pos, child, **kw): ...
    def pane(self, pane, option: Optional[Any] = ..., **kw): ...
    def sashpos(self, index, newpos: Optional[Any] = ...): ...

PanedWindow = Panedwindow

class Progressbar(Widget):
    def __init__(
        self,
        master: Optional[tkinter.Misc] = ...,
        *,
        class_: str = ...,
        cursor: tkinter._Cursor = ...,
        length: tkinter._ScreenUnits = ...,
        maximum: float = ...,
        mode: Literal["determinate", "indeterminate"] = ...,
        name: str = ...,
        orient: Literal["horizontal", "vertical"] = ...,
        phase: int = ...,  # docs say read-only but assigning int to this works
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        value: float = ...,
        variable: Union[tkinter.IntVar, tkinter.DoubleVar] = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: Optional[Dict[str, Any]] = ...,
        *,
        cursor: tkinter._Cursor = ...,
        length: tkinter._ScreenUnits = ...,
        maximum: float = ...,
        mode: Literal["determinate", "indeterminate"] = ...,
        orient: Literal["horizontal", "vertical"] = ...,
        phase: int = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        value: float = ...,
        variable: Union[tkinter.IntVar, tkinter.DoubleVar] = ...,
    ) -> Optional[Dict[str, Tuple[str, str, str, Any, Any]]]: ...
    @overload
    def configure(self, cnf: str) -> Tuple[str, str, str, Any, Any]: ...
    config = configure
    def start(self, interval: Optional[Any] = ...): ...
    def step(self, amount: Optional[Any] = ...): ...
    def stop(self): ...

class Radiobutton(Widget):
    def __init__(
        self,
        master: Optional[tkinter.Misc] = ...,
        *,
        class_: str = ...,
        command: tkinter._ButtonCommand = ...,
        compound: _TtkCompound = ...,
        cursor: tkinter._Cursor = ...,
        image: tkinter._ImageSpec = ...,
        name: str = ...,
        padding: Any = ...,  # undocumented
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        text: Union[float, str] = ...,
        textvariable: tkinter.Variable = ...,
        underline: int = ...,
        value: Any = ...,
        variable: Union[tkinter.Variable, Literal[""]] = ...,
        width: Union[int, Literal[""]] = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: Optional[Dict[str, Any]] = ...,
        *,
        command: tkinter._ButtonCommand = ...,
        compound: _TtkCompound = ...,
        cursor: tkinter._Cursor = ...,
        image: tkinter._ImageSpec = ...,
        padding: Any = ...,
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        text: Union[float, str] = ...,
        textvariable: tkinter.Variable = ...,
        underline: int = ...,
        value: Any = ...,
        variable: Union[tkinter.Variable, Literal[""]] = ...,
        width: Union[int, Literal[""]] = ...,
    ) -> Optional[Dict[str, Tuple[str, str, str, Any, Any]]]: ...
    @overload
    def configure(self, cnf: str) -> Tuple[str, str, str, Any, Any]: ...
    config = configure
    def invoke(self): ...

class Scale(Widget, tkinter.Scale):
    def __init__(
        self,
        master: Optional[tkinter.Misc] = ...,
        *,
        class_: str = ...,
        command: Union[str, Callable[[str], Any]] = ...,
        cursor: tkinter._Cursor = ...,
        from_: float = ...,
        length: tkinter._ScreenUnits = ...,
        name: str = ...,
        orient: Literal["horizontal", "vertical"] = ...,
        state: str = ...,  # undocumented
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        to: float = ...,
        value: float = ...,
        variable: Union[tkinter.IntVar, tkinter.DoubleVar] = ...,
    ) -> None: ...
    @overload  # type: ignore
    def configure(
        self,
        cnf: Optional[Dict[str, Any]] = ...,
        *,
        command: Union[str, Callable[[str], Any]] = ...,
        cursor: tkinter._Cursor = ...,
        from_: float = ...,
        length: tkinter._ScreenUnits = ...,
        orient: Literal["horizontal", "vertical"] = ...,
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        to: float = ...,
        value: float = ...,
        variable: Union[tkinter.IntVar, tkinter.DoubleVar] = ...,
    ) -> Optional[Dict[str, Tuple[str, str, str, Any, Any]]]: ...
    @overload
    def configure(self, cnf: str) -> Tuple[str, str, str, Any, Any]: ...
    # config must be copy/pasted, otherwise ttk.Scale().config is mypy error (don't know why)
    @overload  # type: ignore
    def config(
        self,
        cnf: Optional[Dict[str, Any]] = ...,
        *,
        command: Union[str, Callable[[str], Any]] = ...,
        cursor: tkinter._Cursor = ...,
        from_: float = ...,
        length: tkinter._ScreenUnits = ...,
        orient: Literal["horizontal", "vertical"] = ...,
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        to: float = ...,
        value: float = ...,
        variable: Union[tkinter.IntVar, tkinter.DoubleVar] = ...,
    ) -> Optional[Dict[str, Tuple[str, str, str, Any, Any]]]: ...
    @overload
    def config(self, cnf: str) -> Tuple[str, str, str, Any, Any]: ...
    def get(self, x: Optional[Any] = ..., y: Optional[Any] = ...): ...

class Scrollbar(Widget, tkinter.Scrollbar):
    def __init__(
        self,
        master: Optional[tkinter.Misc] = ...,
        *,
        class_: str = ...,
        command: Union[Callable[..., Optional[Tuple[float, float]]], str] = ...,
        cursor: tkinter._Cursor = ...,
        name: str = ...,
        orient: Literal["horizontal", "vertical"] = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
    ) -> None: ...
    @overload  # type: ignore
    def configure(
        self,
        cnf: Optional[Dict[str, Any]] = ...,
        *,
        command: Union[Callable[..., Optional[Tuple[float, float]]], str] = ...,
        cursor: tkinter._Cursor = ...,
        orient: Literal["horizontal", "vertical"] = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
    ) -> Optional[Dict[str, Tuple[str, str, str, Any, Any]]]: ...
    @overload
    def configure(self, cnf: str) -> Tuple[str, str, str, Any, Any]: ...
    # config must be copy/pasted, otherwise ttk.Scrollbar().config is mypy error (don't know why)
    @overload  # type: ignore
    def config(
        self,
        cnf: Optional[Dict[str, Any]] = ...,
        *,
        command: Union[Callable[..., Optional[Tuple[float, float]]], str] = ...,
        cursor: tkinter._Cursor = ...,
        orient: Literal["horizontal", "vertical"] = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
    ) -> Optional[Dict[str, Tuple[str, str, str, Any, Any]]]: ...
    @overload
    def config(self, cnf: str) -> Tuple[str, str, str, Any, Any]: ...

class Separator(Widget):
    def __init__(
        self,
        master: Optional[tkinter.Misc] = ...,
        *,
        class_: str = ...,
        cursor: tkinter._Cursor = ...,
        name: str = ...,
        orient: Literal["horizontal", "vertical"] = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: Optional[Dict[str, Any]] = ...,
        *,
        cursor: tkinter._Cursor = ...,
        orient: Literal["horizontal", "vertical"] = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
    ) -> Optional[Dict[str, Tuple[str, str, str, Any, Any]]]: ...
    @overload
    def configure(self, cnf: str) -> Tuple[str, str, str, Any, Any]: ...
    config = configure

class Sizegrip(Widget):
    def __init__(
        self,
        master: Optional[tkinter.Misc] = ...,
        *,
        class_: str = ...,
        cursor: tkinter._Cursor = ...,
        name: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: Optional[Dict[str, Any]] = ...,
        *,
        cursor: tkinter._Cursor = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
    ) -> Optional[Dict[str, Tuple[str, str, str, Any, Any]]]: ...
    @overload
    def configure(self, cnf: str) -> Tuple[str, str, str, Any, Any]: ...
    config = configure

if sys.version_info >= (3, 7):
    class Spinbox(Entry):
        def __init__(
            self,
            master: Optional[tkinter.Misc] = ...,
            *,
            background: tkinter._Color = ...,  # undocumented
            class_: str = ...,
            command: Union[Callable[[], Any], str, tkinter._TkinterSequence[str]] = ...,
            cursor: tkinter._Cursor = ...,
            exportselection: bool = ...,  # undocumented
            font: _FontDescription = ...,  # undocumented
            foreground: tkinter._Color = ...,  # undocumented
            format: str = ...,
            from_: float = ...,
            increment: float = ...,
            invalidcommand: tkinter._EntryValidateCommand = ...,  # undocumented
            justify: Literal["left", "center", "right"] = ...,  # undocumented
            name: str = ...,
            show: Any = ...,  # undocumented
            state: str = ...,
            style: str = ...,
            takefocus: tkinter._TakeFocusValue = ...,
            textvariable: tkinter.Variable = ...,  # undocumented
            to: float = ...,
            validate: Literal["none", "focus", "focusin", "focusout", "key", "all"] = ...,
            validatecommand: tkinter._EntryValidateCommand = ...,
            values: tkinter._TkinterSequence[str] = ...,
            width: int = ...,  # undocumented
            wrap: bool = ...,
            xscrollcommand: tkinter._XYScrollCommand = ...,
        ) -> None: ...
        @overload  # type: ignore
        def configure(
            self,
            cnf: Optional[Dict[str, Any]] = ...,
            *,
            background: tkinter._Color = ...,
            command: Union[Callable[[], Any], str, tkinter._TkinterSequence[str]] = ...,
            cursor: tkinter._Cursor = ...,
            exportselection: bool = ...,
            font: _FontDescription = ...,
            foreground: tkinter._Color = ...,
            format: str = ...,
            from_: float = ...,
            increment: float = ...,
            invalidcommand: tkinter._EntryValidateCommand = ...,
            justify: Literal["left", "center", "right"] = ...,
            show: Any = ...,
            state: str = ...,
            style: str = ...,
            takefocus: tkinter._TakeFocusValue = ...,
            textvariable: tkinter.Variable = ...,
            to: float = ...,
            validate: Literal["none", "focus", "focusin", "focusout", "key", "all"] = ...,
            validatecommand: tkinter._EntryValidateCommand = ...,
            values: tkinter._TkinterSequence[str] = ...,
            width: int = ...,
            wrap: bool = ...,
            xscrollcommand: tkinter._XYScrollCommand = ...,
        ) -> Optional[Dict[str, Tuple[str, str, str, Any, Any]]]: ...
        @overload
        def configure(self, cnf: str) -> Tuple[str, str, str, Any, Any]: ...
        config = configure  # type: ignore
        def set(self, value: Any) -> None: ...

class Treeview(Widget, tkinter.XView, tkinter.YView):
    def __init__(
        self,
        master: Optional[tkinter.Misc] = ...,
        *,
        class_: str = ...,
        columns: Union[str, tkinter._TkinterSequence[str]] = ...,
        cursor: tkinter._Cursor = ...,
        displaycolumns: Union[str, tkinter._TkinterSequence[str], tkinter._TkinterSequence[int], Literal["#all"]] = ...,
        height: int = ...,
        name: str = ...,
        padding: tkinter._Padding = ...,
        selectmode: Literal["extended", "browse", "none"] = ...,
        # _TkinterSequences of Literal don't actually work, using str instead.
        #
        # 'tree headings' is same as ['tree', 'headings'], and I wouldn't be
        # surprised if someone was using it.
        show: Union[Literal["tree", "headings", "tree headings"], tkinter._TkinterSequence[str]] = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        xscrollcommand: tkinter._XYScrollCommand = ...,
        yscrollcommand: tkinter._XYScrollCommand = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: Optional[Dict[str, Any]] = ...,
        *,
        columns: Union[str, tkinter._TkinterSequence[str]] = ...,
        cursor: tkinter._Cursor = ...,
        displaycolumns: Union[str, tkinter._TkinterSequence[str], tkinter._TkinterSequence[int], Literal["#all"]] = ...,
        height: int = ...,
        padding: tkinter._Padding = ...,
        selectmode: Literal["extended", "browse", "none"] = ...,
        show: Union[Literal["tree", "headings", "tree headings"], tkinter._TkinterSequence[str]] = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        xscrollcommand: tkinter._XYScrollCommand = ...,
        yscrollcommand: tkinter._XYScrollCommand = ...,
    ) -> Optional[Dict[str, Tuple[str, str, str, Any, Any]]]: ...
    @overload
    def configure(self, cnf: str) -> Tuple[str, str, str, Any, Any]: ...
    config = configure
    def bbox(self, item, column: Optional[Any] = ...): ...  # type: ignore
    def get_children(self, item: Optional[Any] = ...): ...
    def set_children(self, item, *newchildren): ...
    def column(self, column, option: Optional[Any] = ..., **kw): ...
    def delete(self, *items): ...
    def detach(self, *items): ...
    def exists(self, item): ...
    def focus(self, item: Optional[Any] = ...): ...
    def heading(self, column, option: Optional[Any] = ..., **kw): ...
    def identify(self, component, x, y): ...
    def identify_row(self, y): ...
    def identify_column(self, x): ...
    def identify_region(self, x, y): ...
    def identify_element(self, x, y): ...
    def index(self, item): ...
    def insert(self, parent, index, iid: Optional[Any] = ..., **kw): ...
    def item(self, item, option: Optional[Any] = ..., **kw): ...
    def move(self, item, parent, index): ...
    reattach: Any
    def next(self, item): ...
    def parent(self, item): ...
    def prev(self, item): ...
    def see(self, item): ...
    if sys.version_info >= (3, 8):
        def selection(self) -> Tuple[str, ...]: ...
    else:
        def selection(self, selop: Optional[Any] = ..., items: Optional[Any] = ...) -> Tuple[str, ...]: ...
    def selection_set(self, items): ...
    def selection_add(self, items): ...
    def selection_remove(self, items): ...
    def selection_toggle(self, items): ...
    def set(self, item, column: Optional[Any] = ..., value: Optional[Any] = ...): ...
    # There's no tag_unbind() or 'add' argument for whatever reason.
    # Also, it's 'callback' instead of 'func' here.
    @overload
    def tag_bind(
        self, tagname: str, sequence: Optional[str] = ..., callback: Optional[Callable[[tkinter.Event[Treeview]], Any]] = ...
    ) -> str: ...
    @overload
    def tag_bind(self, tagname: str, sequence: Optional[str], callback: str) -> None: ...
    @overload
    def tag_bind(self, tagname: str, *, callback: str) -> None: ...
    def tag_configure(self, tagname, option: Optional[Any] = ..., **kw): ...
    def tag_has(self, tagname, item: Optional[Any] = ...): ...

class LabeledScale(Frame):
    label: Any
    scale: Any
    # TODO: don't any-type **kw. That goes to Frame.__init__.
    def __init__(
        self,
        master: Optional[tkinter.Misc] = ...,
        variable: Optional[Union[tkinter.IntVar, tkinter.DoubleVar]] = ...,
        from_: float = ...,
        to: float = ...,
        *,
        compound: Union[Literal["top"], Literal["bottom"]] = ...,
        **kw: Any,
    ) -> None: ...
    # destroy is overrided, signature does not change
    value: Any

class OptionMenu(Menubutton):
    def __init__(
        self,
        master,
        variable,
        default: Optional[str] = ...,
        *values: str,
        # rest of these are keyword-only because *args syntax used above
        style: str = ...,
        direction: Union[Literal["above"], Literal["below"], Literal["left"], Literal["right"], Literal["flush"]] = ...,
        command: Optional[Callable[[tkinter.StringVar], Any]] = ...,
    ) -> None: ...
    # configure, config, cget, destroy are inherited from Menubutton
    # destroy and __setitem__ are overrided, signature does not change
    def set_menu(self, default: Optional[Any] = ..., *values): ...
