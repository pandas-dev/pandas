from _typeshed import StrPath, Unused
from typing import Literal, overload
from typing_extensions import TypeAlias
from xml.etree.ElementTree import ElementTree

from .Menu import Menu, MenuEntry, Separator, XMLMenuBuilder

_MenuItem: TypeAlias = Menu | MenuEntry | Separator

class MenuEditor:
    menu: Menu
    filename: str
    tree: ElementTree
    parser: XMLMenuBuilder
    filenames: list[str]
    def __init__(self, menu: Menu | None = None, filename: StrPath | None = None, root: bool = False) -> None: ...
    def parse(self, menu: Menu | None = None, filename: StrPath | None = None, root: bool = False) -> None: ...
    def save(self) -> None: ...

    # All "before" or "after" items can be one, the other, or neither, but not both.
    @overload
    def createMenuEntry(
        self,
        parent: Menu | None,
        name: str,
        command: str | None = None,
        genericname: str | None = None,
        comment: str | None = None,
        icon: str | None = None,
        terminal: bool | None = None,
        after: _MenuItem | None = None,
        before: None = None,
    ) -> MenuEntry: ...
    @overload
    def createMenuEntry(
        self,
        parent: Menu | None,
        name: str,
        command: str | None = None,
        genericname: str | None = None,
        comment: str | None = None,
        icon: str | None = None,
        terminal: bool | None = None,
        after: None = None,
        before: _MenuItem | None = None,
    ) -> MenuEntry: ...
    @overload
    def createMenu(
        self,
        parent: Menu | None,
        name: str,
        genericname: str | None = None,
        comment: str | None = None,
        icon: str | None = None,
        after: _MenuItem | None = None,
        before: None = None,
    ) -> Menu: ...
    @overload
    def createMenu(
        self,
        parent: Menu | None,
        name: str,
        genericname: str | None = None,
        comment: str | None = None,
        icon: str | None = None,
        after: None = None,
        before: _MenuItem | None = None,
    ) -> Menu: ...
    @overload
    def createSeparator(self, parent: Menu, after: _MenuItem | None = None, before: None = None) -> Separator: ...
    @overload
    def createSeparator(self, parent: Menu, after: None = None, before: _MenuItem | None = None) -> Separator: ...
    @overload
    def moveMenuEntry(
        self,
        menuentry: MenuEntry,
        oldparent: Menu | None,
        newparent: Menu | None,
        after: _MenuItem | None = None,
        before: None = None,
    ) -> MenuEntry: ...
    @overload
    def moveMenuEntry(
        self,
        menuentry: MenuEntry,
        oldparent: Menu | None,
        newparent: Menu | None,
        after: None = None,
        before: _MenuItem | None = None,
    ) -> MenuEntry: ...
    @overload
    def moveMenu(
        self, menu: Menu, oldparent: Menu, newparent: Menu, after: _MenuItem | None = None, before: None = None
    ) -> Menu: ...
    @overload
    def moveMenu(
        self, menu: Menu, oldparent: Menu, newparent: Menu, after: None = None, before: _MenuItem | None = None
    ) -> Menu: ...
    @overload
    def moveSeparator(
        self, separator: Separator, parent: Menu, after: _MenuItem | None = None, before: None = None
    ) -> Separator: ...
    @overload
    def moveSeparator(
        self, separator: Separator, parent: Menu, after: None = None, before: _MenuItem | None = None
    ) -> Separator: ...
    @overload
    def copyMenuEntry(
        self, menuentry: MenuEntry, oldparent: Unused, newparent: Menu, after: _MenuItem | None = None, before: None = None
    ) -> MenuEntry: ...
    @overload
    def copyMenuEntry(
        self, menuentry: MenuEntry, oldparent: Unused, newparent: Menu, after: None = None, before: _MenuItem | None = None
    ) -> MenuEntry: ...
    def editMenuEntry(
        self,
        menuentry: MenuEntry,
        name: str | None = None,
        genericname: str | None = None,
        comment: str | None = None,
        command: str | None = None,
        icon: str | None = None,
        terminal: bool | None = None,
        nodisplay: bool | None = None,
        hidden: bool | None = None,
    ) -> MenuEntry: ...
    def editMenu(
        self,
        menu: Menu,
        name: str | None = None,
        genericname: str | None = None,
        comment: str | None = None,
        icon: str | None = None,
        nodisplay: bool | None = None,
        hidden: bool | None = None,
    ) -> Menu: ...
    def hideMenuEntry(self, menuentry: MenuEntry) -> None: ...
    def unhideMenuEntry(self, menuentry: MenuEntry) -> None: ...
    def hideMenu(self, menu: Menu) -> None: ...
    def unhideMenu(self, menu: Menu) -> None: ...
    def deleteMenuEntry(self, menuentry: MenuEntry) -> MenuEntry: ...
    def revertMenuEntry(self, menuentry: MenuEntry) -> MenuEntry: ...
    def deleteMenu(self, menu: Menu) -> Menu: ...
    def revertMenu(self, menu: Menu) -> Menu: ...
    def deleteSeparator(self, separator: Separator) -> Separator: ...
    def getAction(self, entry: _MenuItem) -> Literal["none", "revert", "delete"]: ...
