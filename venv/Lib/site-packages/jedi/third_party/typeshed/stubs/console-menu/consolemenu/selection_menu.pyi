from collections.abc import Iterable

from consolemenu import ConsoleMenu as ConsoleMenu
from consolemenu.items import SelectionItem as SelectionItem
from consolemenu.menu_formatter import MenuFormatBuilder
from consolemenu.screen import Screen

class SelectionMenu(ConsoleMenu):
    def __init__(
        self,
        strings: Iterable[str],
        title: str | None = None,
        subtitle: str | None = None,
        screen: Screen | None = None,
        formatter: MenuFormatBuilder | None = None,
        prologue_text: str | None = None,
        epilogue_text: str | None = None,
        show_exit_option: bool = True,
        exit_option_text: str = "Exit",
        clear_screen: bool = True,
    ) -> None: ...
    @classmethod
    def get_selection(
        cls,
        strings: Iterable[str],
        title: str = "Select an option",
        subtitle: str | None = None,
        show_exit_option: bool = True,
        _menu: ConsoleMenu | None = None,
    ) -> int: ...
    def append_string(self, string: str) -> None: ...
