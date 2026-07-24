from collections.abc import Callable

from consolemenu.console_menu import ConsoleMenu
from consolemenu.items import MenuItem as MenuItem

class SelectionItem(MenuItem):
    index: int
    def __init__(
        self, text: str | Callable[[], str], index: int, menu: ConsoleMenu | None = None, menu_char: str | None = None
    ) -> None: ...
    def get_return(self) -> int: ...
